program yelmo_calving

    use nml 
    use ncio  
    use yelmo 
    use lsf_module
    use yelmo_tools, only : get_region_indices
    use topography, only: calc_ice_fraction_new

    use calving_benchmarks
    
    implicit none

    type(yelmo_class) :: yelmo1
    type(yelmo_class) :: yelmo_ref

    type control_type
        character(len=256) :: outfldr
        character(len=256) :: file2D, file1D
        character(len=256) :: file_restart
        character(len=256) :: file_cmip
        character(len=512) :: path_par
        character(len=56)  :: exp

        real(wp) :: time_init, time_end, time, dtt
        real(wp) :: dt2D_out, dt1D_out

        real(wp) :: dx

        logical  :: calvingmip_out

        ! Internal parameters
        character(len=56)  :: domain
        character(len=56)  :: grid_name
        real(wp) :: x0, x1
        integer  :: nx
        integer  :: ny
    end type

    type(control_type) :: ctl

    ! CalvingMIP profile geometry (radial A-H for circular; CapA-D, HalA-D for Thule).
    ! Sample spacing is set to the native model dx so the profiles trace the
    ! resolved fields rather than introducing a coarser interpolation grid.
    type profile_t
        character(len=8)      :: name      ! e.g. "A", "CapA"
        integer               :: n         ! number of sample points along the line
        real(wp), allocatable :: s(:)      ! arc-length distance from start [m]
        real(wp), allocatable :: x(:)      ! sample x coords [m]
        real(wp), allocatable :: y(:)      ! sample y coords [m]
    end type
    type(profile_t), allocatable :: profiles(:)

    real(wp) :: time
    integer  :: n
    
    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)

    ! Assume program is running from the output folder
    ctl%outfldr = "./"

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(ctl%path_par)
    !path_par   = trim(outfldr)//"yelmo_calving.nml" 

    ! Define input and output locations
    ctl%file1D       = trim(ctl%outfldr)//"yelmo1D.nc"
    ctl%file2D       = trim(ctl%outfldr)//"yelmo2D.nc"
    ctl%file_restart = trim(ctl%outfldr)//"yelmo_restart.nc"


    ! Define the domain, grid and experiment from parameter file
    call nml_read(ctl%path_par,"ctl","exp",         ctl%exp)            ! "exp1", "exp2", "exp3", "exp4", "exp5"
    call nml_read(ctl%path_par,"ctl","dx",          ctl%dx)             ! [km] Grid resolution

    ! Timing parameters
    call nml_read(ctl%path_par,"ctl","time_init",   ctl%time_init)      ! [yr] Starting time
    call nml_read(ctl%path_par,"ctl","time_end",    ctl%time_end)       ! [yr] Ending time
    call nml_read(ctl%path_par,"ctl","dtt",         ctl%dtt)            ! [yr] Main loop time step
    call nml_read(ctl%path_par,"ctl","dt2D_out",    ctl%dt2D_out)       ! [yr] Frequency of 2D output
    ctl%dt1D_out = ctl%dtt  ! Set 1D output to frequency of main loop timestep

    ! CalvingMIP output flag (default off if not found)
    ctl%calvingmip_out = .FALSE.
    call nml_read(ctl%path_par,"ctl","calvingmip_out",ctl%calvingmip_out)

    ! Build CalvingMIP output filename: CalvingMIP_EXP{N}_YELMO_AWI.nc
    ! ctl%exp is expected to be "exp1", "exp2", ..., so uppercase by extracting digit.
    ctl%file_cmip = trim(ctl%outfldr)//"CalvingMIP_EXP"//trim(ctl%exp(4:))//"_YELMO_AWI.nc"

    ! Now set internal parameters ===

    ! Define domain and grid size based on experiment
    select case(trim(ctl%exp))
        case("exp1","exp2")
            ctl%domain = "circular"
            ctl%x0 = -800.0
            ctl%x1 =  800.0
        case("exp3","exp4","exp5")
            ctl%domain = "thule"
            ctl%x0 = -800.0
            ctl%x1 =  800.0
        case("advection")
            ctl%domain = "advection"
            ctl%x0     = -800.0
            ctl%x1     = 800.0
        case DEFAULT
            write(*,*) "ctl.exp = ",trim(ctl%domain), " not recognized."
            stop
    end select

    ! Get grid size
    ctl%nx = (ctl%x1-ctl%x0) / ctl%dx + 1
    ctl%ny = ctl%nx

    ! Get grid name
    write(ctl%grid_name,"(a,i2,a2)") trim(ctl%domain)//"-",int(ctl%dx),"KM"
   
    ! === Initialize ice sheet model =====

    ! First, define grid 
    call yelmo_init_grid(yelmo1%grd,ctl%grid_name,units="km",dx=ctl%dx,nx=ctl%nx,dy=ctl%dx,ny=ctl%nx)

    ! Initialize data objects (without loading topography, which will be defined inline below)
    call yelmo_init(yelmo1,filename=ctl%path_par,grid_def="none",time=ctl%time_init, &
                    load_topo=.FALSE.,domain=ctl%domain,grid_name=ctl%grid_name)

    ! === Define initial topography ===
    call calvmip_init(yelmo1%bnd%z_bed,yelmo1%grd%x,yelmo1%grd%y,yelmo1%par%domain)

    ! advection test
    if (.not. yelmo1%par%use_restart) then
        ! If no restart, set ice thickness to zero
        yelmo1%tpo%now%H_ice = 0.0
        yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed 
        select case(trim(ctl%exp))
            case("advection")
            call CircularDomain(yelmo1%tpo%now%lsf,yelmo1%bnd%z_bed,yelmo1%tpo%par%dx)
        case DEFAULT 
            call LSFinit(yelmo1%tpo%now%lsf,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%tpo%par%dx)
        end select
    end if

    ! === Define additional boundary conditions =====

    yelmo1%bnd%z_sl     = 0.0
    yelmo1%bnd%bmb_shlf = 0.0  
    yelmo1%bnd%T_shlf   = yelmo1%bnd%c%T0  
    yelmo1%bnd%H_sed    = 0.0 

    yelmo1%bnd%T_srf    = 223.15 
    yelmo1%bnd%Q_geo    = 42.0

    select case(trim(ctl%exp))
        case("advection")
            yelmo1%bnd%smb      = 0.0
        case DEFAULT 
            yelmo1%bnd%smb      = 0.3
    end select    

    ! Check boundary values 
    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize state variables (dyn,therm,mat)
    call yelmo_init_state(yelmo1,time=ctl%time_init,thrm_method="robin")

    ! == Write initial state ==

    ! Standard yelmo 2D output
    call yelmo_write_init(yelmo1,ctl%file2D,time_init=ctl%time_init,units="years")
    call yelmo_write_step(yelmo1,ctl%file2D,time=ctl%time_init)

    ! Standard yelmo 1D regional output
    call yelmo_write_reg_init(yelmo1,ctl%file1D,time_init=ctl%time_init,units="years",mask=(yelmo1%bnd%mask_ice /= -1))
    call yelmo_write_reg_step(yelmo1,ctl%file1D,time=ctl%time_init)

    ! Optional CalvingMIP output (wiki-spec file with Time1 / Time100 dims)
    if (ctl%calvingmip_out) then
        call calvingmip_init(yelmo1,ctl)
        ! Initial state at t=time_init. The dispatcher decides per-exp whether
        ! this triggers a Time1 entry, a Time100 entry, both, or neither.
        call calvingmip_write_step(yelmo1,ctl,ctl%time_init)
    end if

    ! Store default parameters
    yelmo_ref = yelmo1

    ! Set calving mask if needed
    if (trim(yelmo1%tpo%par%calv_flt_method) .eq. "kill-pos") then
        call set_calving_mask(yelmo1%bnd%calv_mask,yelmo1%grd%x,yelmo1%grd%y,r_lim=750e3_wp)
    end if

    ! Advance timesteps
    do n = 1, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

        ! Get current time 
        time = ctl%time_init + n*ctl%dtt
        
        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)
        
        ! == MODEL OUTPUT =======================================================

        ! Standard yelmo 2D output at user-specified frequency
        if (mod(nint(time*100),nint(ctl%dt2D_out*100))==0) then
            call yelmo_write_step(yelmo1,ctl%file2D,time=time)
        end if

        ! Standard yelmo 1D regional output at main-loop timestep
        if (mod(nint(time*100),nint(ctl%dt1D_out*100))==0) then
            call yelmo_write_reg_step(yelmo1,ctl%file1D,time=time)
        end if

        ! CalvingMIP output (hardcoded frequencies per wiki experiment spec)
        if (ctl%calvingmip_out) then
            call calvingmip_write_step(yelmo1,ctl,time)
        end if

    end do

    ! Write a restart file too
    call yelmo_restart_write(yelmo1,ctl%file_restart,time=time)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctl%time_end-ctl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
    
contains

    subroutine CircularDomain(LSF,zbed,dx)
        
        implicit none
    
        real(wp), intent(OUT) :: LSF(:,:)      ! LSF mask
        real(wp), intent(IN)  :: zbed(:,:)    
        real(wp), intent(IN)  :: dx            ! Model resolution [m]
        
        ! Internal variables
        real(wp) :: rc
        integer  :: i,j,nx,ny
    
        nx = size(zbed,1)
        ny = size(zbed,2)
        rc = 10.0_wp ! grid points below zero
    
        do j=1,ny
        do i=1,nx
    
        LSF(i,j) = (sqrt((0.5*(nx+1)-i)**2 + (0.5*(ny+1)-j)**2) - rc)*dx*1e-3 
    
        end do
        end do
    
        return
    
    end subroutine CircularDomain

    ! ============================================================================
    ! CalvingMIP output: separate file with wiki-spec variables, dims Time1/Time100
    ! File name: CalvingMIP_EXP{N}_YELMO_AWI.nc
    ! Spec: https://github.com/JRowanJordan/CalvingMIP/wiki
    ! ============================================================================

    subroutine calvingmip_init(ylmo,ctl)
        ! Create the CalvingMIP output file and define dimensions:
        !   X, Y    — native model grid in metres
        !   Time1   — unlimited, annual scalars/profiles (grown by nc_time_index)
        !   Time100 — fixed-size, 100-yearly snapshots; pre-populated time values
        !             matching the per-exp schedule so nc_time_index can match.
        ! Note: NetCDF3 (the default for ncio nc_create) allows only one
        ! unlimited dim per file, so Time100 is fixed-size.

        implicit none

        type(yelmo_class),  intent(IN) :: ylmo
        type(control_type), intent(IN) :: ctl

        integer  :: i, n_time100
        real(wp), allocatable :: t100(:)

        ! Build the Time100 coordinate values based on experiment schedule.
        ! Include t=time_init for non-exp1 cases so the initial state has a slot
        ! (matches wiki exp3 schedule "Time0, Time100, Time200, ...").
        select case(trim(ctl%exp))
            case("exp1")
                ! Single final snapshot
                n_time100 = 1
                allocate(t100(n_time100))
                t100(1) = ctl%time_end
            case default
                n_time100 = max(1, int((ctl%time_end - ctl%time_init) / 100.0_wp) + 1)
                allocate(t100(n_time100))
                do i = 1, n_time100
                    t100(i) = ctl%time_init + (i-1)*100.0_wp
                end do
        end select

        ! Warn if the timestep won't honour the wiki cadence. Yearly Time1
        ! writes need dtt=1; 100-yearly Time100 writes need dtt to divide 100.
        select case(trim(ctl%exp))
            case("exp2","exp3","exp4","exp5")
                if (ctl%dtt > 1.0_wp + 1e-6_wp) then
                    write(*,"(a,f6.2,a)") &
                        " calvingmip_init:: WARNING: ctl.dtt = ", ctl%dtt, &
                        " > 1; wiki Time1 cadence is annual but"
                    write(*,"(a)") &
                        "                   output will only land at multiples of dtt."
                end if
                if (abs(100.0_wp - nint(100.0_wp/ctl%dtt)*ctl%dtt) > 1e-6_wp) then
                    write(*,"(a,f6.2,a)") &
                        " calvingmip_init:: WARNING: ctl.dtt = ", ctl%dtt, &
                        " does not divide 100; some Time100 slots will be skipped."
                end if
        end select

        call nc_create(ctl%file_cmip, overwrite=.TRUE., institution="AWI", &
                       description="CalvingMIP output from YELMO ice-sheet model")

        call nc_write_dim(ctl%file_cmip,"X",     x=ylmo%grd%xc, units="m", &
                          long_name="X coordinate (model native grid)", axis="X")
        call nc_write_dim(ctl%file_cmip,"Y",     x=ylmo%grd%yc, units="m", &
                          long_name="Y coordinate (model native grid)", axis="Y")
        call nc_write_dim(ctl%file_cmip,"Time1", x=0.0_wp, dx=1.0_wp, nx=1, &
                          units="years", unlimited=.TRUE.)
        call nc_write_dim(ctl%file_cmip,"Time100", x=t100, units="years")

        ! Set up profile geometry and write each profile's s{name} dim.
        call profiles_setup(profiles, ctl%exp, ylmo%grd%dx)
        do i = 1, size(profiles)
            call nc_write_dim(ctl%file_cmip,"s"//trim(profiles(i)%name), &
                              x=profiles(i)%s, units="m", &
                              long_name="Distance along profile "//trim(profiles(i)%name))
        end do

        return

    end subroutine calvingmip_init

    subroutine calvingmip_write_step(ylmo,ctl,time)
        ! Hardcoded per-experiment schedule. See wiki:
        !   exp1: single final snapshot (scalars + 2D + profiles) on Time100
        !   exp2: yearly scalars on Time1, 100-yearly 2D on Time100
        !   exp3: 100-yearly scalars + 2D on Time100 (profiles yearly on Time1 — commit 2)
        !   exp4, exp5: yearly scalars on Time1, 100-yearly 2D on Time100

        implicit none

        type(yelmo_class),  intent(IN) :: ylmo
        type(control_type), intent(IN) :: ctl
        real(wp),           intent(IN) :: time

        integer :: hundredths

        hundredths = nint(time*100)

        select case(trim(ctl%exp))

            case("exp1")
                ! Single snapshot at final time (Time100, n=1) — scalars, 2D, profiles
                if (abs(time - ctl%time_end) < 0.5_wp*ctl%dtt) then
                    call calvingmip_write_scalars(ylmo,ctl%file_cmip,time,"Time100")
                    call calvingmip_write_2D(ylmo,ctl%file_cmip,time)
                    call calvingmip_write_profiles(ylmo,ctl%file_cmip,time,"Time100")
                end if

            case("exp3")
                ! 100-yearly scalars + 2D on Time100
                if (mod(hundredths, 10000) == 0) then
                    call calvingmip_write_scalars(ylmo,ctl%file_cmip,time,"Time100")
                    call calvingmip_write_2D(ylmo,ctl%file_cmip,time)
                end if
                ! Yearly profiles on Time1
                if (mod(hundredths, 100) == 0) then
                    call calvingmip_write_profiles(ylmo,ctl%file_cmip,time,"Time1")
                end if

            case("exp2","exp4","exp5")
                ! Yearly scalars + profiles on Time1
                if (mod(hundredths, 100) == 0) then
                    call calvingmip_write_scalars(ylmo,ctl%file_cmip,time,"Time1")
                    call calvingmip_write_profiles(ylmo,ctl%file_cmip,time,"Time1")
                end if
                ! 100-yearly 2D on Time100
                if (mod(hundredths, 10000) == 0) then
                    call calvingmip_write_2D(ylmo,ctl%file_cmip,time)
                end if

        end select

        return

    end subroutine calvingmip_write_step

    subroutine calvingmip_write_scalars(ylmo,filename,time,timedim)
        ! Write domain-aggregate scalars on the requested time dimension
        ! (either "Time1" or "Time100"). Variables follow CalvingMIP wiki names:
        ! iareafl, iareagr, lim, limnsw, tendlicalvf, tendligroundf,
        ! iareatotal{NW,NE,SW,SE}.

        implicit none

        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time
        character(len=*),  intent(IN) :: timedim    ! "Time1" or "Time100"

        ! Local
        type(yregions_class) :: reg
        integer  :: ncid, n, i, j
        real(wp) :: rho_ice
        real(wp) :: dx, dy
        real(wp) :: flux_grl, calv_flt
        real(wp) :: iareatotalNW, iareatotalNE, iareatotalSW, iareatotalSE
        logical, allocatable :: mask_grl(:,:), mask_frnt(:,:)
        logical, allocatable :: mask_NW(:,:), mask_NE(:,:), mask_SW(:,:), mask_SE(:,:)

        rho_ice = 917.0_wp
        dx = ylmo%grd%dx
        dy = ylmo%grd%dy

        allocate(mask_grl(ylmo%grd%nx,ylmo%grd%ny))
        allocate(mask_frnt(ylmo%grd%nx,ylmo%grd%ny))
        allocate(mask_NW(ylmo%grd%nx,ylmo%grd%ny))
        allocate(mask_NE(ylmo%grd%nx,ylmo%grd%ny))
        allocate(mask_SW(ylmo%grd%nx,ylmo%grd%ny))
        allocate(mask_SE(ylmo%grd%nx,ylmo%grd%ny))

        reg = ylmo%reg

        mask_grl  = (ylmo%tpo%now%H_ice .gt. 0.0 .and. ylmo%tpo%now%f_grnd .gt. 0.0 &
                                                    .and. ylmo%tpo%now%mask_grz .eq. 0.0)
        mask_frnt = (ylmo%tpo%now%H_ice .gt. 0.0 .and. ylmo%tpo%now%f_grnd .eq. 0.0 &
                                                    .and. ylmo%tpo%now%mask_frnt .eq. 1.0)
        mask_NW = .FALSE.; mask_NE = .FALSE.; mask_SW = .FALSE.; mask_SE = .FALSE.
        do j = 1, ylmo%grd%ny
        do i = 1, ylmo%grd%nx
            if (ylmo%grd%xc(i) <= 0.0_wp .and. ylmo%grd%yc(j) >= 0.0_wp) mask_NW(i,j) = .TRUE.
            if (ylmo%grd%xc(i) >= 0.0_wp .and. ylmo%grd%yc(j) >= 0.0_wp) mask_NE(i,j) = .TRUE.
            if (ylmo%grd%xc(i) <= 0.0_wp .and. ylmo%grd%yc(j) <= 0.0_wp) mask_SW(i,j) = .TRUE.
            if (ylmo%grd%xc(i) >= 0.0_wp .and. ylmo%grd%yc(j) <= 0.0_wp) mask_SE(i,j) = .TRUE.
        end do
        end do

        if (count(mask_grl) > 0) then
            flux_grl = sum(ylmo%dyn%now%uxy_bar*ylmo%tpo%now%H_ice*rho_ice, mask=mask_grl)*dx
        else
            flux_grl = 0.0_wp
        end if

        if (count(mask_frnt) > 0) then
            calv_flt = sum(ylmo%dyn%now%uxy_bar*ylmo%tpo%now%H_ice*rho_ice, mask=mask_frnt)*dx
        else
            calv_flt = 0.0_wp
        end if

        iareatotalNW = count(ylmo%tpo%now%H_ice .gt. 0.0 .and. mask_NW)*dx*dy
        iareatotalNE = count(ylmo%tpo%now%H_ice .gt. 0.0 .and. mask_NE)*dx*dy
        iareatotalSW = count(ylmo%tpo%now%H_ice .gt. 0.0 .and. mask_SW)*dx*dy
        iareatotalSE = count(ylmo%tpo%now%H_ice .gt. 0.0 .and. mask_SE)*dx*dy

        call nc_open(filename,ncid,writable=.TRUE.)
        n = nc_time_index(filename,timedim,time,ncid)
        call nc_write(filename,timedim,time,dim1=timedim,start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"iareafl",reg%A_ice_f*1e6_wp,units="m^2", &
                long_name="Floating ice area", &
                standard_name="floating_ice_shelf_area",dim1=timedim,start=[n],ncid=ncid)
        call nc_write(filename,"iareagr",reg%A_ice_g*1e6_wp,units="m^2", &
                long_name="Grounded ice area", &
                standard_name="grounded_ice_sheet_area",dim1=timedim,start=[n],ncid=ncid)
        call nc_write(filename,"lim",reg%V_ice*rho_ice*1e9_wp,units="kg", &
                long_name="Total ice mass", &
                standard_name="land_ice_mass",dim1=timedim,start=[n],ncid=ncid)
        call nc_write(filename,"limnsw",reg%V_sl*rho_ice*1e9_wp,units="kg", &
                long_name="Mass above flotation", &
                standard_name="land_ice_mass_not_displacing_sea_water",dim1=timedim,start=[n],ncid=ncid)
        call nc_write(filename,"tendlicalvf",calv_flt,units="kg a-1", &
                long_name="Total calving flux", &
                standard_name="tendency_of_land_ice_mass_due_to_calving",dim1=timedim,start=[n],ncid=ncid)
        call nc_write(filename,"tendligroundf",flux_grl,units="kg a-1", &
                long_name="Total grounding line flux", &
                standard_name="tendency_of_grounded_ice_mass",dim1=timedim,start=[n],ncid=ncid)
        call nc_write(filename,"iareatotalNW",iareatotalNW,units="m^2", &
                long_name="Total ice area NorthWest",dim1=timedim,start=[n],ncid=ncid)
        call nc_write(filename,"iareatotalNE",iareatotalNE,units="m^2", &
                long_name="Total ice area NorthEast",dim1=timedim,start=[n],ncid=ncid)
        call nc_write(filename,"iareatotalSW",iareatotalSW,units="m^2", &
                long_name="Total ice area SouthWest",dim1=timedim,start=[n],ncid=ncid)
        call nc_write(filename,"iareatotalSE",iareatotalSE,units="m^2", &
                long_name="Total ice area SouthEast",dim1=timedim,start=[n],ncid=ncid)

        call nc_close(ncid)

        return

    end subroutine calvingmip_write_scalars

    subroutine calvingmip_write_2D(ylmo,filename,time)
        ! Write 2D snapshot fields on the Time100 dimension: lithk, xvelmean,
        ! yvelmean, mask, topg. Velocities are staggered onto aa-nodes.
        ! Mask convention (CalvingMIP wiki): grounded=1, floating=2, ocean=3.

        implicit none

        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time

        integer :: ncid, n, i, j
        character(len=32) :: dims(3)
        integer,  allocatable :: mask_cmip(:,:)
        real(wp), allocatable :: ux_aa(:,:), uy_aa(:,:)

        allocate(mask_cmip(ylmo%grd%nx,ylmo%grd%ny))
        allocate(ux_aa(ylmo%grd%nx,ylmo%grd%ny))
        allocate(uy_aa(ylmo%grd%nx,ylmo%grd%ny))

        dims(1) = "X"
        dims(2) = "Y"
        dims(3) = "Time100"

        mask_cmip = 3   ! ocean by default
        where(ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp) mask_cmip = 2
        where(ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp) mask_cmip = 1

        ux_aa = 0.0_wp
        uy_aa = 0.0_wp
        do j = 2, ylmo%grd%ny-1
        do i = 2, ylmo%grd%nx-1
            ux_aa(i,j) = 0.5_wp*(ylmo%dyn%now%ux_bar(i,j) + ylmo%dyn%now%ux_bar(i-1,j))
            uy_aa(i,j) = 0.5_wp*(ylmo%dyn%now%uy_bar(i,j) + ylmo%dyn%now%uy_bar(i,j-1))
        end do
        end do

        call nc_open(filename,ncid,writable=.TRUE.)
        n = nc_time_index(filename,"Time100",time,ncid)
        call nc_write(filename,"Time100",time,dim1="Time100",start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"lithk",ylmo%tpo%now%H_ice,start=[1,1,n],units="m", &
                long_name="Ice thickness", &
                standard_name="land_ice_thickness",dims=dims,ncid=ncid)
        call nc_write(filename,"xvelmean",ux_aa,start=[1,1,n],units="m a-1", &
                long_name="Vertical-mean X velocity", &
                standard_name="land_ice_vertical_mean_x_velocity",dims=dims,ncid=ncid)
        call nc_write(filename,"yvelmean",uy_aa,start=[1,1,n],units="m a-1", &
                long_name="Vertical-mean Y velocity", &
                standard_name="land_ice_vertical_mean_y_velocity",dims=dims,ncid=ncid)
        call nc_write(filename,"mask",mask_cmip,start=[1,1,n],units="1", &
                long_name="Ice mask (1=grounded, 2=floating, 3=ocean)", &
                dims=dims,ncid=ncid)
        call nc_write(filename,"topg",ylmo%bnd%z_bed,start=[1,1,n],units="m", &
                long_name="Bedrock elevation", &
                standard_name="bedrock_altimetry",dims=dims,ncid=ncid)

        call nc_close(ncid)

        return

    end subroutine calvingmip_write_2D

    subroutine profiles_setup(profs, exp_name, ds)
        ! Build CalvingMIP profile geometry: 8 radial profiles A-H from (0,0)
        ! for the circular domain, or 8 named profiles (Caprona A-D + Halbrane A-D)
        ! for the Thule domain. Sample spacing ds [m] matches the native dx.
        ! Bearings (circular) per wiki: A=0°(N), B=45°(NE), C=90°(E), D=135°(SE),
        ! E=180°(S), F=225°(SW), G=270°(W), H=315°(NW). Compass convention:
        ! x = r*sin(bearing), y = r*cos(bearing).

        implicit none

        type(profile_t), allocatable, intent(out) :: profs(:)
        character(len=*), intent(IN) :: exp_name
        real(wp),         intent(IN) :: ds

        integer  :: ip
        real(wp) :: r_max
        real(wp) :: bearing(8)
        character(len=4), parameter :: circ_names(8) = &
            [character(len=4) :: "A","B","C","D","E","F","G","H"]
        character(len=4), parameter :: thule_names(8) = &
            [character(len=4) :: "CapA","CapB","CapC","CapD", &
                                 "HalA","HalB","HalC","HalD"]
        real(wp) :: thule_endpts(8,4)   ! [x0, y0, x1, y1] in meters
        real(wp) :: b

        allocate(profs(8))

        select case(trim(exp_name))

            case("exp1","exp2","advection")
                ! Circular domain — radial profiles from (0,0) out to r_max.
                r_max = 800.0e3_wp
                bearing = [0.0_wp, 45.0_wp, 90.0_wp, 135.0_wp, &
                           180.0_wp, 225.0_wp, 270.0_wp, 315.0_wp]
                do ip = 1, 8
                    profs(ip)%name = trim(circ_names(ip))
                    b = bearing(ip) * degrees_to_radians
                    call build_profile_line(profs(ip), 0.0_wp, 0.0_wp, &
                                            r_max*sin(b), r_max*cos(b), ds)
                end do

            case("exp3","exp4","exp5")
                ! Thule domain — Caprona and Halbrane profiles per wiki coords.
                thule_endpts(1,:) = [-390.0e3_wp,    0.0_wp, -590.0e3_wp,  450.0e3_wp]  ! CapA
                thule_endpts(2,:) = [ 390.0e3_wp,    0.0_wp,  590.0e3_wp,  450.0e3_wp]  ! CapB
                thule_endpts(3,:) = [-390.0e3_wp,    0.0_wp, -590.0e3_wp, -450.0e3_wp]  ! CapC
                thule_endpts(4,:) = [ 390.0e3_wp,    0.0_wp,  590.0e3_wp, -450.0e3_wp]  ! CapD
                thule_endpts(5,:) = [-150.0e3_wp,    0.0_wp, -150.0e3_wp,  740.0e3_wp]  ! HalA
                thule_endpts(6,:) = [ 150.0e3_wp,    0.0_wp,  150.0e3_wp,  740.0e3_wp]  ! HalB
                thule_endpts(7,:) = [-150.0e3_wp,    0.0_wp, -150.0e3_wp, -740.0e3_wp]  ! HalC
                thule_endpts(8,:) = [ 150.0e3_wp,    0.0_wp,  150.0e3_wp, -740.0e3_wp]  ! HalD
                do ip = 1, 8
                    profs(ip)%name = trim(thule_names(ip))
                    call build_profile_line(profs(ip), &
                            thule_endpts(ip,1), thule_endpts(ip,2), &
                            thule_endpts(ip,3), thule_endpts(ip,4), ds)
                end do

        end select

        return

    end subroutine profiles_setup

    subroutine build_profile_line(prof, x0, y0, x1, y1, ds)
        ! Fill a profile_t with sample points along the line from (x0,y0) to
        ! (x1,y1) at approximate spacing ds. The first sample is the start and
        ! the last sample is the end, so the actual spacing may be slightly
        ! less than ds (ceil division).

        implicit none

        type(profile_t), intent(INOUT) :: prof
        real(wp),        intent(IN)    :: x0, y0, x1, y1, ds

        integer  :: i
        real(wp) :: L, ux, uy

        L = sqrt((x1-x0)**2 + (y1-y0)**2)
        prof%n = max(2, ceiling(L/ds) + 1)

        if (allocated(prof%s)) deallocate(prof%s)
        if (allocated(prof%x)) deallocate(prof%x)
        if (allocated(prof%y)) deallocate(prof%y)
        allocate(prof%s(prof%n), prof%x(prof%n), prof%y(prof%n))

        if (L > 0.0_wp) then
            ux = (x1-x0)/L; uy = (y1-y0)/L
        else
            ux = 0.0_wp;    uy = 0.0_wp
        end if

        do i = 1, prof%n
            prof%s(i) = (i-1) * L / (prof%n - 1)
            prof%x(i) = x0 + prof%s(i)*ux
            prof%y(i) = y0 + prof%s(i)*uy
        end do

        return

    end subroutine build_profile_line

    function bilinear_sample(field, xc, yc, x, y) result(val)
        ! Bilinear interpolation of a real field at (x,y). Returns 0 if (x,y)
        ! falls outside the grid bounds.

        implicit none

        real(wp), intent(IN) :: field(:,:)
        real(wp), intent(IN) :: xc(:), yc(:)
        real(wp), intent(IN) :: x, y
        real(wp) :: val

        integer  :: i, j, nx, ny
        real(wp) :: tx, ty, dx, dy

        nx = size(xc); ny = size(yc)
        val = 0.0_wp

        if (x < xc(1) .or. x > xc(nx) .or. y < yc(1) .or. y > yc(ny)) return

        dx = xc(2) - xc(1)
        dy = yc(2) - yc(1)
        i = min(nx-1, max(1, int((x - xc(1))/dx) + 1))
        j = min(ny-1, max(1, int((y - yc(1))/dy) + 1))

        tx = (x - xc(i))/dx
        ty = (y - yc(j))/dy

        val = (1.0_wp-tx)*(1.0_wp-ty)*field(i,j)   + tx*(1.0_wp-ty)*field(i+1,j) + &
              (1.0_wp-tx)*ty*field(i,j+1)         + tx*ty*field(i+1,j+1)

        return

    end function bilinear_sample

    function nearest_sample_int(field, xc, yc, x, y, default_val) result(val)
        ! Nearest-neighbour sampling of an integer field at (x,y). Returns
        ! default_val if (x,y) falls outside the grid bounds.

        implicit none

        integer,  intent(IN) :: field(:,:)
        real(wp), intent(IN) :: xc(:), yc(:)
        real(wp), intent(IN) :: x, y
        integer,  intent(IN) :: default_val
        integer :: val

        integer  :: i, j, nx, ny
        real(wp) :: dx, dy

        nx = size(xc); ny = size(yc)
        val = default_val

        if (x < xc(1) .or. x > xc(nx) .or. y < yc(1) .or. y > yc(ny)) return

        dx = xc(2) - xc(1)
        dy = yc(2) - yc(1)
        i = min(nx, max(1, nint((x - xc(1))/dx) + 1))
        j = min(ny, max(1, nint((y - yc(1))/dy) + 1))

        val = field(i,j)

        return

    end function nearest_sample_int

    subroutine calvingmip_write_profiles(ylmo,filename,time,timedim)
        ! Sample fields along every profile and write the per-profile variables
        ! to the CalvingMIP file. Variable names follow the wiki:
        !   lithk{P}, xvelmean{P}, yvelmean{P}, mask{P}    on (s{P}, timedim)
        !   xcf{P}, ycf{P}, lithkcf{P}, xvelmeancf{P}, yvelmeancf{P}  on (timedim)
        ! Calving front is detected as the last ice sample (mask in {1,2}) before
        ! the first ocean sample (mask=3) walking outward from s=0.

        implicit none

        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time
        character(len=*),  intent(IN) :: timedim

        ! Local
        integer  :: ncid, nt, ip, i, j, i_lsf, i_inland
        integer,  allocatable :: mask_cmip(:,:)
        real(wp), allocatable :: ux_aa(:,:), uy_aa(:,:), H_eff(:,:)
        real(wp), allocatable :: lithk_p(:), xvel_p(:), yvel_p(:), lsf_p(:)
        integer,  allocatable :: mask_p(:)
        real(wp) :: xcf, ycf, lithkcf, xvelmeancf, yvelmeancf
        real(wp) :: lsf_a, lsf_b, frac
        character(len=32) :: sdim, dim2(2)
        character(len=8)  :: pname

        ! Build 2D mask and aa-staggered velocity arrays once (shared across profiles).
        allocate(mask_cmip(ylmo%grd%nx,ylmo%grd%ny))
        allocate(ux_aa(ylmo%grd%nx,ylmo%grd%ny))
        allocate(uy_aa(ylmo%grd%nx,ylmo%grd%ny))
        allocate(H_eff(ylmo%grd%nx,ylmo%grd%ny))

        mask_cmip = 3
        where(ylmo%tpo%now%f_ice .eq. 1.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp) mask_cmip = 2
        where(ylmo%tpo%now%f_ice .eq. 1.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp) mask_cmip = 1

        ux_aa = 0.0_wp; uy_aa = 0.0_wp; H_eff = 0.0_wp
        do j = 2, ylmo%grd%ny-1
        do i = 2, ylmo%grd%nx-1
            if (ylmo%tpo%now%f_ice(i,j) .gt. 0.0) then
                ux_aa(i,j) = 0.5_wp*(ylmo%dyn%now%ux_bar(i,j) + ylmo%dyn%now%ux_bar(i-1,j))
                uy_aa(i,j) = 0.5_wp*(ylmo%dyn%now%uy_bar(i,j) + ylmo%dyn%now%uy_bar(i,j-1))
                H_eff(i,j) = ylmo%tpo%now%H_ice(i,j)/ylmo%tpo%now%f_ice(i,j)
            end if
        end do
        end do

        call nc_open(filename,ncid,writable=.TRUE.)
        nt = nc_time_index(filename,timedim,time,ncid)
        call nc_write(filename,timedim,time,dim1=timedim,start=[nt],count=[1],ncid=ncid)

        do ip = 1, size(profiles)

            pname = profiles(ip)%name
            sdim  = "s"//trim(pname)
            dim2(1) = trim(sdim)
            dim2(2) = trim(timedim)

            allocate(lithk_p(profiles(ip)%n))
            allocate(xvel_p(profiles(ip)%n))
            allocate(yvel_p(profiles(ip)%n))
            allocate(lsf_p(profiles(ip)%n))
            allocate(mask_p(profiles(ip)%n))

            do i = 1, profiles(ip)%n
                lithk_p(i) = bilinear_sample(H_eff, ylmo%grd%xc, ylmo%grd%yc, &
                                             profiles(ip)%x(i), profiles(ip)%y(i))
                xvel_p(i)  = bilinear_sample(ux_aa, ylmo%grd%xc, ylmo%grd%yc, &
                                             profiles(ip)%x(i), profiles(ip)%y(i))
                yvel_p(i)  = bilinear_sample(uy_aa, ylmo%grd%xc, ylmo%grd%yc, &
                                             profiles(ip)%x(i), profiles(ip)%y(i))
                lsf_p(i)   = bilinear_sample(ylmo%tpo%now%lsf, ylmo%grd%xc, ylmo%grd%yc, &
                                             profiles(ip)%x(i), profiles(ip)%y(i))
                mask_p(i)  = nearest_sample_int(mask_cmip, ylmo%grd%xc, ylmo%grd%yc, &
                                                profiles(ip)%x(i), profiles(ip)%y(i), 3)
            end do

            ! Calving front via LSF zero crossing along the profile: locate
            ! the first sample with lsf > 0 (ocean side) so the preceding
            ! sample sits at lsf <= 0 (ice side). xcf/ycf are linearly
            ! interpolated to the zero crossing for subgrid front position.
            ! lithkcf and velocities are taken one sample further inland
            ! (i_lsf-2) to skip the transition cell, where calving has
            ! pulled H_ice down before calc_ice_fraction has re-flagged the
            ! cell as fractional, producing a misleading thickness dip.
            i_lsf = 0
            do i = 1, profiles(ip)%n
                if (lsf_p(i) > 0.0_wp) then
                    i_lsf = i
                    exit
                end if
            end do
            if (i_lsf >= 2) then
                lsf_a      = lsf_p(i_lsf-1)
                lsf_b      = lsf_p(i_lsf)
                frac       = -lsf_a / (lsf_b - lsf_a)
                xcf        = profiles(ip)%x(i_lsf-1) + frac*(profiles(ip)%x(i_lsf) - profiles(ip)%x(i_lsf-1))
                ycf        = profiles(ip)%y(i_lsf-1) + frac*(profiles(ip)%y(i_lsf) - profiles(ip)%y(i_lsf-1))
                i_inland   = max(1, i_lsf - 2)
                lithkcf    = lithk_p(i_inland)
                xvelmeancf = xvel_p(i_inland)
                yvelmeancf = yvel_p(i_inland)
            else if (i_lsf == 0) then
                ! Front lies beyond the profile end — fall back to the last
                ! profile sample for position and one inland for thickness.
                xcf        = profiles(ip)%x(profiles(ip)%n)
                ycf        = profiles(ip)%y(profiles(ip)%n)
                i_inland   = max(1, profiles(ip)%n - 1)
                lithkcf    = lithk_p(i_inland)
                xvelmeancf = xvel_p(i_inland)
                yvelmeancf = yvel_p(i_inland)
            else
                ! i_lsf == 1: profile starts on the ocean side, no ice.
                xcf        = 0.0_wp
                ycf        = 0.0_wp
                lithkcf    = 0.0_wp
                xvelmeancf = 0.0_wp
                yvelmeancf = 0.0_wp
            end if

            call nc_write(filename,"lithk"//trim(pname), lithk_p, &
                          dim1=sdim,dim2=timedim,start=[1,nt],units="m", &
                          long_name="Ice thickness along profile "//trim(pname),ncid=ncid)
            call nc_write(filename,"xvelmean"//trim(pname), xvel_p, &
                          dim1=sdim,dim2=timedim,start=[1,nt],units="m a-1", &
                          long_name="X velocity along profile "//trim(pname),ncid=ncid)
            call nc_write(filename,"yvelmean"//trim(pname), yvel_p, &
                          dim1=sdim,dim2=timedim,start=[1,nt],units="m a-1", &
                          long_name="Y velocity along profile "//trim(pname),ncid=ncid)
            call nc_write(filename,"mask"//trim(pname), mask_p, &
                          dim1=sdim,dim2=timedim,start=[1,nt],units="1", &
                          long_name="Ice mask along profile "//trim(pname),ncid=ncid)

            call nc_write(filename,"xcf"//trim(pname), xcf, &
                          dim1=timedim,start=[nt],units="m", &
                          long_name="Calving front x position, profile "//trim(pname),ncid=ncid)
            call nc_write(filename,"ycf"//trim(pname), ycf, &
                          dim1=timedim,start=[nt],units="m", &
                          long_name="Calving front y position, profile "//trim(pname),ncid=ncid)
            call nc_write(filename,"lithkcf"//trim(pname), lithkcf, &
                          dim1=timedim,start=[nt],units="m", &
                          long_name="Calving front ice thickness, profile "//trim(pname),ncid=ncid)
            call nc_write(filename,"xvelmeancf"//trim(pname), xvelmeancf, &
                          dim1=timedim,start=[nt],units="m a-1", &
                          long_name="Calving front x velocity, profile "//trim(pname),ncid=ncid)
            call nc_write(filename,"yvelmeancf"//trim(pname), yvelmeancf, &
                          dim1=timedim,start=[nt],units="m a-1", &
                          long_name="Calving front y velocity, profile "//trim(pname),ncid=ncid)

            deallocate(lithk_p, xvel_p, yvel_p, lsf_p, mask_p)

        end do

        call nc_close(ncid)

        return

    end subroutine calvingmip_write_profiles

end program yelmo_calving
