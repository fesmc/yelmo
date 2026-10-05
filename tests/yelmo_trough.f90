program yelmo_trough
    ! For mimicking Feldmann and Levermann (2017, TC) 
    ! and for running mismip+ etc. 

    use omp_lib

    use nml
    use ncio 
    use yelmo 
    use deformation 
    use lsf_module, only : LSFinit
    use timestepping
    use yelmo_tools, only : integrate_trapezoid1D_pt
    use, intrinsic :: iso_fortran_env, only : int64

    implicit none 

    type(tstep_class) :: ts
    type(yelmo_class) :: yelmo1

    character(len=56)  :: domain, grid_name  
    character(len=256) :: outfldr, file2D, file1D, file_restart, file_ts
    character(len=512) :: path_par 
    character(len=56)  :: experiment, res  
    real(wp) :: time_init, time_end, dt1D_out, dt2D_out
    logical  :: write_ts                            ! Write TROUGH-F17 point/transect time series?
    real(wp) :: dtts_out                            ! [yr] Frequency of time-series output
    real(wp) :: pts_x(3)                            ! [km] Points on the centreline (y=0)
    real(wp) :: sec_x                               ! [km] Position of the cross-section (y-transect)
    real(wp) :: H0_init                             ! [m] Initial ice thickness of the slab (TROUGH-F17)
    character(len=512) :: H0_file                   ! TROUGH-F17: file with the initial H_ice on the model grid ("None": slab)
    character(len=56) :: thrm_init                  ! Initial ice temperature: "robin-cold", "robin", "linear"
    logical  :: full_domain                         ! TROUGH-F17: full symmetric domain x = [-lx,lx] (divide inside)?
    real(wp) :: dtt
    integer  :: n

    ! Control parameters 
    real(wp) :: dx 
    real(wp) :: lx, ly, fc, dc, wc
    real(wp) :: x_cf 
    real(wp) :: Tsrf_const, smb_const, Qgeo_const  

    real(wp) :: s06_alpha, s06_H0, s06_W, s06_m
    real(wp) :: fs_H0, fs_dHdx, fs_zb               ! FRONT-SLAB: divide thickness, thickness gradient, bed elevation
    real(wp) :: col_H0, col_zb                      ! COLUMN-SLAB: initial ice thickness, bed elevation
    real(wp) :: B, L  
    real(wp), allocatable :: ux_ref(:,:) 
    real(wp), allocatable :: tau_c_ref(:,:)
    real(wp), allocatable :: H_init(:,:)            ! [m] Analytic initial ice thickness
    real(wp), allocatable :: z_srf_init(:,:)        ! [m] Analytic initial surface elevation

    real(wp) :: xmax, ymin, ymax, x0, y0 
    integer  :: i, j, nx, ny 
    
    real(8)  :: cpu_start_time, cpu_end_time, cpu_dtime  
    integer  :: perr 

    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)
    
    ! Assume program is running from the output folder
    outfldr = "./"

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)
    !path_par   = trim(outfldr)//"yelmo_TROUGH-F17.nml" 
    
    ! Define input and output locations 
    file2D     = trim(outfldr)//"yelmo.nc"
    file1D     = trim(outfldr)//"yelmo_ts.nc"
    file_restart = trim(outfldr)//"yelmo_restart.nc"
    file_ts    = trim(outfldr)//"yelmo_trough_ts.nc"
    
    ! Define the domain, grid and experiment from parameter file
    call nml_read(path_par,"ctrl","domain",       domain)        ! TROUGH-F17, MISMIP+

    ! Timing parameters 
    call nml_read(path_par,"ctrl","time_init",    time_init)     ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_end",     time_end)      ! [yr] Ending time
    call nml_read(path_par,"ctrl","dtt",          dtt)           ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","dt1D_out",     dt1D_out)      ! [yr] Frequency of 1D output 
    call nml_read(path_par,"ctrl","dt2D_out",     dt2D_out)      ! [yr] Frequency of 2D output 

    ! Domain parameters
    call nml_read(path_par,"ctrl","dx",           dx)            ! [km] Grid resolution ! must be multiple of xmax and ymax!!
    call nml_read(path_par,"ctrl","lx",           lx)            ! [km] Trough parameter
    call nml_read(path_par,"ctrl","ly",           ly)            ! [km] Trough parameter
    call nml_read(path_par,"ctrl","fc",           fc)            ! [km] Trough parameter
    call nml_read(path_par,"ctrl","dc",           dc)            ! [km] Trough parameter
    call nml_read(path_par,"ctrl","wc",           wc)            ! [km] Trough parameter
    call nml_read(path_par,"ctrl","x_cf",         x_cf)          ! [km] Trough parameter

    ! TROUGH-F17 only: initial state and point/transect time series
    write_ts  = (trim(domain) .eq. "TROUGH-F17")
    H0_init   = 50.0_wp
    H0_file   = "None"
    thrm_init = "robin-cold"
    full_domain = .FALSE.
    if (trim(domain) .eq. "TROUGH-F17") then
        call nml_read(path_par,"ctrl","H0_init",  H0_init)           ! [m] Initial ice thickness of the slab
        call nml_read(path_par,"ctrl","H0_file",  H0_file)           ! Initial H_ice from file ("None": slab of H0_init)
        call nml_read(path_par,"ctrl","thrm_init",thrm_init)         ! Initial ice temperature
        call nml_read(path_par,"ctrl","full_domain",full_domain)     ! Full symmetric domain x = [-lx,lx]?
        call nml_read(path_par,"ctrl","dtts_out", dtts_out)          ! [yr] Frequency of time-series output
        call nml_read(path_par,"ctrl","pts_x",    pts_x)             ! [km] Points on the centreline (y=0)
        call nml_read(path_par,"ctrl","sec_x",    sec_x)             ! [km] Position of the cross-section
    end if

    ! Schoof domain parameters (only needed by the slab domains)
    select case(trim(domain))
        case("SLAB-S06","RAYMOND")
            call nml_read(path_par,"ctrl_schoof","alpha",s06_alpha)  ! [m/m] Constant slope
            call nml_read(path_par,"ctrl_schoof","H0",   s06_H0)     ! [m]   Constant ice thickness
            call nml_read(path_par,"ctrl_schoof","W",    s06_W)      ! [m]   Half-width weak till
            call nml_read(path_par,"ctrl_schoof","m",    s06_m)      ! []    Exponent
        case("FRONT-SLAB")
            call nml_read(path_par,"ctrl_front","H0",   fs_H0)       ! [m]    Ice thickness at the divide (x=0)
            call nml_read(path_par,"ctrl_front","dHdx", fs_dHdx)     ! [m/km] Thickness decrease with |x|
            call nml_read(path_par,"ctrl_front","zb",   fs_zb)       ! [m]    Flat bed elevation
            full_domain = .TRUE.
        case("COLUMN-SLAB")
            call nml_read(path_par,"ctrl_column","H0",  col_H0)      ! [m]    Initial ice thickness
            call nml_read(path_par,"ctrl_column","zb",  col_zb)      ! [m]    Flat bed elevation
    end select

    ! Simulation parameters 
    call nml_read(path_par,"ctrl","Tsrf_const",   Tsrf_const)    ! [degC]  Surface temperature
    call nml_read(path_par,"ctrl","smb_const",    smb_const)     ! [m/yr]  Surface mass balance
    call nml_read(path_par,"ctrl","Qgeo_const",   Qgeo_const)    ! [mW/m2] Geothermal heat flux
    
    ! === Initialize timestepping ===
    
    call tstep_init(ts,time_init,time_end,method="const",units="year", &
                                            time_ref=1950.0_wp,const_rel=0.0_wp,const_cal=1950.0_wp)

    
    ! Define default grid name for completeness 
    grid_name = trim(domain)
    
    ! Define the domain and grid
    xmax =  lx 
    ymax =  ly/2.0_wp
    ymin = -ly/2.0_wp

    select case(trim(domain))

        case("TROUGH-F17","MISMIP+","FRONT-SLAB")
            ! Channel periodic in y (true wrap, period ny*dx, no halo): centred
            ! grid y_j = (j-jc)*dx with jc = ny/2+1, so that y=0 is a row and the
            ! period is exactly ly. For even ny the wall y=-ly/2 is a row, for
            ! odd ny the wall lies midway between the first and last rows.

            ! With full_domain (TROUGH-F17), x spans [-lx,lx] and the ice
            ! divide at x=0 is an interior point, as in Feldmann and
            ! Levermann (2017). Otherwise x spans [0,lx] with a symmetry
            ! boundary at x=0.

            if (full_domain) then
                ! Symmetric about x=0 (a grid point) for any dx
                nx = 2*int(xmax/dx)+1
                x0 = -real(int(xmax/dx),wp)*dx
            else
                nx = int(xmax/dx)+1
                x0 = 0.0_wp
            end if
            ny = periodic_npts(ly,dx,"ly")
            y0 = -real(ny/2,wp)*dx

        case("SLAB-S06","RAYMOND","COLUMN-SLAB")
            ! Slab periodic in x and y (experiment "SLAB"): x_i = (i-1)*dx with
            ! period exactly lx (no duplicated end column), and the same centred
            ! y-grid as the channel above, with period exactly ly.

            nx = periodic_npts(lx,dx,"lx")
            ny = periodic_npts(ly,dx,"ly")
            x0 = 0.0_wp
            y0 = -real(ny/2,wp)*dx

        case DEFAULT

            nx = int(xmax/dx)+1
            ny = int((ymax-ymin)/dx)+1
            x0 = 0.0_wp
            y0 = ymin

    end select

    call yelmo_init_grid(yelmo1%grd,grid_name,units="km", &
                            x0=x0,dx=dx,nx=nx, &
                            y0=y0,dy=dx,ny=ny)

    ! === Initialize ice sheet model =====

    ! Initialize data objects
    call yelmo_init(yelmo1,filename=path_par,grid_def="none",time=ts%time,load_topo=.FALSE., &
                        domain=domain,grid_name=grid_name)

    ! Load boundary values

    yelmo1%bnd%z_sl     = 0.0
    yelmo1%bnd%bmb_shlf = 0.0 
    yelmo1%bnd%T_shlf   = yelmo1%bnd%c%T0  
    yelmo1%bnd%H_sed    = 0.0 

    yelmo1%bnd%T_srf    = yelmo1%bnd%c%T0 + Tsrf_const   ! [K] 
    yelmo1%bnd%smb      = smb_const         ! [m/yr]
    yelmo1%bnd%Q_geo    = Qgeo_const        ! [mW/m2] 

    ! Check boundary values 
    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize output file 
    call yelmo_write_init(yelmo1,file2D,time_init=ts%time,units="years")
    
    ! Intialize topography (bed always; the initial ice thickness is
    ! applied below unless it comes from the restart file)
    allocate(H_init(yelmo1%grd%G%nx,yelmo1%grd%G%ny))
    allocate(z_srf_init(yelmo1%grd%G%nx,yelmo1%grd%G%ny))

    select case(trim(domain)) 

        case("RAYMOND")
            ! Raymond (2000) domain - constant slope slab

            ! ===== Intialize topography and set parameters =========
        
            ! The tilted bed is not periodic in x: carry the slope as a uniform
            ! background slope (ytopo slope_bg_x), added to the surface and bed
            ! gradients; z_bed and z_srf only contain the periodic (flat) part.
            yelmo1%tpo%par%slope_bg_x = -s06_alpha
            yelmo1%bnd%z_bed = 10000.0_wp

            H_init = s06_H0

            ! Define surface elevation 
            z_srf_init = yelmo1%bnd%z_bed + H_init

            ! Define reference ice thickness (for prescribing boundary values, potentially)
            yelmo1%bnd%H_ice_ref = s06_H0 

            ! Define basal friction 
            yelmo1%dyn%now%cb_ref = 5.2e3 
            where(abs(yelmo1%grd%y) .gt. s06_W) yelmo1%dyn%now%cb_ref = 70e3 

            write(*,*) "RAYMOND: W      = ", s06_W 
            write(*,*) "RAYMOND: ATT    = ", yelmo1%mat%now%ATT(1,1,1)
            write(*,*) "RAYMOND: cb_ref = ", yelmo1%dyn%now%cb_ref(1,1)

            ! Initialiaze ux values to be safe too 
            yelmo1%dyn%now%ux_b   = 500.0 
            yelmo1%dyn%now%ux_bar = 500.0 
            yelmo1%dyn%now%ux_s   = 500.0 
            
        case("SLAB-S06")
            ! Schoof (2006) domain - constant slope slab

            ! ===== Intialize topography and set parameters =========
        
            ! Tilted bed as a uniform background slope (see RAYMOND above)
            yelmo1%tpo%par%slope_bg_x = -s06_alpha
            yelmo1%bnd%z_bed = 10000.0_wp

            H_init = s06_H0
            yelmo1%bnd%H_ice_ref = s06_H0 

            ! Define surface elevation 
            z_srf_init = yelmo1%bnd%z_bed + H_init

            ! Calculate analytical stream function to get tau_c and ux

            allocate(ux_ref(yelmo1%grd%G%nx,yelmo1%grd%G%ny))
            allocate(tau_c_ref(yelmo1%grd%G%nx,yelmo1%grd%G%ny))
            
            call SSA_Schoof2006_analytical_solution_yelmo(ux_ref, tau_c_ref, yelmo1%grd%y, &
                                    s06_alpha,s06_H0,yelmo1%mat%par%rf_const,s06_W,s06_m, &
                                    yelmo1%mat%par%n_glen,yelmo1%bnd%c%rho_ice, yelmo1%bnd%c%g)
            
            ! Assign analytical values (tau_c as a boundary condition, ux as initial condition)
            yelmo1%dyn%now%cb_ref = tau_c_ref
            
            ! Assign initial velocity values to help achieve quicker convergence...
            yelmo1%dyn%now%ux_b   = ux_ref 
            yelmo1%dyn%now%ux_bar = ux_ref 
            yelmo1%dyn%now%ux_s   = ux_ref 
            yelmo1%dyn%now%uy_b   = 0.0_wp 
            yelmo1%dyn%now%uy_bar = 0.0_wp 
            yelmo1%dyn%now%uy_s   = 0.0_wp 

            ! Determine constant L too, for diagnostic output
            L = s06_W / ((1.0_wp+s06_m)**(1.0_wp/s06_m))

            write(*,*) "SLAB-S06: H0          = ", s06_H0 
            write(*,*) "SLAB-S06: alpha       = ", s06_alpha 
            write(*,*) "SLAB-S06: W           = ", s06_W 
            write(*,*) "SLAB-S06: L           = ", L 
            write(*,*) "SLAB-S06: m           = ", s06_m 
            write(*,*) "SLAB-S06: rho g       = ", yelmo1%bnd%c%rho_ice, yelmo1%bnd%c%g
            write(*,*) "SLAB-S06: f           = ", (yelmo1%bnd%c%rho_ice*yelmo1%bnd%c%g*s06_H0)*s06_alpha
            write(*,*) "SLAB-S06: ATT         = ", yelmo1%mat%par%rf_const
            write(*,*) "SLAB-S06: cb_ref      = ", yelmo1%dyn%now%cb_ref(1,1)
            write(*,*) "SLAB-S06: tau_c_ref   = ", tau_c_ref(1,1)
            write(*,*) "SLAB-S06: max(ux_ref) = ", maxval(ux_ref)

        case("TROUGH-F17")
            ! Feldmann and Levermann (2017) domain 

            call trough_f17_topo_init(yelmo1%bnd%z_bed,H_init,z_srf_init, &
                                    yelmo1%grd%G%x*1e-3,yelmo1%grd%G%y*1e-3,fc,dc,wc,x_cf,H0_init)

            ! Optionally replace the slab by an initial thickness from file (same grid)
            if (trim(H0_file) .ne. "None") then
                call nc_read(H0_file,"H_ice",H_init)
                where (abs(spread(yelmo1%grd%G%x*1e-3,2,size(H_init,2))) .gt. x_cf) H_init = 0.0_wp
                z_srf_init = yelmo1%bnd%z_bed + H_init
                write(*,*) "TROUGH-F17: initial H_ice read from ", trim(H0_file)
            end if
        
        case("FRONT-SLAB")
            ! Grounded marine slab with ice fronts at |x| = x_cf (front-stress benchmark)

            call front_slab_topo_init(yelmo1%bnd%z_bed,H_init,z_srf_init, &
                                    yelmo1%grd%G%x*1e-3,fs_H0,fs_dHdx,fs_zb,x_cf)

        case("COLUMN-SLAB")
            ! Flat, periodic slab of uniform thickness: no flow, so each column
            ! only thickens by smb (thermodynamics benchmark)

            yelmo1%bnd%z_bed = col_zb
            H_init           = col_H0
            z_srf_init       = yelmo1%bnd%z_bed + H_init

        case("MISMIP+") 
            ! MISMIP+ domain 

            call trough_mismipp_topo_init(yelmo1%bnd%z_bed,H_init,z_srf_init, &
                                    yelmo1%grd%G%x*1e-3,yelmo1%grd%G%y*1e-3,fc,dc,wc,x_cf)
        
        case("SLAB-SHELF")
            ! Constant slab slope with an ice shelf

            ! call trough_f17_topo_init(yelmo1%bnd%z_bed,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%z_srf, &
            !                         yelmo1%grd%G%x*1e-3,yelmo1%grd%G%y*1e-3,fc,dc,wc,x_cf)
            
            call slab_topo_init(yelmo1%bnd%z_bed,H_init,z_srf_init, &
                                    yelmo1%grd%G%x*1e-3,yelmo1%grd%G%y*1e-3)


        case DEFAULT 

            write(*,*) "yelmo_trough:: Error: domain not recognized: "//trim(domain)
            error stop 1

    end select 


    ! Initial ice thickness, unless yelmo_init loaded it from the restart file
    if (.not. (yelmo1%par%use_restart .and. yelmo1%par%restart_H_ice)) then
        yelmo1%tpo%now%H_ice = H_init
        yelmo1%tpo%now%z_srf = z_srf_init
    end if

    ! Define calving front 
    call define_calving_front(yelmo1%bnd%calv_mask,yelmo1%grd%x*1e-3,x_cf)

    ! Initialize the LSF mask from the topography, if not restarting
    if (.not. yelmo1%par%use_restart) then
        call LSFinit(yelmo1%tpo%now%lsf,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%tpo%par%dx)
    end if

    ! Initialize the yelmo state (dyn,therm,mat)
    call yelmo_init_state(yelmo1,time=ts%time,thrm_method=thrm_init)

    ! Write initial state 
    call write_step_2D(yelmo1,file2D,time=ts%time) 

    ! 1D file 
    call yelmo_write_reg_init(yelmo1,file1D,time_init=ts%time,units="years",mask=(yelmo1%bnd%mask_ice /= MASK_ICE_NONE))
    call yelmo_write_reg_step(yelmo1,file1D,time=ts%time)  

    ! Point/transect time-series file
    if (write_ts) then
        call trough_ts_init(file_ts,yelmo1,pts_x,sec_x,time_init=ts%time)
        call trough_ts_write(file_ts,yelmo1,pts_x,sec_x,time=ts%time)
    end if

    ! Advance timesteps
    call tstep_print_header(ts)

    do while (.not. ts%is_finished)

        ! == Update timestep ===

        call tstep_update(ts,dtt)
        call tstep_print(ts)
        
if (.FALSE.) then
        !if (trim(domain) .eq. "SLAB-SHELF" .and. ts%time_elapsed .ge. 3e3) then 
        if (trim(domain) .eq. "TROUGH-F17" .and. ts%time_elapsed .ge. 3e3) then 

            ! ! Define calving front 
            ! x_cf = 540.0_wp 
            ! call define_calving_front(yelmo1%bnd%calv_mask,yelmo1%grd%x*1e-3,x_cf)

            ! Kill all floating ice now
            yelmo1%tpo%par%calv_flt_method = "kill"

        end if 
end if 

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,ts%time)

        ! == MODEL OUTPUT =======================================================
        ! int64: a default integer overflows for |time| > ~2.1e7 yr
        if (mod(nint(ts%time_elapsed*100,int64),nint(dt2D_out*100,int64))==0) then  
            call write_step_2D(yelmo1,file2D,time=ts%time)    
        end if 

        if (mod(nint(ts%time_elapsed*100,int64),nint(dt1D_out*100,int64))==0) then 
            call yelmo_write_reg_step(yelmo1,file1D,time=ts%time) 
        end if

        if (write_ts) then
            if (mod(nint(ts%time_elapsed*100,int64),nint(dtts_out*100,int64))==0) then
                call trough_ts_write(file_ts,yelmo1,pts_x,sec_x,time=ts%time)
            end if
        end if

        if (mod(ts%time_elapsed,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", ts%time
        end if  

    end do 

    ! Write summary 
    write(*,*) "====== "//trim(domain)//" ======="
    write(*,*) "nz, H0 = ", yelmo1%par%nz_aa, maxval(yelmo1%tpo%now%H_ice)

    ! Write a restart file
    call yelmo_restart_write(yelmo1,file_restart,ts%time)

    ! Finalize program
    call yelmo_end(yelmo1,time=ts%time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*ts%time_elapsed)/(cpu_dtime/3600.0), " kiloyears / hr"
    
contains
    
    function periodic_npts(l,dx,name) result(n)
        ! Number of points in a periodic direction of length l (true wrap,
        ! period n*dx, no halo): l must be a multiple of dx.

        implicit none

        real(wp),         intent(IN) :: l
        real(wp),         intent(IN) :: dx
        character(len=*), intent(IN) :: name
        integer :: n

        n = nint(l/dx)
        if (abs(n*dx-l) .gt. 1e-6_wp*l) then
            write(*,*) "yelmo_trough:: Error: "//name//" must be a multiple of dx in a periodic direction."
            write(*,*) name//", dx = ", l, dx
            error stop 1
        end if

        return

    end function periodic_npts

    subroutine define_calving_front(calv_mask,xx,x_cf)
        ! Define a calving mask in the x direction where 
        ! beyond the position |x| = x_cf ice will be calved. 

        implicit none 

        logical, intent(OUT) :: calv_mask(:,:) 
        real(dp), intent(IN) :: xx(:,:) 
        real(wp), intent(IN) :: x_cf 

        calv_mask = .FALSE. 
        where (abs(xx) .ge. x_cf) calv_mask = .TRUE. 
    
        return 

    end subroutine define_calving_front

    subroutine slab_topo_init(z_bed,H_ice,z_srf,xc,yc)

        implicit none 

        real(prec), intent(OUT) :: z_bed(:,:) 
        real(prec), intent(OUT) :: H_ice(:,:) 
        real(prec), intent(OUT) :: z_srf(:,:) 
        real(dp), intent(IN)  :: xc(:)
        real(dp), intent(IN)  :: yc(:)

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(z_bed,1)
        ny = size(z_bed,2)
        
        ! Define bedrock as a slope
        do j = 1, ny
            z_bed(:,j) = -100.0 - xc
        end do 
        
        ! Set ice thickness to 500 m everywhere initially 
        H_ice = 500.0

        ! Remove ice from deep bed
        where(z_bed .lt. -500.0) H_ice = 0.0 

        ! Adjust for floating ice later, for now assume fully grounded
        z_srf = z_bed + H_ice 

        return 

    end subroutine slab_topo_init

    subroutine front_slab_topo_init(z_bed,H_ice,z_srf,xc,H0,dHdx,zb,x_cf)
        ! Flat bed and an ice thickness decreasing linearly away from the
        ! divide at x=0, uniform in y, ending at ice fronts |x| < x_cf.

        implicit none

        real(wp), intent(OUT) :: z_bed(:,:)
        real(wp), intent(OUT) :: H_ice(:,:)
        real(wp), intent(OUT) :: z_srf(:,:)
        real(dp), intent(IN)  :: xc(:)          ! [km]
        real(wp), intent(IN)  :: H0             ! [m]    Thickness at x=0
        real(wp), intent(IN)  :: dHdx           ! [m/km] Thickness decrease with |x|
        real(wp), intent(IN)  :: zb             ! [m]    Bed elevation
        real(wp), intent(IN)  :: x_cf           ! [km]   Front position

        ! Local variables
        integer :: i

        z_bed = zb
        do i = 1, size(H_ice,1)
            if (abs(xc(i)) .lt. x_cf) then
                H_ice(i,:) = max(H0 - dHdx*abs(xc(i)), 0.0_wp)
            else
                H_ice(i,:) = 0.0_wp
            end if
        end do

        z_srf = z_bed + H_ice

        return

    end subroutine front_slab_topo_init

    subroutine trough_f17_topo_init(z_bed,H_ice,z_srf,xc,yc,fc,dc,wc,x_cf,H0)

        implicit none 

        real(wp), intent(OUT) :: z_bed(:,:) 
        real(wp), intent(OUT) :: H_ice(:,:) 
        real(wp), intent(OUT) :: z_srf(:,:) 
        real(dp), intent(IN)  :: xc(:) 
        real(dp), intent(IN)  :: yc(:)  
        real(wp), intent(IN)  :: fc 
        real(wp), intent(IN)  :: dc 
        real(wp), intent(IN)  :: wc 
        real(wp), intent(IN)  :: x_cf 
        real(wp), intent(IN)  :: H0                     ! [m] Initial ice thickness of the slab

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp) :: zb_x, zb_y 
        real(wp) :: e1, e2 

        real(wp), parameter :: zb_deep = -720.0_wp 

        nx = size(z_bed,1) 
        ny = size(z_bed,2) 

        write(*,*) "params: ", ly,fc,dc,wc,x_cf

        ! == Bedrock elevation == 
        do j = 1, ny
        do i = 1, nx 
            
            ! x-direction 
            zb_x = -150.0_wp - 0.84*abs(xc(i))

            ! y-direction 
            e1 = -2.0*(yc(j)-wc)/fc 
            e2 =  2.0*(yc(j)+wc)/fc 
            zb_y = ( dc / (1.0+exp(e1)) ) + ( dc / (1.0+exp(e2)) ) 

            ! Convolution 
            z_bed(i,j) = max(zb_x + zb_y, zb_deep)

        end do
        end do  

        ! == Ice thickness == 
        H_ice = H0 
        do j = 1, ny 
            where(abs(xc) .gt. x_cf) H_ice(:,j) = 0.0 
        end do 

        ! == Surface elevation == 
        z_srf = z_bed + H_ice

        where(z_srf .lt. 0.0) z_srf = 0.0 

        return 

    end subroutine trough_f17_topo_init

    subroutine trough_mismipp_topo_init(z_bed,H_ice,z_srf,xc,yc,fc,dc,wc,x_cf)

        implicit none 

        real(wp), intent(OUT) :: z_bed(:,:) 
        real(wp), intent(OUT) :: H_ice(:,:) 
        real(wp), intent(OUT) :: z_srf(:,:) 
        real(dp), intent(IN)  :: xc(:) 
        real(dp), intent(IN)  :: yc(:)  
        real(wp), intent(IN)  :: fc 
        real(wp), intent(IN)  :: dc 
        real(wp), intent(IN)  :: wc 
        real(wp), intent(IN)  :: x_cf 

        ! Local variables 
        integer  :: i, j, nx, ny 
        real(wp) :: zb_x, zb_y 
        real(wp) :: e1, e2 

        real(wp) :: x1 

        real(wp), parameter :: zb_deep = -720.0_wp 
        real(wp), parameter :: xbar    =  300.0_wp          ! [km] Characteristic along-flow length scale of the bedrock
        real(wp), parameter :: b0      = -150.00_wp 
        real(wp), parameter :: b2      = -728.80_wp 
        real(wp), parameter :: b4      =  343.91_wp 
        real(wp), parameter :: b6      =  -50.57_wp 


        nx = size(z_bed,1) 
        ny = size(z_bed,2) 

        write(*,*) "params: ", ly,fc,dc,wc,x_cf

        ! == Bedrock elevation == 
        do j = 1, ny
        do i = 1, nx 
            
            ! x-direction 
            x1 = xc(i) / xbar 
            zb_x = b0 + b2*x1**2 + b4*x1**4 + b6*x1**6 

            ! y-direction 
            e1 = -2.0*(yc(j)-wc)/fc 
            e2 =  2.0*(yc(j)+wc)/fc 
            zb_y = ( dc / (1.0+exp(e1)) ) + ( dc / (1.0+exp(e2)) ) 

            ! Convolution 
            z_bed(i,j) = max(zb_x + zb_y, zb_deep)

        end do
        end do  

        ! == Ice thickness == 
        H_ice = 50.0_wp 
        do j = 1, ny 
            where(abs(xc) .gt. x_cf) H_ice(:,j) = 0.0 
        end do 

        ! == Surface elevation == 
        z_srf = z_bed + H_ice

        where(z_srf .lt. 0.0) z_srf = 0.0 

        return 

    end subroutine trough_mismipp_topo_init

    subroutine write_step_2D(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n, i, j, nx, ny  

        nx = ylmo%tpo%par%nx 
        ny = ylmo%tpo%par%ny 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"mask_frnt",ylmo%tpo%now%mask_frnt,units="",long_name="Ice-front mask", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"taul_int_acx",ylmo%dyn%now%taul_int_acx,units="Pa m",long_name="Vertically integrated lateral stress (x)", &
        !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"taul_int_acy",ylmo%dyn%now%taul_int_acy,units="Pa m",long_name="Vertically integrated lateral stress (y)", &
        !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"H_ice_pred",ylmo%tpo%now%pred%H_ice,units="m",long_name="Ice thickness (predicted)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"H_ice_corr",ylmo%tpo%now%corr%H_ice,units="m",long_name="Ice thickness (corrected)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! call nc_write(filename,"dzsdt",ylmo%tpo%now%dzsdt,units="m/a",long_name="Surface elevation change", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice-covered fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"cmb",ylmo%tpo%now%cmb,units="m/a",long_name="Calving mass balance rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"cmb_flt",ylmo%tpo%now%cmb_flt,units="m/a",long_name="Calving mass balance rate (floating)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_thermodynamics ==
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_prime,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="deg C",long_name="Homologous basal ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Basal water moved to hyd (fasthydrology); output as "hyd_W_til".
        call yelmo_write_var(filename,"hyd_W_til",ylmo,n,ncid)
        
        ! == yelmo_material ==
!         call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"visc",ylmo%mat%now%visc,units="Pa a",long_name="Viscosity", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"ATT_bar",ylmo%mat%now%ATT_bar,units="a^-1 Pa^-3",long_name="Vertically averaged rate factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="mW m-2",long_name="Basal ice heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(ylmo%bnd%c%rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"dQsdt",ylmo%thrm%now%dQsdt/(ylmo%bnd%c%rho_ice*ylmo%thrm%now%cp),units="K a-2",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="mW m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_dynamics ==

if (.FALSE.) then
        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
end if

        call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"beta_acx",ylmo%dyn%now%beta_acx,units="Pa a m-1",long_name="Basal friction coefficient (acx)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"beta_acy",ylmo%dyn%now%beta_acy,units="Pa a m-1",long_name="Basal friction coefficient (acy)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a",long_name="Effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface gradient, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface gradient, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"duxdz",ylmo%dyn%now%duxdz,units="1/a",long_name="Vertical shear (x)", &
!                        dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"duydz",ylmo%dyn%now%duydz,units="1/a",long_name="Vertical shear (y)", &
!                        dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"qq_gl_acx",ylmo%dyn%now%qq_gl_acx,units="m2/a",long_name="Grounding line flux (x)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"qq_gl_acy",ylmo%dyn%now%qq_gl_acy,units="m2/a",long_name="Grounding line flux (y)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

if (.FALSE.) then
        call nc_write(filename,"ux",ylmo%dyn%now%ux,units="m/a",long_name="Horizontal velocity (x)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uy",ylmo%dyn%now%uy,units="m/a",long_name="Horizontal velocity (y)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uxy",ylmo%dyn%now%uxy,units="m/a",long_name="Horizontal velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity", &
                      dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
        ! call nc_write(filename,"uz_star",ylmo%dyn%now%uz_star,units="m/a",long_name="Advective vertical velocity", &
        !               dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
        ! call nc_write(filename,"advecxy",ylmo%thrm%now%advecxy,units="m/a",long_name="Horizontal advection", &
        !               dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
end if

!         call nc_write(filename,"f_vbvs",ylmo%dyn%now%f_vbvs,units="1",long_name="Basal to surface velocity fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_shear_bar",ylmo%mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"de",ylmo%mat%now%strn%de,units="yr^-1",long_name="Effective strain rate", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! call nc_write(filename,"de_jac",ylmo%dyn%now%strn%de,units="yr^-1",long_name="Effective strain rate", &
        !               dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! == yelmo_bound ==

!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/yr",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/yr",long_name="Basal mass balance (shelf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! No time dimension::

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

    subroutine trough_ts_init(filename,ylmo,pts_x,sec_x,time_init)
        ! Initialize the TROUGH-F17 time-series file: points on the centreline
        ! (Feldmann and Levermann, 2017, Figs. A1, A3), a cross-section at sec_x
        ! (their Fig. 1c), and centreline and grounded-area means (their Figs. 3, A2)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(yelmo_class), intent(IN) :: ylmo
        real(wp),          intent(IN) :: pts_x(:)       ! [km] Points on the centreline (y=0)
        real(wp),          intent(IN) :: sec_x          ! [km] Position of the cross-section
        real(wp),          intent(IN) :: time_init

        call nc_create(filename)
        call nc_write_dim(filename,"pt",  x=pts_x,units="kilometers",long_name="Point position on the centreline (y=0)")
        call nc_write_dim(filename,"xsec",x=sec_x,dx=1.0_wp,nx=1,units="kilometers",long_name="Cross-section position")
        call nc_write_dim(filename,"yc",  x=ylmo%grd%G%y*1d-3,units="kilometers")
        call nc_write_dim(filename,"time",x=time_init,dx=1.0_wp,nx=1,units="years",unlimited=.TRUE.)

        return

    end subroutine trough_ts_init

    subroutine trough_ts_write(filename,ylmo,pts_x,sec_x,time)
        ! Write one step of the TROUGH-F17 time series. Each aa-node field is
        ! written at the points (name), along the cross-section (name_sec), as
        ! the grounded-ice mean along the centreline (name_cl) and as the
        ! grounded-area mean (name_gr). The x-momentum balance at acx-faces
        ! (calc_xmom_terms_row) is written at the points and as centreline mean.

        implicit none

        character(len=*),  intent(IN) :: filename
        type(yelmo_class), intent(IN) :: ylmo
        real(wp),          intent(IN) :: pts_x(:)       ! [km] Points on the centreline (y=0)
        real(wp),          intent(IN) :: sec_x          ! [km] Position of the cross-section
        real(wp),          intent(IN) :: time

        ! Local variables
        integer  :: ncid, n, i, j, j0, nx, ny
        real(wp) :: dx, dy, x0, sec_year, rho_ice
        real(wp) :: A_grnd, V_grnd, x_gl, calv, uxy_max_grnd, uxy_max_flt
        real(wp), allocatable :: wt(:,:)                ! [--] Grounded-ice weight, f_grnd*f_ice
        real(wp), allocatable :: flux(:,:)              ! [m2 a-1] Ice flux magnitude
        real(wp), allocatable :: Q_strn_int(:,:)        ! [mW m-2] Column-integrated strain heating
        real(wp), allocatable :: T_ice_bar(:,:)         ! [degC] Vertically averaged ice temperature
        real(wp), allocatable :: taud_x(:)              ! [Pa] Driving stress (flow direction)
        real(wp), allocatable :: taub_x(:)              ! [Pa] Basal stress
        real(wp), allocatable :: taulon_x(:)            ! [Pa] Longitudinal resistance
        real(wp), allocatable :: taulat_x(:)            ! [Pa] Lateral resistance
        real(wp), allocatable :: taures_x(:)            ! [Pa] Remainder of the balance
        real(wp), allocatable :: wt_x(:)                ! [--] Grounded-face weight, f_grnd_acx

        nx = ylmo%tpo%par%nx
        ny = ylmo%tpo%par%ny
        dx = ylmo%tpo%par%dx
        dy = ylmo%tpo%par%dy
        x0 = ylmo%grd%G%x(1)
        j0 = minloc(abs(ylmo%grd%G%y),1)                ! Centreline row (y=0)

        sec_year = ylmo%bnd%c%sec_year
        rho_ice  = ylmo%bnd%c%rho_ice

        allocate(wt(nx,ny),flux(nx,ny),Q_strn_int(nx,ny),T_ice_bar(nx,ny))
        allocate(taud_x(nx),taub_x(nx),taulon_x(nx),taulat_x(nx),taures_x(nx),wt_x(nx))

        ! Grounded-ice weight, flux and column-integrated strain heating
        wt = 0.0_wp
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp) wt = ylmo%tpo%now%f_grnd*ylmo%tpo%now%f_ice

        flux = ylmo%tpo%now%H_ice*ylmo%dyn%now%uxy_bar

        do j = 1, ny
        do i = 1, nx
            ! [J a-1 m-3] => [mW m-2]
            Q_strn_int(i,j) = ylmo%tpo%now%H_ice(i,j) &
                    * integrate_trapezoid1D_pt(ylmo%thrm%now%Q_strn(i,j,:),ylmo%par%zeta_aa) * 1e3_wp/sec_year
            T_ice_bar(i,j)  = 0.0_wp
            if (ylmo%tpo%now%H_ice(i,j) .gt. 0.0_wp) T_ice_bar(i,j) = &
                    integrate_trapezoid1D_pt(ylmo%thrm%now%T_ice(i,j,:),ylmo%par%zeta_aa) - ylmo%bnd%c%T0
        end do
        end do

        ! x-momentum balance along the centreline
        call calc_xmom_terms_row(taud_x,taub_x,taulon_x,taulat_x,taures_x, &
                                 ylmo%dyn%now%ux_bar,ylmo%dyn%now%uy_bar,ylmo%dyn%now%visc_eff_int, &
                                 ylmo%dyn%now%taud_acx,ylmo%dyn%now%taub_acx,j0,dx,dy)
        wt_x = ylmo%tpo%now%f_grnd_acx(:,j0)
        wt_x(1)       = 0.0_wp
        wt_x(nx-1:nx) = 0.0_wp

        ! Scalars
        A_grnd = sum(wt)*dx*dy*1e-6_wp                              ! [km2]
        V_grnd = sum(wt*ylmo%tpo%now%H_ice)*dx*dy*1e-9_wp           ! [km3]
        calv   = -sum(ylmo%tpo%now%cmb)*dx*dy*rho_ice*1e-12_wp      ! [Gt a-1]

        uxy_max_grnd = 0.0_wp
        if (any(wt .gt. 0.0_wp)) uxy_max_grnd = maxval(ylmo%dyn%now%uxy_bar,mask=wt .gt. 0.0_wp)
        uxy_max_flt  = 0.0_wp
        if (any(ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp)) &
            uxy_max_flt = maxval(ylmo%dyn%now%uxy_bar, &
                                 mask=ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp)

        ! Grounding line on the centreline: first crossing of f_grnd = 0.5 from upstream
        x_gl = 0.0_wp
        do i = 1, nx-1
            if (ylmo%tpo%now%f_grnd(i,j0) .ge. 0.5_wp .and. ylmo%tpo%now%f_grnd(i+1,j0) .lt. 0.5_wp) then
                x_gl = (ylmo%grd%G%x(i) + dx*(ylmo%tpo%now%f_grnd(i,j0)-0.5_wp) &
                            / (ylmo%tpo%now%f_grnd(i,j0)-ylmo%tpo%now%f_grnd(i+1,j0)))*1e-3_wp
                exit
            end if
        end do

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"A_grnd",A_grnd,units="km2",long_name="Grounded ice area", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_grnd",V_grnd,units="km3",long_name="Grounded ice volume", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"x_gl",x_gl,units="km",long_name="Grounding-line position on the centreline", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"calv",calv,units="Gt a-1",long_name="Calving flux", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_max_grnd",uxy_max_grnd,units="m/a",long_name="Maximum speed, grounded ice", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_max_flt",uxy_max_flt,units="m/a",long_name="Maximum speed, floating ice", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"ssa_lim_n",ylmo%dyn%par%ssa_lim_n,units="1",long_name="Faces at the velocity limit", &
                      dim1="time",start=[n],ncid=ncid)

        ! aa-node fields
        call ts_write_aa(filename,ncid,n,"H_ice",ylmo%tpo%now%H_ice,"m","Ice thickness", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"z_srf",ylmo%tpo%now%z_srf,"m","Surface elevation", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"f_grnd",ylmo%tpo%now%f_grnd,"1","Grounded fraction", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"uxy_bar",ylmo%dyn%now%uxy_bar,"m/a","Vertically averaged velocity magnitude", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"uxy_b",ylmo%dyn%now%uxy_b,"m/a","Basal sliding velocity magnitude", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"uxy_s",ylmo%dyn%now%uxy_s,"m/a","Surface velocity magnitude", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"flux",flux,"m2/a","Ice flux magnitude", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"taud",ylmo%dyn%now%taud,"Pa","Driving stress magnitude", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"taub",ylmo%dyn%now%taub,"Pa","Basal stress magnitude", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"N_eff",ylmo%dyn%now%N_eff,"Pa","Effective pressure", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"c_bed",ylmo%dyn%now%c_bed,"Pa","Bed friction coefficient (till yield stress)", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"W_til",ylmo%hyd%now%W_til,"m","Till water thickness", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"f_pmp",ylmo%thrm%now%f_pmp,"1","Fraction of grid point at pmp", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"T_prime_b",ylmo%thrm%now%T_prime_b,"deg C","Homologous basal ice temperature", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"bmb_grnd",ylmo%thrm%now%bmb_grnd,"m/a","Grounded basal mass balance", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"Q_b",ylmo%thrm%now%Q_b,"mW m-2","Basal frictional heating", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"Q_strn_int",Q_strn_int,"mW m-2","Column-integrated strain heating", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"Q_ice_b",ylmo%thrm%now%Q_ice_b,"mW m-2","Basal ice heat flux (positive up)", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"T_ice_bar",T_ice_bar,"degC","Vertically averaged ice temperature", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"ATT_bar",ylmo%mat%now%ATT_bar,"a^-1 Pa^-3","Vertically averaged rate factor", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"visc_eff_int",ylmo%dyn%now%visc_eff_int,"Pa a m","Depth-integrated effective viscosity", &
                         pts_x,sec_x,wt,j0,x0,dx)
        call ts_write_aa(filename,ncid,n,"de",ylmo%dyn%now%strn2D%de,"a^-1","Vertically averaged effective strain rate", &
                         pts_x,sec_x,wt,j0,x0,dx)

        ! x-momentum balance (acx-faces)
        call ts_write_acx(filename,ncid,n,"ux_bar_x",ylmo%dyn%now%ux_bar(:,j0),"m/a","Vertically averaged velocity (x)", &
                          pts_x,wt_x,x0,dx)
        call ts_write_acx(filename,ncid,n,"taud_x",taud_x,"Pa","Driving stress (x, flow direction)", &
                          pts_x,wt_x,x0,dx)
        call ts_write_acx(filename,ncid,n,"taub_x",taub_x,"Pa","Basal stress (x)", &
                          pts_x,wt_x,x0,dx)
        call ts_write_acx(filename,ncid,n,"taulon_x",taulon_x,"Pa","Longitudinal resistance (x)", &
                          pts_x,wt_x,x0,dx)
        call ts_write_acx(filename,ncid,n,"taulat_x",taulat_x,"Pa","Lateral resistance (x)", &
                          pts_x,wt_x,x0,dx)
        call ts_write_acx(filename,ncid,n,"taures_x",taures_x,"Pa","Remainder of the x-momentum balance", &
                          pts_x,wt_x,x0,dx)

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine trough_ts_write

    subroutine ts_write_aa(filename,ncid,n,name,var,units,long_name,pts_x,sec_x,wt,j0,x0,dx)
        ! aa-node field: points, cross-section, centreline and grounded-area means

        implicit none

        character(len=*), intent(IN) :: filename
        integer,          intent(IN) :: ncid
        integer,          intent(IN) :: n                ! Time index
        character(len=*), intent(IN) :: name
        real(wp),         intent(IN) :: var(:,:)
        character(len=*), intent(IN) :: units
        character(len=*), intent(IN) :: long_name
        real(wp),         intent(IN) :: pts_x(:)         ! [km] Points on the centreline (y=0)
        real(wp),         intent(IN) :: sec_x            ! [km] Position of the cross-section
        real(wp),         intent(IN) :: wt(:,:)          ! [--] Grounded-ice weight
        integer,          intent(IN) :: j0               ! Centreline row
        real(wp),         intent(IN) :: x0               ! [m] x of the first aa-node
        real(wp),         intent(IN) :: dx               ! [m]

        ! Local variables
        integer  :: k, jj, ny
        real(wp) :: v_pt(size(pts_x))
        real(wp) :: v_sec(size(var,2))
        real(wp) :: v_cl, v_gr

        ny = size(var,2)

        do k = 1, size(pts_x)
            v_pt(k) = interp_x(var(:,j0),pts_x(k)*1e3_wp,x0,dx)
        end do

        do jj = 1, ny
            v_sec(jj) = interp_x(var(:,jj),sec_x*1e3_wp,x0,dx)
        end do

        v_cl = 0.0_wp
        if (sum(wt(:,j0)) .gt. 0.0_wp) v_cl = sum(wt(:,j0)*var(:,j0)) / sum(wt(:,j0))

        v_gr = 0.0_wp
        if (sum(wt) .gt. 0.0_wp) v_gr = sum(wt*var) / sum(wt)

        call nc_write(filename,name,v_pt,units=units,long_name=long_name//" (points)", &
                      dim1="pt",dim2="time",start=[1,n],count=[size(pts_x),1],ncid=ncid)
        call nc_write(filename,name//"_sec",v_sec,units=units,long_name=long_name//" (cross-section)", &
                      dim1="yc",dim2="time",start=[1,n],count=[ny,1],ncid=ncid)
        call nc_write(filename,name//"_cl",v_cl,units=units,long_name=long_name//" (centreline grounded mean)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,name//"_gr",v_gr,units=units,long_name=long_name//" (grounded-area mean)", &
                      dim1="time",start=[n],ncid=ncid)

        return

    end subroutine ts_write_aa

    subroutine ts_write_acx(filename,ncid,n,name,var,units,long_name,pts_x,wt_x,x0,dx)
        ! acx-face values along the centreline: points and grounded mean

        implicit none

        character(len=*), intent(IN) :: filename
        integer,          intent(IN) :: ncid
        integer,          intent(IN) :: n                ! Time index
        character(len=*), intent(IN) :: name
        real(wp),         intent(IN) :: var(:)
        character(len=*), intent(IN) :: units
        character(len=*), intent(IN) :: long_name
        real(wp),         intent(IN) :: pts_x(:)         ! [km] Points on the centreline (y=0)
        real(wp),         intent(IN) :: wt_x(:)          ! [--] Grounded-face weight
        real(wp),         intent(IN) :: x0               ! [m] x of the first aa-node
        real(wp),         intent(IN) :: dx               ! [m]

        ! Local variables
        integer  :: k
        real(wp) :: v_pt(size(pts_x))
        real(wp) :: v_cl

        do k = 1, size(pts_x)
            v_pt(k) = interp_x(var,pts_x(k)*1e3_wp,x0+0.5_wp*dx,dx)
        end do

        v_cl = 0.0_wp
        if (sum(wt_x) .gt. 0.0_wp) v_cl = sum(wt_x*var) / sum(wt_x)

        call nc_write(filename,name,v_pt,units=units,long_name=long_name//" (points)", &
                      dim1="pt",dim2="time",start=[1,n],count=[size(pts_x),1],ncid=ncid)
        call nc_write(filename,name//"_cl",v_cl,units=units,long_name=long_name//" (centreline grounded mean)", &
                      dim1="time",start=[n],ncid=ncid)

        return

    end subroutine ts_write_acx

    subroutine calc_xmom_terms_row(taud_x,taub_x,taulon_x,taulat_x,taures_x, &
                                   ux,uy,visc_int,taud_acx,taub_acx,j,dx,dy)
        ! x-momentum balance of the depth-integrated stress balance at the
        ! acx-faces (i+1/2,j) of row j, with positive values in the +x
        ! (flow) direction:
        !
        !     taud = taub + taulon + taulat + taures
        !
        ! with taud = -taud_acx (driving stress), taub = taub_acx (basal),
        ! taulon = -d/dx[2*nuH*(2*dudx+dvdy)] (longitudinal resistance) and
        ! taulat = -d/dy[nuH*(dudy+dvdx)] (lateral resistance), nuH = visc_int.
        ! The membrane stresses use the C-grid stencil of the SSA solver
        ! (normal stresses on aa-nodes, shear stresses on ab-nodes). taures
        ! holds the rest, e.g. the velocity-limit drag and the difference
        ! from the viscosity staggering used by the solver.
        ! Faces at i = 1 and i >= nx-1 are set to zero.

        implicit none

        real(wp), intent(OUT) :: taud_x(:)              ! [Pa] Driving stress
        real(wp), intent(OUT) :: taub_x(:)              ! [Pa] Basal stress
        real(wp), intent(OUT) :: taulon_x(:)            ! [Pa] Longitudinal resistance
        real(wp), intent(OUT) :: taulat_x(:)            ! [Pa] Lateral resistance
        real(wp), intent(OUT) :: taures_x(:)            ! [Pa] Remainder
        real(wp), intent(IN)  :: ux(:,:)                ! [m a-1] Vertically averaged velocity (acx-nodes)
        real(wp), intent(IN)  :: uy(:,:)                ! [m a-1] Vertically averaged velocity (acy-nodes)
        real(wp), intent(IN)  :: visc_int(:,:)          ! [Pa a m] Depth-integrated effective viscosity (aa-nodes)
        real(wp), intent(IN)  :: taud_acx(:,:)          ! [Pa] Driving stress (acx-nodes)
        real(wp), intent(IN)  :: taub_acx(:,:)          ! [Pa] Basal stress (acx-nodes)
        integer,  intent(IN)  :: j                      ! Row index
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy

        ! Local variables
        integer  :: i, ii, jj, k, nx
        real(wp) :: Txx(2)                              ! [Pa m] Normal stress at aa-nodes i, i+1
        real(wp) :: Txy(2)                              ! [Pa m] Shear stress at ab-nodes (i+1/2,j-1/2), (i+1/2,j+1/2)
        real(wp) :: visc_ab

        nx = size(ux,1)

        taud_x   = 0.0_wp
        taub_x   = 0.0_wp
        taulon_x = 0.0_wp
        taulat_x = 0.0_wp
        taures_x = 0.0_wp

        do i = 2, nx-2

            do k = 1, 2
                ii = i+k-1
                Txx(k) = 2.0_wp*visc_int(ii,j)*( 2.0_wp*(ux(ii,j)-ux(ii-1,j))/dx + (uy(ii,j)-uy(ii,j-1))/dy )
            end do

            do k = 1, 2
                jj = j+k-2
                visc_ab = 0.25_wp*(visc_int(i,jj)+visc_int(i+1,jj)+visc_int(i,jj+1)+visc_int(i+1,jj+1))
                Txy(k)  = visc_ab*( (ux(i,jj+1)-ux(i,jj))/dy + (uy(i+1,jj)-uy(i,jj))/dx )
            end do

            taud_x(i)   = -taud_acx(i,j)
            taub_x(i)   =  taub_acx(i,j)
            taulon_x(i) = -(Txx(2)-Txx(1))/dx
            taulat_x(i) = -(Txy(2)-Txy(1))/dy
            taures_x(i) = taud_x(i) - taub_x(i) - taulon_x(i) - taulat_x(i)

        end do

        return

    end subroutine calc_xmom_terms_row

    function interp_x(var,x,x0,dx) result(var_x)
        ! Linear interpolation of var, defined at x0 + (i-1)*dx, to position x
        ! (held within the end points)

        implicit none

        real(wp), intent(IN) :: var(:)
        real(wp), intent(IN) :: x
        real(wp), intent(IN) :: x0
        real(wp), intent(IN) :: dx
        real(wp) :: var_x

        ! Local variables
        integer  :: i
        real(wp) :: w

        i = floor((x-x0)/dx) + 1
        i = max(1,min(size(var)-1,i))
        w = (x - (x0+real(i-1,wp)*dx)) / dx
        w = max(0.0_wp,min(1.0_wp,w))

        var_x = (1.0_wp-w)*var(i) + w*var(i+1)

        return

    end function interp_x

    ! == Analytical solution by Schoof 2006 for the "SSA_icestream" benchmark experiment
  elemental subroutine SSA_Schoof2006_analytical_solution_yelmo(u, tauc, y, tantheta, h0, A_flow, W, m,  &
                            n_glen,ice_density, grav)
      
    implicit none
    
    ! In/output variables:
    real(wp),                            intent(OUT)   :: u             ! Ice velocity in the x-direction
    real(wp),                            intent(OUT)   :: tauc          ! Till yield stress
    real(dp),                            intent(IN)    :: y             ! y-coordinate
    real(wp),                            intent(IN)    :: tantheta      ! Surface slope in the x-direction
    real(wp),                            intent(IN)    :: h0            ! Ice thickness
    real(wp),                            intent(IN)    :: A_flow        ! Ice flow factor
    real(wp),                            intent(IN)    :: W             ! Ice-stream half width (m)
    real(wp),                            intent(IN)    :: m             ! Ice stream exponent
    real(wp),                            intent(IN)    :: n_glen
    real(wp),                            intent(IN)    :: ice_density
    real(wp),                            intent(IN)    :: grav
    
    ! Local variables:
    real(dp) :: B, f, L, ua, ub, uc, ud, ue   
    real(dp) :: ux, taud, H, yy 

    ! Calculate the gravitational driving stress f
    f = ice_density * grav * h0 * tantheta
    
    ! Calculate the ice hardness factor B
    B = A_flow**(-1._dp/n_glen)
    
    ! Determine constant L (ice-stream width)
    L = W / ((1.0_dp+m)**(1.0_dp/m))

    ! Calculate the till yield stress across the stream
    tauc = f * ABS(y/L)**m
    
    taud = f 
    H    = h0
    yy   = y 
    ux = -2.0 * taud**3 * L**4 / (B**3 * H**3) * ( ((yy/L)**4 - (m+1.0)**(4.0/m))/4.0 - 3.0*( abs(yy/L)**(m+4.0) &
    - (m+1.0)**(1.0+4.0/m) )/((m+1.0)*(m+4.0)) + 3.0*( abs(yy/L)**(2.0*m+4.0) - (m+1.0)**(2.0+4.0/m) )/((m+1.0)**2*(2.0*m+4.0)) &
    - ( abs(yy/L)**(3.0*m+4.0) - (m+1.0)**(3.0+4.0/m) )/ ( (m+1.0)**3*(3.0*m+4.0)) )

if (.TRUE.) then

    u = ux 

else 

    ! Calculate the analytical solution for u
    ua = -2._dp * f**3 * L**4 / (B**3 * h0**3)
    ub = ( 1._dp / 4._dp                           ) * (   (y/L)**     4._dp  - (m+1._dp)**(       4._dp/m) )
    uc = (-3._dp / ((m+1._dp)    * (      m+4._dp))) * (ABS(y/L)**(  m+4._dp) - (m+1._dp)**(1._dp+(4._dp/m)))
    ud = ( 3._dp / ((m+1._dp)**2 * (2._dp*m+4._dp))) * (ABS(y/L)**(2*m+4._dp) - (m+1._dp)**(2._dp+(4._dp/m)))
    ue = (-1._dp / ((m+1._dp)**3 * (3._dp*m+4._dp))) * (ABS(y/L)**(3*m+4._dp) - (m+1._dp)**(3._dp+(4._dp/m)))
    u = ua * (ub + uc + ud + ue)
    
end if

    ! Outside the ice-stream, velocity is zero
    IF (ABS(y) > W) u = 0._dp
    
    end subroutine SSA_Schoof2006_analytical_solution_yelmo

  ! == Analytical solution by Schoof 2006 for the "SSA_icestream" benchmark experiment
  ELEMENTAL SUBROUTINE SSA_Schoof2006_analytical_solution(U, tauc, y, tantheta, h0, A_flow, L, m,  &
                            n_glen,ice_density, grav)
      
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(wp),                            INTENT(OUT)   :: U             ! Ice velocity in the x-direction
    REAL(wp),                            INTENT(OUT)   :: tauc          ! Till yield stress
    REAL(wp),                            INTENT(IN)    :: y             ! y-coordinate
    REAL(wp),                            INTENT(IN)    :: tantheta      ! Surface slope in the x-direction
    REAL(wp),                            INTENT(IN)    :: h0            ! Ice thickness
    REAL(wp),                            INTENT(IN)    :: A_flow        ! Ice flow factor
    REAL(wp),                            INTENT(IN)    :: L             ! Ice-stream width (m), default 40e3
    REAL(wp),                            INTENT(IN)    :: m             ! Ice stream exponent
    REAL(wp),                            INTENT(IN)    :: n_glen
    REAL(wp),                            INTENT(IN)    :: ice_density
    REAL(wp),                            INTENT(IN)    :: grav
    
    ! Local variables:
    REAL(wp)                                           :: B, f, W, ua, ub, uc, ud, ue   
    
    ! Calculate the gravitational driving stress f
    f = ice_density * grav * h0 * tantheta
    
    ! Calculate the ice hardness factor B
    B = A_flow**(-1._dp/n_glen)
    
    ! Calculate the "ice stream half-width" W
    W = L * (m+1._dp)**(1._dp/m)
    
    ! Calculate the till yield stress across the stream
    tauc = f * ABS(y/L)**m
    
    ! Calculate the analytical solution for u
    ua = -2._dp * f**3 * L**4 / (B**3 * h0**3)
    ub = ( 1._dp / 4._dp                           ) * (   (y/L)**     4._dp  - (m+1._dp)**(       4._dp/m) )
    uc = (-3._dp / ((m+1._dp)    * (      m+4._dp))) * (ABS(y/L)**(  m+4._dp) - (m+1._dp)**(1._dp+(4._dp/m)))
    ud = ( 3._dp / ((m+1._dp)**2 * (2._dp*m+4._dp))) * (ABS(y/L)**(2*m+4._dp) - (m+1._dp)**(2._dp+(4._dp/m)))
    ue = (-1._dp / ((m+1._dp)**3 * (3._dp*m+4._dp))) * (ABS(y/L)**(3*m+4._dp) - (m+1._dp)**(3._dp+(4._dp/m)))
    u = ua * (ub + uc + ud + ue)
    
    ! Outside the ice-stream, velocity is zero
    IF (ABS(y) > w) U = 0._dp
    
  END SUBROUTINE SSA_Schoof2006_analytical_solution

end program yelmo_trough



