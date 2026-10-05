program yelmo_bench
    ! Driver for the benchmark protocol (docs/dev/benchmark-protocol).
    !
    ! The driver holds no geometry. It reads a fixture generated in Julia
    ! (tests/bench, e.g. make_fixture.jl island4) with the grid axes, the
    ! boundary fields and the initial ice thickness, and runs Yelmo with the
    ! parameters of the parameter file. Forcing that depends on the evolving
    ! state is computed online (tests/bench_forcing.f90), selected in &bench.
    ! A run with time_end = time_init is diagnostic: it writes the initial
    ! state, including the initial velocity solution, and stops.

    use nml
    use ncio
    use yelmo
    use lsf_module, only : LSFinit
    use yelmo_grid, only : yelmo_init_grid_fromaxes
    use bench_forcing
    use, intrinsic :: iso_fortran_env, only : int64

    implicit none

    type ctrl_type
        character(len=512) :: path_par
        character(len=256) :: file2D, file1D, file_restart
        real(wp) :: time_init, time_end, dtt
        real(wp) :: dt2D_out, dt1D_out
    end type

    type bench_type
        character(len=512) :: fixture       ! Path to the fixture (NetCDF)
        character(len=56)  :: tsrf_method   ! "fixed" (from fixture) | "lapse"
        character(len=56)  :: bmb_method    ! "fixed" (from fixture) | "mismip+"
        real(wp) :: T_sl                    ! [degC] lapse: temperature at sea level
        real(wp) :: lapse                   ! [K m-1] lapse: lapse rate
        real(wp) :: Omega                   ! [yr-1] mismip+: melt-rate factor
        real(wp) :: Hc0                     ! [m] mismip+: water-column scale
        real(wp) :: z0                      ! [m] mismip+: depth above which melt is zero
    end type

    type(yelmo_class) :: yelmo1
    type(ctrl_type)   :: ctl
    type(bench_type)  :: bch

    character(len=56)  :: domain, grid_name
    real(wp), allocatable :: xc(:), yc(:)
    integer,  allocatable :: mask_ice(:,:)
    real(wp) :: time
    integer  :: n, nx, ny

    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime

    call yelmo_cpu_time(cpu_start_time)

    ! === Parameters =====================================================

    call yelmo_load_command_line_args(ctl%path_par)

    ctl%file1D       = "yelmo_ts.nc"
    ctl%file2D       = "yelmo.nc"
    ctl%file_restart = "yelmo_restart.nc"

    call nml_read(ctl%path_par,"ctrl","time_init",  ctl%time_init)     ! [yr] Starting time
    call nml_read(ctl%path_par,"ctrl","time_end",   ctl%time_end)      ! [yr] Ending time
    call nml_read(ctl%path_par,"ctrl","dtt",        ctl%dtt)           ! [yr] Main loop time step
    call nml_read(ctl%path_par,"ctrl","dt2D_out",   ctl%dt2D_out)      ! [yr] Frequency of 2D output
    ctl%dt1D_out = ctl%dtt

    call nml_read(ctl%path_par,"bench","fixture",     bch%fixture)
    call nml_read(ctl%path_par,"bench","tsrf_method", bch%tsrf_method)
    call nml_read(ctl%path_par,"bench","bmb_method",  bch%bmb_method)
    call nml_read(ctl%path_par,"bench","T_sl",        bch%T_sl)
    call nml_read(ctl%path_par,"bench","lapse",       bch%lapse)
    call nml_read(ctl%path_par,"bench","Omega",       bch%Omega)
    call nml_read(ctl%path_par,"bench","Hc0",         bch%Hc0)
    call nml_read(ctl%path_par,"bench","z0",          bch%z0)

    call yelmo_check_enum("bench","tsrf_method",bch%tsrf_method,"fixed|lapse")
    call yelmo_check_enum("bench","bmb_method", bch%bmb_method, "fixed|mismip+")

    call nml_read(ctl%path_par,"yelmo","domain",    domain)
    call nml_read(ctl%path_par,"yelmo","grid_name", grid_name)

    ! === Grid and initialization ========================================

    ! Grid axes from the fixture [km]
    nx = nc_size(bch%fixture,"xc")
    ny = nc_size(bch%fixture,"yc")
    allocate(xc(nx),yc(ny),mask_ice(nx,ny))
    call nc_read(bch%fixture,"xc",xc)
    call nc_read(bch%fixture,"yc",yc)
    call yelmo_init_grid_fromaxes(yelmo1%grd,grid_name,xc*1e3_wp,yc*1e3_wp)

    call nc_read(bch%fixture,"mask_ice",mask_ice)

    call yelmo_init(yelmo1,filename=ctl%path_par,grid_def="none",time=ctl%time_init, &
                    load_topo=.FALSE.,domain=domain,grid_name=grid_name,mask_ice=mask_ice)

    ! Boundary fields
    call nc_read(bch%fixture,"z_bed",    yelmo1%bnd%z_bed)
    call nc_read(bch%fixture,"z_sl",     yelmo1%bnd%z_sl)
    call nc_read(bch%fixture,"smb_ref",  yelmo1%bnd%smb)
    call nc_read(bch%fixture,"T_srf",    yelmo1%bnd%T_srf)
    call nc_read(bch%fixture,"Q_geo",    yelmo1%bnd%Q_geo)
    call nc_read(bch%fixture,"bmb_shlf", yelmo1%bnd%bmb_shlf)
    call nc_read(bch%fixture,"T_shlf",   yelmo1%bnd%T_shlf)
    call nc_read(bch%fixture,"H_sed",    yelmo1%bnd%H_sed)

    if (.not. yelmo1%par%use_restart) then
        call nc_read(bch%fixture,"H_ice",yelmo1%tpo%now%H_ice)
        call LSFinit(yelmo1%tpo%now%lsf,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%tpo%par%dx)
    end if

    call yelmo_print_bound(yelmo1%bnd)

    ! The fixture holds the online forcing evaluated on the initial geometry
    ! (Julia counterparts), which initializes the thermodynamics.
    call yelmo_init_state(yelmo1,time=ctl%time_init,thrm_method="robin")

    ! Compare the online forcing with the fixture values (fresh start only)
    if (.not. yelmo1%par%use_restart) call bench_check_forcing(yelmo1,bch)

    ! === Initial output =================================================

    time = ctl%time_init

    call yelmo_write_init(yelmo1,ctl%file2D,time_init=ctl%time_init,units="years")
    call yelmo_write_step(yelmo1,ctl%file2D,time=ctl%time_init)

    call yelmo_write_reg_init(yelmo1,ctl%file1D,time_init=ctl%time_init,units="years", &
                                                mask=(yelmo1%bnd%mask_ice /= MASK_ICE_NONE))
    call yelmo_write_reg_step(yelmo1,ctl%file1D,time=ctl%time_init)

    ! === Time loop ======================================================

    do n = 1, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

        time = ctl%time_init + n*ctl%dtt

        call bench_update_forcing(yelmo1,bch)

        call yelmo_update(yelmo1,time)

        ! int64: a default integer overflows for |time| > ~2.1e7 yr
        if (mod(nint(time*100,int64),nint(ctl%dt2D_out*100,int64))==0) then
            call yelmo_write_step(yelmo1,ctl%file2D,time=time)
        end if

        if (mod(nint(time*100,int64),nint(ctl%dt1D_out*100,int64))==0) then
            call yelmo_write_reg_step(yelmo1,ctl%file1D,time=time)
        end if

    end do

    if (ctl%time_end .gt. ctl%time_init) then
        call yelmo_restart_write(yelmo1,ctl%file_restart,time=time)
    end if

    call yelmo_end(yelmo1,time=time)

    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    if (ctl%time_end .gt. ctl%time_init) then
        write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctl%time_end-ctl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
    end if

contains

    subroutine bench_update_forcing(ylmo,bch)
        ! Update the online forcing from the current state.

        implicit none

        type(yelmo_class), intent(INOUT) :: ylmo
        type(bench_type),  intent(IN)    :: bch

        select case(trim(bch%tsrf_method))
            case("fixed")
                ! Keep T_srf from the fixture
            case("lapse")
                call bench_tsrf_lapse(ylmo%bnd%T_srf,ylmo%tpo%now%z_srf,bch%T_sl,bch%lapse,ylmo%bnd%c%T0)
        end select

        select case(trim(bch%bmb_method))
            case("fixed")
                ! Keep bmb_shlf from the fixture
            case("mismip+")
                call bench_bmb_mismipplus(ylmo%bnd%bmb_shlf,ylmo%tpo%now%H_ice,ylmo%bnd%z_bed,ylmo%bnd%z_sl, &
                                          ylmo%bnd%c%rho_ice,ylmo%bnd%c%rho_sw,bch%Omega,bch%Hc0,bch%z0)
        end select

        return

    end subroutine bench_update_forcing

    subroutine bench_check_forcing(ylmo,bch)
        ! Print the maximum difference between the online forcing on the
        ! initial state and the fixture values, which the Julia counterparts
        ! computed on the same geometry. Partially ice-covered cells
        ! (0 < f_ice < 1) are excluded, since Yelmo defines their surface
        ! elevation with the effective front thickness.

        implicit none

        type(yelmo_class), intent(IN) :: ylmo
        type(bench_type),  intent(IN) :: bch

        real(wp), allocatable :: var(:,:)
        logical,  allocatable :: mask(:,:)

        allocate(var(size(ylmo%bnd%z_bed,1),size(ylmo%bnd%z_bed,2)))
        mask = ylmo%tpo%now%f_ice .eq. 0.0_wp .or. ylmo%tpo%now%f_ice .eq. 1.0_wp

        write(*,"(a,i0,a)") "bench:: forcing check excludes ", count(.not. mask), " partially ice-covered cells"

        if (trim(bch%tsrf_method) .eq. "lapse") then
            call bench_tsrf_lapse(var,ylmo%tpo%now%z_srf,bch%T_sl,bch%lapse,ylmo%bnd%c%T0)
            write(*,"(a,g12.4,a)") "bench:: T_srf online - fixture, max abs diff    = ", &
                                    maxval(abs(var-ylmo%bnd%T_srf),mask=mask), " K"
        end if

        if (trim(bch%bmb_method) .eq. "mismip+") then
            call bench_bmb_mismipplus(var,ylmo%tpo%now%H_ice,ylmo%bnd%z_bed,ylmo%bnd%z_sl, &
                                      ylmo%bnd%c%rho_ice,ylmo%bnd%c%rho_sw,bch%Omega,bch%Hc0,bch%z0)
            write(*,"(a,g12.4,a)") "bench:: bmb_shlf online - fixture, max abs diff = ", &
                                    maxval(abs(var-ylmo%bnd%bmb_shlf),mask=mask), " m/yr"
        end if

        return

    end subroutine bench_check_forcing

end program yelmo_bench
