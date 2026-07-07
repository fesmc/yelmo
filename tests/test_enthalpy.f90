program test_enthalpy
    ! Standalone single-column thermodynamics validation driver for Yelmo.
    !
    ! Purpose: exercise and validate the enthalpy solver (calc_enth_column)
    ! against the temperature solver (calc_temp_column) and the Kleiner et al.
    ! (2015) enthalpy benchmarks. See docs/physics/enthalpy-transition-plan.md.
    !
    ! Experiments (1st cmdline arg, default "cold-limit"):
    !   cold-limit : identical cold column solved by both solvers; asserts that
    !                enth reduces to temp in the cold limit (cr=1, no temperate ice).
    !   kleiner-a  : Kleiner (2015) Exp A - transient basal melt under a
    !                time-varying surface temperature, no flow. Reference:
    !                tests/data/Kleiner2015/Kleiner2015_EXPA_Fig2-IIIa-melt.txt
    !
    ! Usage: test_enthalpy.x [experiment] [solver] [nz]
    !   solver: "temp" | "enth" | "both" (default: experiment-dependent)
    !   nz    : number of aa-nodes (default 51)
    !
    ! Output: output/test_enthalpy_<experiment>_<solver>.nc

    use yelmo_defs,     only : wp, prec, ybound_const_class
    use ncio
    use yelmo_grid,     only : calc_zeta
    use thermodynamics, only : convert_to_enthalpy, calc_T_pmp, &
                               calc_specific_heat_capacity, calc_thermal_conductivity
    use ice_enthalpy,   only : calc_temp_column, calc_enth_column, calc_dzeta_terms

    implicit none

    ! Column state
    type column_class
        integer :: nz_aa, nz_ac
        real(wp), allocatable :: zeta_aa(:), zeta_ac(:), dzeta_a(:), dzeta_b(:)
        real(wp), allocatable :: enth(:), T_ice(:), omega(:), T_pmp(:)
        real(wp), allocatable :: cp(:), kt(:), advecxy(:), uz(:), Q_strn(:)
        real(wp) :: H_ice, T_srf, T_shlf, smb, Q_b, Q_rock, W_til, f_grnd
        real(wp) :: bmb, Q_ice_b, H_cts
    end type

    type(ybound_const_class) :: c
    character(len=256) :: experiment, solver, arg
    integer :: nz, narg
    real(wp) :: cr_arg

    ! Load physical constants (Kleiner 2015 / EISMINT values, set directly)
    call set_constants(c)

    ! Parse arguments
    experiment = "cold-limit"
    solver     = ""
    nz         = 51
    cr_arg     = -1.0_wp        ! <0 => use solver default
    narg = command_argument_count()
    if (narg .ge. 1) call get_command_argument(1,experiment)
    if (narg .ge. 2) call get_command_argument(2,solver)
    if (narg .ge. 3) then
        call get_command_argument(3,arg); read(arg,*) nz
    end if
    if (narg .ge. 4) then
        call get_command_argument(4,arg); read(arg,*) cr_arg
    end if

    ! Experiment-specific physical constants (Kleiner 2015 Table A1)
    if (trim(experiment) .eq. "kleiner-b") then
        c%T_pmp_beta = 0.0_wp          ! no pressure-melting dependence
        c%L_ice      = 3.35e5_wp       ! Exp B latent heat
    end if

    select case(trim(experiment))
        case("cold-limit")
            if (trim(solver) .eq. "") solver = "both"
            call run_experiment(c,"cold-limit","temp",nz,cr_arg)
            call run_experiment(c,"cold-limit","enth",nz,cr_arg)
            call compare_cold_limit(c,nz)
        case("kleiner-a")
            if (trim(solver) .eq. "") solver = "enth"
            call run_experiment(c,"kleiner-a",trim(solver),nz,cr_arg)
        case("kleiner-b")
            if (trim(solver) .eq. "") solver = "enth"
            call run_experiment(c,"kleiner-b",trim(solver),nz,cr_arg)
        case DEFAULT
            write(*,*) "test_enthalpy:: unknown experiment: ", trim(experiment)
            write(*,*) "  choose one of: cold-limit, kleiner-a, kleiner-b"
            stop 1
    end select

contains

    subroutine set_constants(c)
        ! Kleiner et al. (2015) constants (also consistent with EISMINT)
        type(ybound_const_class), intent(OUT) :: c
        c%sec_year   = 31556926.0_wp
        c%g          = 9.81_wp
        c%T0         = 273.15_wp
        c%rho_ice    = 910.0_wp
        c%rho_w      = 1000.0_wp
        c%rho_sw     = 1028.0_wp
        c%rho_rock   = 2000.0_wp
        c%L_ice      = 3.34e5_wp
        c%T_pmp_beta = 7.9e-8_wp       ! Kleiner (2015) Table A1 Clausius-Clapeyron
        return
    end subroutine set_constants

    subroutine column_alloc(col,nz)
        type(column_class), intent(INOUT) :: col
        integer, intent(IN) :: nz
        col%nz_aa = nz
        col%nz_ac = nz+1        ! calc_zeta convention: nz_ac = nz_aa + 1
        allocate(col%zeta_aa(nz), col%dzeta_a(nz), col%dzeta_b(nz))
        allocate(col%zeta_ac(nz+1))
        allocate(col%enth(nz), col%T_ice(nz), col%omega(nz), col%T_pmp(nz))
        allocate(col%cp(nz), col%kt(nz), col%advecxy(nz), col%Q_strn(nz))
        allocate(col%uz(nz+1))
        return
    end subroutine column_alloc

    subroutine setup_experiment(col,c,experiment,nz)
        ! Define the fixed column geometry, forcing and initial state.
        type(column_class),       intent(INOUT) :: col
        type(ybound_const_class), intent(IN)    :: c
        character(len=*),         intent(IN)    :: experiment
        integer,                  intent(IN)    :: nz

        integer :: k, nz_ac
        real(wp), allocatable :: zeta_aa(:), zeta_ac(:)
        real(wp) :: T_init, z, sin_gam, tau_gam
        real(wp), parameter :: A_glen = 5.3e-24_wp             ! [Pa-3 s-1] Exp B rate factor
        real(wp), parameter :: deg2rad = 3.14159265358979_wp/180.0_wp

        call column_alloc(col,nz)

        ! Vertical grid (linear, zeta==height, k=1 base .. k=nz surface)
        call calc_zeta(zeta_aa,zeta_ac,nz_ac,nz,zeta_scale="linear",zeta_exp=1.0_wp)
        col%zeta_aa = zeta_aa
        col%zeta_ac = zeta_ac
        call calc_dzeta_terms(col%dzeta_a,col%dzeta_b,col%zeta_aa,col%zeta_ac)

        ! Common defaults
        col%f_grnd  = 1.0_wp
        col%W_til   = 0.0_wp
        col%Q_b     = 0.0_wp        ! no basal frictional heating
        col%advecxy = 0.0_wp        ! no horizontal advection
        col%uz      = 0.0_wp        ! no vertical advection
        col%Q_strn  = 0.0_wp        ! no strain heating
        col%bmb     = 0.0_wp
        col%Q_ice_b = 0.0_wp
        col%H_cts   = 0.0_wp

        ! Constant material properties (Kleiner: cp, kt fixed)
        col%cp = 2009.0_wp                        ! [J kg-1 K-1]
        col%kt = 2.1_wp * c%sec_year              ! [W m-1 K-1] => [J a-1 m-1 K-1]

        select case(trim(experiment))
            case("cold-limit")
                ! Deliberately cold column: base stays well below pmp so no
                ! temperate ice forms and enth must reduce exactly to temp.
                col%H_ice  = 1500.0_wp
                col%T_srf  = c%T0 - 40.0_wp       ! -40 C
                col%smb    = 0.1_wp               ! [m a-1]
                col%Q_rock = 15.0_wp              ! [mW m-2] low geothermal -> cold base
            case("kleiner-a")
                col%H_ice  = 1000.0_wp
                col%T_srf  = surf_temp_kleiner_a(0.0_wp,c)
                col%smb    = 0.0_wp
                col%Q_rock = 42.0_wp              ! [mW m-2] (0.042 W m-2)
            case("kleiner-b")
                ! Steady polythermal column: 200 m slab inclined 4 deg, driven by
                ! strain heating with constant downward advection (Kleiner Exp B).
                col%H_ice  = 200.0_wp
                col%T_srf  = c%T0 - 3.0_wp
                col%smb    = 0.2_wp
                col%Q_rock = 0.0_wp               ! no geothermal flux
                ! Constant downward velocity vz = -a_s = -0.2 m/a (ac-nodes)
                col%uz = -0.2_wp
                ! Strain heating Psi(z) = 2 A (rho g sin gamma)^4 (H-z)^4 [W m-3]
                ! (from Eq 16-17); passed to the solver in [J a-1 m-3].
                sin_gam = sin(4.0_wp*deg2rad)
                tau_gam = c%rho_ice*c%g*sin_gam                ! [Pa m-1]
                do k = 1, nz
                    z = col%zeta_aa(k)*col%H_ice
                    col%Q_strn(k) = 2.0_wp*A_glen*tau_gam**4*(col%H_ice-z)**4 * c%sec_year
                end do
            case DEFAULT
                stop "setup_experiment:: unknown experiment"
        end select

        ! Pressure melting point profile
        do k = 1, nz
            col%T_pmp(k) = calc_T_pmp(col%H_ice,col%zeta_aa(k),c%T0,c%T_pmp_beta,c%rho_ice,c%g)
        end do

        ! Initial state: uniform temperature (Exp B starts at -1.5 C; others cold)
        if (trim(experiment) .eq. "kleiner-b") then
            T_init = c%T0 - 1.5_wp
        else
            T_init = col%T_srf
        end if
        do k = 1, nz
            col%T_ice(k) = T_init
        end do
        col%omega = 0.0_wp
        call convert_to_enthalpy(col%enth,col%T_ice,col%omega,col%T_pmp,col%cp,c%L_ice)

        col%T_shlf = col%T_pmp(1)

        return
    end subroutine setup_experiment

    function surf_temp_kleiner_a(time,c) result(T_srf)
        ! Kleiner (2015) Exp A surface-temperature forcing (their phases I-III):
        !   Ts = -30 C for 0    <= t < 100 ka   (initial, cold)
        !   Ts = -10 C for 100  <= t < 150 ka   (warming)
        !   Ts = -30 C for t   >= 150 ka        (cooling)
        real(wp), intent(IN) :: time      ! [a]
        type(ybound_const_class), intent(IN) :: c
        real(wp) :: T_srf
        real(wp) :: tka
        tka = time / 1000.0_wp
        if (tka .lt. 100.0_wp) then
            T_srf = c%T0 - 30.0_wp
        else if (tka .lt. 150.0_wp) then
            T_srf = c%T0 - 10.0_wp
        else
            T_srf = c%T0 - 30.0_wp
        end if
        return
    end function surf_temp_kleiner_a

    subroutine run_experiment(c,experiment,solver,nz,cr_override)
        type(ybound_const_class), intent(IN) :: c
        character(len=*),         intent(IN) :: experiment, solver
        integer,                  intent(IN) :: nz
        real(wp),                 intent(IN) :: cr_override   ! <0 => use solver default

        type(column_class) :: col
        character(len=256) :: filename
        real(wp) :: time, time_end, dt, dt_out, time_out
        real(wp) :: enth_cr, omega_max, Q_lith
        real(wp) :: ab_warm, ab_cold, ab_warm_an, ab_cold_an
        integer  :: n

        ab_warm = -999.0_wp
        ab_cold = -999.0_wp


        ! Enthalpy solver parameters
        omega_max = 0.01_wp
        if (trim(solver) .eq. "temp") then
            enth_cr = 1.0_wp        ! temp <=> enth with cr=1, omega_max=0
        else
            enth_cr = 1.0e-3_wp
        end if
        ! Exp B holds a temperate layer with omega up to ~2%; do not clip it
        ! (Kleiner applies no water-content restriction). The analytic solution is
        ! the K0 -> 0 limit, so use a small conductivity ratio by default; the
        ! solution converges to the analytic as cr decreases (paper: CR 1e-1..1e-5).
        if (trim(experiment) .eq. "kleiner-b") then
            omega_max = 1.0_wp
            enth_cr   = 1.0e-4_wp
        end if
        if (cr_override .ge. 0.0_wp) enth_cr = cr_override

        ! Time control
        select case(trim(experiment))
            case("cold-limit")
                time_end = 50000.0_wp;  dt = 5.0_wp;  dt_out = 1000.0_wp
            case("kleiner-a")
                time_end = 300000.0_wp; dt = 1.0_wp;  dt_out = 100.0_wp
            case("kleiner-b")
                time_end = 50000.0_wp;  dt = 2.0_wp;  dt_out = 500.0_wp
            case DEFAULT
                time_end = 50000.0_wp;  dt = 5.0_wp;  dt_out = 1000.0_wp
        end select

        call setup_experiment(col,c,experiment,nz)

        write(filename,"(a)") "output/test_enthalpy_"//trim(experiment)//"_"//trim(solver)//".nc"
        call write_init(col,filename)

        time     = 0.0_wp
        time_out = 0.0_wp
        call write_step(col,filename,time)

        do while (time .lt. time_end - 1e-6_wp)

            ! Update time-varying surface forcing
            if (trim(experiment) .eq. "kleiner-a") then
                col%T_srf = surf_temp_kleiner_a(time,c)
            end if

            ! Solve the column
            select case(trim(solver))
                case("temp")
                    call calc_temp_column(col%enth,col%T_ice,col%omega,col%bmb,col%Q_ice_b, &
                            col%H_cts,col%T_pmp,col%cp,col%kt,col%advecxy,col%uz,col%Q_strn, &
                            col%Q_b,col%Q_rock,col%T_srf,col%T_shlf,col%H_ice,col%W_til,col%f_grnd, &
                            col%zeta_aa,col%zeta_ac,col%dzeta_a,col%dzeta_b,omega_max,c%T0, &
                            c%rho_ice,c%rho_w,c%L_ice,c%sec_year,dt)
                case("enth")
                    ! enth solver takes basal fluxes pre-converted to [J a-1 m-2]
                    Q_lith = col%Q_rock*1e-3_wp*c%sec_year
                    call calc_enth_column(col%enth,col%T_ice,col%omega,col%bmb,col%Q_ice_b, &
                            col%H_cts,col%T_pmp,col%cp,col%kt,col%advecxy,col%uz,col%Q_strn, &
                            col%Q_b,Q_lith,col%T_srf,col%T_shlf,col%H_ice,col%W_til,col%f_grnd, &
                            col%zeta_aa,col%zeta_ac,col%dzeta_a,col%dzeta_b,enth_cr,omega_max,c%T0, &
                            c%rho_ice,c%rho_w,c%L_ice,c%sec_year,dt)
                case DEFAULT
                    write(*,*) "run_experiment:: unknown solver: ", trim(solver); stop 1
            end select

            ! Basal water: accumulate freely with no drainage or cap, following
            ! Kleiner (2015) Exp A. Basal melt (bmb<0) adds water, freeze-on
            ! (bmb>0) removes it; floored at zero. Same water-equivalent
            ! conversion as the solver's W_til_predicted.
            col%W_til = max(0.0_wp, col%W_til - col%bmb*(c%rho_w/c%rho_ice)*dt)

            ! Sample steady basal melt (as +melt, mm/a) in the warm and cold phases
            if (trim(experiment) .eq. "kleiner-a") then
                if (abs(time-148000.0_wp) .lt. 0.5_wp*dt) ab_warm = -col%bmb*1000.0_wp
                if (abs(time-168000.0_wp) .lt. 0.5_wp*dt) ab_cold = -col%bmb*1000.0_wp
            end if

            time = time + dt

            if (time - time_out .ge. dt_out - 1e-6_wp) then
                call write_step(col,filename,time)
                time_out = time
            end if

        end do

        write(*,"(a,a,a,a,a,f10.2,a,es12.4)") "  [",trim(experiment),"/",trim(solver), &
                "] T_base = ", col%T_ice(1), " K, bmb = ", col%bmb

        ! T4: compare steady polythermal structure to the Kleiner (2015) Exp B
        ! analytic solution (CTS at 19 m, base water content 2.07%).
        if (trim(experiment) .eq. "kleiner-b" .and. trim(solver) .eq. "enth") then
            write(*,*) ""
            write(*,*) "=== Kleiner Exp B polythermal structure vs analytic (T4) ==="
            write(*,"(a,f7.2,a)")   "  CTS height = ", col%H_cts,   " m     (analytic 19.0 m)"
            write(*,"(a,f8.4,a)")   "  base omega = ", col%omega(1), "      (analytic 0.0207)"
            if (abs(col%H_cts-19.0_wp) .lt. 2.0_wp .and. &
                abs(col%omega(1)-0.0207_wp) .lt. 0.10_wp*0.0207_wp) then
                write(*,*) "  PASS: CTS within 2 m and base omega within 10% of analytic."
            else
                write(*,*) "  FAIL: exceeds tolerance (try larger nz / smaller cr)."
            end if
            write(*,*) ""
        end if

        ! T3: compare warm/cold steady melt to the Kleiner (2015) analytic solution
        if (trim(experiment) .eq. "kleiner-a" .and. trim(solver) .eq. "enth") then
            ab_warm_an = analytic_ab_steady(c%T0-10.0_wp,c,col%H_ice)
            ab_cold_an = analytic_ab_steady(c%T0-30.0_wp,c,col%H_ice)
            write(*,*) ""
            write(*,*) "=== Kleiner Exp A steady melt vs analytic (T3) ==="
            write(*,"(a,f8.3,a,f8.3,a)") "  warm a_b = ", ab_warm, " mm/a  (analytic ", ab_warm_an, ")"
            write(*,"(a,f8.3,a,f8.3,a)") "  cold a_b = ", ab_cold, " mm/a  (analytic ", ab_cold_an, ")"
            if (abs(ab_warm-ab_warm_an) .lt. 0.05_wp*abs(ab_warm_an) .and. &
                abs(ab_cold-ab_cold_an) .lt. 0.05_wp*abs(ab_cold_an)) then
                write(*,*) "  PASS: within 5% of analytic (increase nz to tighten)."
            else
                write(*,*) "  FAIL: exceeds 5% of analytic."
            end if
            write(*,*) ""
        end if

        return
    end subroutine run_experiment

    function analytic_ab_steady(T_srf,c,H_ice) result(ab)
        ! Kleiner (2015) Exp A analytic steady basal melt rate (Eq. A14) for a
        ! cold slab with the base held at the pressure melting point. Returns
        ! mm/a, positive = melt.
        real(wp), intent(IN) :: T_srf
        type(ybound_const_class), intent(IN) :: c
        real(wp), intent(IN) :: H_ice
        real(wp) :: ab
        real(wp) :: T_pmp_b, q_i, kt_si
        real(wp), parameter :: q_geo = 0.042_wp   ! [W m-2]
        kt_si   = 2.1_wp                            ! [W m-1 K-1]
        T_pmp_b = c%T0 - c%T_pmp_beta*c%rho_ice*c%g*H_ice
        q_i     = kt_si*(T_pmp_b-T_srf)/H_ice       ! cold-ice conductive flux up [W m-2]
        ab      = (q_geo - q_i)/(c%rho_ice*c%L_ice)*c%sec_year*1000.0_wp
        return
    end function analytic_ab_steady

    subroutine compare_cold_limit(c,nz)
        ! Re-run both solvers and compare final temperature profiles.
        type(ybound_const_class), intent(IN) :: c
        integer, intent(IN) :: nz
        type(column_class) :: cold_t, cold_e
        real(wp) :: dt, time, time_end, omega_max, Q_lith
        real(wp) :: dmax
        integer :: k

        omega_max = 0.0_wp
        time_end  = 50000.0_wp;  dt = 5.0_wp

        call setup_experiment(cold_t,c,"cold-limit",nz)
        call setup_experiment(cold_e,c,"cold-limit",nz)

        time = 0.0_wp
        do while (time .lt. time_end - 1e-6_wp)
            call calc_temp_column(cold_t%enth,cold_t%T_ice,cold_t%omega,cold_t%bmb,cold_t%Q_ice_b, &
                    cold_t%H_cts,cold_t%T_pmp,cold_t%cp,cold_t%kt,cold_t%advecxy,cold_t%uz,cold_t%Q_strn, &
                    cold_t%Q_b,cold_t%Q_rock,cold_t%T_srf,cold_t%T_shlf,cold_t%H_ice,cold_t%W_til,cold_t%f_grnd, &
                    cold_t%zeta_aa,cold_t%zeta_ac,cold_t%dzeta_a,cold_t%dzeta_b,omega_max,c%T0, &
                    c%rho_ice,c%rho_w,c%L_ice,c%sec_year,dt)

            Q_lith = cold_e%Q_rock*1e-3_wp*c%sec_year
            call calc_enth_column(cold_e%enth,cold_e%T_ice,cold_e%omega,cold_e%bmb,cold_e%Q_ice_b, &
                    cold_e%H_cts,cold_e%T_pmp,cold_e%cp,cold_e%kt,cold_e%advecxy,cold_e%uz,cold_e%Q_strn, &
                    cold_e%Q_b,Q_lith,cold_e%T_srf,cold_e%T_shlf,cold_e%H_ice,cold_e%W_til,cold_e%f_grnd, &
                    cold_e%zeta_aa,cold_e%zeta_ac,cold_e%dzeta_a,cold_e%dzeta_b,1.0_wp,omega_max,c%T0, &
                    c%rho_ice,c%rho_w,c%L_ice,c%sec_year,dt)
            time = time + dt
        end do

        dmax = 0.0_wp
        do k = 1, nz
            dmax = max(dmax, abs(cold_t%T_ice(k)-cold_e%T_ice(k)))
        end do

        write(*,*) ""
        write(*,*) "=== cold-limit equivalence (T2) ==="
        write(*,"(a,es12.4,a)") "  max|T_temp - T_enth| = ", dmax, " K"
        if (dmax .lt. 1e-2_wp) then
            write(*,*) "  PASS: enth reduces to temp in the cold limit."
        else
            write(*,*) "  FAIL: solvers disagree in the cold limit."
        end if
        write(*,*) ""
        return
    end subroutine compare_cold_limit

    subroutine write_init(col,filename)
        type(column_class), intent(IN) :: col
        character(len=*),   intent(IN) :: filename
        call nc_create(filename)
        call nc_write_dim(filename,"zeta",   x=col%zeta_aa, units="1")
        call nc_write_dim(filename,"zeta_ac",x=col%zeta_ac, units="1")
        call nc_write_dim(filename,"time",   x=0.0_wp,dx=1.0_wp,nx=1,units="years",unlimited=.TRUE.)
        return
    end subroutine write_init

    subroutine write_step(col,filename,time)
        type(column_class), intent(IN) :: col
        character(len=*),   intent(IN) :: filename
        real(wp),           intent(IN) :: time
        integer    :: ncid, n
        real(wp)   :: time_prev
        character(len=12), parameter :: vd = "zeta"

        call nc_open(filename,ncid,writable=.TRUE.)
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid)
        if (abs(time-time_prev) .gt. 1e-5_wp .or. n .eq. 0) n = n+1

        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"enth",   col%enth,  units="J kg-1",long_name="Ice enthalpy",           dim1=vd,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_ice",  col%T_ice, units="K",     long_name="Ice temperature",        dim1=vd,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"omega",  col%omega, units="1",     long_name="Water content fraction", dim1=vd,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_pmp",  col%T_pmp, units="K",     long_name="Pressure melting point", dim1=vd,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_prime",col%T_ice-col%T_pmp,units="K",long_name="Homologous temperature",dim1=vd,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"bmb",    col%bmb,   units="m a-1", long_name="Basal mass balance",     dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"H_cts",  col%H_cts, units="m",     long_name="CTS height",             dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"W_til",  col%W_til, units="m",     long_name="Basal till water",       dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"T_srf",  col%T_srf, units="K",     long_name="Surface temperature",    dim1="time",start=[n],ncid=ncid)
        call nc_close(ncid)
        return
    end subroutine write_step

end program test_enthalpy
