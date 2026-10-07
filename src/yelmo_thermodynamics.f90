
module yelmo_thermodynamics

    use nml 
    use yelmo_defs 
    use fast_hydrology, only : TRANSPORT_NONE
    use yelmo_grid, only : calc_zeta
    use yelmo_tools, only : smooth_gauss_2D, smooth_gauss_3D, gauss_values, fill_borders_2D, fill_borders_3D, &
            boundary_code, get_neighbor_indices_bc_codes, get_periodic_directions
    
    use thermodynamics 
    use ice_enthalpy
    use solver_advection, only : calc_advec2D

    implicit none
    
    private
    public :: calc_ytherm 
    public :: ytherm_par_load, ytherm_alloc, ytherm_dealloc 

contains

    subroutine calc_ytherm(thrm,tpo,dyn,mat,bnd,hyd,time)

        implicit none

        type(ytherm_class), intent(INOUT) :: thrm
        type(ytopo_class),  intent(IN)    :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ybound_class), intent(IN)    :: bnd
        type(hydro_class),  intent(IN)    :: hyd
        real(wp),         intent(IN)    :: time

        ! Local variables
        integer :: i, j, k, nx, ny
        real(wp) :: dt
        real(wp), allocatable :: dTdz_b_now(:,:)
        real(wp), allocatable :: C_cap(:,:)      ! [m/a ice equiv.] freeze-on capacity for the basal BC
        real(wp), allocatable :: Q_wat(:,:)      ! [mW m-2] water-side basal heat (Q_diss + Q_sens)
        character(len=56)     :: cap_source      ! thrm%par%cap_source with "auto" resolved

        logical, parameter :: calculate_Q_strn_derivative = .FALSE.

        nx = thrm%par%nx
        ny = thrm%par%ny

        ! Initialize time if necessary 
        if (thrm%par%time .gt. dble(time)) then 
            thrm%par%time = dble(time)
        end if 

        ! Get time step and advance current time 
        dt            = dble(time) - thrm%par%time 
        thrm%par%time = dble(time) 
        

        ! === Determine some thermal properties === 

        ! Calculate the specific heat capacity of the ice
        if (thrm%par%use_const_cp) then 
            thrm%now%cp  = thrm%par%const_cp
        else  
            thrm%now%cp  = calc_specific_heat_capacity(thrm%now%T_ice)
        end if 
        
        ! Calculate the heat conductivity of the ice
        if (thrm%par%use_const_kt) then 
            thrm%now%kt  = thrm%par%const_kt
        else  
            thrm%now%kt  = calc_thermal_conductivity(thrm%now%T_ice,bnd%c%sec_year)
        end if 

        ! The thermodynamics uses the column of the dynamics: thickness
        ! tpo%now%H_ice_dyn (H_eff in partial front cells and on the H_eff floor)

        ! Calculate the pressure-corrected melting point (in Kelvin)
        do k = 1, thrm%par%nz_aa  
            thrm%now%T_pmp(:,:,k) = calc_T_pmp(tpo%now%H_ice_dyn,thrm%par%z%zeta_aa(k), &
                                        bnd%c%T0,bnd%c%T_pmp_beta,bnd%c%rho_ice,bnd%c%g)
        end do 

        ! === Calculate heat source terms (Yelmo vertical grid) === 

        select case(thrm%par%qb_method)
            case(1)   ! "faces" == taub*u formed on the acx/acy faces, averaged to aa-nodes (energy-consistent)
                ! Calculate the basal frictional heating (from face products)
                call calc_basal_heating_faces(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy,tpo%now%f_ice, &
                                beta1=thrm%par%dt_beta(1),beta2=thrm%par%dt_beta(2),sec_year=bnd%c%sec_year,boundaries=thrm%par%boundaries)
            case(2)   ! "faces-nodes" == taub*u formed on the acx/acy faces, to quadrature nodes, averaged to aa-nodes (energy-consistent), default
                ! Calculate the basal frictional heating (from face products at quadrature-nodes)
                call calc_basal_heating_faces_nodes(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy,tpo%now%f_ice, &
                                beta1=thrm%par%dt_beta(1),beta2=thrm%par%dt_beta(2),sec_year=bnd%c%sec_year,boundaries=thrm%par%boundaries)
            case(3)     ! "aa" == simple stagger to aa-nodes directly
                ! Calculate the basal frictional heating (from aa-nodes)
                call calc_basal_heating_simplestagger(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy, &
                                                    beta1=thrm%par%dt_beta(1),beta2=thrm%par%dt_beta(2),sec_year=bnd%c%sec_year, &
                                                    boundaries=thrm%par%boundaries)
            case(4)   ! "nodes" == Gaussian quadrature to aa-node
                ! Calculate the basal frictional heating (from quadrature-nodes)
                call calc_basal_heating_nodes(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy,tpo%now%f_ice, &
                                beta1=thrm%par%dt_beta(1),beta2=thrm%par%dt_beta(2),sec_year=bnd%c%sec_year,boundaries=thrm%par%boundaries)
            case DEFAULT

        end select
        
        ! Calculate internal strain heating

        select case(trim(thrm%par%strain_heating))

            case("full")
                ! Calculate strain heating from strain rate tensor and viscosity (general approach)

                call calc_strain_heating(thrm%now%Q_strn,mat%now%strn%de,mat%now%visc,thrm%now%cp,bnd%c%rho_ice, &
                                                                        thrm%par%dt_beta(1),thrm%par%dt_beta(2))

            case("sia")
                ! Calculate strain heating from SIA approximation

                call calc_strain_heating_sia(thrm%now%Q_strn,dyn%now%ux,dyn%now%uy,tpo%now%dzsdx,tpo%now%dzsdy, &
                                      thrm%now%cp,tpo%now%H_ice_dyn,bnd%c%rho_ice,bnd%c%g,thrm%par%z%zeta_aa,thrm%par%z%zeta_ac, &
                                      thrm%par%dt_beta(1),thrm%par%dt_beta(2))

            case("none")
                ! No internal strain heating

                thrm%now%Q_strn = 0.0_wp

        end select
        
        ! Diagnose rate of change of strain heating w.r.t. temperature (dQsdT)
        if (calculate_Q_strn_derivative) then
            call calc_strain_heating_temp_derivative(thrm%now%dQsdT,thrm%now%Q_strn,thrm%now%T_ice,thrm%now%cp,tpo%now%H_ice_dyn,tpo%now%f_ice, &
                                                        thrm%par%z%zeta_aa,bnd%c%rho_ice,thrm%par%dx,thrm%par%dy,thrm%par%boundaries)
        else
            thrm%now%dQsdT = 0.0
        end if

        ! Ensure that Q_rock is defined. At initialization, 
        ! it may have a value of zero. In this case, set equal 
        ! to Q_geo to be consistent with equilibrium bedrock conditions. 
        if (maxval(thrm%now%Q_rock) .eq. 0.0) then 
            thrm%now%Q_rock = bnd%Q_geo 
        end if

        if ( dt .gt. 0.0 ) then
            ! Ice thermodynamics should evolve, perform calculations.
            ! Basal water (W_til) is owned by hyd and updated separately
            ! after thermodynamics. For the basal BC decision inside the
            ! enthalpy solver, a thin one-line forward-Euler predictor
            ! W_til_predicted = W_til - bmb*dt is computed locally inside
            ! calc_enth_column from the current hyd%now%W_til.

            select case(trim(thrm%par%method))

                case("enth","temp") 
                    ! Perform enthalpy/temperature solving via advection-diffusion equation
                    ! Note: method==temp solves the columns with calc_temp_column (no water
                    ! content, centred vertical advection, "wtil" basal BC), with
                    ! enth_cr=1.0 and omega_max=0.0 prescribed in par_load(). 

                    if (trim(thrm%par%method) .eq. "enth") then 

                        ! Calculate the explicit horizontal advection term using enthalpy from previous timestep
                        call calc_advec_horizontal_3D(thrm%now%advecxy,thrm%now%enth,tpo%now%H_ice, &
                                            dyn%now%ux,dyn%now%uy,thrm%par%dx,dt,thrm%par%advecxy_order, &
                                            thrm%par%advecxy_cfl,thrm%par%advecxy_nmax,thrm%par%boundaries)

                    else

                        ! Calculate the explicit horizontal advection term using temperature from previous timestep
                        call calc_advec_horizontal_3D(thrm%now%advecxy,thrm%now%T_ice,tpo%now%H_ice, &
                                            dyn%now%ux,dyn%now%uy,thrm%par%dx,dt,thrm%par%advecxy_order, &
                                            thrm%par%advecxy_cfl,thrm%par%advecxy_nmax,thrm%par%boundaries)
                    
                    end if 

                    ! Freeze-on capacity (used only when basal_bc_method="capacity") and
                    ! water-side basal heat (used under either basal BC rule).
                    ! hyd stores them in SI: C_frz [m/s ice equiv.], Q_diss/Q_sens [W m-2].
                    allocate(C_cap(nx,ny), Q_wat(nx,ny))
                    cap_source = thrm%par%cap_source
                    if (trim(cap_source) .eq. "auto") then
                        ! The transport model's own C if there is one, else the bucket's stock.
                        if (hyd%par%method_transport .ne. TRANSPORT_NONE) then
                            cap_source = "hyd"
                        else
                            cap_source = "till"
                        end if
                    end if
                    select case(trim(cap_source))
                        case("hyd")
                            C_cap = hyd%now%C_frz * bnd%c%sec_year
                        case("till")
                            ! Stock estimate from the bucket: all till water above the
                            ! floor refrozen over this step, converted to ice equivalent.
                            C_cap = (bnd%c%rho_w/bnd%c%rho_ice) * max(hyd%now%W_til - thrm%par%cap_W_floor, 0.0_wp) / dt
                        case("water")
                            ! Stock estimate from the water thickness: all water above
                            ! the floor refrozen over this step, converted to ice equivalent.
                            C_cap = (bnd%c%rho_w/bnd%c%rho_ice) * max(hyd%now%W - thrm%par%cap_W_floor, 0.0_wp) / dt
                        case DEFAULT    ! "none"
                            C_cap = 0.0_wp
                    end select
                    Q_wat = (hyd%now%Q_diss + hyd%now%Q_sens) * 1e3_wp

                    ! Now calculate the thermodynamics:

                    call calc_ytherm_enthalpy_3D(thrm%now%enth,thrm%now%T_ice,thrm%now%omega,thrm%now%bmb_grnd, &
                                thrm%now%Q_ice_b,thrm%now%H_cts,thrm%now%T_pmp,thrm%now%cp,thrm%now%kt,thrm%now%advecxy, &
                                dyn%now%ux,dyn%now%uy,dyn%now%uz_star,thrm%now%Q_strn,thrm%now%Q_b,thrm%now%Q_rock,bnd%T_srf, &
                                tpo%now%H_ice_dyn,tpo%now%f_ice,tpo%now%z_srf,hyd%now%W_til,tpo%now%H_grnd, &
                                tpo%now%f_grnd,thrm%par%z%zeta_aa,thrm%par%z%zeta_ac,thrm%par%z%dzeta_a,thrm%par%z%dzeta_b, &
                                thrm%par%enth_cr,thrm%par%omega_max,thrm%par%H_ice_thin,bnd%c%rho_ice,bnd%c%rho_sw,bnd%c%rho_w,bnd%c%L_ice,bnd%c%T0, &
                                bnd%c%sec_year,dt,thrm%par%method,thrm%par%solver_advec,thrm%par%enth_integral, &
                                thrm%par%boundaries,C_cap,Q_wat,thrm%par%basal_bc_method,thrm%par%cap_eps, &
                                thrm%now%bmb_grnd_star,thrm%now%bc_b,thrm%now%bmb_clamp,thrm%now%melt_int, &
                                thrm%par%gl_temperate)

                    deallocate(C_cap, Q_wat)

                case("robin")
                    ! Use Robin solution for ice temperature

                    call define_temp_robin_3D(thrm%now%enth,thrm%now%T_ice,thrm%now%omega,thrm%now%T_pmp,thrm%par%const_cp,thrm%par%const_kt, &
                                       thrm%now%Q_rock,bnd%T_srf,tpo%now%H_ice_dyn,hyd%now%W_til,bnd%smb, &
                                       thrm%now%bmb_grnd,tpo%now%f_grnd,thrm%par%z%zeta_aa, &
                                       bnd%c%rho_ice,bnd%c%L_ice,bnd%c%sec_year,cold=.FALSE.,enth_integral=thrm%par%enth_integral)

                case("robin-cold")
                    ! Use Robin solution for ice temperature averaged with cold linear profile
                    ! to ensure cold ice at the base

                    call define_temp_robin_3D(thrm%now%enth,thrm%now%T_ice,thrm%now%omega,thrm%now%T_pmp,thrm%par%const_cp,thrm%par%const_kt, &
                                       thrm%now%Q_rock,bnd%T_srf,tpo%now%H_ice_dyn,hyd%now%W_til,bnd%smb, &
                                       thrm%now%bmb_grnd,tpo%now%f_grnd,thrm%par%z%zeta_aa, &
                                       bnd%c%rho_ice,bnd%c%L_ice,bnd%c%sec_year,cold=.TRUE.,enth_integral=thrm%par%enth_integral)

                case("linear")
                    ! Use linear solution for ice temperature

                    ! Calculate the ice temperature (eventually water content and enthalpy too)
                    call define_temp_linear_3D(thrm%now%enth,thrm%now%T_ice,thrm%now%omega,thrm%now%cp,tpo%now%H_ice_dyn,bnd%T_srf,thrm%par%z%zeta_aa, &
                                        bnd%c%T0,bnd%c%rho_ice,bnd%c%L_ice,bnd%c%T_pmp_beta,bnd%c%g,enth_integral=thrm%par%enth_integral)

                case("fixed") 
                    ! Pass - do nothing, use the enth/temp/omega fields as they are defined

                case("prescribed")
                    ! T_ice has been set externally (yelmo_init_state only):
                    ! cap at the melting point, omega = 0, consistent enthalpy

                    call define_temp_prescribed_3D(thrm%now%enth,thrm%now%T_ice,thrm%now%omega,thrm%now%T_pmp, &
                                                   bnd%c%L_ice,enth_integral=thrm%par%enth_integral)

                case DEFAULT 

                    write(*,*) "ytherm:: Error: thermodynamics option not recognized: method = ", trim(thrm%par%method)
                    error stop 1

            end select 

            ! (No basal-water bucket update here - hyd owns W_til and is
            ! advanced separately after calc_ytherm by calc_yhyd.)


            ! ==== Bedrock ======================================

            ! Update the bedrock temperature profile 
            ! (using basal ice temperature from previous timestep)
            select case(trim(thrm%par%rock_method))

                case("equil")
                    ! Prescribe bedrock temperature profile assuming 
                    ! equilibrium with the bed surface temperature 
                    ! (ie, no active bedrock) 

                    call define_temp_bedrock_3D(thrm%now%T_rock,thrm%now%Q_rock, &
                                             thrm%par%kt_rock,bnd%Q_geo,thrm%now%T_ice(:,:,1), &
                                             thrm%par%H_rock,thrm%par%zr%zeta_aa,bnd%c%sec_year)

                case("active")
                    ! Solve thermodynamic equation for the bedrock 

                    call calc_ytherm_temp_bedrock_3D(thrm%now%T_rock,thrm%now%Q_rock, &
                                    thrm%now%T_ice(:,:,1),thrm%now%T_pmp(:,:,1),thrm%par%rhoc_rock,thrm%par%kt_rock, &
                                    thrm%par%H_rock,tpo%now%H_ice_dyn,tpo%now%H_grnd,thrm%now%Q_ice_b,bnd%Q_geo, &
                                    thrm%par%zr%zeta_aa,thrm%par%zr%zeta_ac,thrm%par%zr%dzeta_a,thrm%par%zr%dzeta_b, &
                                    bnd%c%rho_ice,bnd%c%rho_sw,bnd%c%T0,bnd%c%sec_year,dt)

                case("fixed") 
                    ! Pass - do nothing, use the enth/temp/omega fields as they are defined

                case DEFAULT 

                    write(*,*) "calc_ytherm:: Error: rock_method not recognized."
                    write(*,*) "rock_method = ", trim(thrm%par%rock_method)

            end select 

            ! =======================================================

        end if 

        ! Calculate homologous temperature everywhere and at the base
        ! (T_prime = T_ice - T_pmp; 0 == temperate, negative below the melting point)
        thrm%now%T_prime   = thrm%now%T_ice - thrm%now%T_pmp
        thrm%now%T_prime_b = thrm%now%T_prime(:,:,1)
        
        ! Calculate gridpoint fraction at the pressure melting point
        call calc_f_pmp(thrm%now%f_pmp,thrm%now%T_ice(:,:,1),thrm%now%T_pmp(:,:,1), &
                                                        tpo%now%f_grnd,thrm%par%gamma)

!         if (yelmo_log) then 
!             if (count(tpo%now%H_ice.gt.0.0) .gt. 0) then 
!                 write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytherm:: time = ", thrm%par%time, dt, &
!                     sum(thrm%now%T_ice(:,:,thrm%par%nz_aa),mask=tpo%now%H_ice.gt.0.0)/real(count(tpo%now%H_ice.gt.0.0))
!             else 
!                 write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytherm:: time = ", thrm%par%time, dt, 0.0 
!             end if 
!         end if 

        return

    end subroutine calc_ytherm

    subroutine calc_ytherm_enthalpy_3D(enth,T_ice,omega,bmb_grnd,Q_ice_b,H_cts,T_pmp,cp,kt,advecxy,ux,uy,uz,Q_strn,Q_b,Q_rock, &
                                        T_srf,H_ice_dyn,f_ice,z_srf,W_til,H_grnd,f_grnd,zeta_aa,zeta_ac,dzeta_a,dzeta_b, &
                                        cr,omega_max,H_ice_thin,rho_ice,rho_sw,rho_w,L_ice,T0,sec_year,dt,solver,solver_advec,enth_integral, &
                                        boundaries,C_cap,Q_wat,basal_bc_method,cap_eps,bmb_grnd_star,bc_b,bmb_clamp,melt_int, &
                                        gl_temperate)
        ! This wrapper subroutine breaks the thermodynamics problem into individual columns,
        ! which are solved independently by calling calc_enth_column.
        ! The column is that of the dynamics (thickness H_ice_dyn, paired with
        ! uz_star). It is solved in fully ice-covered cells (f_ice == 1) thicker
        ! than H_ice_thin; partial front cells get a linear profile and are then
        ! filled from their fully ice-covered neighbours.

        ! Note zeta=height, k=1 base, k=nz surface 
        
        !$ use omp_lib

        implicit none 

        real(wp), intent(INOUT) :: enth(:,:,:)    ! [J kg-1] Ice enthalpy
        real(wp), intent(INOUT) :: T_ice(:,:,:)   ! [K] Ice column temperature
        real(wp), intent(INOUT) :: omega(:,:,:)   ! [--] Ice water content
        real(wp), intent(INOUT) :: bmb_grnd(:,:)  ! [m a-1] Basal mass balance (melting is negative)
        real(wp), intent(OUT)   :: Q_ice_b(:,:)   ! [mW m-2] Basal ice heat flux 
        real(wp), intent(OUT)   :: H_cts(:,:)     ! [m] Height of the cold-temperate transition surface (CTS)
        real(wp), intent(INOUT) :: T_pmp(:,:,:)   ! [K] Pressure melting point temp.
        real(wp), intent(IN)    :: cp(:,:,:)      ! [J kg-1 K-1] Specific heat capacity
        real(wp), intent(IN)    :: kt(:,:,:)      ! [J a-1 m-1 K-1] Heat conductivity 
        real(wp), intent(IN)    :: advecxy(:,:,:) ! [J kg-1 a-1] (enth) or [K a-1] (temp) Horizontal advection 
        real(wp), intent(IN)    :: ux(:,:,:)      ! [m a-1] Horizontal x-velocity 
        real(wp), intent(IN)    :: uy(:,:,:)      ! [m a-1] Horizontal y-velocity 
        real(wp), intent(IN)    :: uz(:,:,:)      ! [m a-1] Vertical velocity 
        real(wp), intent(IN)    :: Q_strn(:,:,:)  ! [J a-1 m-3] Internal strain heat production in ice
        real(wp), intent(IN)    :: Q_b(:,:)       ! [J a-1 m-2] Basal frictional heat production 
        real(wp), intent(IN)    :: Q_rock(:,:)    ! [mW m-2] Heat flux at bed surface from bedrock (like Q_geo)
        real(wp), intent(IN)    :: T_srf(:,:)     ! [K] Surface temperature 
        real(wp), intent(IN)    :: H_ice_dyn(:,:) ! [m] Active column thickness (tpo%now%H_ice_dyn)
        real(wp), intent(IN)    :: f_ice(:,:)     ! [--] Area fraction ice cover
        real(wp), intent(IN)    :: z_srf(:,:)     ! [m] Surface elevation 
        real(wp), intent(IN)    :: W_til(:,:)     ! [m] Basal till water thickness (from hyd)
        real(wp), intent(IN)    :: H_grnd(:,:)    ! [--] Ice thickness above flotation
        real(wp), intent(IN)    :: f_grnd(:,:)    ! [--] Grounded fraction
        real(wp), intent(IN)    :: zeta_aa(:)     ! [--] Vertical sigma coordinates (zeta==height), aa-nodes
        real(wp), intent(IN)    :: zeta_ac(:)     ! [--] Vertical sigma coordinates (zeta==height), ac-nodes
        real(wp), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(wp), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(wp), intent(IN)    :: cr             ! [--] Conductivity ratio for temperate ice (kappa_temp = enth_cr*kappa_cold)
        real(wp), intent(IN)    :: omega_max      ! [--] Maximum allowed water content fraction
        real(wp), intent(IN)    :: H_ice_thin     ! [m] Thickness threshold below which the column solver is skipped
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw
        real(wp), intent(IN)    :: rho_w
        real(wp), intent(IN)    :: L_ice
        real(wp), intent(IN)    :: T0
        real(wp), intent(IN)    :: sec_year 
        real(wp), intent(IN)    :: dt             ! [a] Time step 
        character(len=*), intent(IN) :: solver      ! "enth" or "temp"
        character(len=*), intent(IN) :: solver_advec    ! "expl" or "impl-upwind"
        logical,          intent(IN) :: enth_integral   ! use integral (A2) enthalpy definition?
        character(len=*), intent(IN) :: boundaries      ! Boundary treatment
        real(wp),         intent(IN) :: C_cap(:,:)      ! [m/a ice equiv.] Freeze-on capacity (basal_bc_method="capacity")
        real(wp),         intent(IN) :: Q_wat(:,:)      ! [mW m-2] Water-side basal heat, Q_diss + Q_sens
        character(len=*), intent(IN) :: basal_bc_method ! "wtil" or "capacity"
        real(wp),         intent(IN) :: cap_eps         ! [m/a ice equiv.] Capacity below which the bed counts as dry
        real(wp),         intent(OUT) :: bmb_grnd_star(:,:) ! [m/a] bmb of a base held at T_pmp (capacity rule)
        real(wp),         intent(OUT) :: bc_b(:,:)          ! [--] basal BC used: 0 not grounded/solved, 1 held at T_pmp, 2 flux
        logical,          intent(IN)  :: gl_temperate       ! Hold grounded bases next to the ocean (f_grnd=0 neighbour) at T_pmp
        real(wp),         intent(OUT) :: bmb_clamp(:,:)     ! [m/a] freeze-on removed by the capacity safety clamp
        real(wp),         intent(OUT) :: melt_int(:,:)      ! [m/a ice equiv.] englacial water drained to the bed

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa, nz_ac  
        integer :: i1, i2, j1, j2 
        integer :: im1, ip1, jm1, jp1 
        integer :: ii(3), jj(3) 
        integer :: BC 
        logical :: per_x, per_y 
        real(wp) :: T_shlf, H_grnd_lim, f_scalar, T_base  
        real(wp) :: H_ice_now 
        real(wp) :: wt_neighb(3,3) 
        real(wp) :: wt_tot
        logical  :: gl_temp 

        ! ajr symtest
        logical :: is_symmetric 

        nx    = size(T_ice,1)
        ny    = size(T_ice,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        ! Solve all points in periodic directions (true wrap: every point is
        ! an interior point), otherwise only the interior; non-periodic
        ! borders are filled from their interior neighbors below.
        BC = boundary_code(boundaries)
        call get_periodic_directions(per_x,per_y,BC)

        i1 = 2
        i2 = nx-1
        if (per_x) then
            i1 = 1
            i2 = nx
        end if

        j1 = 2
        j2 = ny-1
        if (per_y) then
            j1 = 1
            j2 = ny
        end if

        ! ===================================================

        !$omp parallel do collapse(2) schedule(dynamic,64) private(i,j,im1,ip1,jm1,jp1,H_ice_now,T_shlf,T_base,gl_temp)
        do j = j1, j2
        do i = i1, i2 
            
            H_ice_now = H_ice_dyn(i,j)

            ! Fully grounded cell next to floating ice or open ocean: the bed is
            ! wetted by the ocean, so the base is held temperate (ytherm.gl_temperate)
            gl_temp = .FALSE.
            if (gl_temperate .and. f_grnd(i,j) .eq. 1.0_wp) then
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
                gl_temp = f_grnd(im1,j) .eq. 0.0_wp .or. f_grnd(ip1,j) .eq. 0.0_wp .or. &
                          f_grnd(i,jm1) .eq. 0.0_wp .or. f_grnd(i,jp1) .eq. 0.0_wp
            end if

            ! For floating points, calculate the approximate marine-shelf temperature 
            ! ajr, later this should come from an external model, and T_shlf would
            ! be the boundary variable directly
            if (f_grnd(i,j) .lt. 1.0) then 

                ! Calculate approximate marine freezing temp, limited to pressure melting point 
                T_shlf = calc_T_base_shlf_approx(H_ice_now,T_pmp(i,j,1),H_grnd(i,j),T0,rho_ice,rho_sw)

            else 
                ! Assigned for safety 

                T_shlf   = T_pmp(i,j,1)

            end if 

            if (f_ice(i,j) .eq. 1.0 .and. H_ice_now .gt. H_ice_thin) then 
                ! Thick ice exists, call thermodynamic solver for the column

                if (trim(solver) .eq. "enth") then 

                    call calc_enth_column(enth(i,j,:),T_ice(i,j,:),omega(i,j,:),bmb_grnd(i,j),Q_ice_b(i,j), &
                            H_cts(i,j),T_pmp(i,j,:),cp(i,j,:),kt(i,j,:),advecxy(i,j,:),uz(i,j,:),Q_strn(i,j,:), &
                            Q_b(i,j),Q_rock(i,j),T_srf(i,j),T_shlf,H_ice_now,W_til(i,j),f_grnd(i,j),zeta_aa, &
                            zeta_ac,dzeta_a,dzeta_b,cr,omega_max,T0,rho_ice,rho_w,L_ice,sec_year,dt,enth_integral, &
                            basal_bc_method,C_cap(i,j),Q_wat(i,j),cap_eps, &
                            bmb_grnd_star(i,j),bc_b(i,j),bmb_clamp(i,j),melt_int_out=melt_int(i,j),gl_temp=gl_temp)

                else

                    call calc_temp_column(enth(i,j,:),T_ice(i,j,:),omega(i,j,:),bmb_grnd(i,j),Q_ice_b(i,j), &
                            H_cts(i,j),T_pmp(i,j,:),cp(i,j,:),kt(i,j,:),advecxy(i,j,:),uz(i,j,:),Q_strn(i,j,:), &
                            Q_b(i,j),Q_rock(i,j),T_srf(i,j),T_shlf,H_ice_now,W_til(i,j),f_grnd(i,j),zeta_aa, &
                            zeta_ac,dzeta_a,dzeta_b,omega_max,T0,rho_ice,rho_w,L_ice,sec_year,dt,enth_integral, &
                            gl_temp=gl_temp)
                    bmb_grnd_star(i,j) = 0.0_wp
                    bc_b(i,j)          = 0.0_wp
                    bmb_clamp(i,j)     = 0.0_wp
                    melt_int(i,j)      = 0.0_wp

                end if

            else 
                ! Ice is at margin, too thin or zero: prescribe linear temperature profile
                ! between temperate ice at base and surface temperature 
                ! (accounting for floating/grounded nature via T_base)

                if (f_grnd(i,j) .lt. 1.0) then 
                    ! Impose T_shlf for the basal temperature
                    T_base = T_shlf 
                else
                    ! Impose temperature at the pressure melting point of grounded ice 
                    T_base = T_pmp(i,j,1) 
                end if 

                T_ice(i,j,:)  = define_temp_linear_column(T_srf(i,j),T_base,T_pmp(i,j,nz_aa),zeta_aa)
                omega(i,j,:)  = 0.0_wp
                call convert_to_enthalpy_ice(enth(i,j,:),T_ice(i,j,:),omega(i,j,:),T_pmp(i,j,:),L_ice,enth_integral)
                bmb_grnd(i,j) = 0.0_wp
                Q_ice_b(i,j)  = 0.0_wp
                H_cts(i,j)    = 0.0_wp
                bmb_grnd_star(i,j) = 0.0_wp
                bc_b(i,j)          = 0.0_wp
                bmb_clamp(i,j)     = 0.0_wp
                melt_int(i,j)      = 0.0_wp

            end if 

        end do 
        end do 
        !$omp end parallel do

! ajr symtest: check BCs for symmetry
if (.FALSE.) then

        ! diva (moving)
        ! i = 25
        ! j = 18 

        ! sia (moving)
        i = 20
        j = 25 

        write(*,*)
        call check_symmetry_2D(T_ice(:,:,1),"T_ice_b",i,j,"x",is_symmetric)
        
        if (.not. is_symmetric) then
            call check_symmetry_2D(H_ice_dyn,"H_ice_dyn",i,j,"x")
            call check_symmetry_2D(f_ice,"f_ice",i,j,"x")
            call check_symmetry_2D(bmb_grnd,"bmb_grnd",i,j,"x")
            call check_symmetry_2D(Q_strn(:,:,1),"Q_strn_b",i,j,"x")
            call check_symmetry_2D(Q_b,"Q_b",i,j,"x")
            call check_symmetry_2D(T_srf,"T_srf",i,j,"x")
            call check_symmetry_2D(uz(:,:,1),"uz_b",i,j,"x")
            call check_symmetry_2D(W_til,"W_til",i,j,"x")
            call check_symmetry_2D(Q_rock,"Q_rock",i,j,"x")
            call check_symmetry_2D(Q_ice_b,"Q_ice_b",i,j,"x")
            call check_symmetry_2D(advecxy(:,:,1),"advecxy_b",i,j,"x")

            stop "Symmetry!"
        end if

        write(*,*) 

end if

if (.TRUE.) then
        ! Extrapolate thermodynamics to ice-free and partially ice-covered 
        ! neighbors to the ice margin.
        ! (Helps with stability to give good values of ATT to newly advected points)
        !$omp parallel do collapse(2) private(i,j,k,im1,ip1,jm1,jp1,ii,jj,wt_neighb,wt_tot)
        do j = j1, j2
        do i = i1, i2 
            
            if (f_ice(i,j) .lt. 1.0) then 

                ! 3x3 neighborhood with BC-aware neighbor indices
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
                ii = [im1,i,ip1]
                jj = [jm1,j,jp1]

                wt_neighb = 0.0 
                where (f_ice(ii,jj) .eq. 1.0) wt_neighb = 1.0 
                wt_tot = sum(wt_neighb)

                if (wt_tot .gt. 0.0) then 
                    ! Ice covered neighbor(s) found, assign average of neighbors 

                    ! Normalize weights 
                    wt_neighb = wt_neighb / wt_tot 

                    do k = 1, nz_aa 
                        enth(i,j,k)  = sum(enth(ii,jj,k) *wt_neighb)
                        T_ice(i,j,k) = sum(T_ice(ii,jj,k)*wt_neighb)
                        omega(i,j,k) = sum(omega(ii,jj,k)*wt_neighb)
                        T_pmp(i,j,k) = sum(T_pmp(ii,jj,k)*wt_neighb)
                    end do 

                end if 

            end if 
    
        end do 
        end do 
        !$omp end parallel do
end if 

        ! Fill in non-periodic borders from interior neighbors
        call fill_borders_3D(enth,    nfill=1,fill_x=.not.per_x,fill_y=.not.per_y)
        call fill_borders_3D(T_ice,   nfill=1,fill_x=.not.per_x,fill_y=.not.per_y)
        call fill_borders_3D(omega,   nfill=1,fill_x=.not.per_x,fill_y=.not.per_y)
        call fill_borders_2D(bmb_grnd,nfill=1,fill_x=.not.per_x,fill_y=.not.per_y)
        call fill_borders_2D(Q_ice_b, nfill=1,fill_x=.not.per_x,fill_y=.not.per_y)
        call fill_borders_2D(H_cts,   nfill=1,fill_x=.not.per_x,fill_y=.not.per_y)
        call fill_borders_2D(bmb_grnd_star,nfill=1,fill_x=.not.per_x,fill_y=.not.per_y)
        call fill_borders_2D(bc_b,         nfill=1,fill_x=.not.per_x,fill_y=.not.per_y)
        call fill_borders_2D(bmb_clamp,    nfill=1,fill_x=.not.per_x,fill_y=.not.per_y)
        call fill_borders_2D(melt_int,     nfill=1,fill_x=.not.per_x,fill_y=.not.per_y)
        
        return 

    end subroutine calc_ytherm_enthalpy_3D

    subroutine check_symmetry_2D(var,varnm,i,j,dir,is_symmetric)

        implicit none

        real(wp), intent(IN) :: var(:,:)
        character(len=*), intent(IN) :: varnm
        integer,  intent(IN) :: i, j  
        character(len=*), intent(IN) :: dir     ! Direction to check "x" or "y"
        logical, optional, intent(OUT) :: is_symmetric 

        ! Local variables
        integer :: imid, jmid
        integer :: is, js 
        
        imid = (size(var,1)-1)/2 + 1 
        jmid = (size(var,2)-1)/2 + 1 
        
        ! Get symmetric counterparts
        if (dir .eq. "x") then 
            js = j 
            is = imid - (i-imid)
        else if (dir .eq. "y") then 
            is = i 
            js = jmid - (j-jmid)
        else
            write(error_unit,*) "check_symmetry_2D:: Error: argument 'dir' must be 'x' or 'y'."
            error stop 1
        end if

        write(*,"(a4,a12,2f15.3,g18.6)") "sym: ", trim(varnm), var(i,j), var(is,js), abs(var(is,js)-var(i,j))

        if (present(is_symmetric)) then
            if (abs(var(is,js)-var(i,j)) .lt. 1e-3) then
                is_symmetric = .TRUE.
            else 
                is_symmetric = .FALSE.
            end if
        end if

        return

    end subroutine check_symmetry_2D

    subroutine calc_ytherm_temp_bedrock_3D(T_rock,Q_rock,T_ice_b,T_pmp_b,rhoc_rock,kt_rock,H_rock, &
                                                H_ice,H_grnd,Q_ice_b,Q_geo,zeta_aa,zeta_ac,dzeta_a,dzeta_b, &
                                                rho_ice,rho_sw,T0,sec_year,dt)
        ! This wrapper subroutine breaks the thermodynamics problem into individual columns,
        ! which are solved independently by calling calc_temp_bedrock_column

        ! Note zeta=height, k=1 base, k=nz surface 
        
        !$ use omp_lib

        implicit none 

        real(wp), intent(INOUT) :: T_rock(:,:,:)      ! [K] Bedrock temperature
        real(wp), intent(OUT)   :: Q_rock(:,:)        ! [mW m-2] Bed surface heat flux 
        real(wp), intent(IN)    :: T_ice_b(:,:)       ! [K] Ice temperature at ice base
        real(wp), intent(IN)    :: T_pmp_b(:,:)       ! [K] Pressure melting point temp at ice base.
        real(wp), intent(IN)    :: rhoc_rock          ! [J m-3 K-1] Volumetric heat capacity
        real(wp), intent(IN)    :: kt_rock            ! [J a-1 m-1 K-1] Heat conductivity
        real(wp), intent(IN)    :: H_rock             ! [m] Bedrock thickness 
        real(wp), intent(IN)    :: H_ice(:,:)         ! [m] Ice column thickness (tpo%now%H_ice_dyn)
        real(wp), intent(IN)    :: H_grnd(:,:)        ! [--] Ice thickness above flotation 
        real(wp), intent(IN)    :: Q_ice_b(:,:)       ! [mW m-2] Ice base heat flux
        real(wp), intent(IN)    :: Q_geo(:,:)         ! [mW m-2] Geothermal heat flux deep in bedrock
        real(wp), intent(IN)    :: zeta_aa(:)         ! [--] Vertical sigma coordinates (zeta==height), aa-nodes
        real(wp), intent(IN)    :: zeta_ac(:)         ! [--] Vertical sigma coordinates (zeta==height), ac-nodes
        real(wp), intent(IN)    :: dzeta_a(:)         ! nz_aa [--] Solver discretization helper variable ak
        real(wp), intent(IN)    :: dzeta_b(:)         ! nz_aa [--] Solver discretization helper variable bk
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw
        real(wp), intent(IN)    :: T0 
        real(wp), intent(IN)    :: sec_year 
        real(wp), intent(IN)    :: dt                 ! [a] Time step 

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa, nz_ac  
        real(wp) :: T_base

        nx    = size(T_rock,1)
        ny    = size(T_rock,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        ! ===================================================

        ! ajr: openmp problematic here - leads to NaNs
        !$omp parallel do collapse(2) private(i,j,T_base)
        do j = 1, ny
        do i = 1, nx 

            ! For floating points, calculate the approximate marine-shelf temperature 
            ! although really the temperature at the bottom of the ocean is needed
            if (H_grnd(i,j) .lt. 0.0) then 

                ! Calculate approximate marine freezing temp, limited to pressure melting point 
                T_base = calc_T_base_shlf_approx(H_ice(i,j),T_pmp_b(i,j),H_grnd(i,j),T0,rho_ice,rho_sw)

            else 
                ! Assign ice basal temperature
                T_base   = T_ice_b(i,j)

            end if 

            if (H_ice(i,j) .gt. 0.0) then 
                ! Call thermodynamic solver for the column

                call calc_temp_bedrock_column(T_rock(i,j,:),Q_rock(i,j),  &
                        rhoc_rock,kt_rock,Q_ice_b(i,j),Q_geo(i,j),T_base,H_rock,zeta_aa, &
                        zeta_ac,dzeta_a,dzeta_b,sec_year,dt)
            
            else 
                ! Assume equilibrium conditions: impose linear temperature 
                ! profile following Q_geo and T_base

                call define_temp_bedrock_column(T_rock(i,j,:),kt_rock,H_rock, &
                                                                    T_base,Q_geo(i,j),zeta_aa,sec_year)

            end if 

        end do 
        end do 
        !$omp end parallel do

        return 

    end subroutine calc_ytherm_temp_bedrock_3D
    
    subroutine ytherm_par_load(par,filename,group,zeta_aa,zeta_ac,nx,ny,dx,init)

        type(ytherm_param_class), intent(OUT) :: par
        character(len=*),         intent(IN)  :: filename
        character(len=*),         intent(IN)  :: group          ! Usually "ytherm"
        real(wp),                 intent(IN)  :: zeta_aa(:)  
        real(wp),                 intent(IN)  :: zeta_ac(:)  
        integer,                  intent(IN)  :: nx, ny 
        real(wp),                 intent(IN)  :: dx 
        logical, optional,        intent(IN)  :: init

        ! Local variables
        logical :: init_pars
        integer :: k

        character(len=*), parameter :: def_file   = "input/yelmo_defaults.nml"
        character(len=*), parameter :: def_ytherm = "ytherm"

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE.

        call nml_validate(filename,def_file,group,defaults_group=def_ytherm)

        ! Store local parameter values in output object
        call nml_read(filename,group,"method",         par%method,           init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"qb_method",      par%qb_method,        init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"dt_method",      par%dt_method,        init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"solver_advec",   par%solver_advec,     init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"advecxy_order",   par%advecxy_order,    init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"advecxy_cfl",     par%advecxy_cfl,      init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"advecxy_nmax",    par%advecxy_nmax,     init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"gamma",          par%gamma,            init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"strain_heating", par%strain_heating,   init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"use_const_cp",   par%use_const_cp,     init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"const_cp",       par%const_cp,         init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"use_const_kt",   par%use_const_kt,     init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"const_kt",       par%const_kt,         init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"enth_cr",        par%enth_cr,          init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"omega_max",      par%omega_max,        init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"H_ice_thin",     par%H_ice_thin,       init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"enth_cp_method",  par%enth_cp_method,  init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        par%enth_integral = (trim(par%enth_cp_method) .eq. "integral")

        call nml_read(filename,group,"basal_bc_method",par%basal_bc_method,  init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"cap_source",     par%cap_source,       init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"cap_W_floor",    par%cap_W_floor,      init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"cap_eps",        par%cap_eps,          init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"gl_temperate",   par%gl_temperate,     init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)

        call yelmo_check_enum(group,"basal_bc_method", par%basal_bc_method, "capacity|wtil")
        call yelmo_check_enum(group,"cap_source",      par%cap_source,      "auto|hyd|till|water|none")

        if (trim(par%basal_bc_method) .eq. "wtil") then
            write(io_unit_err,*) "ytherm_par_load:: warning: basal_bc_method='wtil' is deprecated; use 'capacity'."
        end if

        if (trim(par%basal_bc_method) .eq. "capacity" .and. trim(par%method) .ne. "enth") then
            ! Only the enthalpy column has the capacity rule; the other
            ! solvers keep their own basal treatment.
            write(io_unit_err,*) "ytherm_par_load:: note: basal_bc_method='capacity' applies to method='enth' only; ", &
                                 "using 'wtil' with method=", trim(par%method)
            par%basal_bc_method = "wtil"
        end if

        if (par%cap_W_floor .lt. 0.0_wp .or. par%cap_eps .lt. 0.0_wp) then
            write(io_unit_err,*) "ytherm_par_load:: error: cap_W_floor and cap_eps must be >= 0; got ", par%cap_W_floor, par%cap_eps
            stop
        end if

        ! Note: till_rate and H_w_max moved to &fhyd (par%bucket%till_rate
        ! and par%W_til_max in fasthydrology). They are no longer read here.

        call nml_read(filename,group,"rock_method",    par%rock_method,      init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"nzr_aa",         par%nzr_aa,           init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"zeta_scale_rock",par%zeta_scale_rock,  init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"zeta_exp_rock",  par%zeta_exp_rock,    init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"H_rock",         par%H_rock,           init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"rhoc_rock",      par%rhoc_rock,        init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)
        call nml_read(filename,group,"kt_rock",        par%kt_rock,          init=init_pars,defaults_file=def_file,defaults_group=def_ytherm)

        ! Validate parameter values
        call yelmo_check_enum(group,"method",          par%method,          "enth|temp|robin|robin-cold|linear|fixed")
        call yelmo_check_enum(group,"dt_method",       par%dt_method,       "FE|AB|SAM")
        call yelmo_check_enum(group,"solver_advec",    par%solver_advec,    "expl|impl-upwind")
        call yelmo_check_enum(group,"rock_method",     par%rock_method,     "equil|active|fixed")
        call yelmo_check_enum(group,"zeta_scale_rock", par%zeta_scale_rock, "linear|exp-inv")
        call yelmo_check_enum(group,"enth_cp_method",  par%enth_cp_method,  "const|integral")
        call yelmo_check_enum(group,"strain_heating",  par%strain_heating,  "full|sia|none")

        if (par%nzr_aa .lt. 2) then
            write(io_unit_err,*) "ytherm_par_load:: error: nzr_aa must be >= 2; got ", par%nzr_aa
            error stop 1
        end if

        if (par%qb_method .lt. 1 .or. par%qb_method .gt. 4) then
            write(io_unit_err,*) "ytherm_par_load:: error: qb_method must be 1, 2, 3 or 4; got ", par%qb_method
            error stop 1
        end if

        if (par%advecxy_order .ne. 1 .and. par%advecxy_order .ne. 2) then
            write(io_unit_err,*) "ytherm_par_load:: error: advecxy_order must be 1 or 2; got ", par%advecxy_order
            error stop 1
        end if

        if (par%advecxy_cfl .le. 0.0_wp) then
            write(io_unit_err,*) "ytherm_par_load:: error: advecxy_cfl must be > 0; got ", par%advecxy_cfl
            error stop 1
        end if

        if (par%advecxy_nmax .lt. 1) then
            write(io_unit_err,*) "ytherm_par_load:: error: advecxy_nmax must be >= 1; got ", par%advecxy_nmax
            error stop 1
        end if

        ! In case of method=="temp", prescribe some parameters
        if (trim(par%method) .eq. "temp") then
            par%enth_cr   = 1.0_wp
            par%omega_max = 0.0_wp
        end if

        ! Set internal parameters
        par%nx  = nx
        par%ny  = ny 
        par%dx  = dx
        par%dy  = dx  
        par%nz_aa = size(zeta_aa,1)     ! bottom, layer centers, top 
        par%nz_ac = size(zeta_ac,1)     ! layer borders

        if (allocated(par%z%zeta_aa)) deallocate(par%z%zeta_aa)
        allocate(par%z%zeta_aa(par%nz_aa))
        par%z%zeta_aa = zeta_aa 
        
        if (allocated(par%z%zeta_ac)) deallocate(par%z%zeta_ac)
        allocate(par%z%zeta_ac(par%nz_ac))
        par%z%zeta_ac = zeta_ac 
        
        ! Calculate ice column dzeta terms 
        call calc_dzeta_terms(par%z%dzeta_a,par%z%dzeta_b,par%z%zeta_aa,par%z%zeta_ac)

        ! == Bedrock == 

        ! Calculate zeta_aa and zeta_ac 
        call calc_zeta(par%zr%zeta_aa,par%zr%zeta_ac,par%nzr_ac,par%nzr_aa, &
                                        par%zeta_scale_rock,par%zeta_exp_rock)

        ! Calculate bedrock dzeta terms 
        call calc_dzeta_terms(par%zr%dzeta_a,par%zr%dzeta_b, &
                                    par%zr%zeta_aa,par%zr%zeta_ac)

        ! Define current time as unrealistic value
        par%time = 1000000000   ! [a] 1 billion years in the future

        ! Intialize timestepping parameters to Forward Euler (beta2=0: no contribution from previous timestep)
        par%dt_zeta     = 1.0 
        par%dt_beta(1)  = 1.0 
        par%dt_beta(2)  = 0.0 

        ! Define how boundaries of grid should be treated 
        ! This should only be modified by the dom%par%experiment variable
        ! in yelmo_init. By default set boundaries to zero 
        par%boundaries = "zeros" 
        
        return

    end subroutine ytherm_par_load

    subroutine ytherm_alloc(now,nx,ny,nz_aa,nz_ac,nzr_aa)

        implicit none 

        type(ytherm_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny, nz_aa, nz_ac, nzr_aa

        call ytherm_dealloc(now)

        allocate(now%enth(nx,ny,nz_aa))
        allocate(now%T_ice(nx,ny,nz_aa))
        allocate(now%omega(nx,ny,nz_aa))
        allocate(now%T_pmp(nx,ny,nz_aa))
        allocate(now%T_prime(nx,ny,nz_aa))
        allocate(now%bmb_grnd(nx,ny))
        allocate(now%f_pmp(nx,ny))
        allocate(now%Q_strn(nx,ny,nz_aa))
        allocate(now%dQsdT(nx,ny,nz_aa))
        allocate(now%Q_b(nx,ny))
        allocate(now%Q_ice_b(nx,ny))
        allocate(now%cp(nx,ny,nz_aa))
        allocate(now%kt(nx,ny,nz_aa))
        allocate(now%H_cts(nx,ny))
        allocate(now%bmb_grnd_star(nx,ny))
        allocate(now%bc_b(nx,ny))
        allocate(now%bmb_clamp(nx,ny))
        allocate(now%melt_int(nx,ny))
        allocate(now%T_prime_b(nx,ny))
        allocate(now%advecxy(nx,ny,nz_aa))

        allocate(now%Q_rock(nx,ny))
        allocate(now%T_rock(nx,ny,nzr_aa))

        now%enth        = 0.0
        now%T_ice       = 0.0
        now%omega       = 0.0  
        now%T_pmp       = 0.0
        now%T_prime     = 0.0
        now%bmb_grnd    = 0.0 
        now%f_pmp       = 0.0 
        now%Q_strn      = 0.0 
        now%dQsdT       = 0.0 
        now%Q_b         = 0.0 
        now%Q_ice_b     = 0.0 
        now%cp          = 0.0 
        now%kt          = 0.0 
        now%H_cts       = 0.0
        now%bmb_grnd_star = 0.0
        now%bc_b          = 0.0
        now%bmb_clamp     = 0.0
        now%melt_int      = 0.0
        now%T_prime_b   = 0.0

        now%advecxy     = 0.0

        now%Q_rock      = 0.0 
        now%T_rock      = 0.0 

        return

    end subroutine ytherm_alloc

    subroutine ytherm_dealloc(now)

        implicit none 

        type(ytherm_state_class), intent(INOUT) :: now

        if (allocated(now%enth))        deallocate(now%enth)
        if (allocated(now%T_ice))       deallocate(now%T_ice)
        if (allocated(now%omega))       deallocate(now%omega)
        if (allocated(now%T_pmp))       deallocate(now%T_pmp)
        if (allocated(now%T_prime))     deallocate(now%T_prime)
        if (allocated(now%bmb_grnd))    deallocate(now%bmb_grnd)
        if (allocated(now%f_pmp))       deallocate(now%f_pmp)
        if (allocated(now%Q_strn))      deallocate(now%Q_strn)
        if (allocated(now%dQsdT))       deallocate(now%dQsdT)
        if (allocated(now%Q_b))         deallocate(now%Q_b)
        if (allocated(now%Q_ice_b))     deallocate(now%Q_ice_b)
        if (allocated(now%cp))          deallocate(now%cp)
        if (allocated(now%kt))          deallocate(now%kt)
        if (allocated(now%H_cts))       deallocate(now%H_cts)
        if (allocated(now%bmb_grnd_star)) deallocate(now%bmb_grnd_star)
        if (allocated(now%bc_b))          deallocate(now%bc_b)
        if (allocated(now%bmb_clamp))     deallocate(now%bmb_clamp)
        if (allocated(now%melt_int))      deallocate(now%melt_int)
        if (allocated(now%T_prime_b))   deallocate(now%T_prime_b)

        if (allocated(now%advecxy))     deallocate(now%advecxy)

        if (allocated(now%Q_rock))      deallocate(now%Q_rock)
        if (allocated(now%T_rock))      deallocate(now%T_rock)

        return 

    end subroutine ytherm_dealloc

end module yelmo_thermodynamics
