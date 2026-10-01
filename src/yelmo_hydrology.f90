module yelmo_hydrology
    ! Thin yelmo-side wrapper around fasthydrology. Owns the
    ! init / init_state / step calls; builds the input arrays
    ! (mask, bmb_w, A_glen) expected by fasthydrology from the
    ! corresponding yelmo state.
    !
    ! fasthydrology owns the till water W_til (hyd%now%W_til, read by
    ! ytherm) and the effective pressure hyd%now%N, which calc_ydyn_neff
    ! copies (or subgrid-averages, ydyn.neff_nxi) into dyn%now%N_eff.

    use yelmo_defs
    use yelmo_tools, only : boundary_code, get_periodic_directions
    use fast_hydrology, only : hydro_init, hydro_init_state, hydro_update

    implicit none

    private
    public :: yhyd_par_load
    public :: yhyd_init_state
    public :: calc_yhyd

contains

    subroutine yhyd_par_load(hyd, filename, group, nx, ny, dx, dy, c, cnst, boundaries)
        ! Load fasthydrology parameters from the yelmo namelist. The grid
        ! spacing is passed to hydro_init so the user does not have to
        ! keep dx / dy in sync in the namelist (they are no longer
        ! namelist-loaded by fasthydrology), and so are the periodic
        ! directions of the domain, so that fasthydrology's neighbour
        ! stencils wrap there and its border BC is not applied.
        !
        ! cnst hands fasthydrology the domain's physical constants, and
        ! c%sec_year its calendar year, so that N and p_w are consistent with
        ! yelmo's overburden and flotation criterion (marine N vanishes where
        ! yelmo's H_grnd does), the K24 melt and opening terms use the same
        ! L_ice as ytherm, and the per-year rates below convert with the year
        ! the domain actually chose. fasthydrology resolves all of it inside
        ! hydro_init, before anything derived from the densities is computed.

        type(hydro_class),        intent(INOUT) :: hyd
        character(len=*),         intent(IN)    :: filename
        character(len=*),         intent(IN)    :: group
        integer,                  intent(IN)    :: nx, ny
        real(wp),                 intent(IN)    :: dx, dy
        type(ybound_const_class), intent(IN)    :: c        ! Physical constants of the domain
        type(phys_const_class),   intent(IN)    :: cnst     ! The same set, as the shared record
        character(len=*),         intent(IN)    :: boundaries

        ! Local variables
        logical :: per_x, per_y

        call get_periodic_directions(per_x,per_y,boundary_code(boundaries))

        ! Defaults overlay and typo validation now happen inside
        ! hydro_init (fast_hydrology) using input/yelmo_defaults.nml.
        call hydro_init(hyd, filename, nx, ny, dx, dy, group=group, &
                        cnst=cnst, sec_year=c%sec_year, &
                        periodic_x=per_x, periodic_y=per_y)

        return

    end subroutine yhyd_par_load

    subroutine yhyd_init_state(hyd, bnd, tpo, time)
        ! Seed fasthydrology's state from the freshly-initialized
        ! yelmo topo / boundary fields. Must be called after f_ice /
        ! f_grnd are populated (i.e. after yelmo_init_topo) and before
        ! the first calc_yhyd.

        type(hydro_class),  intent(INOUT) :: hyd
        type(ybound_class), intent(IN)    :: bnd
        type(ytopo_class),  intent(IN)    :: tpo
        real(wp),           intent(IN)    :: time

        call hydro_init_state(hyd, tpo%now%H_ice, bnd%z_bed, tpo%now%f_ice, tpo%now%f_grnd, time)

        return

    end subroutine yhyd_init_state

    subroutine calc_yhyd(hyd, tpo, dyn, mat, thrm, bnd, time)
        ! Advance fasthydrology by one yelmo timestep using the current
        ! state of all upstream yelmo components.
        !
        ! Argument mapping (yelmo  ->  fasthydrology):
        !   tpo%now%H_ice           ->  H_ice
        !   bnd%z_bed               ->  z_bed
        !   bnd%z_sl                ->  z_sl
        !   tpo%now%f_ice           ->  f_ice
        !   tpo%now%f_grnd          ->  f_grnd
        !   {f_ice >= 0.5 .and.
        !    f_grnd > 0.0}          ->  mask    (1.0 = active hydrology cell)
        !   -thrm%now%bmb_grnd *
        !       rho_ice / rho_w     ->  mdot    (water-equivalent, +ve = source)
        !   dyn%now%uxy_b           ->  uxy_b
        !   bnd%Q_geo               ->  G      (geothermal heat into the bed)
        !   thrm%now%Q_ice_b        ->  q_T    (conductive heat into the ice)
        !   dyn%now%ux_b, uy_b      ->  ux_b, uy_b (C-grid, staggered friction)
        !   dyn%now%taub            ->  taub   (prescribed-field sliding law)
        !   dyn%now%cb_ref          ->  c_till (regularized-Coulomb-field law)
        !   i_eb = 0 for now: Yelmo's drained englacial water (melt_internal)
        !   is not yet kept as a field.
        !
        ! K24 builds its water source from these terms (and the friction and
        ! dissipation heat it computes itself); mdot (from bmb_grnd) drives the
        ! bucket only. G and q_T are mW m-2 in Yelmo, W m-2 in fasthydrology.
        !   mat%now%ATT(:,:,1)      ->  A_glen  (basal layer; zeta_aa(1) = 0)
        !   time                    ->  time    [a]
        !
        ! fasthydrology's rate inputs (mdot, uxy_b, A_glen) are SI (per second),
        ! while yelmo's are per year. All three are divided by the domain's
        ! sec_year, which yhyd_par_load also handed to fasthydrology as the
        ! constant it converts the time step to dt_sec with, so that e.g.
        ! mdot*dt_sec is exactly bmb_w*dt [m] and the solution does not depend
        ! on the choice of seconds per year.
        !
        ! fasthydrology's hydro_update internally skips work when
        ! dt = time - hyd%now%time is non-positive, so we can call
        ! unconditionally and let it handle the no-step case.

        type(hydro_class),  intent(INOUT) :: hyd
        type(ytopo_class),  intent(IN)    :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd
        real(wp),           intent(IN)    :: time

        ! Local scratch arrays sized to the yelmo grid
        real(wp), allocatable :: mask(:,:)
        real(wp), allocatable :: bmb_w(:,:)
        real(wp), allocatable :: uxy_b(:,:)
        real(wp), allocatable :: A_glen_b(:,:)
        real(wp), allocatable :: G(:,:), q_T(:,:), i_eb(:,:)
        real(wp), allocatable :: ux_b(:,:), uy_b(:,:)
        integer :: nx, ny

        nx = size(tpo%now%H_ice, 1)
        ny = size(tpo%now%H_ice, 2)

        allocate(mask(nx,ny))
        allocate(bmb_w(nx,ny))
        allocate(uxy_b(nx,ny))
        allocate(A_glen_b(nx,ny))
        allocate(G(nx,ny), q_T(nx,ny), i_eb(nx,ny), ux_b(nx,ny), uy_b(nx,ny))

        ! Active-hydrology mask: grounded ice cells.
        where (tpo%now%f_ice >= 0.5_wp .and. tpo%now%f_grnd > 0.0_wp)
            mask = 1.0_wp
        elsewhere
            mask = 0.0_wp
        end where

        ! Convert grounded basal mass balance (ice-equivalent m/a,
        ! +ve = accumulation) to water-equivalent melt (+ve = water source).
        ! fasthydrology's mdot contract is SI [m/s], so divide by sec_year;
        ! omitting this overscales the source by ~3.16e7 and pegs any melting
        ! cell straight to W_til_max.
        bmb_w = -thrm%now%bmb_grnd * (bnd%c%rho_ice / bnd%c%rho_w) / bnd%c%sec_year

        ! Basal sliding speed [m/a] -> [m/s]. K24's sliding laws compare it
        ! with u0 given in m/s, and its N_inf balances sliding opening
        ! against melt opening (Q*|grad phi|, already per second).
        uxy_b = dyn%now%uxy_b / bnd%c%sec_year

        ! Basal Glen-A [Pa^-3 a^-1] -> [Pa^-3 s^-1]. zeta_aa(1) = 0 in
        ! yelmo, so index 1 is the base.
        A_glen_b = mat%now%ATT(:,:,1) / bnd%c%sec_year

        ! Terms of the K24 water source [W m-2]; i_eb [m/s water-equivalent].
        G    = bnd%Q_geo        * 1e-3_wp
        q_T  = thrm%now%Q_ice_b * 1e-3_wp
        i_eb = 0.0_wp

        ! C-grid basal velocities [m/a] -> [m/s] (staggered friction).
        ux_b = dyn%now%ux_b / bnd%c%sec_year
        uy_b = dyn%now%uy_b / bnd%c%sec_year

        call hydro_update(hyd, tpo%now%H_ice, bnd%z_bed, bnd%z_sl,        &
                          tpo%now%f_ice, tpo%now%f_grnd, mask,            &
                          bmb_w, G, q_T, i_eb, uxy_b, A_glen_b, time,     &
                          ux_b=ux_b, uy_b=uy_b, taub=dyn%now%taub,        &
                          c_till=dyn%now%cb_ref)

        deallocate(mask, bmb_w, uxy_b, A_glen_b)
        deallocate(G, q_T, i_eb, ux_b, uy_b)

        return

    end subroutine calc_yhyd

end module yelmo_hydrology
