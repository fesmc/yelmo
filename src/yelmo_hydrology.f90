module yelmo_hydrology
    ! Thin yelmo-side wrapper around fasthydrology. Owns the
    ! init / init_state / step calls; builds the input arrays
    ! (mask, bmb_w, A_glen) expected by fasthydrology from the
    ! corresponding yelmo state.
    !
    ! fasthydrology writes hyd%now%N back to yelmo, and dyn%now%N_eff
    ! is read from hyd%now%N when ydyn neff_method=6 (see calc_ydyn_neff).
    ! Other paths (thrm%now%H_w, the legacy ytherm bucket, neff_method=1..5)
    ! still co-exist during the migration and are untouched here.

    use yelmo_defs
    use fast_hydrology, only : hydro_init, hydro_init_state, hydro_update, SEC_PER_YEAR

    implicit none

    private
    public :: yhyd_par_load
    public :: yhyd_init_state
    public :: calc_yhyd

contains

    subroutine yhyd_par_load(hyd, filename, group, nx, ny, dx, dy)
        ! Load fasthydrology parameters from the yelmo namelist. The grid
        ! spacing is passed to hydro_init so the user does not have to
        ! keep dx / dy in sync in the namelist (they are no longer
        ! namelist-loaded by fasthydrology).

        type(hydro_class), intent(INOUT) :: hyd
        character(len=*),  intent(IN)    :: filename
        character(len=*),  intent(IN)    :: group
        integer,           intent(IN)    :: nx, ny
        real(wp),          intent(IN)    :: dx, dy

        ! Defaults overlay and typo validation now happen inside
        ! hydro_init (fast_hydrology) using input/yelmo_defaults.nml.
        call hydro_init(hyd, filename, nx, ny, dx, dy, group=group)

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

        call hydro_init_state(hyd, bnd%z_bed, tpo%now%f_ice, tpo%now%f_grnd, time)

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
        !       rho_ice / rho_w     ->  bmb_w   (water-equivalent m/a, +ve = source)
        !   dyn%now%uxy_b           ->  uxy_b
        !   mat%now%ATT(:,:,1)      ->  A_glen  (basal layer; zeta_aa(1) = 0)
        !   time                    ->  time
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
        real(wp), allocatable :: A_glen_b(:,:)
        integer :: nx, ny

        nx = size(tpo%now%H_ice, 1)
        ny = size(tpo%now%H_ice, 2)

        allocate(mask(nx,ny))
        allocate(bmb_w(nx,ny))
        allocate(A_glen_b(nx,ny))

        ! Active-hydrology mask: grounded ice cells.
        where (tpo%now%f_ice >= 0.5_wp .and. tpo%now%f_grnd > 0.0_wp)
            mask = 1.0_wp
        elsewhere
            mask = 0.0_wp
        end where

        ! Convert grounded basal mass balance (ice-equivalent m/a,
        ! +ve = accumulation) to water-equivalent melt (+ve = water source).
        ! fasthydrology's mdot contract is SI [m/s], so divide by SEC_PER_YEAR;
        ! omitting this overscales the source by ~3.16e7 and pegs any melting
        ! cell straight to W_til_max.
        bmb_w = -thrm%now%bmb_grnd * (bnd%c%rho_ice / bnd%c%rho_w) / SEC_PER_YEAR

        ! Basal Glen-A. zeta_aa(1) = 0 in yelmo, so index 1 is the base.
        A_glen_b = mat%now%ATT(:,:,1)

        call hydro_update(hyd, tpo%now%H_ice, bnd%z_bed, bnd%z_sl,        &
                          tpo%now%f_ice, tpo%now%f_grnd, mask,            &
                          bmb_w, dyn%now%uxy_b, A_glen_b, time)

        deallocate(mask, bmb_w, A_glen_b)

        return

    end subroutine calc_yhyd

end module yelmo_hydrology
