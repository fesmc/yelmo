module bench_forcing
    ! Online forcing of the benchmark protocol (docs/dev/benchmark-protocol),
    ! used by tests/yelmo_bench.f90. Each routine has a Julia counterpart in
    ! tests/bench/src (YelmoBench), which generates the fixed fields and is
    ! used to check the forcing written by the model.

    use yelmo_defs, only : wp

    implicit none

    private
    public :: bench_tsrf_lapse
    public :: bench_bmb_mismipplus

contains

    elemental subroutine bench_tsrf_lapse(T_srf,z_srf,T_sl,lapse,T0)
        ! Surface temperature from a lapse rate on the surface elevation.
        ! Julia counterpart: island4_tsrf

        implicit none

        real(wp), intent(OUT) :: T_srf      ! [K]
        real(wp), intent(IN)  :: z_srf      ! [m] Surface elevation
        real(wp), intent(IN)  :: T_sl       ! [degC] Temperature at sea level
        real(wp), intent(IN)  :: lapse      ! [K m-1] Lapse rate
        real(wp), intent(IN)  :: T0         ! [K] Freezing temperature

        T_srf = T0 + T_sl - lapse*z_srf

        return

    end subroutine bench_tsrf_lapse

    elemental subroutine bench_bmb_mismipplus(bmb_shlf,H_ice,z_bed,z_sl,rho_ice,rho_sw,Omega,Hc0,z0)
        ! Basal mass balance below floating ice, MISMIP+ Ice1 parameterization
        ! (Asay-Davis et al., 2016): m = Omega tanh(H_c/Hc0) max(z0 - z_d, 0),
        ! where z_d is the depth of the ice base and H_c the water-column
        ! thickness. Zero for grounded ice and ice-free points.
        ! Julia counterpart: island4_bmb

        implicit none

        real(wp), intent(OUT) :: bmb_shlf   ! [m/yr] Negative for melt
        real(wp), intent(IN)  :: H_ice      ! [m]
        real(wp), intent(IN)  :: z_bed      ! [m]
        real(wp), intent(IN)  :: z_sl       ! [m]
        real(wp), intent(IN)  :: rho_ice    ! [kg m-3]
        real(wp), intent(IN)  :: rho_sw     ! [kg m-3]
        real(wp), intent(IN)  :: Omega      ! [yr-1]
        real(wp), intent(IN)  :: Hc0        ! [m]
        real(wp), intent(IN)  :: z0         ! [m]

        ! Local variables
        real(wp) :: z_d, H_c

        z_d = max(z_bed, z_sl - H_ice*rho_ice/rho_sw)
        H_c = z_d - z_bed

        bmb_shlf = -Omega*tanh(H_c/Hc0)*max(z0 - z_d, 0.0_wp)

        return

    end subroutine bench_bmb_mismipplus

end module bench_forcing
