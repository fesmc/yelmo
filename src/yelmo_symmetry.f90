module yelmo_symmetry
    ! Reusable symmetry-regression check for 2D fields.
    !
    ! Computes the normalized asymmetry of a 2D field under four reflections
    ! (left<->right, top<->bottom, 180deg rotation, and transpose for square
    ! grids). Intended for verifying that symmetric benchmark configurations
    ! (e.g. EISMINT moving/expa/expf, and later CalvingMIP exp1) preserve the
    ! symmetry of the imposed forcing over the course of a simulation.
    !
    ! Two usage modes:
    !   * ONLINE  - called at the end of a benchmark run (see tests/yelmo_benchmarks.f90)
    !   * OFFLINE - via the standalone driver tests/test_symmetry.f90 pointed at
    !               a yelmo2D.nc (or any NetCDF file with a 2D+time field).

    use yelmo_defs, only : wp

    implicit none

    private

    ! Values with |v| > SYM_VALID_MAX are treated as invalid/missing (sentinel
    ! guard, e.g. -9999 missing values or other fill flags).
    real(wp), parameter :: SYM_VALID_MAX = 9000.0_wp

    type sym_metrics_class
        ! Per-reflection metrics: L1 = mean(|reflected-H|)/Hmax,
        !                         Linf = max(|reflected-H|)/Hmax,
        ! evaluated over cells where BOTH sides are valid.
        real(wp) :: l1_lr,  linf_lr      ! left  <-> right
        real(wp) :: l1_tb,  linf_tb      ! top   <-> bottom
        real(wp) :: l1_rot, linf_rot     ! 180-degree rotation
        real(wp) :: l1_tr,  linf_tr      ! transpose (square grids only)
        logical  :: has_tr               ! .TRUE. only when nx==ny
        real(wp) :: hmax                 ! max(|H|) over valid cells
        integer  :: n_valid              ! number of valid cells in H
    end type sym_metrics_class

    public :: sym_metrics_class
    public :: calc_symmetry_metrics
    public :: report_symmetry

contains

    subroutine calc_symmetry_metrics(H, tol, met, all_pass)
        ! Compute normalized asymmetry metrics of field H under four reflections.

        implicit none

        real(wp),                intent(IN)  :: H(:,:)
        real(wp),                intent(IN)  :: tol
        type(sym_metrics_class), intent(OUT) :: met
        logical,                 intent(OUT) :: all_pass

        ! Local variables
        integer  :: nx, ny, i, j, cnt
        real(wp) :: hmax, av
        real(wp), allocatable :: R(:,:)

        nx = size(H,1)
        ny = size(H,2)

        ! Hmax and valid-cell count over cells where |H| <= sentinel
        hmax = 0.0_wp
        cnt  = 0
        do j = 1, ny
        do i = 1, nx
            if (abs(H(i,j)) <= SYM_VALID_MAX) then
                cnt = cnt + 1
                av  = abs(H(i,j))
                if (av > hmax) hmax = av
            end if
        end do
        end do

        met%hmax    = hmax
        met%n_valid = cnt
        met%has_tr  = (nx == ny)

        ! Initialize all metrics to zero
        met%l1_lr  = 0.0_wp ; met%linf_lr  = 0.0_wp
        met%l1_tb  = 0.0_wp ; met%linf_tb  = 0.0_wp
        met%l1_rot = 0.0_wp ; met%linf_rot = 0.0_wp
        met%l1_tr  = 0.0_wp ; met%linf_tr  = 0.0_wp

        ! Guard: no ice (or entirely invalid) => trivially symmetric
        if (hmax == 0.0_wp) then
            all_pass = .TRUE.
            return
        end if

        allocate(R(nx,ny))

        ! Left <-> Right
        R = H(nx:1:-1,:)
        call refl_metric(H, R, hmax, met%l1_lr, met%linf_lr)

        ! Top <-> Bottom
        R = H(:,ny:1:-1)
        call refl_metric(H, R, hmax, met%l1_tb, met%linf_tb)

        ! 180-degree rotation
        R = H(nx:1:-1,ny:1:-1)
        call refl_metric(H, R, hmax, met%l1_rot, met%linf_rot)

        ! Transpose (only meaningful for square grids)
        if (met%has_tr) then
            R = transpose(H)
            call refl_metric(H, R, hmax, met%l1_tr, met%linf_tr)
        end if

        deallocate(R)

        ! Overall pass = all COMPUTED Linf metrics below tolerance
        all_pass = (met%linf_lr  < tol) .and. &
                   (met%linf_tb  < tol) .and. &
                   (met%linf_rot < tol)
        if (met%has_tr) all_pass = all_pass .and. (met%linf_tr < tol)

        return

    end subroutine calc_symmetry_metrics

    subroutine refl_metric(H, R, hmax, l1, linf)
        ! Normalized L1 and Linf asymmetry of reflected field R vs original H,
        ! over cells where both sides are valid (|value| <= sentinel).

        implicit none

        real(wp), intent(IN)  :: H(:,:)
        real(wp), intent(IN)  :: R(:,:)
        real(wp), intent(IN)  :: hmax          ! assumed > 0 by caller
        real(wp), intent(OUT) :: l1
        real(wp), intent(OUT) :: linf

        ! Local variables
        integer  :: nx, ny, i, j, cnt
        real(wp) :: d, s, mx

        nx = size(H,1)
        ny = size(H,2)

        s   = 0.0_wp
        mx  = 0.0_wp
        cnt = 0

        do j = 1, ny
        do i = 1, nx
            if (abs(H(i,j)) <= SYM_VALID_MAX .and. abs(R(i,j)) <= SYM_VALID_MAX) then
                d   = abs(R(i,j) - H(i,j))
                s   = s + d
                if (d > mx) mx = d
                cnt = cnt + 1
            end if
        end do
        end do

        if (cnt > 0) then
            l1   = (s / real(cnt,wp)) / hmax
            linf = mx / hmax
        else
            l1   = 0.0_wp
            linf = 0.0_wp
        end if

        return

    end subroutine refl_metric

    subroutine report_symmetry(label, met, tol, all_pass)
        ! Print a compact one-row-per-reflection symmetry table plus overall result.

        implicit none

        character(len=*),        intent(IN) :: label
        type(sym_metrics_class), intent(IN) :: met
        real(wp),                intent(IN) :: tol
        logical,                 intent(IN) :: all_pass

        write(*,*) ""
        write(*,*) "=== Symmetry check: "//trim(label)//" ==="
        write(*,"(a,g13.4,a,i0)") "  Hmax = ", met%hmax, "   valid cells = ", met%n_valid
        write(*,"(a,g13.4)")      "  tolerance (Linf/Hmax) < ", tol
        write(*,"(a)") "  reflection      L1/Hmax       Linf/Hmax     result"
        call print_row("L<->R    ", met%l1_lr,  met%linf_lr,  tol, .TRUE.)
        call print_row("T<->B    ", met%l1_tb,  met%linf_tb,  tol, .TRUE.)
        call print_row("rot180   ", met%l1_rot, met%linf_rot, tol, .TRUE.)
        call print_row("transpose", met%l1_tr,  met%linf_tr,  tol, met%has_tr)
        write(*,*) ""
        if (all_pass) then
            write(*,"(a)") "  Symmetry check: overall PASS"
        else
            write(*,"(a)") "  Symmetry check: overall FAIL"
        end if
        write(*,*) ""

        return

    end subroutine report_symmetry

    subroutine print_row(name, l1, linf, tol, computed)
        ! Print a single reflection row (or N/A when the metric was not computed).

        implicit none

        character(len=*), intent(IN) :: name
        real(wp),         intent(IN) :: l1
        real(wp),         intent(IN) :: linf
        real(wp),         intent(IN) :: tol
        logical,          intent(IN) :: computed

        character(len=4) :: res

        if (.not. computed) then
            write(*,"(a,a14,a14,a)") "  "//name//"  ", "N/A", "N/A", "    N/A"
        else
            if (linf < tol) then
                res = "PASS"
            else
                res = "FAIL"
            end if
            write(*,"(a,2g14.4,a)") "  "//name//"  ", l1, linf, "    "//trim(res)
        end if

        return

    end subroutine print_row

end module yelmo_symmetry
