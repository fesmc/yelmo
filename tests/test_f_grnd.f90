program test_f_grnd
    ! Grounded fraction of the bilinear H_grnd interpolant (gl_sep = 3,
    ! calc_f_grnd_subgrid_area / bilinear_grounded_fraction).
    !
    ! 1. Unit square, exact values: linear fields (straight grounding line,
    !    cross term d = 0), the hyperbola x*y = 1/4 (0.75 - ln(4)/4), an
    !    exact saddle (two straight zero lines crossing: 1/2), H_grnd = 0
    !    (grounded), and invariance under the 8 symmetries of the square.
    ! 2. Cells: random fields (smooth, isolated grounded cells, small
    !    integers with ties and zeros) against dense sampling of the same
    !    interpolant (fesm-utils calc_subgrid_array_quad).
    ! 3. A cell grounded at its centre next to deep ocean has f_grnd > 0.

    use yelmo_defs, only : wp
    use topography, only : calc_f_grnd_subgrid_area, bilinear_grounded_fraction
    use subgrid,    only : calc_subgrid_array_quad

    implicit none

    integer,  parameter :: n = 16, nxi = 300, n_fields = 30
    real(wp), parameter :: tol_exact  = 1.0e-5_wp
    real(wp), parameter :: tol_sample = 5.0e-3_wp     ! dense sampling with nxi*nxi points

    real(wp) :: v(4), f0, err, err_max
    real(wp) :: h(n,n), g(n,n), gx(n,n), gy(n,n), gab(n,n), vint(nxi,nxi), fref
    integer  :: k, i, j, n_fail
    integer, allocatable :: seed(:)

    n_fail = 0

    ! Fixed seed, so that the random fields are the same in every run
    call random_seed(size=k)
    allocate(seed(k))
    seed = 20261008
    call random_seed(put=seed)

    ! === 1. Unit square =====================================================

    call check_exact("linear in x, root at 0.3",   [-0.3_wp,0.7_wp,-0.3_wp,0.7_wp],  0.7_wp)
    call check_exact("linear, x + y = 1/2",        [-0.5_wp,0.5_wp,0.5_wp,1.5_wp],   0.875_wp)
    call check_exact("x*y >= 1/4",                 [-0.25_wp,-0.25_wp,-0.25_wp,0.75_wp], &
                                                   0.75_wp - 0.25_wp*log(4.0_wp))
    call check_exact("exact saddle",               [1.0_wp,-1.0_wp,-1.0_wp,1.0_wp],  0.5_wp)
    call check_exact("saddle, zero lines at edge", [0.0_wp,0.0_wp,1.0_wp,-1.0_wp],   0.5_wp)
    call check_exact("all zero (grounded)",        [0.0_wp,0.0_wp,0.0_wp,0.0_wp],    1.0_wp)
    call check_exact("all floating",               [-1.0_wp,-2.0_wp,-3.0_wp,-4.0_wp],0.0_wp)

    ! Symmetries of the square: (h00,h10,h01,h11) rotated and reflected
    err_max = 0.0_wp
    do k = 1, 2000
        call random_number(v)
        v  = 100.0_wp*(v-0.5_wp)
        f0 = bilinear_grounded_fraction(v(1),v(2),v(3),v(4))
        err_max = max(err_max, abs(bilinear_grounded_fraction(v(2),v(4),v(1),v(3)) - f0))   ! rotation
        err_max = max(err_max, abs(bilinear_grounded_fraction(v(4),v(3),v(2),v(1)) - f0))   ! rotation by pi
        err_max = max(err_max, abs(bilinear_grounded_fraction(v(2),v(1),v(4),v(3)) - f0))   ! reflection in x
        err_max = max(err_max, abs(bilinear_grounded_fraction(v(1),v(3),v(2),v(4)) - f0))   ! transpose
    end do
    call report("symmetries of the square (2000 random squares)", err_max, tol_exact)

    ! === 2. Cells against dense sampling ===================================

    err_max = 0.0_wp
    do k = 1, n_fields
        call random_number(h)
        select case(mod(k,3))
            case(0)
                h = 200.0_wp*(h-0.5_wp)
            case(1)
                h = merge(-300.0_wp, 1.0_wp+50.0_wp*h, h .lt. 0.7_wp)
            case(2)
                h = real(nint(6.0_wp*(h-0.5_wp)),wp)
        end select

        call calc_f_grnd_subgrid_area(g,gx,gy,gab,h,"zeros")

        do j = 2, n-1
        do i = 2, n-1
            call calc_subgrid_array_quad(vint,h,nxi,i,j,i-1,i+1,j-1,j+1)
            fref    = real(count(vint .ge. 0.0_wp),wp) / real(nxi*nxi,wp)
            err_max = max(err_max, abs(g(i,j)-fref))
        end do
        end do

        if (any(g .lt. 0.0_wp) .or. any(g .gt. 1.0_wp) .or. any(g .ne. g)) then
            write(*,*) "FAIL: f_grnd outside [0,1] or NaN in field ", k
            n_fail = n_fail + 1
        end if
    end do
    call report("cells vs dense sampling (random fields)", err_max, tol_sample)

    ! === 3. Grounded centre next to deep ocean =============================

    h = -300.0_wp
    h(n/2,n/2) = 20.0_wp
    call calc_f_grnd_subgrid_area(g,gx,gy,gab,h,"zeros")
    write(*,'(a,f8.4)') "  isolated grounded cell: f_grnd = ", g(n/2,n/2)
    if (.not. (g(n/2,n/2) .gt. 0.0_wp)) then
        write(*,*) "FAIL: isolated grounded cell has f_grnd = 0"
        n_fail = n_fail + 1
    end if

    if (n_fail .eq. 0) then
        write(*,*) "test_f_grnd: PASSED"
    else
        write(*,*) "test_f_grnd: FAILED (", n_fail, " checks)"
        error stop 1
    end if

contains

    subroutine check_exact(label,c,f_exact)
        character(len=*), intent(IN) :: label
        real(wp),         intent(IN) :: c(4)            ! h00, h10, h01, h11
        real(wp),         intent(IN) :: f_exact
        call report(label, abs(bilinear_grounded_fraction(c(1),c(2),c(3),c(4)) - f_exact), tol_exact)
    end subroutine check_exact

    subroutine report(label,e,tol)
        character(len=*), intent(IN) :: label
        real(wp),         intent(IN) :: e, tol
        if (e .le. tol) then
            write(*,'(a,a,es10.2)') "  ok    ", label//": err = ", e
        else
            write(*,'(a,a,es10.2,a,es10.2)') "  FAIL  ", label//": err = ", e, " > ", tol
            n_fail = n_fail + 1
        end if
    end subroutine report

end program test_f_grnd
