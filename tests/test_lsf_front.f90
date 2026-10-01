program test_lsf_front
    ! Subgrid front following the level set (front_subgrid with use_lsf).
    !
    ! Cases (dx = 8 km, z_sl = 0, dt = 1 yr):
    !   1. a_lsf of a 1-cell-wide tongue is continuous in the tongue lsf
    !      (0.667 at lsf = -0.5, near 1 at lsf = -0.99)
    !   2. a_lsf of a half-covered cell (lsf = 0 between -1 and +1) is 0.5
    !   3. marine-grounded cliff, stationary level set, 5 m/yr refill:
    !      repeated trimming holds H = a_lsf*H_nb (no compounding)
    !   4. isolated 2-cell-wide floating tongue, stationary level set:
    !      no interior neighbour, so the thickness is unchanged
    !   5. outflow from a partial cell into an ice-free cell: closed without
    !      a prescribed front, open where a_front >= A_FRONT_MIN
    !   6. marine-grounded front cell on a deeper bed, entirely behind the
    !      level set (a_lsf = 1): full; without the level set, partial from
    !      the neighbour reference

    use yelmo_defs,        only : wp, A_FRONT_MIN
    use topography,        only : calc_lsf_area_fraction, calc_ice_fraction
    use mass_conservation, only : calc_G_lsf_front
    use velocity_general,  only : set_inactive_margins

    implicit none

    real(wp), parameter :: dx = 8000.0_wp, dt = 1.0_wp
    real(wp), parameter :: rho_ice = 917.0_wp, rho_sw = 1028.0_wp
    real(wp), parameter :: H_eff_min = 50.0_wp, dHdx = 0.0_wp
    real(wp), parameter :: tol = 1.0e-3_wp
    integer,  parameter :: n_steps = 40

    integer :: n_fail

    n_fail = 0

    write(*,*) "==============================================="
    write(*,*) " Subgrid front following the level set"
    write(*,*) "==============================================="

    call test_tongue_area()
    call test_half_cell()
    call test_cliff()
    call test_tongue_trim()
    call test_gate()
    call test_behind_front()

    write(*,*)
    if (n_fail .eq. 0) then
        write(*,*) " PASS"
    else
        write(*,'(a,i0,a)') "  FAIL (", n_fail, " cases)"
        error stop 1
    end if

contains

    subroutine check(name,ok,val,ref)
        character(len=*), intent(IN) :: name
        logical,  intent(IN) :: ok
        real(wp), intent(IN) :: val, ref
        write(*,'(2x,a40,2f12.4,3x,a)') name, val, ref, merge("ok  ","FAIL",ok)
        if (.not. ok) n_fail = n_fail + 1
    end subroutine check

    subroutine test_tongue_area()
        ! 1-cell-wide tongue along x (row j = 4) in the ocean
        integer,  parameter :: nx = 9, ny = 7
        real(wp) :: lsf(nx,ny), a(nx,ny)

        lsf = 1.0_wp
        lsf(3:7,4) = -0.5_wp
        call calc_lsf_area_fraction(a,lsf,"infinite")
        call check("1 1-wide tongue a_lsf, lsf=-0.5",abs(a(5,4)-2.0_wp/3.0_wp) .lt. tol,a(5,4),2.0_wp/3.0_wp)

        lsf(3:7,4) = -0.99_wp
        call calc_lsf_area_fraction(a,lsf,"infinite")
        call check("1 1-wide tongue a_lsf, lsf=-0.99",a(5,4) .gt. 0.99_wp,a(5,4),0.995_wp)
    end subroutine test_tongue_area

    subroutine test_half_cell()
        ! Straight front through the centre of column i = 4
        integer,  parameter :: nx = 8, ny = 5
        real(wp) :: lsf(nx,ny), a(nx,ny)

        lsf = 1.0_wp
        lsf(1:3,:) = -1.0_wp
        lsf(4,:)   =  0.0_wp
        call calc_lsf_area_fraction(a,lsf,"infinite")
        call check("2 half-covered cell a_lsf",abs(a(4,3)-0.5_wp) .lt. tol,a(4,3),0.5_wp)
    end subroutine test_half_cell

    subroutine test_cliff()
        ! Marine-grounded cliff, uniform in y: ice in i <= 6, front in cell 6
        integer,  parameter :: nx = 12, ny = 4
        real(wp) :: H(nx,ny), lsf(nx,ny), a(nx,ny), cmb(nx,ny), z_bed(nx,ny), z_sl(nx,ny)
        real(wp) :: H_tgt, H_flot
        integer  :: n

        z_bed = -300.0_wp
        z_sl  = 0.0_wp
        H     = 0.0_wp
        H(1:6,:) = 1000.0_wp
        lsf   = 1.0_wp
        lsf(1:5,:) = -1.0_wp
        lsf(6,:)   = -0.4_wp
        call calc_lsf_area_fraction(a,lsf,"infinite")

        do n = 1, n_steps
            H(6,:) = H(6,:) + 5.0_wp*dt
            call calc_G_lsf_front(cmb,H,a,z_bed,z_sl,rho_ice,rho_sw,"marine",H_eff_min,dHdx,dx,dt,"infinite")
            H = H + cmb*dt
        end do

        H_tgt  = a(6,2)*1000.0_wp
        H_flot = 300.0_wp*rho_sw/rho_ice
        call check("3 marine cliff H after trims",abs(H(6,2)-H_tgt) .lt. tol*H_tgt,H(6,2),H_tgt)
        call check("3 marine cliff stays grounded",H(6,2) .gt. H_flot,H(6,2),H_flot)
        call check("3 interior cell unchanged",H(5,2) .eq. 1000.0_wp,H(5,2),1000.0_wp)
    end subroutine test_cliff

    subroutine test_tongue_trim()
        ! Isolated 2-cell-wide floating tongue, every cell a front cell
        integer,  parameter :: nx = 12, ny = 12
        real(wp) :: H(nx,ny), lsf(nx,ny), a(nx,ny), cmb(nx,ny), z_bed(nx,ny), z_sl(nx,ny)
        integer  :: n

        z_bed = -1500.0_wp
        z_sl  = 0.0_wp
        H     = 0.0_wp
        H(4:7,6:7) = 400.0_wp
        lsf   = 1.0_wp
        lsf(4:7,6:7) = -0.9_wp
        call calc_lsf_area_fraction(a,lsf,"infinite")

        do n = 1, n_steps
            call calc_G_lsf_front(cmb,H,a,z_bed,z_sl,rho_ice,rho_sw,"marine",H_eff_min,dHdx,dx,dt,"infinite")
            H = H + cmb*dt
        end do

        call check("4 2-wide tongue min a_lsf",minval(a(4:7,6:7)) .ge. A_FRONT_MIN,minval(a(4:7,6:7)),A_FRONT_MIN)
        call check("4 2-wide tongue min H after trims",minval(H(4:7,6:7)) .eq. 400.0_wp,minval(H(4:7,6:7)),400.0_wp)
    end subroutine test_tongue_trim

    subroutine test_gate()
        ! Faces along x: full | partial | ice-free | ice-free
        integer,  parameter :: nx = 4, ny = 3
        real(wp) :: f_ice(nx,ny), ux(nx,ny), uy(nx,ny), a_front(nx,ny)

        f_ice = 0.0_wp
        f_ice(1,:) = 1.0_wp
        f_ice(2,:) = 0.5_wp

        ux = 100.0_wp; uy = 0.0_wp
        call set_inactive_margins(ux,uy,f_ice,"infinite")
        call check("5 partial->free, no front: closed",ux(2,2) .eq. 0.0_wp,ux(2,2),0.0_wp)

        a_front = 0.0_wp
        a_front(1:2,:) = 1.0_wp
        a_front(3,:)   = 0.5_wp
        ux = 100.0_wp; uy = 0.0_wp
        call set_inactive_margins(ux,uy,f_ice,"infinite",a_front)
        call check("5 partial->free, a_front=0.5: open",ux(2,2) .eq. 100.0_wp,ux(2,2),100.0_wp)
        call check("5 free->free beyond front: closed",ux(3,2) .eq. 0.0_wp,ux(3,2),0.0_wp)

        a_front(3,:) = 0.5_wp*A_FRONT_MIN
        ux = 100.0_wp; uy = 0.0_wp
        call set_inactive_margins(ux,uy,f_ice,"infinite",a_front)
        call check("5 partial->free, a_front<min: closed",ux(2,2) .eq. 0.0_wp,ux(2,2),0.0_wp)
    end subroutine test_gate

    subroutine test_behind_front()
        ! Grounded ice in i <= 6, front cell 6 on a 400 m deeper bed, level-set
        ! front in the ice-free cell 7
        integer,  parameter :: nx = 12, ny = 4
        real(wp) :: H(nx,ny), lsf(nx,ny), a(nx,ny), z_bed(nx,ny), z_sl(nx,ny)
        real(wp) :: f_ice(nx,ny), H_eff(nx,ny)

        z_bed = -300.0_wp
        z_bed(6,:) = -700.0_wp
        z_sl  = 0.0_wp
        H     = 0.0_wp
        H(1:5,:) = 1000.0_wp
        H(6,:)   = 800.0_wp
        lsf   = 1.0_wp
        lsf(1:7,:) = -1.0_wp
        call calc_lsf_area_fraction(a,lsf,"infinite")

        call calc_ice_fraction(f_ice,H_eff,H,z_bed,z_sl,rho_ice,rho_sw,"marine",H_eff_min,dHdx,dx,"infinite",a)
        call check("6 behind level set: a_lsf",a(6,2) .eq. 1.0_wp,a(6,2),1.0_wp)
        call check("6 behind level set: full",f_ice(6,2) .eq. 1.0_wp,f_ice(6,2),1.0_wp)

        call calc_ice_fraction(f_ice,H_eff,H,z_bed,z_sl,rho_ice,rho_sw,"marine",H_eff_min,dHdx,dx,"infinite")
        call check("6 no level set: partial",f_ice(6,2) .lt. 1.0_wp,f_ice(6,2),1.0_wp)
    end subroutine test_behind_front

end program test_lsf_front
