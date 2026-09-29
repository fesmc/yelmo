program test_front_strain
    ! Strain rate at calving fronts for a uniformly translating slab.
    !
    ! A rigid translation has zero strain rate everywhere, including at the
    ! front. Faces with no adjacent ice hold zero velocity (as in the SSA
    ! solution), so a discretisation that differences or averages across
    ! them produces a spurious front strain rate of order u/dx.
    !
    ! Cases (dx = 8 km, uniform thickness, flat geometry):
    !   1. straight x-front, flow normal to the front     (U = 1000 m/yr)
    !   2. straight x-front, flow parallel to the front   (V =  500 m/yr)
    !   3. 45-degree staircase front, flow normal to it   (U = V = 707 m/yr)
    !   4. no front, linear shear u = a*y (exact dxy = a/2), interior accuracy
    !
    ! For cases 1-3 the effective strain rate de at ice cells should be ~0;
    ! the test reports its maximum over front cells for the 3D Jacobian path
    ! (calc_strain_rate_tensor_jac_quad3D, used for the calving stress) and
    ! the 2D path (calc_strain_rate_tensor_2D, same scheme as the DIVA
    ! viscosity). Case 4 checks that interior values stay exact.

    use yelmo_defs,  only : wp, jacobian_3D_class, strain_2D_class, strain_3D_class
    use deformation, only : calc_jacobian_vel_3D_uxyterms, calc_strain_rate_tensor_jac_quad3D, &
                            calc_strain_rate_tensor_2D

    implicit none

    integer,  parameter :: nx = 24, ny = 24, nz_aa = 5
    real(wp), parameter :: dx = 8000.0_wp, dy = 8000.0_wp
    real(wp), parameter :: de_max = 1.0e3_wp
    real(wp), parameter :: tol_rel = 1.0e-3_wp     ! max front de relative to u/dx
    character(len=*), parameter :: bnd = "periodic"

    real(wp) :: zeta_aa(nz_aa), zeta_ac(nz_aa+1)
    real(wp) :: f_ice(nx,ny), H_ice(nx,ny), f_grnd(nx,ny), zero2D(nx,ny)
    real(wp) :: ux(nx,ny,nz_aa), uy(nx,ny,nz_aa), uz(nx,ny,nz_aa+1)
    logical  :: is_front(nx,ny)
    type(jacobian_3D_class) :: jvel
    type(strain_3D_class)   :: strn
    type(strain_2D_class)   :: strn2D, strn2D_b

    real(wp) :: de3, de2, scale, err4
    integer  :: k, icase, n_fail
    logical  :: passed

    do k = 1, nz_aa
        zeta_aa(k) = real(k-1,wp)/real(nz_aa-1,wp)
    end do
    zeta_ac(1) = 0.0_wp
    do k = 2, nz_aa
        zeta_ac(k) = 0.5_wp*(zeta_aa(k-1)+zeta_aa(k))
    end do
    zeta_ac(nz_aa+1) = 1.0_wp

    call alloc_all()

    H_ice  = 500.0_wp
    f_grnd = 0.0_wp
    zero2D = 0.0_wp
    uz     = 0.0_wp

    write(*,*) "==============================================="
    write(*,*) " Front strain rate for a translating slab"
    write(*,*) "==============================================="
    write(*,'(a)') "  case                          max de 3D [1/yr]   max de 2D [1/yr]   u/dx [1/yr]"

    n_fail = 0

    do icase = 1, 3
        call set_case(icase,scale)
        call run_paths()
        de3 = maxval(strn2D%de,   mask=is_front)
        de2 = maxval(strn2D_b%de, mask=is_front)
        write(*,'(2x,a28,3es19.4)') case_name(icase), de3, de2, scale
        if (de3 .gt. tol_rel*scale .or. de2 .gt. tol_rel*scale) n_fail = n_fail + 1
    end do

    ! Case 4: no front, linear shear u = a*y, exact dxy = a/2
    call set_case(4,scale)
    call run_paths()
    err4 = max( maxval(abs(strn2D%dxy   - 0.5_wp*scale)), &
                maxval(abs(strn2D_b%dxy - 0.5_wp*scale)) )
    write(*,'(2x,a28,es19.4,a)') case_name(4), err4, "  (max |dxy - a/2|)"
    if (err4 .gt. tol_rel*scale) n_fail = n_fail + 1

    passed = (n_fail .eq. 0)
    write(*,*)
    if (passed) then
        write(*,*) " PASS"
    else
        write(*,'(a,i0,a)') "  FAIL (", n_fail, " cases)"
        stop 1
    end if

contains

    function case_name(icase) result(name)
        integer, intent(IN) :: icase
        character(len=28) :: name
        select case(icase)
            case(1); name = "1 x-front, normal flow"
            case(2); name = "2 x-front, parallel flow"
            case(3); name = "3 45-deg front, normal flow"
            case(4); name = "4 no front, linear shear"
        end select
    end function case_name

    subroutine set_case(icase,scale)
        ! Ice mask, face velocities (zero on faces with no adjacent ice)
        integer,  intent(IN)  :: icase
        real(wp), intent(OUT) :: scale
        integer  :: i, j, ip1, jp1
        real(wp) :: U, V, a

        f_ice = 0.0_wp
        do j = 1, ny
        do i = 1, nx
            select case(icase)
                case(1,2)
                    if (i .le. nx/2) f_ice(i,j) = 1.0_wp
                case(3)
                    if (i + j .le. nx) f_ice(i,j) = 1.0_wp
                case(4)
                    f_ice(i,j) = 1.0_wp
            end select
        end do
        end do

        select case(icase)
            case(1); U = 1000.0_wp; V = 0.0_wp
            case(2); U = 0.0_wp;    V = 500.0_wp
            case(3); U = 707.0_wp;  V = 707.0_wp
            case(4); U = 0.0_wp;    V = 0.0_wp
        end select
        a = 1.0e-2_wp       ! [1/yr] shear rate for case 4

        ux = 0.0_wp
        uy = 0.0_wp
        do j = 1, ny
        do i = 1, nx
            ip1 = i+1; if (ip1 .gt. nx) ip1 = 1
            jp1 = j+1; if (jp1 .gt. ny) jp1 = 1
            if (f_ice(i,j) .eq. 1.0_wp .or. f_ice(ip1,j) .eq. 1.0_wp) then
                if (icase .eq. 4) then
                    ux(i,j,:) = a*real(j,wp)*dy
                else
                    ux(i,j,:) = U
                end if
            end if
            if (f_ice(i,j) .eq. 1.0_wp .or. f_ice(i,jp1) .eq. 1.0_wp) uy(i,j,:) = V
        end do
        end do

        ! Front cells: ice cells with an ice-free edge neighbour, away from
        ! the periodic seam of the staircase case
        is_front = .FALSE.
        do j = 3, ny-2
        do i = 3, nx-2
            if (f_ice(i,j) .eq. 1.0_wp) then
                if (f_ice(i+1,j) .lt. 1.0_wp .or. f_ice(i-1,j) .lt. 1.0_wp .or. &
                    f_ice(i,j+1) .lt. 1.0_wp .or. f_ice(i,j-1) .lt. 1.0_wp) is_front(i,j) = .TRUE.
            end if
        end do
        end do
        if (icase .eq. 4) is_front = .TRUE.

        if (icase .eq. 4) then
            scale = a
        else
            scale = sqrt(U**2+V**2)/dx
        end if

    end subroutine set_case

    subroutine run_paths()
        ! 3D Jacobian path (calving stress) and 2D path (DIVA viscosity scheme)
        integer :: i, j

        jvel%dzx = 0.0_wp
        jvel%dzy = 0.0_wp
        jvel%dzz = 0.0_wp
        call calc_jacobian_vel_3D_uxyterms(jvel,ux,uy,uz,H_ice,f_ice,f_grnd,zero2D,zero2D, &
                                            zero2D,zero2D,zeta_aa,zeta_ac,dx,dy,bnd)
        call calc_strain_rate_tensor_jac_quad3D(strn,strn2D,jvel,H_ice,f_ice,f_grnd, &
                                            zeta_aa,zeta_ac,dx,dy,de_max,bnd)
        call calc_strain_rate_tensor_2D(strn2D_b,ux(:,:,1),uy(:,:,1),H_ice,f_ice,f_grnd,dx,dy,de_max,bnd)

        ! Case 4 reports dxy; exclude the periodic seam rows there
        do j = 1, ny
        do i = 1, nx
            if (j .le. 2 .or. j .ge. ny-1) then
                strn2D%dxy(i,j)   = 0.5_wp*1.0e-2_wp
                strn2D_b%dxy(i,j) = 0.5_wp*1.0e-2_wp
            end if
        end do
        end do

    end subroutine run_paths

    subroutine alloc_all()
        allocate(jvel%dxx(nx,ny,nz_aa),jvel%dxy(nx,ny,nz_aa),jvel%dxz(nx,ny,nz_aa))
        allocate(jvel%dyx(nx,ny,nz_aa),jvel%dyy(nx,ny,nz_aa),jvel%dyz(nx,ny,nz_aa))
        allocate(jvel%dzx(nx,ny,nz_aa+1),jvel%dzy(nx,ny,nz_aa+1),jvel%dzz(nx,ny,nz_aa+1))
        allocate(strn%dxx(nx,ny,nz_aa),strn%dyy(nx,ny,nz_aa),strn%dxy(nx,ny,nz_aa))
        allocate(strn%dxz(nx,ny,nz_aa),strn%dyz(nx,ny,nz_aa),strn%de(nx,ny,nz_aa))
        allocate(strn%div(nx,ny,nz_aa),strn%f_shear(nx,ny,nz_aa))
        call alloc_2D(strn2D)
        call alloc_2D(strn2D_b)
    end subroutine alloc_all

    subroutine alloc_2D(s)
        type(strain_2D_class), intent(INOUT) :: s
        allocate(s%dxx(nx,ny),s%dyy(nx,ny),s%dxy(nx,ny),s%dxz(nx,ny),s%dyz(nx,ny))
        allocate(s%de(nx,ny),s%div(nx,ny),s%f_shear(nx,ny),s%eps_eig_1(nx,ny),s%eps_eig_2(nx,ny))
        s%dxz = 0.0_wp
        s%dyz = 0.0_wp
        s%f_shear = 0.0_wp
    end subroutine alloc_2D

end program test_front_strain
