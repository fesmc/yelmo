module yelmo_tools
    ! Generic functions and subroutines that could be used in many contexts:
    ! math, vectors, sorting, etc. 

    use yelmo_defs, only : sp, dp, wp, missing_value, TOL_UNDERFLOW, pi, &
                            io_unit_err

    use, intrinsic :: iso_fortran_env, only : int32, int64

    !$ use omp_lib
    
    implicit none 

    integer, parameter :: BND_ZEROS    = 0
    integer, parameter :: BND_INFINITE = 1
    integer, parameter :: BND_MISMIP3D = 2
    integer, parameter :: BND_TROUGH   = 3
    integer, parameter :: BND_PERIODIC = 4
    integer, parameter :: BND_PERIODIC_X = 5
    integer, parameter :: BND_PERIODIC_Y = 6
    
    interface is_finite
        module procedure is_finite_sp
        module procedure is_finite_dp
    end interface

    private 
    public :: is_finite
    public :: get_region_indices
    public :: get_neighbor_indices
    public :: get_neighbor_indices_bc_codes
    public :: get_periodic_directions
    public :: calc_magnitude 
    public :: calc_magnitude_from_staggered

    public :: calc_gradient_acx
    public :: calc_gradient_acy
    public :: calc_gradient_column_ac

    public :: mean_mask
    public :: minmax

    public :: set_boundaries_2D_aa
    public :: set_boundaries_3D_aa

    public :: fill_borders_2D
    public :: fill_borders_3D 
    
    public :: smooth_gauss_2D
    public :: smooth_gauss_3D
    public :: gauss_values

    public :: adjust_topography_gradients 

    ! Integration functions
    public :: test_integration
    public :: integrate_trapezoid1D_pt
    public :: integrate_trapezoid1D_1D
    public :: calc_vertical_integrated_2D
    public :: calc_vertical_integrated_3D
    
    ! Boundary constants (for converting string definitions to integers for faster computations)
    public :: boundary_code
    public :: BND_ZEROS, BND_INFINITE, BND_MISMIP3D, BND_TROUGH, BND_PERIODIC, BND_PERIODIC_X, BND_PERIODIC_Y

contains 

    elemental function is_finite_sp(x) result(ok)
        ! True unless x is NaN or Inf (all exponent bits set). An integer
        ! test on the bit pattern, so it is not folded away under -Ofast or
        ! -ffast-math, which may do that to ieee_is_nan and x /= x.
        real(sp), intent(IN) :: x
        logical :: ok
        integer(int32), parameter :: exp_mask = int(z'7F800000',int32)
        ok = iand(transfer(x,0_int32),exp_mask) .ne. exp_mask
    end function is_finite_sp

    elemental function is_finite_dp(x) result(ok)
        ! Double-precision version of is_finite_sp
        real(dp), intent(IN) :: x
        logical :: ok
        integer(int64), parameter :: exp_mask = int(z'7FF0000000000000',int64)
        ok = iand(transfer(x,0_int64),exp_mask) .ne. exp_mask
    end function is_finite_dp


    subroutine get_region_indices(i1,i2,j1,j2,nx,ny,irange,jrange)
        ! Get indices for a region based on bounds. 
        ! If no bounds provided, use whole domain.

        implicit none

        integer, intent(OUT) :: i1
        integer, intent(OUT) :: i2
        integer, intent(OUT) :: j1
        integer, intent(OUT) :: j2
        integer, intent(IN)  :: nx
        integer, intent(IN)  :: ny
        integer, intent(IN), optional :: irange(2)
        integer, intent(IN), optional :: jrange(2)

        
        if (present(irange)) then
            i1 = irange(1)
            i2 = irange(2)
        else
            i1 = 1
            i2 = nx
        end if

        if (present(jrange)) then
            j1 = jrange(1)
            j2 = jrange(2)
        else
            j1 = 1
            j2 = ny
        end if

        return

    end subroutine get_region_indices

    subroutine get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
        ! String-based wrapper of get_neighbor_indices_bc_codes, so that
        ! both interfaces share one definition of each boundary treatment.

        implicit none

        integer, intent(OUT) :: im1 
        integer, intent(OUT) :: ip1 
        integer, intent(OUT) :: jm1 
        integer, intent(OUT) :: jp1 
        integer, intent(IN)  :: i 
        integer, intent(IN)  :: j
        integer, intent(IN)  :: nx 
        integer, intent(IN)  :: ny
        
        character(len=*), intent(IN) :: boundaries

        call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,boundary_code(boundaries))

        return

    end subroutine get_neighbor_indices

    subroutine get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

        implicit none

        integer, intent(OUT) :: im1 
        integer, intent(OUT) :: ip1 
        integer, intent(OUT) :: jm1 
        integer, intent(OUT) :: jp1 
        integer, intent(IN)  :: i 
        integer, intent(IN)  :: j
        integer, intent(IN)  :: nx 
        integer, intent(IN)  :: ny
        integer, intent(IN)  :: BC

        select case(BC)

            case(BND_ZEROS,BND_INFINITE)
                ! Clamp neighbour indices to the domain edge so that ghost
                ! reads return the edge value (Neumann-zero / zero-gradient
                ! extension). Previously BND_ZEROS fell into the periodic
                ! DEFAULT branch, silently wrapping the field around the
                ! domain — see issue #34 follow-up.
                im1 = max(i-1,1)
                ip1 = min(i+1,nx)
                jm1 = max(j-1,1)
                jp1 = min(j+1,ny)

            case(BND_MISMIP3D,BND_TROUGH)
                im1 = max(i-1,1)
                ip1 = min(i+1,nx)
                jm1 = j-1
                if (jm1 .eq. 0)    jm1 = ny
                jp1 = j+1
                if (jp1 .eq. ny+1) jp1 = 1

            case(BND_PERIODIC)
                ! Periodic in x and y: true wrap with period nx (ny), i.e.,
                ! cells 1 and nx (1 and ny) are neighbors and every grid
                ! point is a regular interior point (no halo/ghost cells).

                im1 = i-1
                if (im1 .eq. 0)    im1 = nx
                ip1 = i+1
                if (ip1 .eq. nx+1) ip1 = 1

                jm1 = j-1
                if (jm1 .eq. 0)    jm1 = ny
                jp1 = j+1
                if (jp1 .eq. ny+1) jp1 = 1

            case(BND_PERIODIC_X)
                ! Periodic in x (true wrap, period nx),
                ! infinite (clamped) in y

                im1 = i-1
                if (im1 .eq. 0)    im1 = nx
                ip1 = i+1
                if (ip1 .eq. nx+1) ip1 = 1

                jm1 = max(j-1,1)
                jp1 = min(j+1,ny)

            case(BND_PERIODIC_Y)
                ! Periodic in y (true wrap, period ny),
                ! infinite (clamped) in x

                im1 = max(i-1,1)
                ip1 = min(i+1,nx)

                jm1 = j-1
                if (jm1 .eq. 0)    jm1 = ny
                jp1 = j+1
                if (jp1 .eq. ny+1) jp1 = 1

            case DEFAULT

                write(io_unit_err,*) "get_neighbor_indices_bc_codes:: Error: boundary code not recognized: ", BC
                error stop 1

        end select

        return

    end subroutine get_neighbor_indices_bc_codes

    subroutine get_periodic_directions(per_x,per_y,BC)
        ! Which directions wrap (true wrap, no halo cells) for a boundary
        ! code, consistent with get_neighbor_indices_bc_codes. In a periodic
        ! direction every point is an interior point, so loops must cover
        ! the full index range and no border values may be overwritten.

        implicit none

        logical, intent(OUT) :: per_x
        logical, intent(OUT) :: per_y
        integer, intent(IN)  :: BC

        select case(BC)

            case(BND_ZEROS,BND_INFINITE)
                per_x = .FALSE.
                per_y = .FALSE.

            case(BND_MISMIP3D,BND_TROUGH)
                per_x = .FALSE.
                per_y = .TRUE.

            case(BND_PERIODIC)
                per_x = .TRUE.
                per_y = .TRUE.

            case(BND_PERIODIC_X)
                per_x = .TRUE.
                per_y = .FALSE.

            case(BND_PERIODIC_Y)
                per_x = .FALSE.
                per_y = .TRUE.

            case DEFAULT

                write(io_unit_err,*) "get_periodic_directions:: Error: boundary code not recognized: ", BC
                error stop 1

        end select

        return

    end subroutine get_periodic_directions

    function boundary_code(boundaries) result(code)

        implicit none
        
        character(len=*), intent(in) :: boundaries
        integer :: code

        select case(trim(boundaries))
            case("zeros");      code = BND_ZEROS
            case("infinite");   code = BND_INFINITE
            case("MISMIP3D");   code = BND_MISMIP3D
            case("TROUGH");     code = BND_TROUGH
            case("periodic");   code = BND_PERIODIC
            case("periodic-x"); code = BND_PERIODIC_X
            case("periodic-y"); code = BND_PERIODIC_Y
            case("mask");       code = BND_INFINITE
            case default
                write(io_unit_err,*) "boundary_code:: Error: Boundary string not recognized: "//trim(boundaries)
                error stop 1
        end select

        return

    end function boundary_code

    elemental function calc_magnitude(u,v) result(umag)
        ! Get the vector magnitude from two components at the same grid location

        implicit none 

        real(wp), intent(IN)  :: u, v 
        real(wp) :: umag 

        umag = sqrt(u*u+v*v)

        return

    end function calc_magnitude
    
    function calc_magnitude_from_staggered(u,v,f_ice,boundaries) result(umag)
        ! Calculate the centered (aa-nodes) magnitude of a vector 
        ! from the staggered (ac-nodes) components

        implicit none 
        
        real(wp), intent(IN)  :: u(:,:), v(:,:)
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp) :: umag(size(u,1),size(u,2)) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1
        real(wp) :: unow, vnow 
        real(wp) :: f1, f2, H1, H2 
        integer  :: BC

        nx = size(u,1)
        ny = size(u,2) 

        umag = 0.0_wp 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            
            if (f_ice(i,j) .eq. 1.0) then 
                unow = 0.5*(u(im1,j)+u(i,j))
                vnow = 0.5*(v(i,jm1)+v(i,j))
                
                if (abs(unow) .lt. TOL_UNDERFLOW) unow = 0.0_wp 
                if (abs(vnow) .lt. TOL_UNDERFLOW) vnow = 0.0_wp 
            else 
                unow = 0.0 
                vnow = 0.0 
            end if 

            umag(i,j) = sqrt(unow*unow+vnow*vnow)

            if (abs(umag(i,j)) .lt. TOL_UNDERFLOW) umag(i,j) = 0.0_wp

        end do 
        end do 

        return

    end function calc_magnitude_from_staggered

    subroutine calc_gradient_acx(dvardx,var,f_ice,dx,grad_lim,zero_outside,boundaries,slope_bg)
        ! Calculate gradient on ac-nodes, accounting for ice margin if needed

        implicit none 

        real(wp), intent(OUT) :: dvardx(:,:) 
        real(wp), intent(IN)  :: var(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: grad_lim 
        logical,  intent(IN)  :: zero_outside 
        character(len=*), intent(IN) :: boundaries  ! Boundary conditions to apply 
        real(wp), intent(IN), optional :: slope_bg  ! Uniform background slope not contained in var
        
        ! Local variables 
        integer  :: i, j, nx, ny 
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: V0, V1 
        integer  :: BC

        nx = size(var,1)
        ny = size(var,2)

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,V0,V1)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            
            V0 = var(i,j) 
            V1 = var(ip1,j) 

            if (zero_outside) then 

                if (f_ice(i,j)   .lt. 1.0) V0 = 0.0 
                if (f_ice(ip1,j) .lt. 1.0) V1 = 0.0 
                
            end if 

            dvardx(i,j) = (V1-V0)/dx 

        end do 
        end do
        !$omp end parallel do

        ! Special case for infinite boundary conditions - ensure that slope 
        ! is the same, not the variable itself.
        select case(trim(boundaries))

            case("infinite","MISMIP3D","TROUGH","mask","periodic-y")
                dvardx(1,:)  = dvardx(2,:)
                dvardx(nx,:) = dvardx(nx-1,:)

        end select

        ! Add the background slope, so that the limit below bounds the total slope
        if (present(slope_bg)) dvardx = dvardx + slope_bg

        ! Finally, ensure that gradient is beneath desired limit 
        call minmax(dvardx,grad_lim)

        return 

    end subroutine calc_gradient_acx
    
subroutine calc_gradient_acy(dvardy,var,f_ice,dy,grad_lim,zero_outside,boundaries,slope_bg)
        ! Calculate gradient on ac-nodes, accounting for ice margin if needed

        implicit none 

        real(wp), intent(OUT) :: dvardy(:,:) 
        real(wp), intent(IN)  :: var(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: dy 
        real(wp), intent(IN)  :: grad_lim 
        logical,  intent(IN)  :: zero_outside 
        character(len=*), intent(IN) :: boundaries  ! Boundary conditions to apply 
        real(wp), intent(IN), optional :: slope_bg  ! Uniform background slope not contained in var
        
        ! Local variables 
        integer  :: i, j, nx, ny 
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: V0, V1 
        integer  :: BC

        nx = size(var,1)
        ny = size(var,2)

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,V0,V1)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            
            V0 = var(i,j) 
            V1 = var(i,jp1) 

            if (zero_outside) then 

                if (f_ice(i,j)   .lt. 1.0) V0 = 0.0 
                if (f_ice(i,jp1) .lt. 1.0) V1 = 0.0 
                
            end if 

            dvardy(i,j) = (V1-V0)/dy

        end do 
        end do
        !$omp end parallel do
        
        ! Special case for infinite boundary conditions - ensure that slope 
        ! is the same, not the variable itself.
        select case(trim(boundaries))

            case("infinite","mask","periodic-x")
                dvardy(:,1)  = dvardy(:,2)
                dvardy(:,ny) = dvardy(:,ny-1)

        end select

        ! Add the background slope, so that the limit below bounds the total slope
        if (present(slope_bg)) dvardy = dvardy + slope_bg

        ! Finally, ensure that gradient is beneath desired limit 
        call minmax(dvardy,grad_lim)

        return 

    end subroutine calc_gradient_acy

    subroutine calc_gradient_column_ac(dvdx_c,dvdy_c,dvdx,dvdy,f_ice,mask_ocn,boundaries)
        ! Gradient of the ice-column geometry (surface or base elevation) on
        ! ac-nodes, for the sigma-coordinate transform. A face between an
        ! ice-covered cell and ice-free ocean holds the jump to sea level (a
        ! cliff, the calving front), not a slope of the column: there the
        ! gradient of the adjacent face on the ice side is used, if that face
        ! lies between two ice-covered cells (otherwise zero). Faces to ice-free
        ! land keep their gradient (the margin slope). Faces without ice on
        ! either side are zero.

        implicit none 

        real(wp), intent(OUT) :: dvdx_c(:,:)        ! acx-nodes
        real(wp), intent(OUT) :: dvdy_c(:,:)        ! acy-nodes
        real(wp), intent(IN)  :: dvdx(:,:)          ! acx-nodes
        real(wp), intent(IN)  :: dvdy(:,:)          ! acy-nodes
        real(wp), intent(IN)  :: f_ice(:,:)         ! aa-nodes, ice-covered where f_ice == 1
        logical,  intent(IN)  :: mask_ocn(:,:)      ! aa-nodes, bed below sea level
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1
        integer :: im2, ip2, jm2, jp2
        integer :: BC
        logical :: ice_0, ice_1

        nx = size(f_ice,1)
        ny = size(f_ice,2)

        BC = boundary_code(boundaries)

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,im2,ip2,jm2,jp2,ice_0,ice_1)
        do j = 1, ny 
        do i = 1, nx 

            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! acx-node between (i,j) and (ip1,j)
            ice_0 = f_ice(i,j)   .eq. 1.0_wp
            ice_1 = f_ice(ip1,j) .eq. 1.0_wp

            dvdx_c(i,j) = 0.0_wp
            if (ice_0 .and. ice_1) then 
                dvdx_c(i,j) = dvdx(i,j)
            else if ( (ice_0 .and. .not. mask_ocn(ip1,j)) .or. &
                      (ice_1 .and. .not. mask_ocn(i,j)) ) then 
                ! Margin to ice-free land
                dvdx_c(i,j) = dvdx(i,j)
            else if (ice_0) then 
                ! Ice on the left: face between (im1,j) and (i,j)
                if (im1 .ne. i) then 
                    if (f_ice(im1,j) .eq. 1.0_wp) dvdx_c(i,j) = dvdx(im1,j)
                end if 
            else if (ice_1) then 
                ! Ice on the right: face between (ip1,j) and (ip2,j)
                call get_neighbor_indices_bc_codes(im2,ip2,jm2,jp2,ip1,j,nx,ny,BC)
                if (ip2 .ne. ip1) then 
                    if (f_ice(ip2,j) .eq. 1.0_wp) dvdx_c(i,j) = dvdx(ip1,j)
                end if 
            end if 

            ! acy-node between (i,j) and (i,jp1)
            ice_0 = f_ice(i,j)   .eq. 1.0_wp
            ice_1 = f_ice(i,jp1) .eq. 1.0_wp

            dvdy_c(i,j) = 0.0_wp
            if (ice_0 .and. ice_1) then 
                dvdy_c(i,j) = dvdy(i,j)
            else if ( (ice_0 .and. .not. mask_ocn(i,jp1)) .or. &
                      (ice_1 .and. .not. mask_ocn(i,j)) ) then 
                ! Margin to ice-free land
                dvdy_c(i,j) = dvdy(i,j)
            else if (ice_0) then 
                ! Ice below: face between (i,jm1) and (i,j)
                if (jm1 .ne. j) then 
                    if (f_ice(i,jm1) .eq. 1.0_wp) dvdy_c(i,j) = dvdy(i,jm1)
                end if 
            else if (ice_1) then 
                ! Ice above: face between (i,jp1) and (i,jp2)
                call get_neighbor_indices_bc_codes(im2,ip2,jm2,jp2,i,jp1,nx,ny,BC)
                if (jp2 .ne. jp1) then 
                    if (f_ice(i,jp2) .eq. 1.0_wp) dvdy_c(i,j) = dvdy(i,jp1)
                end if 
            end if 

        end do 
        end do 
        !$omp end parallel do

        return 

    end subroutine calc_gradient_column_ac
    
    function mean_mask(var,mask) result(ave)

        implicit none 

        real(wp), intent(IN) :: var(:,:) 
        logical,    intent(IN) :: mask(:,:) 
        real(wp) :: ave 
        integer :: n 

        n = count(mask)
        
        if (n .gt. 0) then 
            ave = sum(var,mask=mask) / real(n,wp)
        else 
            ave = 0.0 
        end if 

        return 

    end function mean_mask
    
    elemental subroutine minmax(var,var_lim)

        implicit none 

        real(wp), intent(INOUT) :: var 
        real(wp), intent(IN)    :: var_lim 

        if (var .lt. -var_lim) then 
            var = -var_lim 
        else if (var .gt. var_lim) then 
            var =  var_lim 
        end if 

        return 

    end subroutine minmax

    subroutine set_boundaries_2D_aa(var,boundaries,var_ref)

        implicit none 

        real(wp), intent(INOUT) :: var(:,:) 
        character(len=*), intent(IN) :: boundaries 
        real(wp), intent(IN), optional :: var_ref(:,:) 

        ! Local variables 
        integer :: nx, ny  

        nx = size(var,1) 
        ny = size(var,2) 

        select case(trim(boundaries))

            case("zeros","EISMINT")

                ! Set border values to zero
                var(1,:)  = 0.0
                var(nx,:) = 0.0

                var(:,1)  = 0.0
                var(:,ny) = 0.0

            case("periodic","periodic-xy") 

                ! Periodic x and y: true wrap (period nx, ny), all points
                ! are interior points, so there are no halo cells to set.

            case("periodic-x")

                ! Periodic x: true wrap (period nx), nothing to set.

                ! Infinite y (free-slip too)
                var(:,1)  = var(:,2)
                var(:,ny) = var(:,ny-1)

            case("periodic-y")

                ! Periodic y: true wrap (period ny), nothing to set.

                ! Infinite x (free-slip too)
                var(1,:)  = var(2,:)
                var(nx,:) = var(nx-1,:)

            case("MISMIP3D")

                ! === MISMIP3D =====
                var(1,:)    = var(2,:)          ! x=0, Symmetry 
                var(nx,:)   = 0.0               ! x=800km, no ice
                
!                var(:,1)    = var(:,2)          ! y=-50km, Free-slip condition
!                var(:,ny)   = var(:,ny-1)       ! y= 50km, Free-slip condition

            case("TROUGH")

                ! === MISMIP3D =====
                var(1,:)    = var(2,:)          ! x=0, Symmetry 
                var(nx,:)   = 0.0               ! x=800km, no ice
                
            case("infinite","mask")
                ! Set border points equal to inner neighbors

                call fill_borders_2D(var,nfill=1)

            case("fixed")
                ! Set border points equal to prescribed values from array

                call fill_borders_2D(var,nfill=1,fill=var_ref)

            case DEFAULT

                write(io_unit_err,*) "set_boundaries_2D_aa:: error: boundary method not recognized."
                write(io_unit_err,*) "boundaries = ", trim(boundaries)
                error stop 1

        end select 

        return

    end subroutine set_boundaries_2D_aa

    subroutine set_boundaries_3D_aa(var,boundaries,var_ref)

        implicit none 

        real(wp), intent(INOUT) :: var(:,:,:) 
        character(len=*), intent(IN) :: boundaries 
        real(wp), intent(IN), optional :: var_ref(:,:,:) 

        ! Local variables 
        integer :: nx, ny, nz  
        integer :: k 

        nx = size(var,1) 
        ny = size(var,2) 
        nz = size(var,3) 

        if (present(var_ref)) then

            do k = 1, nz 
                call set_boundaries_2D_aa(var(:,:,k),boundaries,var_ref(:,:,k))
            end do 
        
        else 
        
            do k = 1, nz 
                call set_boundaries_2D_aa(var(:,:,k),boundaries)
            end do 
        
        end if 
        
        return

    end subroutine set_boundaries_3D_aa

    subroutine fill_borders_2D(var,nfill,fill,fill_x,fill_y)

        implicit none

        real(wp), intent(INOUT) :: var(:,:)
        integer,    intent(IN)    :: nfill        ! How many neighbors to fill in
        real(wp), intent(IN), optional :: fill(:,:) ! Values to impose
        logical,  intent(IN), optional :: fill_x    ! Fill the x-borders? (default: true)
        logical,  intent(IN), optional :: fill_y    ! Fill the y-borders? (default: true)

        ! Local variables
        integer :: i, j, nx, ny, q
        logical :: do_x, do_y

        nx = size(var,1)
        ny = size(var,2)

        do_x = .TRUE.
        if (present(fill_x)) do_x = fill_x
        do_y = .TRUE.
        if (present(fill_y)) do_y = fill_y

        if (present(fill)) then
            ! Fill with prescribed values from array 'fill'

            do q = 1, nfill
                if (do_x) then
                    var(q,:)      = fill(nfill+1,:)
                    var(nx-q+1,:) = fill(nx-nfill,:)
                end if
                if (do_y) then
                    var(:,q)      = fill(:,nfill+1)
                    var(:,ny-q+1) = fill(:,ny-nfill)
                end if
            end do

        else
            ! Fill with interior neighbor values

            do q = 1, nfill
                if (do_x) then
                    var(q,:)      = var(nfill+1,:)
                    var(nx-q+1,:) = var(nx-nfill,:)
                end if
                if (do_y) then
                    var(:,q)      = var(:,nfill+1)
                    var(:,ny-q+1) = var(:,ny-nfill)
                end if
            end do

        end if

        return

    end subroutine fill_borders_2D

    subroutine fill_borders_3D(var,nfill,fill_x,fill_y)
        ! 3rd dimension is not filled (should be vertical dimension)

        implicit none

        real(wp), intent(INOUT) :: var(:,:,:)
        integer,    intent(IN)    :: nfill        ! How many neighbors to fill in
        logical,  intent(IN), optional :: fill_x    ! Fill the x-borders? (default: true)
        logical,  intent(IN), optional :: fill_y    ! Fill the y-borders? (default: true)

        ! Local variables
        integer :: i, j, nx, ny, q
        logical :: do_x, do_y

        nx = size(var,1)
        ny = size(var,2)

        do_x = .TRUE.
        if (present(fill_x)) do_x = fill_x
        do_y = .TRUE.
        if (present(fill_y)) do_y = fill_y

        do q = 1, nfill
            if (do_x) then
                var(q,:,:)      = var(nfill+1,:,:)
                var(nx-q+1,:,:) = var(nx-nfill,:,:)
            end if
            if (do_y) then
                var(:,q,:)      = var(:,nfill+1,:)
                var(:,ny-q+1,:) = var(:,ny-nfill,:)
            end if
        end do

        return

    end subroutine fill_borders_3D

    subroutine smooth_gauss_3D(var,dx,f_sigma,mask_apply,mask_use)

        ! Smooth out strain heating to avoid noise 

        implicit none

        real(wp),   intent(INOUT) :: var(:,:,:)      ! nx,ny,nz_aa: 3D variable
        real(wp),   intent(IN)    :: dx 
        real(wp),   intent(IN)    :: f_sigma  
        logical,    intent(IN), optional :: mask_apply(:,:) 
        logical,    intent(IN), optional :: mask_use(:,:) 

        ! Local variables
        integer    :: k, nz_aa  

        nz_aa = size(var,3)

        do k = 1, nz_aa 
             call smooth_gauss_2D(var(:,:,k),dx,f_sigma,mask_apply,mask_use)
        end do 

        return 

    end subroutine smooth_gauss_3D
    
    subroutine smooth_gauss_2D(var,dx,f_sigma,mask_apply,mask_use)
        ! Smooth out a field to avoid noise 
        ! mask_apply designates where smoothing should be applied 
        ! mask_use   designates which points can be considered in the smoothing filter 

        implicit none

        real(wp),   intent(INOUT) :: var(:,:)      ! [nx,ny] 2D variable
        real(wp),   intent(IN)    :: dx 
        real(wp),   intent(IN)    :: f_sigma  
        logical,    intent(IN), optional :: mask_apply(:,:) 
        logical,    intent(IN), optional :: mask_use(:,:) 

        ! Local variables
        integer  :: i, j, nx, ny, n, n2
        real(wp) :: sigma    
        real(wp), allocatable :: filter0(:,:), filter(:,:) 
        real(wp), allocatable :: var_old(:,:) 
        logical,  allocatable :: mask_apply_local(:,:) 
        logical,  allocatable :: mask_use_local(:,:)

        nx    = size(var,1)
        ny    = size(var,2)

        ! Safety check
        if (f_sigma .lt. 1.0_wp) then 
            write(io_unit_err,*) ""
            write(io_unit_err,*) "smooth_gauss_2D:: Error: f_sigma must be >= 1."
            write(io_unit_err,*) "f_sigma: ", f_sigma 
            write(io_unit_err,*) "dx:      ", dx 
            error stop 1
        end if 

        ! Get smoothing radius as standard devation of Gaussian function
        sigma = dx*f_sigma 

        ! Determine half-width of filter as 3-sigma
        n2 = 3*ceiling(f_sigma)

        ! Get total number of points for filter window in each direction
        n = 2*n2+1
        
        allocate(var_old(nx+2*n2,ny+2*n2))
        allocate(mask_apply_local(nx+2*n2,ny+2*n2))
        allocate(mask_use_local(nx+2*n2,ny+2*n2))
        allocate(filter0(n,n))
        allocate(filter(n,n))

        ! Check whether mask_apply is available 
        if (present(mask_apply)) then 
            ! use mask_use to define neighborhood points
            
            mask_apply_local = .FALSE. 
            mask_apply_local(n2+1:n2+nx,n2+1:n2+ny) = mask_apply 

        else
            ! Assume that everywhere should be smoothed

            mask_apply_local = .FALSE. 
            mask_apply_local(n2+1:n2+nx,n2+1:n2+ny) = .TRUE.
        
        end if

        ! Check whether mask_use is available 
        if (present(mask_use)) then 
            ! use mask_use to define neighborhood points
            
            mask_use_local = .TRUE.
            mask_use_local(n2+1:n2+nx,n2+1:n2+ny) = mask_use 

        else
            ! Assume that mask_apply also gives the points to use for smoothing 

            mask_use_local = mask_apply_local
        
        end if

        ! Calculate default 2D Gaussian smoothing kernel
        filter0 = gauss_values(dx,dx,sigma=sigma,n=n)

        var_old = 0.0
        var_old(n2+1:n2+nx,n2+1:n2+ny) = var

        ! Fill edge halos by mirror reflection about each border
        var_old(1:n2,n2+1:n2+ny)            = var(n2:1:-1,:)
        var_old(n2+nx+1:nx+2*n2,n2+1:n2+ny) = var(nx:nx-n2+1:-1,:)
        var_old(n2+1:n2+nx,1:n2)            = var(:,n2:1:-1)
        var_old(n2+1:n2+nx,n2+ny+1:ny+2*n2) = var(:,ny:ny-n2+1:-1)

        ! Fill corner halos by double mirror reflection
        var_old(1:n2,1:n2)                       = var(n2:1:-1,n2:1:-1)
        var_old(n2+nx+1:nx+2*n2,1:n2)            = var(nx:nx-n2+1:-1,n2:1:-1)
        var_old(1:n2,n2+ny+1:ny+2*n2)            = var(n2:1:-1,ny:ny-n2+1:-1)
        var_old(n2+nx+1:nx+2*n2,n2+ny+1:ny+2*n2) = var(nx:nx-n2+1:-1,ny:ny-n2+1:-1)
        
        !$omp parallel do collapse(2) private(i,j,filter)
        do j = n2+1, n2+ny 
        do i = n2+1, n2+nx 

            if (mask_apply_local(i,j)) then 
                ! Apply smoothing to this point 

                ! Limit filter input to neighbors of interest
                filter = filter0 
                where(.not. mask_use_local(i-n2:i+n2,j-n2:j+n2)) filter = 0.0

                ! If neighbors are available, normalize and perform smoothing  
                if (sum(filter) .gt. 0.0) then 
                    filter = filter/sum(filter)
                    var(i-n2,j-n2) = sum(var_old(i-n2:i+n2,j-n2:j+n2)*filter) 
                end if  

            end if 

        end do 
        end do 
        !$omp end parallel do

        return 

    end subroutine smooth_gauss_2D

    function gauss_values(dx,dy,sigma,n) result(filt)
        ! Calculate 2D Gaussian smoothing kernel
        ! https://en.wikipedia.org/wiki/Gaussian_blur

        implicit none 

        real(wp), intent(IN) :: dx 
        real(wp), intent(IN) :: dy 
        real(wp), intent(IN) :: sigma 
        integer,    intent(IN) :: n 
        real(wp) :: filt(n,n) 

        ! Local variables 
        real(wp) :: x, y  
        integer    :: n2, i, j, i1, j1  

        if (mod(n,2) .ne. 1) then 
            write(*,*) "gauss_values:: error: n can only be odd."
            write(*,*) "n = ", n 
        end if 

        n2 = (n-1)/2 

        do j = -n2, n2 
        do i = -n2, n2 
            x = i*dx 
            y = j*dy 

            i1 = i+1+n2 
            j1 = j+1+n2 
            filt(i1,j1) = 1.0/(2.0*pi*sigma**2)*exp(-(x**2+y**2)/(2*sigma**2))

        end do 
        end do 
        
        ! Normalize to ensure sum to 1
        filt = filt / sum(filt)

        return 

    end function gauss_values

    ! ================================================================================
    !
    ! Regularizing/smoothing functions 
    !
    ! ================================================================================

    subroutine adjust_topography_gradients(z_bed,H_ice,grad_lim,dx,boundaries)
        ! Smooth the bedrock topography and corresponding ice thickness,
        ! so that specified limit on gradients is not exceeded. Only
        ! apply smoothing directly in places where gradient is too large.

        implicit none 

        real(wp), intent(INOUT) :: z_bed(:,:) 
        real(wp), intent(INOUT) :: H_ice(:,:) 
        real(wp), intent(IN)    :: grad_lim 
        real(wp), intent(IN)    :: dx 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, q, nx, ny 
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: dy
        real(wp), allocatable :: dzbdx(:,:)
        real(wp), allocatable :: dzbdy(:,:)
        real(wp), allocatable :: f_ice(:,:)
        logical,  allocatable :: mask_apply(:,:) 
        logical,  allocatable :: mask_use(:,:) 
        integer  :: BC

        integer, parameter :: iter_max = 50 

        nx = size(z_bed,1)
        ny = size(z_bed,2) 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        dy = dx 

        allocate(dzbdx(nx,ny))
        allocate(dzbdy(nx,ny))
        allocate(f_ice(nx,ny))
        allocate(mask_apply(nx,ny))
        allocate(mask_use(nx,ny))

        ! Smooth z_bed in specific locations if gradients are exceeded.

        mask_use = .TRUE. 
        f_ice    = 1.0 

        do q = 1, iter_max

            ! Calculate bedrock gradients (f_ice and grad_lim are not used)
            call calc_gradient_acx(dzbdx,z_bed,f_ice,dx,grad_lim=100.0_wp, &
                                        zero_outside=.FALSE.,boundaries=boundaries)
            call calc_gradient_acy(dzbdy,z_bed,f_ice,dy,grad_lim=100.0_wp, &
                                        zero_outside=.FALSE.,boundaries=boundaries)

            ! Determine where gradients are too large
            mask_apply = .FALSE.
            do j = 3, ny-3
            do i = 3, nx-3
                
                ! Get neighbor indices
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

                if (abs(dzbdx(i,j)) .ge. grad_lim) then 
                    mask_apply(i,j)   = .TRUE. 
                    mask_apply(ip1,j) = .TRUE. 
                end if

                if (abs(dzbdy(i,j)) .ge. grad_lim) then 
                    mask_apply(i,j)   = .TRUE. 
                    mask_apply(i,jp1) = .TRUE. 
                end if
            end do 
            end do 

            write(*,*) "z_bed smoothing: ", q, count(mask_apply),  &
                                        maxval(abs(dzbdx(3:nx-3,3:ny-3))), &
                                        maxval(abs(dzbdy(3:nx-3,3:ny-3)))

            if (count(mask_apply) .eq. 0) exit 

            ! Smooth z_bed at desired locations, and H_ice so that H_ice avoids spurious patterns
            call smooth_gauss_2D(z_bed,dx=dx,f_sigma=2.0_wp,mask_apply=mask_apply,mask_use=mask_use)
            call smooth_gauss_2D(H_ice,dx=dx,f_sigma=2.0_wp,mask_apply=mask_apply,mask_use=mask_use)
            
        end do

        return

    end subroutine adjust_topography_gradients

    ! === Generic integration functions ============

    function calc_vertical_integrated_3D(var,zeta) result(var_int)
        ! Vertically integrate a field 3D field (nx,ny,nz)
        ! layer by layer (in the z-direction), return a 3D array

        implicit none

        real(wp), intent(IN) :: var(:,:,:)
        real(wp), intent(IN) :: zeta(:)
        real(wp) :: var_int(size(var,1),size(var,2),size(var,3))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(var,1)
        ny = size(var,2)

        !$omp parallel do collapse(2) private(i,j)
        do j = 1, ny
        do i = 1, nx
            var_int(i,j,:) = integrate_trapezoid1D_1D(var(i,j,:),zeta)
        end do
        end do
        !$omp end parallel do

        return

    end function calc_vertical_integrated_3D

    function calc_vertical_integrated_2D(var,zeta) result(var_int)
        ! Vertically integrate a field 3D field (nx,ny,nz) 
        ! to the surface, return a 2D array (nx,ny)
        
        implicit none

        real(wp), intent(IN) :: var(:,:,:)
        real(wp), intent(IN) :: zeta(:)
        real(wp) :: var_int(size(var,1),size(var,2))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(var,1)
        ny = size(var,2)

        !$omp parallel do collapse(2) private(i,j)
        do j = 1, ny
        do i = 1, nx
            var_int(i,j) = integrate_trapezoid1D_pt(var(i,j,:),zeta)
        end do
        end do
        !$omp end parallel do 

        return

    end function calc_vertical_integrated_2D
    
    function integrate_trapezoid1D_pt(var,zeta) result(var_int)
        ! Integrate a variable from the base to height zeta(nk) in the ice column.
        ! The value of the integral using the trapezium rule can be found using
        ! integral = (b - a)*((f(a) +f(b))/2 + Σ_1_n-1(f(k)) )/n 
        ! Returns a point of integrated value of var at level zeta(nk).

        implicit none

        real(wp), intent(IN) :: var(:)
        real(wp), intent(IN) :: zeta(:)
        real(wp) :: var_int

        ! Local variables 
        integer :: k, nk
        real(wp) :: var_mid 
        
        nk = size(var,1)

        ! Initial value is zero
        var_int = 0.0_wp 

        ! Intermediate values include sum of all previous values 
        ! Take current value as average between points
        do k = 2, nk
            var_mid = 0.5_wp*(var(k)+var(k-1))
            var_int = var_int + var_mid*(zeta(k) - zeta(k-1))
        end do

        return

    end function integrate_trapezoid1D_pt

    function integrate_trapezoid1D_1D(var,zeta) result(var_int)
        ! Integrate a variable from the base to each layer zeta of the ice column.
        ! Note this is designed assuming indices 1 = base, nk = surface 
        ! The value of the integral using the trapezium rule can be found using
        ! integral = (b - a)*((f(a) +f(b))/2 + Σ_1_n-1(f(k)) )/n 
        ! Returns a 1D array with integrated value at each level 

        implicit none

        real(wp), intent(IN) :: var(:)
        real(wp), intent(IN) :: zeta(:)
        real(wp) :: var_int(size(var,1))

        ! Local variables 
        integer    :: k, nk
        real(wp) :: var_mid 

        nk = size(var,1)

        ! Initial value is zero
        var_int(1:nk) = 0.0_wp 

        ! Intermediate values include sum of all previous values 
        ! Take current value as average between points
        do k = 2, nk
            var_mid = 0.5_wp*(var(k)+var(k-1))
            var_int(k:nk) = var_int(k:nk) + var_mid*(zeta(k) - zeta(k-1))
        end do
        
        return

    end function integrate_trapezoid1D_1D

    subroutine simpne(x,y,result)
        !*****************************************************************************80
        !
        !! SIMPNE approximates the integral of unevenly spaced data.
        !
        !  Discussion:
        !
        !    The routine repeatedly interpolates a 3-point Lagrangian polynomial 
        !    to the data and integrates that exactly.
        !
        !  Modified:
        !
        !    10 February 2006
        !
        !  Reference:
        !
        !    Philip Davis, Philip Rabinowitz,
        !    Methods of Numerical Integration,
        !    Second Edition,
        !    Dover, 2007,
        !    ISBN: 0486453391,
        !    LC: QA299.3.D28.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) NTAB, number of data points.  
        !    NTAB must be at least 3.
        !
        !    Input, real ( kind = 8 ) X(NTAB), contains the X values of the data,
        !    in order.
        !
        !    Input, real ( kind = 8 ) Y(NTAB), contains the Y values of the data.
        !
        !    Output, real ( kind = 8 ) RESULT.
        !    RESULT is the approximate value of the integral.
        
        implicit none

        real(wp) :: x(:)
        real(wp) :: y(:)
        real(wp) :: result

        integer :: ntab

        real(wp) :: del(3)
        real(wp) :: e
        real(wp) :: f
        real(wp) :: feints
        real(wp) :: g(3)
        integer    :: i
        integer    :: n
        real(wp) :: pi(3)
        real(wp) :: sum1

        real(wp) :: x1
        real(wp) :: x2
        real(wp) :: x3

        ntab = size(x,1) 

        result = 0.0D+00

        if ( ntab <= 2 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SIMPNE - Fatal error!'
            write ( *, '(a)' ) '  NTAB <= 2.'
            error stop 1
        end if
     
        n = 1
     
        do
     
            x1 = x(n)
            x2 = x(n+1)
            x3 = x(n+2)
            e = x3 * x3- x1 * x1
            f = x3 * x3 * x3 - x1 * x1 * x1
            feints = x3 - x1

            del(1) = x3 - x2
            del(2) = x1 - x3
            del(3) = x2 - x1

            g(1) = x2 + x3
            g(2) = x1 + x3
            g(3) = x1 + x2

            pi(1) = x2 * x3
            pi(2) = x1 * x3
            pi(3) = x1 * x2
     
            sum1 = 0.0D+00
            do i = 1, 3
                sum1 = sum1 + y(n-1+i) * del(i) &
                    * ( f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints )
            end do
            result = result - sum1 / ( del(1) * del(2) * del(3) )
     
            n = n + 2

            if ( ntab <= n + 1 ) then
            exit
            end if

        end do
     
        if ( mod ( ntab, 2 ) /= 0 ) then
            return
        end if

        n = ntab - 2
        x3 = x(ntab)
        x2 = x(ntab-1)
        x1 = x(ntab-2)
        e = x3 * x3 - x2 * x2
        f = x3 * x3 * x3 - x2 * x2 * x2
        feints = x3 - x2

        del(1) = x3 - x2
        del(2) = x1 - x3
        del(3) = x2 - x1

        g(1) = x2 + x3
        g(2) = x1 + x3
        g(3) = x1 + x2

        pi(1) = x2 * x3
        pi(2) = x1 * x3
        pi(3) = x1 * x2
     
        sum1 = 0.0D+00
        do i = 1, 3
            sum1 = sum1 + y(n-1+i) * del(i) * &
                ( f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints )
        end do
     
        result = result - sum1 / ( del(1) * del(2) * del(3) )
     
        return

    end subroutine simpne 

    subroutine test_integration()

        implicit none 

        ! Local variables
        integer :: i, j, n, k, t 
        integer :: nn(11) 
        real(wp), allocatable :: zeta0(:) 
        real(wp), allocatable :: zeta(:)
        real(wp), allocatable :: var0(:) 
        real(wp), allocatable :: var(:) 
        real(wp), allocatable :: var_ints(:)
        real(wp) :: var_int 
        real(wp) :: var_int_00

        write(*,*) "=== test_integration ======"
        
        nn = [11,21,31,41,51,61,71,81,91,101,1001]

        do k = 1, size(nn)

            n = nn(k)

            allocate(zeta0(n))
            allocate(zeta(n))
            allocate(var0(n))
            allocate(var(n))
            allocate(var_ints(n))

            do i = 1, n 
                zeta0(i) = real(i-1)/real(n-1)
!                 var0(i) = real(i-1)
                var0(i)  = (n-1)-real(i-1)
            end do 

            ! Linear zeta 
            zeta = zeta0
            var = var0 

            ! Non-linear zeta 
            zeta = zeta0*zeta0 
            do i = 1, n 
                do j = 1, n 
                    if (zeta0(j) .ge. zeta(i)) exit 
                end do 

                if (zeta0(j) .eq. zeta(i)) then 
                    var(i) = var0(j) 
                else 
                    var(i) = var0(j-1) + (var0(j)-var0(j-1))*(zeta(i)-zeta0(j-1))/(zeta0(j)-zeta0(j-1))
                end if 
            end do 

!             do i = 1, n 
!                 write(*,*) zeta0(i), var0(i), zeta(i), var(i) 
!             end do 
!             stop 
            
            ! Analytical average value 
!             var_int_00 = real(n-1)/2.0

            do t = 1, 10000
                
                ! Test trapezoid1D solver 
                var_int  = integrate_trapezoid1D_pt(var,zeta)

                ! Determine "analytical" value from simpson approximation solver 
                call simpne(zeta,var,var_int_00)
            end do 

            ! Test trapezoid 1D_1D solver, check last value for full average over column
!             var_ints = integrate_trapezoid1D_1D(var,zeta)
!             var_int  = var_ints(n) 

            write(*,*) "mean (0:",n ,") = ", var_int_00, var_int, var_int-var_int_00, 100.0*(var_int-var_int_00)/var_int_00 

            deallocate(zeta0,var0,zeta,var,var_ints)

        end do 

        stop 

        return 

    end subroutine test_integration
    
end module yelmo_tools
