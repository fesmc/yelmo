module lsf_module
    ! ----------------------------------------------------------------------
    ! Level-set function (LSF) for flux-form calving.
    !
    ! Convention: lsf < 0 => ice domain, lsf > 0 => ocean. The zero
    ! level set is the calving front.
    !
    ! Public surface:
    !   - LSFinit         : initialise phi to +-1 from H_ice / z_bed / z_sl
    !   - LSFupdate       : advect phi at w = u_bar + cr (extrapolated into
    !                       the ocean) with a dedicated advective-form
    !                       upwind solver, then saturate to [-1, 1]
    !   - LSFredistance   : Sussman/Osher Hamilton-Jacobi redistancing to
    !                       restore |grad phi| ~= 1 without moving the
    !                       zero level set. Replaces the older ad-hoc
    !                       neighbour-snap and periodic ±1 re-flag.
    !
    ! Mirrors the design in Yelmo.jl/src/topo/lsf.jl.
    ! ----------------------------------------------------------------------

    use yelmo_defs,        only : sp, dp, wp, prec, TOL, TOL_UNDERFLOW, MISSING_VALUE, io_unit_err
    use yelmo_tools,       only : boundary_code, get_neighbor_indices_bc_codes, get_periodic_directions
    use topography,        only : calc_H_eff
    use, intrinsic :: iso_fortran_env, only : int64

    implicit none

    private

    ! === LSF routines ===
    public :: LSFinit
    public :: LSFupdate
    public :: LSFredistance
    public :: LSFsnap

    ! === Ocean extrapolation routines ===
    public :: extrapolate_ocn_acx
    public :: extrapolate_ocn_acy

contains
    ! ===================================================================
    !
    !                        LSF functions
    !
    ! ===================================================================

    subroutine LSFinit(LSF,H_ice,z_bed,z_sl,dx)

        implicit none

        real(wp), intent(OUT) :: LSF(:,:)       ! LSF mask
        real(wp), intent(IN)  :: H_ice(:,:)     ! Ice thickness
        real(wp), intent(IN)  :: z_bed(:,:)     ! Bedrock elevation
        real(wp), intent(IN)  :: z_sl(:,:)      ! Sea level
        real(wp), intent(IN)  :: dx             ! Model resolution

        ! Initialize LSF value at ocean value
        LSF = 1.0_wp

        ! Assign values
        where(H_ice .gt. 0.0_wp) LSF = -1.0_wp
        where(z_bed .gt. z_sl)   LSF = -1.0_wp

        return

    end subroutine LSFinit

    subroutine LSFupdate(dlsf,lsf,cr_acx,cr_acy,u_acx,v_acy,dx,dy,dt,boundaries)
        ! Advect the LSF with the front velocity w = u_bar + cr:
        !
        !   d phi / dt + w . grad phi = 0
        !
        ! i.e. the level-set equation in advective form (phi is not a
        ! conserved quantity, so no phi*div(w) term as in the flux-form
        ! ice-thickness solvers). w lives on the ac-nodes, where the
        ! calving rates are defined; see calc_lsf_advec_rate for the
        ! upwind discretisation. The explicit update is sub-cycled so that
        ! it is stable for any model timestep dt.

        implicit none

        real(wp),       intent(INOUT) :: dlsf(:,:)               ! [1/yr] LSF rate of change
        real(wp),       intent(INOUT) :: lsf(:,:)                ! LSF to be advected (aa-nodes)
        real(wp),       intent(INOUT) :: cr_acx(:,:),cr_acy(:,:) ! [m/yr] calving rate (vertical)
        real(wp),       intent(IN)    :: u_acx(:,:)              ! [m/a] 2D velocity, x-direction (ac-nodes)
        real(wp),       intent(IN)    :: v_acy(:,:)              ! [m/a] 2D velocity, y-direction (ac-nodes)
        real(wp),       intent(IN)    :: dx                      ! [m] Horizontal resolution, x-direction
        real(wp),       intent(IN)    :: dy                      ! [m] Horizontal resolution, y-direction
        real(wp),       intent(IN)    :: dt                      ! [a]   Timestep
        character(len=*), intent(IN)  :: boundaries              ! Boundary condition string (neighbour indices)

        ! Local variables
        integer  :: i, j, n, nx, ny, n_sub
        integer  :: im1, ip1, jm1, jp1
        integer  :: BC
        real(wp) :: w_max, rate_max, dt_max, dt_sub
        real(wp), allocatable :: wx(:,:), wy(:,:), lsf_n(:,:), lsf_dot(:,:)

        real(wp), parameter :: cfl = 0.5_wp                      ! CFL number on max|w|

        nx = size(lsf,1)
        ny = size(lsf,2)
        allocate(wx(nx,ny))
        allocate(wy(nx,ny))
        allocate(lsf_n(nx,ny))
        allocate(lsf_dot(nx,ny))

        BC = boundary_code(boundaries)

        dlsf = 0.0_wp  ! LSF change in a time dt

        ! Only advect if dt > 0
        if (dt .le. 0.0_wp) return

        ! Net LSF velocity: dynamic velocity + calving retreat rate
        wx = u_acx + cr_acx
        wy = v_acy + cr_acy

        ! Extrapolate the front velocity from the faces adjacent to ice
        ! into the ocean, so that upwind advection near the front sees it
        ! (also where u = 0, i.e. a stagnant front retreating at rate cr).
        call extrapolate_ocn_acx(wx,lsf,boundaries)
        call extrapolate_ocn_acy(wy,lsf,boundaries)

        ! Sub-step size: CFL = 0.5 on max|w| over the domain, and at most the
        ! positivity limit of the face-upwind update, dt*(a+b+c+d) <= 1, with
        ! a,b,c,d the inflow speeds/dx through the four cell faces (only
        ! binding where the flow converges on a cell from several sides).
        w_max    = max(maxval(abs(wx)),maxval(abs(wy)))
        rate_max = 0.0_wp
        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1) reduction(max:rate_max)
        do j = 1, ny
        do i = 1, nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            rate_max = max(rate_max, (max(wx(im1,j),0.0_wp) - min(wx(i,j),0.0_wp)) / dx &
                                   + (max(wy(i,jm1),0.0_wp) - min(wy(i,j),0.0_wp)) / dy)
        end do
        end do
        !$omp end parallel do

        dt_max = dt
        if (w_max    .gt. 0.0_wp) dt_max = min(dt_max, cfl*min(dx,dy)/w_max)
        if (rate_max .gt. 0.0_wp) dt_max = min(dt_max, 1.0_wp/rate_max)
        n_sub  = ceiling(dt/dt_max)
        dt_sub = dt / real(n_sub,wp)

        ! Advect with explicit sub-steps
        lsf_n = lsf
        do n = 1, n_sub
            call calc_lsf_advec_rate(lsf_dot,lsf,wx,wy,dx,dy,BC)
            lsf = lsf + dt_sub*lsf_dot
        end do

        ! Rate of change over the whole timestep (before saturation)
        dlsf = (lsf - lsf_n) / dt

        ! Saturate to [-1, 1] as a guardrail against upwind diffusion.
        ! LSFredistance is what restores |grad phi| ~= 1; this just keeps
        ! the field bounded.
        where(lsf .gt.  1.0_wp) lsf =  1.0_wp
        where(lsf .lt. -1.0_wp) lsf = -1.0_wp

        return

    end subroutine LSFupdate

    subroutine LSFredistance(lsf,dx,dy,n_iter,boundaries)
        ! Sussman/Osher Hamilton-Jacobi redistancing.
        !
        ! Solves   d phi / d tau + sgn(phi0) (|grad phi| - 1) = 0
        !
        ! discretised with the Godunov upwind scheme for |grad phi| and a
        ! smoothed sign function sgn(phi0) = phi0 / sqrt(phi0^2 + eps^2),
        ! eps = max(dx, dy). phi0 is the LSF at the start of redistancing
        ! and is held fixed for the duration of the iteration; this is
        ! what keeps the zero level set from drifting.
        !
        ! Pseudo-timestep dtau = 0.5 * min(dx, dy) satisfies the CFL limit
        ! for the explicit Godunov scheme. n_iter ~ 5 is enough for
        ! near-saturated input fields.
        !
        ! Port of lsf_redistance! in Yelmo.jl/src/topo/lsf.jl.

        implicit none

        real(wp),         intent(INOUT) :: lsf(:,:)
        real(wp),         intent(IN)    :: dx
        real(wp),         intent(IN)    :: dy
        integer,          intent(IN)    :: n_iter
        character(len=*), intent(IN)    :: boundaries

        ! Local variables
        integer  :: i, j, n, nx, ny
        integer  :: im1, ip1, jm1, jp1
        integer  :: BC
        real(wp) :: eps, dtau
        real(wp) :: a, b, c, d
        real(wp) :: s0, sgn
        real(wp) :: gx2, gy2, grad
        real(wp), allocatable :: phi0(:,:), new_phi(:,:)

        if (n_iter .le. 0) return

        nx = size(lsf,1)
        ny = size(lsf,2)
        allocate(phi0(nx,ny))
        allocate(new_phi(nx,ny))

        BC   = boundary_code(boundaries)
        eps  = max(dx,dy)
        dtau = 0.5_wp * min(dx,dy)

        ! Freeze the sign field at the start of redistancing. Using the
        ! evolving lsf here would let the zero level set drift.
        phi0 = lsf

        do n = 1, n_iter

            !$omp parallel do collapse(2) default(shared) &
            !$omp private(i,j,im1,ip1,jm1,jp1,a,b,c,d,s0,sgn,gx2,gy2,grad)
            do j = 1, ny
            do i = 1, nx

                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

                ! One-sided differences
                a = (lsf(i,j)   - lsf(im1,j)) / dx
                b = (lsf(ip1,j) - lsf(i,j)  ) / dx
                c = (lsf(i,j)   - lsf(i,jm1)) / dy
                d = (lsf(i,jp1) - lsf(i,j)  ) / dy

                ! Smoothed sign of phi0
                s0  = phi0(i,j)
                sgn = s0 / sqrt(s0*s0 + eps*eps)

                ! Godunov upwind |grad phi|^2
                if (sgn .gt. 0.0_wp) then
                    gx2 = max(max(a, 0.0_wp)**2, min(b, 0.0_wp)**2)
                    gy2 = max(max(c, 0.0_wp)**2, min(d, 0.0_wp)**2)
                else
                    gx2 = max(min(a, 0.0_wp)**2, max(b, 0.0_wp)**2)
                    gy2 = max(min(c, 0.0_wp)**2, max(d, 0.0_wp)**2)
                end if
                grad = sqrt(gx2 + gy2)

                new_phi(i,j) = lsf(i,j) - dtau * sgn * (grad - 1.0_wp)

            end do
            end do
            !$omp end parallel do

            lsf = new_phi

        end do

        deallocate(phi0)
        deallocate(new_phi)

        return

    end subroutine LSFredistance

    subroutine LSFsnap(lsf,time_now,dt_lsf,boundaries)
        ! Legacy LSF discipline (alternative to LSFredistance).
        !
        !   - Neighbour-snap: cells further than n_band edge steps from a
        !     sign change are snapped to +-1. The band of free cells around
        !     the front is n_band cells wide on each side. With one cell
        !     (the original rule) the cell ahead of the front is held at +1
        !     until the front cell changes sign, and the front moves at only
        !     0.87 w for small Courant numbers (0.92 at C=0.13); with two
        !     cells it moves at 0.96-0.98 w.
        !
        !   - Periodic full-field reflag: every dt_lsf years, reset phi to
        !     exact +-1 by sign. Disabled when dt_lsf <= 0.
        !
        ! This is the pre-#34 front-cleanup scheme. Kept selectable via
        ! ytopo_par%lsf_method = "snap" for runs that need the old
        ! behaviour or want to compare algorithms.

        implicit none

        real(wp),         intent(INOUT) :: lsf(:,:)
        real(wp),         intent(IN)    :: time_now      ! [yr] current model time
        real(wp),         intent(IN)    :: dt_lsf        ! [yr] reflag interval (<= 0 disables)
        character(len=*), intent(IN)    :: boundaries

        integer :: i, j, n, nx, ny, im1, ip1, jm1, jp1, BC
        logical :: is_pos
        logical, allocatable :: band(:,:), band_new(:,:)

        integer, parameter :: n_band = 2                 ! Free cells on each side of the front

        nx = size(lsf,1)
        ny = size(lsf,2)
        BC = boundary_code(boundaries)

        allocate(band(nx,ny),band_new(nx,ny))

        ! Front band: cells with an edge neighbour of opposite sign (lsf <= 0: ice side)
        do j = 1, ny
        do i = 1, nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            is_pos = lsf(i,j) .gt. 0.0_wp
            band(i,j) = ((lsf(im1,j) .gt. 0.0_wp) .neqv. is_pos) .or. ((lsf(ip1,j) .gt. 0.0_wp) .neqv. is_pos) .or. &
                        ((lsf(i,jm1) .gt. 0.0_wp) .neqv. is_pos) .or. ((lsf(i,jp1) .gt. 0.0_wp) .neqv. is_pos)
        end do
        end do

        ! Widen the band to n_band edge steps from the front
        do n = 2, n_band
            band_new = band
            do j = 1, ny
            do i = 1, nx
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
                if (band(im1,j) .or. band(ip1,j) .or. band(i,jm1) .or. band(i,jp1)) band_new(i,j) = .TRUE.
            end do
            end do
            band = band_new
        end do

        ! Snap the cells outside the band
        where (.not. band .and. lsf .gt. 0.0_wp) lsf =  1.0_wp
        where (.not. band .and. lsf .le. 0.0_wp) lsf = -1.0_wp

        if (dt_lsf .gt. 0.0_wp) then
            ! int64: a default integer overflows for |time_now| > ~2.1e7 yr
            if (mod(nint(time_now*100,int64),nint(dt_lsf*100,int64)) == 0) then
                where(lsf .gt. 0.0_wp) lsf =  1.0_wp
                where(lsf .le. 0.0_wp) lsf = -1.0_wp
            end if
        end if

        return

    end subroutine LSFsnap

    ! ===================================================================
    !
    ! Internal functions
    !
    ! ===================================================================

    subroutine calc_lsf_advec_rate(lsf_dot,lsf,wx,wy,dx,dy,BC)
        ! Rate of change of the LSF, lsf_dot = -w . grad(lsf), from
        ! first-order upwinding with the face (ac-node) velocities.
        ! In x, for cell i with faces i-1/2 (wx(im1,j)) and i+1/2 (wx(i,j)):
        !
        !   (w phi_x)_i = max(w_{i-1/2},0) * (phi_i     - phi_{i-1}) / dx
        !               + min(w_{i+1/2},0) * (phi_{i+1} - phi_i    ) / dx
        !
        ! i.e. each face that carries flow into the cell contributes the
        ! one-sided gradient on its side (and y likewise). This is the
        ! donor-cell flux form minus phi_i*div(w), so it reduces to the
        ! standard upwind scheme where w is uniform, and it uses the front
        ! velocity exactly on the face where the calving rate is defined.
        ! Stable and monotone for dt*(a+b+c+d) <= 1 (see LSFupdate).

        implicit none

        real(wp), intent(OUT) :: lsf_dot(:,:)
        real(wp), intent(IN)  :: lsf(:,:)
        real(wp), intent(IN)  :: wx(:,:)                        ! [m/yr] LSF velocity (acx-nodes)
        real(wp), intent(IN)  :: wy(:,:)                        ! [m/yr] LSF velocity (acy-nodes)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy
        integer,  intent(IN)  :: BC

        ! Local variables
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1

        nx = size(lsf,1)
        ny = size(lsf,2)

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
        do j = 1, ny
        do i = 1, nx

            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            lsf_dot(i,j) = -( max(wx(im1,j),0.0_wp) * (lsf(i,j)   - lsf(im1,j)) / dx &
                            + min(wx(i,j),  0.0_wp) * (lsf(ip1,j) - lsf(i,j)  ) / dx &
                            + max(wy(i,jm1),0.0_wp) * (lsf(i,j)   - lsf(i,jm1)) / dy &
                            + min(wy(i,j),  0.0_wp) * (lsf(i,jp1) - lsf(i,j)  ) / dy )

        end do
        end do
        !$omp end parallel do

        return

    end subroutine calc_lsf_advec_rate

    ! ===================================================================
    !
    ! Oceanic extrapolation routines.
    !
    ! ===================================================================

    subroutine extrapolate_ocn_acx(wx,lsf,boundaries)
        ! Extrapolate the LSF velocity on acx-nodes into the ocean along
        ! each row. Source faces are those adjacent to ice (lsf <= 0 in at
        ! least one of the two neighbouring cells, consistent with calving
        ! where lsf > 0), i.e. interior and front faces. They keep their
        ! value u + cr, so the retreat rate is extended from a stagnant
        ! (u = 0) front too. See extrapolate_ocn_1D for the fill rule.

        implicit none

        real(wp),         intent(INOUT) :: wx(:,:)
        real(wp),         intent(IN)    :: lsf(:,:)
        character(len=*), intent(IN)    :: boundaries

        ! Local variables
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1
        integer :: BC
        logical :: per_x, per_y
        logical, allocatable :: src(:)

        nx = size(wx,1)
        ny = size(wx,2)
        allocate(src(nx))

        BC = boundary_code(boundaries)
        call get_periodic_directions(per_x,per_y,BC)

        do j = 1, ny
            do i = 1, nx
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
                src(i) = (lsf(i,j) .le. 0.0_wp) .or. (lsf(ip1,j) .le. 0.0_wp)
            end do
            call extrapolate_ocn_1D(wx(:,j),src,per_x)
        end do

        deallocate(src)

        return

    end subroutine extrapolate_ocn_acx

    subroutine extrapolate_ocn_acy(wy,lsf,boundaries)
        ! Same as extrapolate_ocn_acx but on acy-nodes, along each column.

        implicit none

        real(wp),         intent(INOUT) :: wy(:,:)
        real(wp),         intent(IN)    :: lsf(:,:)
        character(len=*), intent(IN)    :: boundaries

        ! Local variables
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1
        integer :: BC
        logical :: per_x, per_y
        logical, allocatable :: src(:)

        nx = size(wy,1)
        ny = size(wy,2)
        allocate(src(ny))

        BC = boundary_code(boundaries)
        call get_periodic_directions(per_x,per_y,BC)

        do i = 1, nx
            do j = 1, ny
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
                src(j) = (lsf(i,j) .le. 0.0_wp) .or. (lsf(i,jp1) .le. 0.0_wp)
            end do
            call extrapolate_ocn_1D(wy(i,:),src,per_y)
        end do

        deallocate(src)

        return

    end subroutine extrapolate_ocn_acy

    subroutine extrapolate_ocn_1D(var,src,periodic)
        ! Fill the non-source points of a line with the value of the
        ! nearest source point, or with the mean of the two nearest ones
        ! (left and right) if they are equally far. Unlike a sequential
        ! sweep, the result does not depend on the sweep direction, so
        ! mirror symmetries of the domain are preserved. The line wraps
        ! around if periodic; otherwise points beyond the outermost source
        ! point take its value. A line without source points is unchanged.
        ! Two sweeps: O(n).

        implicit none

        real(wp), intent(INOUT) :: var(:)
        logical,  intent(IN)    :: src(:)
        logical,  intent(IN)    :: periodic

        ! Local variables
        integer :: k, m, n, k0, k_src, dl, dr
        integer, allocatable :: kl(:), kr(:)

        n = size(var)

        if (.not. any(src)) return

        allocate(kl(n))
        allocate(kr(n))

        ! Nearest source point to the left of each point (0: none).
        ! Sweep rightward; if periodic, start just after the last source
        ! point so that the wrap-around is included.
        k0    = 0
        k_src = 0
        if (periodic) k0    = findloc(src,.TRUE.,dim=1,back=.TRUE.)
        if (periodic) k_src = k0
        do m = 1, n
            k = modulo(k0+m-1,n) + 1
            if (src(k)) k_src = k
            kl(k) = k_src
        end do

        ! Nearest source point to the right of each point (0: none).
        ! Sweep leftward, likewise.
        k0    = n+1
        k_src = 0
        if (periodic) k0    = findloc(src,.TRUE.,dim=1)
        if (periodic) k_src = k0
        do m = 1, n
            k = modulo(k0-m-1,n) + 1
            if (src(k)) k_src = k
            kr(k) = k_src
        end do

        do k = 1, n
            if (src(k)) cycle
            if (kl(k) .eq. 0) then
                var(k) = var(kr(k))
            else if (kr(k) .eq. 0) then
                var(k) = var(kl(k))
            else
                dl = modulo(k-kl(k),n)
                dr = modulo(kr(k)-k,n)
                if (dl .lt. dr) then
                    var(k) = var(kl(k))
                else if (dr .lt. dl) then
                    var(k) = var(kr(k))
                else
                    var(k) = 0.5_wp*(var(kl(k))+var(kr(k)))
                end if
            end if
        end do

        deallocate(kl)
        deallocate(kr)

        return

    end subroutine extrapolate_ocn_1D

end module lsf_module
