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
    use yelmo_tools,       only : boundary_code, get_neighbor_indices_bc_codes
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

        ! Extrapolate LSF velocities outside of the ice domain so that
        ! upwind advection near the front sees a non-zero front velocity.
        call extrapolate_ocn_acx(wx,wx,u_acx)
        call extrapolate_ocn_acy(wy,wy,v_acy)

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
        !   - Neighbour-snap: at each cell, if all four neighbours share
        !     sign with phi, snap phi to +-1. Gauss-Seidel pass in (i,j)
        !     order; in-place updates.
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

        integer :: i, j, nx, ny, im1, ip1, jm1, jp1, BC

        nx = size(lsf,1)
        ny = size(lsf,2)
        BC = boundary_code(boundaries)

        do j = 1, ny
        do i = 1, nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            if (lsf(i,j) .gt. 0.0_wp) then
                if ((lsf(im1,j) .gt. 0.0_wp) .and. (lsf(ip1,j) .gt. 0.0_wp) .and. &
                    (lsf(i,jm1) .gt. 0.0_wp) .and. (lsf(i,jp1) .gt. 0.0_wp)) then
                    lsf(i,j) =  1.0_wp
                end if
            else
                if ((lsf(im1,j) .le. 0.0_wp) .and. (lsf(ip1,j) .le. 0.0_wp) .and. &
                    (lsf(i,jm1) .le. 0.0_wp) .and. (lsf(i,jp1) .le. 0.0_wp)) then
                    lsf(i,j) = -1.0_wp
                end if
            end if
        end do
        end do

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

    subroutine extrapolate_ocn_acx(mask_fill,mask_orig,mask_ac)
        ! Fill ocean cells along the x-axis by nearest-filled-neighbour
        ! sweep. A cell is treated as "ocean" if mask_ac == 0 there.
        ! Single forward+backward pass per row: O(nx*ny) total work.

        implicit none

        real(wp), intent(INOUT) :: mask_fill(:,:)
        real(wp), intent(IN)    :: mask_orig(:,:)
        real(wp), intent(IN)    :: mask_ac(:,:)

        ! Local variables
        integer :: i, j, nx, ny
        logical, allocatable :: filled(:,:)

        nx = size(mask_orig,1)
        ny = size(mask_orig,2)
        allocate(filled(nx,ny))

        filled    = mask_ac .ne. 0.0_wp
        mask_fill = mask_orig

        if (.not. any(filled)) then
            ! No filled cells to extrapolate from
            deallocate(filled)
            return
        end if

        do j = 1, ny
            ! Forward sweep: rightward fill.
            do i = 2, nx
                if (.not. filled(i,j) .and. filled(i-1,j)) then
                    mask_fill(i,j) = mask_fill(i-1,j)
                    filled(i,j)    = .true.
                end if
            end do
            ! Backward sweep: leftward fill (catches any unfilled left tails).
            do i = nx-1, 1, -1
                if (.not. filled(i,j) .and. filled(i+1,j)) then
                    mask_fill(i,j) = mask_fill(i+1,j)
                    filled(i,j)    = .true.
                end if
            end do
        end do

        deallocate(filled)

        return

    end subroutine extrapolate_ocn_acx

    subroutine extrapolate_ocn_acy(mask_fill,mask_orig,mask_ac)
        ! Same as extrapolate_ocn_acx but along the y-axis.

        implicit none

        real(wp), intent(INOUT) :: mask_fill(:,:)
        real(wp), intent(IN)    :: mask_orig(:,:)
        real(wp), intent(IN)    :: mask_ac(:,:)

        ! Local variables
        integer :: i, j, nx, ny
        logical, allocatable :: filled(:,:)

        nx = size(mask_orig,1)
        ny = size(mask_orig,2)
        allocate(filled(nx,ny))

        filled    = mask_ac .ne. 0.0_wp
        mask_fill = mask_orig

        if (.not. any(filled)) then
            ! No filled cells to extrapolate from
            deallocate(filled)
            return
        end if

        do i = 1, nx
            ! Forward sweep: upward fill.
            do j = 2, ny
                if (.not. filled(i,j) .and. filled(i,j-1)) then
                    mask_fill(i,j) = mask_fill(i,j-1)
                    filled(i,j)    = .true.
                end if
            end do
            ! Backward sweep: downward fill.
            do j = ny-1, 1, -1
                if (.not. filled(i,j) .and. filled(i,j+1)) then
                    mask_fill(i,j) = mask_fill(i,j+1)
                    filled(i,j)    = .true.
                end if
            end do
        end do

        deallocate(filled)

        return

    end subroutine extrapolate_ocn_acy

end module lsf_module
