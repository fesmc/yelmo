module solver_ssa_ac

    use yelmo_defs, only : sp, dp, wp, io_unit_err, TOL, TOL_UNDERFLOW, is_equal, &
                           MASK_FRNT_FLOAT, MASK_FRNT_MARINE, MASK_FRNT_GRND, MASK_FRNT_ICE_FREE_LAND
    use yelmo_tools, only : boundary_code, get_neighbor_indices_bc_codes

    use solver_linear
    use ncio        ! For diagnostic outputting only 

    implicit none 

    private 
    public :: set_ssa_masks
    public :: ssa_diagnostics_write_init
    public :: ssa_diagnostics_write_step

    ! Routines that make use of the linear_solver_class object defined in the module solver_linear.F90:
    public :: linear_solver_save_velocity
    public :: linear_solver_matrix_ssa_ac_csr_2D

    ! Helper used by alternative SSA assemblers (e.g. solver_ssa_ac_energy):
    public :: stagger_visc_aa_ab

contains
    
    subroutine linear_solver_save_velocity(ux,uy,lgs,ulim)
        ! Extract velocity solution from lgs object. 

        implicit none 

        real(wp), intent(OUT) :: ux(:,:)                ! [m yr^-1] Horizontal velocity x
        real(wp), intent(OUT) :: uy(:,:)                ! [m yr^-1] Horizontal velocity y
        type(linear_solver_class), intent(IN) :: lgs 
        real(wp), intent(IN)    :: ulim 

        ! Local variables 
        integer :: i, j, n, nr 

        do n = 1, lgs%nmax-1, 2

            i = lgs%n2i((n+1)/2)
            j = lgs%n2j((n+1)/2)

            nr = n
            ux(i,j) = lgs%x_value(nr)

            nr = n+1
            uy(i,j) = lgs%x_value(nr)

        end do

        ! Limit the velocity generally =====================
        call limit_vel(ux,ulim)
        call limit_vel(uy,ulim)

        return

    end subroutine linear_solver_save_velocity

    subroutine linear_solver_matrix_ssa_ac_csr_2D(lgs,ux,uy,beta_acx,beta_acy, &
                            N_aa,ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice,taud_acx, &
                            taud_acy,taul_int_acx,taul_int_acy,dx,dy,beta_min,boundaries,lateral_bc)
        ! Define sparse matrices A*x=b in format 'compressed sparse row' (csr)
        ! for the SSA momentum balance equations with velocity components
        ! ux and uy defined on ac-nodes (right and top borders of i,j grid cell)
        ! Store sparse matrices in linear_solver_class object 'lgs' for later use.

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs
        real(wp), intent(IN) :: ux(:,:)                 ! [m yr^-1] Horizontal velocity x (acx-nodes)
        real(wp), intent(IN) :: uy(:,:)                 ! [m yr^-1] Horizontal velocity y (acy-nodes)
        real(wp), intent(IN) :: beta_acx(:,:)           ! [Pa yr m^-1] Basal friction (acx-nodes)
        real(wp), intent(IN) :: beta_acy(:,:)           ! [Pa yr m^-1] Basal friction (acy-nodes)
        real(wp), intent(IN) :: N_aa(:,:)               ! [Pa yr m] Vertically integrated viscosity (aa-nodes)
        integer,  intent(IN) :: ssa_mask_acx(:,:)       ! [--] Mask to determine ssa solver actions (acx-nodes)
        integer,  intent(IN) :: ssa_mask_acy(:,:)       ! [--] Mask to determine ssa solver actions (acy-nodes)
        integer,  intent(IN) :: mask_frnt(:,:)          ! [--] Ice-front mask 
        real(wp), intent(IN) :: H_ice(:,:)              ! [m]  Ice thickness (aa-nodes)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: taud_acx(:,:)           ! [Pa] Driving stress (acx nodes)
        real(wp), intent(IN) :: taud_acy(:,:)           ! [Pa] Driving stress (acy nodes)
        real(wp), intent(IN) :: taul_int_acx(:,:)       ! [Pa m] Vertically integrated lateral stress (acx nodes)
        real(wp), intent(IN) :: taul_int_acy(:,:)       ! [Pa m] Vertically integrated lateral stress (acy nodes) 
        real(wp), intent(IN) :: dx, dy
        real(wp), intent(IN) :: beta_min                ! [Pa yr m^-1] Minimum allowed basal friction for grounded ice

        character(len=*), intent(IN) :: boundaries 
        character(len=*), intent(IN) :: lateral_bc

        ! Local variables
        integer  :: nx, ny
        integer  :: i, j, k, n, m 
        integer  :: nc, nr
        real(wp) :: inv_dx, inv_dxdx 
        real(wp) :: inv_dy, inv_dydy 
        real(wp) :: inv_dxdy, inv_2dxdy, inv_4dxdy 
        real(wp) :: beta_now 
        
        real(wp), allocatable :: N_ab(:,:)

        ! Boundary conditions (bcs) counterclockwise unit circle 
        ! 1: x, right-border
        ! 2: y, upper-border 
        ! 3: x, left--border 
        ! 4: y, lower-border 
        character(len=56) :: bcs(4)
        logical :: bc_per(4), bc_free(4)          ! bcs(k) is "periodic" / "free-slip"

        integer :: im1, ip1, jm1, jp1 
        real(wp) :: N_aa_now
        integer  :: n_grnd_x, n_grnd_y, n_beta_x, n_beta_y

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Safety check for initialization
        if (.not. allocated(lgs%x_value)) then 
            ! Object 'lgs' has not been initialized yet, do so now.

            call linear_solver_init(lgs,nx,ny,nvar=2,n_terms=9)

        end if 

        ! Define border conditions (only choices are: no-slip, free-slip, periodic)
        select case(trim(boundaries)) 

            case("MISMIP3D")

                bcs(1) = "free-slip"
                bcs(2) = "periodic"
                bcs(3) = "no-slip"
                bcs(4) = "periodic" 

                ! Note: the following BCs might be more correct for
                ! the experimental design, but they can lead more
                ! easily to instabilities. Keeping the y-axis periodic
                ! helps to avoid that, and the results should be 
                ! effectively the same.
                ! bcs(1) = "free-slip"
                ! bcs(2) = "free-slip"
                ! bcs(3) = "no-slip"
                ! bcs(4) = "free-slip" 

            case("TROUGH")

                bcs(1) = "free-slip"
                bcs(2) = "periodic"
                bcs(3) = "no-slip"      ! "free-slip"?
                bcs(4) = "periodic" 

            case("periodic")

                bcs(1:4) = "periodic" 

            case("periodic-x")

                bcs(1) = "periodic"
                bcs(2) = "free-slip"
                bcs(3) = "periodic"
                bcs(4) = "free-slip"

            case("periodic-y")

                bcs(1) = "free-slip"
                bcs(2) = "periodic"
                bcs(3) = "free-slip"
                bcs(4) = "periodic"
            
            case("infinite","mask")

                bcs(1:4) = "free-slip"

            case("zeros")

                bcs(1:4) = "no-slip" 
                
            case DEFAULT 

                bcs(1:4) = "no-slip" 
                
        end select 

        ! Evaluate the border types once, not per cell in the assembly loops
        bc_per  = bcs .eq. "periodic"
        bc_free = bcs .eq. "free-slip"

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        
        allocate(N_ab(nx,ny))

        ! Define some factors

        inv_dx      = 1.0_wp / dx 
        inv_dxdx    = 1.0_wp / (dx*dx)
        inv_dy      = 1.0_wp / dy 
        inv_dydy    = 1.0_wp / (dy*dy)
        inv_dxdy    = 1.0_wp / (dx*dy)
        inv_2dxdy   = 1.0_wp / (2.0_wp*dx*dy)
        inv_4dxdy   = 1.0_wp / (4.0_wp*dx*dy)
        

        ! Calculate the staggered depth-integrated viscosity 
        ! at the grid-cell corners (ab-nodes). 
        call stagger_visc_aa_ab(N_ab,N_aa,H_ice,f_ice,boundaries)
        

        !-------- Assembly of the system of linear equations
        !             (matrix storage: compressed sparse row CSR) --------

        lgs%a_ptr(1) = 1

        k = 0

        ! Counters of inner grounded rows (and those with beta > 0) for the
        ! beta consistency check after assembly
        n_grnd_x = 0
        n_beta_x = 0
        n_grnd_y = 0
        n_beta_y = 0

        do n=1, lgs%nmax-1, 2

            i = lgs%n2i((n+1)/2)
            j = lgs%n2j((n+1)/2)

            ! Get neighbor indices assuming periodic domain
            ! (all other boundary conditions are treated with special cases below)
            im1 = i-1
            if (im1 .eq. 0)    im1 = nx 
            ip1 = i+1
            if (ip1 .eq. nx+1) ip1 = 1 

            jm1 = j-1
            if (jm1 .eq. 0)    jm1 = ny
            jp1 = j+1
            if (jp1 .eq. ny+1) jp1 = 1 
                
            ! ------ Equations for ux ---------------------------

            nr = n   ! row counter

            ! == Treat special cases first ==

            if (ssa_mask_acx(i,j) .eq. 0) then
                ! SSA set to zero velocity for this point

                k = k+1
                lgs%a_value(k) = 1.0_wp   ! diagonal element only
                lgs%a_index(k) = nr

                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = 0.0_wp
            
            else if (ssa_mask_acx(i,j) .eq. -1) then 
                ! Assign prescribed boundary velocity to this point
                ! (eg for prescribed velocity corresponding to 
                ! analytical grounding line flux, or for a regional domain)

                k = k+1
                lgs%a_value(k)  = 1.0   ! diagonal element only
                lgs%a_index(k)  = nr

                lgs%b_value(nr) = ux(i,j)
                lgs%x_value(nr) = ux(i,j)
            
            else if (i .eq. 1 .and. .not. bc_per(3)) then 
                ! Left boundary 

                if (bc_free(3)) then 
                
                    nc = 2*lgs%ij2n(i,j)-1          ! column counter for ux(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(ip1,j)-1        ! column counter for ux(ip1,j)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp
                    lgs%a_index(k) = nc

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = ux(i,j)

                else ! bcs(3) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0   ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0
                    lgs%x_value(nr) = 0.0

                end if 

            else if (i .eq. nx .and. .not. bc_per(1)) then 
                ! Right boundary 
                
                if (bc_free(1)) then 
                
                    nc = 2*lgs%ij2n(i,j)-1          ! column counter for ux(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(nx-1,j)-1       ! column counter for ux(nx-1,j)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp
                    lgs%a_index(k) = nc

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = ux(i,j)

                else ! bcs(3) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0   ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0
                    lgs%x_value(nr) = 0.0

                end if 

            else if (j .eq. 1 .and. .not. bc_per(4)) then 
                ! Lower boundary 

                if (bc_free(4)) then 

                    nc = 2*lgs%ij2n(i,j)-1          ! column counter for ux(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,jp1)-1       ! column counter for ux(i,jp1)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp
                    lgs%a_index(k) = nc

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = ux(i,j)


                else ! bcs(4) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp   ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

            else if (j .eq. ny .and. .not. bc_per(2)) then 
                ! Upper boundary 

                if (bc_free(2)) then 
                    
                    nc = 2*lgs%ij2n(i,j)-1          ! column counter for ux(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,ny-1)-1       ! column counter for ux(i,ny-1)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp
                    lgs%a_index(k) = nc

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = ux(i,j)

                else ! bcs(2) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp   ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

            else if (ssa_mask_acx(i,j) .eq. 3) then 
                ! Lateral boundary condition should be applied here 

                if (is_equal(f_ice(i,j),1.0_wp) .and. f_ice(ip1,j) .lt. 1.0) then 
                    ! === Case 1: ice-free to the right ===

                    N_aa_now = N_aa(i,j)
                    
                    nc = 2*lgs%ij2n(im1,j)-1
                        ! smallest nc (column counter), for ux(im1,j)
                    k = k+1
                    lgs%a_value(k) = -4.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc 

                    nc = 2*lgs%ij2n(i,jm1)
                        ! next nc (column counter), for uy(i,jm1)
                    k = k+1
                    lgs%a_value(k) = -2.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,j)-1
                        ! next nc (column counter), for ux(i,j)
                    k = k+1
                    lgs%a_value(k) = 4.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,j)
                        ! next nc (column counter), for uy(i,j)
                    k = k+1
                    lgs%a_value(k) = 2.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    ! Assign matrix values
                    lgs%b_value(nr) = taul_int_acx(i,j) 
                    lgs%x_value(nr) = ux(i,j)
                    
                else 
                    ! Case 2: ice-free to the left
                    
                    N_aa_now = N_aa(ip1,j)
                    
                    nc = 2*lgs%ij2n(i,j)-1
                        ! next nc (column counter), for ux(i,j)
                    k = k+1
                    lgs%a_value(k) = -4.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(ip1,jm1)
                        ! next nc (column counter), for uy(ip1,jm1)
                    k  = k+1
                    lgs%a_value(k) = -2.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(ip1,j)-1
                        ! next nc (column counter), for ux(ip1,j)
                    k = k+1
                    lgs%a_value(k) = 4.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(ip1,j)
                        ! largest nc (column counter), for uy(ip1,j)
                    k  = k+1
                    lgs%a_value(k) = 2.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    ! Assign matrix values
                    lgs%b_value(nr) = taul_int_acx(i,j) 
                    lgs%x_value(nr) = ux(i,j)
                
                end if 

            else
                ! === Inner SSA solution === 

                beta_now = beta_acx(i,j)
                if (ssa_mask_acx(i,j) .eq. 1 .and. beta_acx(i,j) .eq. 0.0) beta_now = beta_min

                if (ssa_mask_acx(i,j) .eq. 1) then
                    n_grnd_x = n_grnd_x + 1
                    if (beta_acx(i,j) .gt. 0.0) n_beta_x = n_beta_x + 1
                end if

                ! -- vx terms -- 

                nc = 2*lgs%ij2n(i,j)-1          ! column counter for ux(i,j)
                k = k+1
                lgs%a_value(k) = -4.0_wp*inv_dxdx*(N_aa(ip1,j)+N_aa(i,j)) &
                                 -1.0_wp*inv_dydy*(N_ab(i,j)+N_ab(i,jm1)) &
                                 -beta_now
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(ip1,j)-1        ! column counter for ux(ip1,j)
                k = k+1
                lgs%a_value(k) =  4.0_wp*inv_dxdx*N_aa(ip1,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(im1,j)-1        ! column counter for ux(im1,j)
                k = k+1
                lgs%a_value(k) =  4.0_wp*inv_dxdx*N_aa(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jp1)-1        ! column counter for ux(i,jp1)
                k = k+1
                lgs%a_value(k) =  1.0_wp*inv_dydy*N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jm1)-1        ! column counter for ux(i,jm1)
                k = k+1
                lgs%a_value(k) =  1.0_wp*inv_dydy*N_ab(i,jm1)
                lgs%a_index(k) = nc

                ! -- vy terms -- 
                
                nc = 2*lgs%ij2n(i,j)            ! column counter for uy(i,j)
                k = k+1
                lgs%a_value(k) = -2.0_wp*inv_dxdy*N_aa(i,j)     &
                                 -1.0_wp*inv_dxdy*N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(ip1,j)          ! column counter for uy(ip1,j)
                k = k+1
                lgs%a_value(k) =  2.0_wp*inv_dxdy*N_aa(ip1,j)   &
                                 +1.0_wp*inv_dxdy*N_ab(i,j)
                lgs%a_index(k) = nc
                
                nc = 2*lgs%ij2n(ip1,jm1)        ! column counter for uy(ip1,jm1)
                k = k+1
                lgs%a_value(k) = -2.0_wp*inv_dxdy*N_aa(ip1,j)   &
                                 -1.0_wp*inv_dxdy*N_ab(i,jm1)
                lgs%a_index(k) = nc
                
                nc = 2*lgs%ij2n(i,jm1)          ! column counter for uy(i,jm1)
                k = k+1
                lgs%a_value(k) =  2.0_wp*inv_dxdy*N_aa(i,j)   &
                                 +1.0_wp*inv_dxdy*N_ab(i,jm1)
                lgs%a_index(k) = nc
                

                lgs%b_value(nr) = taud_acx(i,j)
                lgs%x_value(nr) = ux(i,j)

            end if

            lgs%a_ptr(nr+1) = k+1   ! row is completed, store index to next row
            
            ! ------ Equations for uy ---------------------------

            nr = n+1   ! row counter

            ! == Treat special cases first ==
            
            if (ssa_mask_acy(i,j) .eq. 0) then
                ! SSA not active here, velocity set to zero

                k = k+1
                lgs%a_value(k)  = 1.0_wp   ! diagonal element only
                lgs%a_index(k)  = nr

                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = 0.0_wp

            else if (ssa_mask_acy(i,j) .eq. -1) then 
                ! Assign prescribed boundary velocity to this point
                ! (eg for prescribed velocity corresponding to analytical grounding line flux)

                k = k+1
                lgs%a_value(k)  = 1.0   ! diagonal element only
                lgs%a_index(k)  = nr

                lgs%b_value(nr) = uy(i,j)
                lgs%x_value(nr) = uy(i,j)
            
            else if (j .eq. 1 .and. .not. bc_per(4)) then 
                ! lower boundary 

                if (bc_free(4)) then 

                    nc = 2*lgs%ij2n(i,j)            ! column counter for uy(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,jp1)            ! column counter for uy(i,jp1)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp
                    lgs%a_index(k) = nc

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = uy(i,j)

                else ! bcs(4) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp        ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

            else if (j .eq. ny .and. .not. bc_per(2)) then 
                ! Upper boundary 

                if (bc_free(2)) then 

                    nc = 2*lgs%ij2n(i,j)            ! column counter for uy(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,ny-1)         ! column counter for uy(i,ny-1)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp
                    lgs%a_index(k) = nc

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = uy(i,j)

                else ! bcs(2) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp        ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

            else if (i .eq. 1 .and. .not. bc_per(3)) then 
                ! Left boundary 

                if (bc_free(3)) then 

                    nc = 2*lgs%ij2n(i,j)            ! column counter for uy(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(ip1,j)          ! column counter for uy(ip1,j)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp
                    lgs%a_index(k) = nc

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = uy(i,j)

                else ! bcs(2) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp        ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

            else if (i .eq. nx .and. .not. bc_per(1)) then 
                ! Right boundary 

                if (bc_free(1)) then 

                    nc = 2*lgs%ij2n(i,j)            ! column counter for uy(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(nx-1,j)         ! column counter for uy(nx-1,j)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp
                    lgs%a_index(k) = nc

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = uy(i,j)

                else ! bcs(2) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp        ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

            else if (ssa_mask_acy(i,j) .eq. 3) then 
                ! Lateral boundary condition should be applied here 

                if (is_equal(f_ice(i,j),1.0) .and. f_ice(i,jp1) .lt. 1.0) then 
                    ! === Case 1: ice-free to the top ===

                    N_aa_now = N_aa(i,j)

                    nc = 2*lgs%ij2n(im1,j)-1
                        ! smallest nc (column counter), for ux(im1,j)
                    k = k+1
                    lgs%a_value(k) = -2.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,jm1)
                        ! next nc (column counter), for uy(i,jm1)
                    k = k+1
                    lgs%a_value(k) = -4.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,j)-1
                        ! next nc (column counter), for ux(i,j)
                    k = k+1
                    lgs%a_value(k) = 2.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,j)
                        ! next nc (column counter), for uy(i,j)
                    k = k+1
                    lgs%a_value(k) = 4.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    ! Assign matrix values
                    lgs%b_value(nr) = taul_int_acy(i,j)
                    lgs%x_value(nr) = uy(i,j)
                    
                else
                    ! === Case 2: ice-free to the bottom ===
                    
                    N_aa_now = N_aa(i,jp1)

                    nc = 2*lgs%ij2n(im1,jp1)-1
                        ! next nc (column counter), for ux(im1,jp1)
                    k = k+1
                    lgs%a_value(k) = -2.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,j)
                        ! next nc (column counter), for uy(i,j)
                    k = k+1
                    lgs%a_value(k) = -4.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,jp1)-1
                        ! next nc (column counter), for ux(i,jp1)
                    k = k+1
                    lgs%a_value(k) = 2.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,jp1)
                        ! next nc (column counter), for uy(i,jp1)
                    k = k+1
                    lgs%a_value(k) = 4.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    ! Assign matrix values
                    lgs%b_value(nr) = taul_int_acy(i,j)
                    lgs%x_value(nr) = uy(i,j)
         
                end if 

            else
                ! === Inner SSA solution === 

                beta_now = beta_acy(i,j)
                if (ssa_mask_acy(i,j) .eq. 1 .and. beta_acy(i,j) .eq. 0.0) beta_now = beta_min

                if (ssa_mask_acy(i,j) .eq. 1) then
                    n_grnd_y = n_grnd_y + 1
                    if (beta_acy(i,j) .gt. 0.0) n_beta_y = n_beta_y + 1
                end if

                ! -- vy terms -- 

                nc = 2*lgs%ij2n(i,j)        ! column counter for uy(i,j)
                k = k+1
                lgs%a_value(k) = -4.0_wp*inv_dydy*(N_aa(i,jp1)+N_aa(i,j))   &
                                 -1.0_wp*inv_dxdx*(N_ab(i,j)+N_ab(im1,j))   &
                                 -beta_now
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jp1)      ! column counter for uy(i,jp1)
                k = k+1
                lgs%a_value(k) =  4.0_wp*inv_dydy*N_aa(i,jp1)
                lgs%a_index(k) = nc 

                nc = 2*lgs%ij2n(i,jm1)      ! column counter for uy(i,jm1)
                k = k+1
                lgs%a_value(k) =  4.0_wp*inv_dydy*N_aa(i,j)
                lgs%a_index(k) = nc
                
                nc = 2*lgs%ij2n(ip1,j)      ! column counter for uy(ip1,j)
                k = k+1
                lgs%a_value(k) =  1.0_wp*inv_dxdx*N_ab(i,j)
                lgs%a_index(k) = nc
                
                nc = 2*lgs%ij2n(im1,j)      ! column counter for uy(im1,j)
                k = k+1
                lgs%a_value(k) =  1.0_wp*inv_dxdx*N_ab(im1,j)
                lgs%a_index(k) = nc
                
                ! -- vx terms -- 

                nc = 2*lgs%ij2n(i,j)-1      ! column counter for ux(i,j)
                k = k+1
                lgs%a_value(k) = -2.0_wp*inv_dxdy*N_aa(i,j)     &
                                 -1.0_wp*inv_dxdy*N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jp1)-1    ! column counter for ux(i,jp1)
                k = k+1
                lgs%a_value(k) =  2.0_wp*inv_dxdy*N_aa(i,jp1)     &
                                 +1.0_wp*inv_dxdy*N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(im1,jp1)-1  ! column counter for ux(im1,jp1)
                k = k+1
                lgs%a_value(k) = -2.0_wp*inv_dxdy*N_aa(i,jp1)     &
                                 -1.0_wp*inv_dxdy*N_ab(im1,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(im1,j)-1  ! column counter for ux(im1,j)
                k = k+1
                lgs%a_value(k) =  2.0_wp*inv_dxdy*N_aa(i,j)     &
                                 +1.0_wp*inv_dxdy*N_ab(im1,j)
                lgs%a_index(k) = nc

                lgs%b_value(nr) = taud_acy(i,j)
                lgs%x_value(nr) = uy(i,j)

            end if

            lgs%a_ptr(nr+1) = k+1   ! row is completed, store index to next row

        end do

        ! Consistency check: ensure beta is defined well for grounded ice.
        ! Only inner rows (momentum equations with a friction term) are counted;
        ! border and lateral-bc rows do not use beta.
        if ( (n_grnd_x .gt. 0 .and. n_beta_x .eq. 0) .or. &
             (n_grnd_y .gt. 0 .and. n_beta_y .eq. 0) ) then
            ! No inner grounded points found with a non-zero beta,
            ! something was not well-defined/well-initialized, give a warning
            ! with some statistics. In the assembly above, beta=beta_min
            ! was used for these points.

            write(*,*)
            write(*,"(a)") "linear_solver_matrix_ssa_ac_csr_2D:: Warning: beta appears to be zero everywhere for grounded ice."
            write(*,*) "inner grounded acx rows: ", n_grnd_x, ", with beta_acx > 0: ", n_beta_x
            write(*,*) "inner grounded acy rows: ", n_grnd_y, ", with beta_acy > 0: ", n_beta_y
            write(*,*)

        end if

        ! Done: A, x and b matrices in Ax=b have been populated 
        ! and stored in lgs object. 

        return

    end subroutine linear_solver_matrix_ssa_ac_csr_2D

    subroutine check_base_slope(is_steep,zb0,zb1,dx,lim)

        logical,  intent(OUT) :: is_steep
        real(wp), intent(IN) :: zb0         ! [m]
        real(wp), intent(IN) :: zb1         ! [m]
        real(wp), intent(IN) :: dx          ! [m]
        real(wp), intent(IN) :: lim         ! [dx/dx] = [unitless]

        if ( abs(zb1-zb0) / dx .gt. lim ) then 
            is_steep = .TRUE. 
        else 
            is_steep = .FALSE. 
        end if 

        return

    end subroutine check_base_slope


    subroutine set_ssa_masks(ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice, &
                                        f_grnd,z_base,z_sl,dx,use_ssa,lateral_bc,boundaries)
        ! Define where ssa calculations should be performed
        ! Note: could be binary, but perhaps also distinguish 
        ! grounding line/zone to use this mask for later gl flux corrections
        ! mask = -1: no ssa calculated, velocity imposed/unchanged
        ! mask = 0: no ssa calculated, velocity set to zero
        ! mask = 1: shelfy-stream ssa calculated 
        ! mask = 2: shelf ssa calculated 
        ! mask = 3: ssa lateral boundary condition applied
        ! mask = 4: ssa lateral boundary, but treated as inner ssa

        ! Note: the parameter gradbase_max is used to check slope of ice base. 
        ! If at a given point, it is greater than this limit, the ssa solver
        ! will be disabled in this direction. gradbase_max=0.1 is a relatively
        ! high value, but is reached for points next to deep troughs in Antarctica,
        ! and next to some fjords in Greenland. Steeper slopes are present
        ! in higher-resolution topographies typically.
        
        implicit none 
        
        integer,  intent(OUT) :: ssa_mask_acx(:,:) 
        integer,  intent(OUT) :: ssa_mask_acy(:,:)
        integer,  intent(IN)  :: mask_frnt(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: z_base(:,:)
        real(wp), intent(IN)  :: z_sl(:,:)
        real(wp), intent(IN)  :: dx 
        logical,  intent(IN)  :: use_ssa       ! SSA is actually active now? 
        character(len=*), intent(IN) :: lateral_bc 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1
        integer  :: BC
        real(wp) :: H_acx, H_acy
        logical  :: is_steep 
        logical  :: is_convergent 
        
        integer  :: mask_lat

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        
        ! Set boundary condition code
        BC = boundary_code(boundaries)

        select case(trim(lateral_bc))
            case("none","floating","float","marine","all")
                ! ok
            case DEFAULT
                write(io_unit_err,*) "set_ssa_masks:: error: ssa_lat_bc parameter value not recognized."
                write(io_unit_err,*) "ydyn.ssa_lat_bc = ", lateral_bc
                stop 
        end select

        ! Initially no active ssa points, all velocities set to zero
        ssa_mask_acx = 0
        ssa_mask_acy = 0
        
        if (use_ssa) then 

            ! Step 2: define ssa solver masks

            do j = 1, ny
            do i = 1, nx

                ! Get neighbor indices
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)


                ! == x-direction ===

                if (is_equal(f_ice(i,j),1.0_wp) .or. is_equal(f_ice(ip1,j),1.0_wp)) then
                
                    ! Current ac-node is border of an ice covered cell in x-direction
                    
                    if (f_grnd(i,j) .gt. 0.0_wp .or. f_grnd(ip1,j) .gt. 0.0_wp) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_acx(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_acx(i,j) = 2
                    end if 

                    ! SPECIAL CASE: floating ice next to ice-free land,
                    ! then set ssa mask to zero (ie, set velocity to zero)
                    if (ssa_mask_acx(i,j) .eq. 2) then 

                        if ( is_equal(f_grnd(i,j),0.0_wp) .and. &
                               (f_grnd(ip1,j) .gt. 0.0 .and. H_ice(ip1,j) .eq. 0.0) ) then 

                            ssa_mask_acx(i,j) = 0

                        else if ( (f_grnd(i,j) .gt. 0.0 .and. H_ice(i,j) .eq. 0.0) .and. &
                                    f_grnd(ip1,j) .eq. 0.0 ) then

                            ssa_mask_acx(i,j) = 0 

                        end if 

                    end if 

                end if
                
                ! Overwrite above if this face is an ice front (lateral bc, or deactivated)
                mask_lat = front_face_mask(mask_frnt(i,j),mask_frnt(ip1,j),lateral_bc)
                if (mask_lat .gt. 0) ssa_mask_acx(i,j) = mask_lat

                ! == y-direction ===

                if (is_equal(f_ice(i,j),1.0_wp) .or. is_equal(f_ice(i,jp1),1.0_wp)) then
                
                    ! Current ac-node is border of an ice covered cell in x-direction
                    
                    if (f_grnd(i,j) .gt. 0.0 .or. f_grnd(i,jp1) .gt. 0.0) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_acy(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_acy(i,j) = 2
                    end if 

                    ! SPECIAL CASE: floating ice next to ice-free land,
                    ! then set ssa mask to zero (ie, set velocity to zero)
                    if (ssa_mask_acy(i,j) .eq. 2) then 

                        if ( f_grnd(i,j) .eq. 0.0 .and. &
                               (f_grnd(i,jp1) .gt. 0.0 .and. H_ice(i,jp1) .eq. 0.0) ) then 

                            ssa_mask_acy(i,j) = 0

                        else if ( (f_grnd(i,j) .gt. 0.0 .and. H_ice(i,j) .eq. 0.0) .and. &
                                    f_grnd(i,jp1) .eq. 0.0 ) then

                            ssa_mask_acy(i,j) = 0 

                        end if 

                    end if 

                end if

                ! Overwrite above if this face is an ice front (lateral bc, or deactivated)
                mask_lat = front_face_mask(mask_frnt(i,j),mask_frnt(i,jp1),lateral_bc)
                if (mask_lat .gt. 0) ssa_mask_acy(i,j) = mask_lat

            end do 
            end do

        end if 

        return
        
    end subroutine set_ssa_masks

    integer function front_face_mask(code_a,code_b,lateral_bc) result(mask_lat)
        ! SSA mask value for the face between two points with ice-front codes
        ! code_a and code_b (MASK_FRNT_* of yelmo_defs):
        ! 0: not a front face, 3: lateral bc applied, 4: front treated as inner ssa.
        ! A front is only treated as floating or marine across faces whose
        ! ice-free side is ocean. Across ice-free land it is a front grounded
        ! above sea level, whatever the bed of the ice-covered point.

        implicit none

        integer,          intent(IN) :: code_a
        integer,          intent(IN) :: code_b
        character(len=*), intent(IN) :: lateral_bc

        ! Local variables
        integer :: code_ice, code_free

        mask_lat = 0

        if (code_a .gt. 0 .and. code_b .lt. 0) then
            code_ice  = code_a
            code_free = code_b
        else if (code_a .lt. 0 .and. code_b .gt. 0) then
            code_ice  = code_b
            code_free = code_a
        else
            return
        end if

        if (code_free .eq. MASK_FRNT_ICE_FREE_LAND) code_ice = MASK_FRNT_GRND

        select case(trim(lateral_bc))
            case("none")
                mask_lat = 4
            case("floating","float")
                mask_lat = merge(3,4,code_ice .eq. MASK_FRNT_FLOAT)
            case("marine")
                mask_lat = merge(3,4,code_ice .eq. MASK_FRNT_FLOAT .or. code_ice .eq. MASK_FRNT_MARINE)
            case("all")
                mask_lat = 3
        end select

    end function front_face_mask
    
! === INTERNAL ROUTINES ==== 

    subroutine stagger_visc_aa_ab(visc_ab,visc,H_ice,f_ice,boundaries)

        implicit none 

        real(wp), intent(OUT) :: visc_ab(:,:) 
        real(wp), intent(IN)  :: visc(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, k
        integer :: im1, ip1, jm1, jp1 
        integer :: nx, ny 
        integer :: BC

        nx = size(visc,1)
        ny = size(visc,2)

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Initialisation
        visc_ab = 0.0_wp 

        ! Stagger viscosity only using contributions from neighbors that have ice  
        do i = 1, nx 
        do j = 1, ny 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            visc_ab(i,j) = 0.0_wp
            k=0

            if (is_equal(f_ice(i,j),1.0_wp)) then
                k = k+1                              ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(i,j)
            end if

            if (is_equal(f_ice(ip1,j),1.0_wp)) then
                k = k+1                                  ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(ip1,j)
            end if

            if (is_equal(f_ice(i,jp1),1.0_wp)) then
                k = k+1                                  ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(i,jp1)
            end if

            if (is_equal(f_ice(ip1,jp1),1.0_wp)) then
                k = k+1                                      ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(ip1,jp1)
            end if

            if (k .gt. 0) visc_ab(i,j) = visc_ab(i,j)/real(k,wp)

        end do
        end do

        return 

    end subroutine stagger_visc_aa_ab

    elemental subroutine limit_vel(u,u_lim)
        ! Apply a velocity limit (for stability)

        implicit none 

        real(wp), intent(INOUT) :: u  
        real(wp), intent(IN)    :: u_lim

        real(wp), parameter :: tol = TOL_UNDERFLOW

        ! Explicit comparisons (not min/max) so that a NaN passes through
        ! unchanged and is caught by yelmo_check_kill, instead of being
        ! silently mapped to +u_lim.
        if (u .gt.  u_lim) u =  u_lim
        if (u .lt. -u_lim) u = -u_lim

        ! Also avoid underflow errors 
        if (abs(u) .lt. tol) u = 0.0 

        return 

    end subroutine limit_vel

    ! === DIAGNOSTIC OUTPUT ROUTINES ===

    subroutine ssa_diagnostics_write_init(filename,nx,ny,time_init)

        implicit none 

        character(len=*),  intent(IN) :: filename 
        integer,           intent(IN) :: nx 
        integer,           intent(IN) :: ny
        real(wp),          intent(IN) :: time_init

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",     x=0.0_wp,dx=1.0_wp,nx=nx,units="gridpoints")
        call nc_write_dim(filename,"yc",     x=0.0_wp,dx=1.0_wp,nx=ny,units="gridpoints")
        call nc_write_dim(filename,"time",   x=time_init,dx=1.0_wp,nx=1,units="iter",unlimited=.TRUE.)

        return

    end subroutine ssa_diagnostics_write_init

    subroutine ssa_diagnostics_write_step(filename,ux,uy,L2_norm,beta_acx,beta_acy,visc_int, &
                                        ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,H_ice,f_ice,taud_acx,taud_acy, &
                                     taul_int_acx,taul_int_acy,H_grnd,z_sl,z_bed,z_srf,ux_prev,uy_prev,time)

        implicit none 
        
        character(len=*),  intent(IN) :: filename
        real(wp), intent(IN) :: ux(:,:) 
        real(wp), intent(IN) :: uy(:,:) 
        real(wp), intent(IN) :: L2_norm
        real(wp), intent(IN) :: beta_acx(:,:) 
        real(wp), intent(IN) :: beta_acy(:,:) 
        real(wp), intent(IN) :: visc_int(:,:) 
        integer,  intent(IN) :: ssa_mask_acx(:,:) 
        integer,  intent(IN) :: ssa_mask_acy(:,:) 
        real(wp), intent(IN) :: ssa_err_acx(:,:) 
        real(wp), intent(IN) :: ssa_err_acy(:,:) 
        real(wp), intent(IN) :: H_ice(:,:) 
        real(wp), intent(IN) :: f_ice(:,:) 
        real(wp), intent(IN) :: taud_acx(:,:) 
        real(wp), intent(IN) :: taud_acy(:,:) 
        real(wp), intent(IN) :: taul_int_acx(:,:) 
        real(wp), intent(IN) :: taul_int_acy(:,:) 
        real(wp), intent(IN) :: H_grnd(:,:) 
        real(wp), intent(IN) :: z_sl(:,:) 
        real(wp), intent(IN) :: z_bed(:,:) 
        real(wp), intent(IN) :: z_srf(:,:)
        real(wp), intent(IN) :: ux_prev(:,:) 
        real(wp), intent(IN) :: uy_prev(:,:) 
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n, i, j, nx, ny  
        real(wp) :: time_prev 

        nx = size(ux,1)
        ny = size(ux,2) 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write the variables 

        call nc_write(filename,"ux",ux,units="m/yr",long_name="ssa velocity (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy",uy,units="m/yr",long_name="ssa velocity (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_diff",ux-ux_prev,units="m/yr",long_name="ssa velocity difference (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_diff",uy-uy_prev,units="m/yr",long_name="ssa velocity difference (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"L2_norm",L2_norm,dim1="time",start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"beta_acx",beta_acx,units="Pa yr m^-1",long_name="Dragging coefficient (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_acy",beta_acy,units="Pa yr m^-1",long_name="Dragging coefficient (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"visc_int",visc_int,units="Pa yr",long_name="Vertically integrated effective viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ssa_mask_acx",ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ssa_err_acx",ssa_err_acx,units="1",long_name="SSA err (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_err_acy",ssa_err_acy,units="1",long_name="SSA err (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_ice",H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",f_ice,units="1",long_name="Ice-covered fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taud_acx",taud_acx,units="Pa",long_name="Driving stress (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",taud_acy,units="Pa",long_name="Driving stress (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taul_int_acx",taul_int_acx,units="Pa",long_name="Vertically integrated lateral stress (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taul_int_acy",taul_int_acy,units="Pa",long_name="Vertically integrated lateral stress (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_grnd",H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",z_sl,units="m",long_name="Sea level", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_bed",z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine ssa_diagnostics_write_step

end module solver_ssa_ac
