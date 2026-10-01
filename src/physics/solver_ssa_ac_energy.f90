module solver_ssa_ac_energy
    ! SSA momentum solver, energy formulation.
    !
    ! Assembles the symmetric positive-(semi)definite Hessian K of the
    ! discretised SSA energy
    !
    !     E = sum_cells   N_aa (2 u_x^2 + 2 v_y^2 + 2 u_x v_y) dx dy      ! membrane
    !       + sum_corners 1/2 N_ab (u_y + v_x)^2 dx dy                   ! shear
    !       + sum_faces   1/2 beta u^2 dx dy                             ! basal drag
    !       - sum_faces   f u                                            ! driving / front work
    !
    ! with (N, beta) frozen per Picard iteration, so E is quadratic in (u, v)
    ! and K is symmetric; CG can be used in place of BiCGStab.
    !
    ! Staggering convention: yelmo-fortran C-grid. ux(i,j) lives at the RIGHT
    ! face of aa-cell (i,j); uy(i,j) lives at the TOP face; corner (i,j) is
    ! the top-right corner of cell (i,j). (Yelmo.jl uses ux on the LEFT face.)
    !
    ! Assembly is element by element: each cell and corner term is a 4x4
    ! local Hessian over its four velocities, and each velocity is mapped to
    ! a matrix unknown before scattering:
    !
    !   - free:      its own row (inner and front faces, ssa_mask = 1, 3, 4)
    !   - Dirichlet: ssa_mask = 0 (u = 0), -1 (u prescribed) and no-slip
    !                domain edges; the column is lifted to the RHS
    !   - tied:      free-slip domain edges (u_edge = u_inner, as in the
    !                residual solver); the unknown is folded into its inner
    !                root (K_red = T^T K T), and takes the root's solution
    !                after the solve (lgs%copy_from)
    !   - ghosts beyond a domain edge: periodic wrap, free-slip copy of the
    !                edge velocity, or zero (no-slip)
    !
    ! Because every term enters through a local Hessian and a linear map,
    ! K is symmetric for any mask and boundary type. At inner faces the
    ! result equals the residual formulation:
    !
    !     K_inner = - A_residual_inner * dx * dy,   b_inner = - taud_inner * dx * dy
    !
    ! At a calving front (ssa_mask = 3) the RHS is the boundary work
    ! +-taul_int*dy only: the front-face driving stress (taken across the
    ! ice front) is the same front force.
    !
    ! Each matrix row gathers its terms from the elements around its
    ! unknown(s), so rows are independent and assembled in parallel.

    use yelmo_defs, only : sp, dp, wp, io_unit_err, TOL, TOL_UNDERFLOW
    use solver_linear
    use solver_ssa_ac, only : stagger_visc_aa_ab

    implicit none

    private
    public :: linear_solver_matrix_ssa_ac_csr_2D_energy

    ! Kinds of velocity unknowns
    integer, parameter :: DOF_FREE = 0
    integer, parameter :: DOF_DIR  = 1
    integer, parameter :: DOF_TIED = 2

    ! Upper limit of entries in one matrix row (inner rows have 9)
    integer, parameter :: NNZ_ROW_MAX = 48

contains

    subroutine linear_solver_matrix_ssa_ac_csr_2D_energy(lgs,ux,uy,beta_acx,beta_acy, &
                            N_aa,ssa_mask_acx,ssa_mask_acy,H_ice,f_ice,taud_acx, &
                            taud_acy,taul_int_acx,taul_int_acy,dx,dy,beta_min,boundaries)
        ! Energy-formulation analogue of linear_solver_matrix_ssa_ac_csr_2D.
        ! Same argument list so the two assemblers are interchangeable from the
        ! Picard loop. Assembles K (symmetric) and b such that K * [u; v] = b.

        implicit none

        type(linear_solver_class), intent(INOUT) :: lgs
        real(wp), intent(IN) :: ux(:,:)                 ! [m yr^-1] horizontal velocity x (acx-nodes)
        real(wp), intent(IN) :: uy(:,:)                 ! [m yr^-1] horizontal velocity y (acy-nodes)
        real(wp), intent(IN) :: beta_acx(:,:)           ! [Pa yr m^-1] basal friction (acx-nodes)
        real(wp), intent(IN) :: beta_acy(:,:)           ! [Pa yr m^-1] basal friction (acy-nodes)
        real(wp), intent(IN) :: N_aa(:,:)               ! [Pa yr m] vertically integrated viscosity (aa-nodes)
        integer,  intent(IN) :: ssa_mask_acx(:,:)       ! [--] ssa solver action mask (acx-nodes)
        integer,  intent(IN) :: ssa_mask_acy(:,:)       ! [--] ssa solver action mask (acy-nodes)
        real(wp), intent(IN) :: H_ice(:,:)              ! [m]  ice thickness (aa-nodes)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: taud_acx(:,:)           ! [Pa] driving stress (acx-nodes)
        real(wp), intent(IN) :: taud_acy(:,:)           ! [Pa] driving stress (acy-nodes)
        real(wp), intent(IN) :: taul_int_acx(:,:)       ! [Pa m] vertically integrated lateral stress (acx-nodes)
        real(wp), intent(IN) :: taul_int_acy(:,:)       ! [Pa m] vertically integrated lateral stress (acy-nodes)
        real(wp), intent(IN) :: dx, dy
        real(wp), intent(IN) :: beta_min                ! [Pa yr m^-1] minimum allowed basal friction
        character(len=*), intent(IN) :: boundaries

        ! Local variables
        integer  :: nx, ny, nmax
        integer  :: i, j, n, m, q, depth, nb, nnz
        real(dp) :: dxdy, dyodx, dxody
        real(dp) :: bval
        integer  :: cols(NNZ_ROW_MAX)
        real(dp) :: vals(NNZ_ROW_MAX)
        logical  :: per_x, per_y

        ! Boundary conditions counterclockwise unit circle:
        ! 1: x, right border; 2: y, upper; 3: x, left; 4: y, lower
        character(len=56) :: bcs(4)
        logical :: bc_per(4), bc_free(4)          ! bcs(k) is "periodic" / "free-slip"

        real(wp), allocatable :: N_ab(:,:)
        integer,  allocatable :: dkind(:)         ! DOF_FREE, DOF_DIR or DOF_TIED
        integer,  allocatable :: dpart(:)         ! tied: partner unknown (before resolving)
        integer,  allocatable :: droot(:)         ! free: itself; tied: free root unknown
        real(dp), allocatable :: dval(:)          ! Dirichlet value
        integer,  allocatable :: row_nnz(:)
        integer,  allocatable :: mem_ptr(:)       ! Tied unknowns folded into each root
        integer,  allocatable :: mem_list(:)
        integer,  allocatable :: mem_fill(:)

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Initialise the lgs object if needed (n_terms=9 matches the residual assembler)
        if (.not. allocated(lgs%x_value)) then
            call linear_solver_init(lgs,nx,ny,nvar=2,n_terms=9)
        end if

        nmax = lgs%nmax

        ! Define border conditions (only choices: no-slip, free-slip, periodic)
        select case(trim(boundaries))
            case("MISMIP3D")
                bcs(1) = "free-slip"; bcs(2) = "periodic"
                bcs(3) = "no-slip";   bcs(4) = "periodic"
            case("TROUGH")
                bcs(1) = "free-slip"; bcs(2) = "periodic"
                bcs(3) = "no-slip";   bcs(4) = "periodic"
            case("periodic")
                bcs(1:4) = "periodic"
            case("periodic-x")
                bcs(1) = "periodic";  bcs(2) = "free-slip"
                bcs(3) = "periodic";  bcs(4) = "free-slip"
            case("periodic-y")
                bcs(1) = "free-slip"; bcs(2) = "periodic"
                bcs(3) = "free-slip"; bcs(4) = "periodic"
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
        per_x   = bc_per(1) .and. bc_per(3)
        per_y   = bc_per(2) .and. bc_per(4)

        ! Stencil prefactors
        dxdy  = real(dx,dp)*real(dy,dp)
        dyodx = real(dy,dp)/real(dx,dp)
        dxody = real(dx,dp)/real(dy,dp)

        ! Stagger depth-integrated viscosity to ab-nodes
        allocate(N_ab(nx,ny))
        call stagger_visc_aa_ab(N_ab,N_aa,f_ice,boundaries)

        ! ================================================================
        ! 1. Classify the unknowns and resolve tied chains to free roots
        ! ================================================================

        allocate(dkind(nmax),dpart(nmax),droot(nmax),dval(nmax),row_nnz(nmax))

        !$omp parallel do collapse(2) private(i,j,n)
        do j = 1, ny
        do i = 1, nx
            n = 2*lgs%ij2n(i,j)-1
            call classify_ux(i,j,dkind(n),dpart(n),dval(n))
            n = 2*lgs%ij2n(i,j)
            call classify_uy(i,j,dkind(n),dpart(n),dval(n))
        end do
        end do
        !$omp end parallel do

        !$omp parallel do private(n,m,depth)
        do n = 1, nmax
            m = n
            depth = 0
            do while (dkind(m) .eq. DOF_TIED)
                m = dpart(m)
                depth = depth + 1
                if (depth .gt. 8) then
                    write(io_unit_err,*) "linear_solver_matrix_ssa_ac_csr_2D_energy:: Error: &
                        &free-slip chain does not end at a free unknown, row ", n
                    error stop 1
                end if
            end do
            droot(n) = m
        end do
        !$omp end parallel do

        ! A tied unknown whose chain ends at a Dirichlet unknown is Dirichlet
        do n = 1, nmax
            if (dkind(n) .eq. DOF_TIED) then
                if (dkind(droot(n)) .eq. DOF_DIR) then
                    dkind(n) = DOF_DIR
                    dval(n)  = dval(droot(n))
                    droot(n) = n
                end if
            end if
        end do

        ! Tied unknowns folded into each root (edge rows only; built serially)
        allocate(mem_ptr(nmax+1),mem_fill(nmax))
        mem_fill = 0
        do n = 1, nmax
            if (dkind(n) .eq. DOF_TIED) mem_fill(droot(n)) = mem_fill(droot(n)) + 1
        end do
        mem_ptr(1) = 1
        do n = 1, nmax
            mem_ptr(n+1) = mem_ptr(n) + mem_fill(n)
        end do
        allocate(mem_list(max(mem_ptr(nmax+1)-1,1)))
        mem_fill = 0
        do n = 1, nmax
            if (dkind(n) .eq. DOF_TIED) then
                m = droot(n)
                mem_list(mem_ptr(m)+mem_fill(m)) = n
                mem_fill(m) = mem_fill(m) + 1
            end if
        end do

        ! ================================================================
        ! 2. Count the entries of each row, then set the row pointers
        ! ================================================================

        !$omp parallel do schedule(dynamic,1024) private(n,nb,cols,vals,bval)
        do n = 1, nmax
            if (dkind(n) .eq. DOF_FREE) then
                call assemble_row(n,nb,cols,vals,bval)
                row_nnz(n) = nb
            else
                row_nnz(n) = 1
            end if
        end do
        !$omp end parallel do

        lgs%a_ptr(1) = 1
        do n = 1, nmax
            lgs%a_ptr(n+1) = lgs%a_ptr(n) + row_nnz(n)
        end do

        nnz = lgs%a_ptr(nmax+1)-1
        if (nnz .gt. size(lgs%a_value)) then
            deallocate(lgs%a_value,lgs%a_index)
            allocate(lgs%a_value(nnz),lgs%a_index(nnz))
            lgs%n_sprs = nnz
        end if

        ! ================================================================
        ! 3. Fill the rows (CSR, columns ascending)
        ! ================================================================

        !$omp parallel do schedule(dynamic,1024) private(n,nb,cols,vals,bval,i,j,q)
        do n = 1, nmax

            q = (n+1)/2
            i = lgs%n2i(q)
            j = lgs%n2j(q)

            select case(dkind(n))

                case(DOF_FREE)
                    call assemble_row(n,nb,cols,vals,bval)
                    call sort_row(nb,cols,vals)
                    lgs%a_index(lgs%a_ptr(n):lgs%a_ptr(n+1)-1) = cols(1:nb)
                    lgs%a_value(lgs%a_ptr(n):lgs%a_ptr(n+1)-1) = vals(1:nb)
                    lgs%b_value(n)   = bval
                    lgs%copy_from(n) = 0
                    if (mod(n,2) .eq. 1) then
                        lgs%x_value(n) = ux(i,j)
                    else
                        lgs%x_value(n) = uy(i,j)
                    end if

                case(DOF_DIR)
                    lgs%a_index(lgs%a_ptr(n)) = n
                    lgs%a_value(lgs%a_ptr(n)) = 1.0_dp
                    lgs%b_value(n)   = dval(n)
                    lgs%x_value(n)   = dval(n)
                    lgs%copy_from(n) = 0

                case(DOF_TIED)
                    ! Decoupled identity row; the unknown takes its root's solution
                    lgs%a_index(lgs%a_ptr(n)) = n
                    lgs%a_value(lgs%a_ptr(n)) = 1.0_dp
                    lgs%b_value(n)   = 0.0_dp
                    lgs%x_value(n)   = 0.0_dp
                    lgs%copy_from(n) = droot(n)

            end select

        end do
        !$omp end parallel do

        return

    contains

        ! ---- Unknown classification ----------------------------------

        subroutine classify_ux(i,j,kind,part,val)
            ! Same border order as the residual assembler: left, right, lower, upper.
            integer,  intent(IN)  :: i, j
            integer,  intent(OUT) :: kind, part
            real(dp), intent(OUT) :: val

            kind = DOF_FREE
            part = 0
            val  = 0.0_dp

            if (ssa_mask_acx(i,j) .eq. 0) then
                kind = DOF_DIR
            else if (ssa_mask_acx(i,j) .eq. -1) then
                kind = DOF_DIR
                val  = ux(i,j)
            else if (i .eq. 1 .and. .not. bc_per(3)) then
                call edge(bc_free(3),2*lgs%ij2n(2,j)-1,kind,part)
            else if (i .eq. nx .and. .not. bc_per(1)) then
                call edge(bc_free(1),2*lgs%ij2n(nx-1,j)-1,kind,part)
            else if (j .eq. 1 .and. .not. bc_per(4)) then
                call edge(bc_free(4),2*lgs%ij2n(i,2)-1,kind,part)
            else if (j .eq. ny .and. .not. bc_per(2)) then
                call edge(bc_free(2),2*lgs%ij2n(i,ny-1)-1,kind,part)
            end if

        end subroutine classify_ux

        subroutine classify_uy(i,j,kind,part,val)
            ! Same border order as the residual assembler: lower, upper, left, right.
            integer,  intent(IN)  :: i, j
            integer,  intent(OUT) :: kind, part
            real(dp), intent(OUT) :: val

            kind = DOF_FREE
            part = 0
            val  = 0.0_dp

            if (ssa_mask_acy(i,j) .eq. 0) then
                kind = DOF_DIR
            else if (ssa_mask_acy(i,j) .eq. -1) then
                kind = DOF_DIR
                val  = uy(i,j)
            else if (j .eq. 1 .and. .not. bc_per(4)) then
                call edge(bc_free(4),2*lgs%ij2n(i,2),kind,part)
            else if (j .eq. ny .and. .not. bc_per(2)) then
                call edge(bc_free(2),2*lgs%ij2n(i,ny-1),kind,part)
            else if (i .eq. 1 .and. .not. bc_per(3)) then
                call edge(bc_free(3),2*lgs%ij2n(2,j),kind,part)
            else if (i .eq. nx .and. .not. bc_per(1)) then
                call edge(bc_free(1),2*lgs%ij2n(nx-1,j),kind,part)
            end if

        end subroutine classify_uy

        subroutine edge(free_slip,inner,kind,part)
            ! Domain-edge unknown: tied to its inner neighbour (free-slip) or zero (no-slip)
            logical, intent(IN)  :: free_slip
            integer, intent(IN)  :: inner
            integer, intent(OUT) :: kind, part
            if (free_slip) then
                kind = DOF_TIED
                part = inner
            else
                kind = DOF_DIR
                part = 0
            end if
        end subroutine edge

        ! ---- Velocities of the elements, with ghosts beyond the edges --------

        integer function dof_ux(ii,jj) result(d)
            ! Unknown of ux(ii,jj); ii, jj may lie one beyond the domain.
            ! Returns 0 for a zero ghost (no-slip side).
            integer, intent(IN) :: ii, jj
            integer :: ic, jc
            ic = ii
            jc = jj
            d  = 0
            if (.not. wrap_index(ic,nx,per_x,bc_free(3),bc_free(1))) return
            if (.not. wrap_index(jc,ny,per_y,bc_free(4),bc_free(2))) return
            d = 2*lgs%ij2n(ic,jc)-1
        end function dof_ux

        integer function dof_uy(ii,jj) result(d)
            integer, intent(IN) :: ii, jj
            integer :: ic, jc
            ic = ii
            jc = jj
            d  = 0
            if (.not. wrap_index(ic,nx,per_x,bc_free(3),bc_free(1))) return
            if (.not. wrap_index(jc,ny,per_y,bc_free(4),bc_free(2))) return
            d = 2*lgs%ij2n(ic,jc)
        end function dof_uy

        logical function wrap_index(k,nk,periodic,free_lo,free_hi) result(ok)
            ! Map an index one beyond the domain: periodic wrap, free-slip
            ! copy of the edge value, or no value (zero ghost, no-slip).
            integer, intent(INOUT) :: k
            integer, intent(IN)    :: nk
            logical, intent(IN)    :: periodic, free_lo, free_hi
            ok = .TRUE.
            if (k .lt. 1) then
                if (periodic) then
                    k = k + nk
                else if (free_lo) then
                    k = 1
                else
                    ok = .FALSE.
                end if
            else if (k .gt. nk) then
                if (periodic) then
                    k = k - nk
                else if (free_hi) then
                    k = nk
                else
                    ok = .FALSE.
                end if
            end if
        end function wrap_index

        logical function element_index(k,nk,periodic) result(ok)
            ! Cells and corners exist inside the domain (wrapped if periodic)
            integer, intent(INOUT) :: k
            integer, intent(IN)    :: nk
            logical, intent(IN)    :: periodic
            ok = .TRUE.
            if (k .lt. 1 .or. k .gt. nk) then
                if (periodic) then
                    k = modulo(k-1,nk) + 1
                else
                    ok = .FALSE.
                end if
            end if
        end function element_index

        ! ---- Row assembly ------------------------------------------------

        subroutine assemble_row(r,nb,cols,vals,bval)
            ! Row of the root unknown r: all terms of the energy that involve r
            ! or an unknown tied to r, differentiated with respect to r.
            integer,  intent(IN)  :: r
            integer,  intent(OUT) :: nb
            integer,  intent(OUT) :: cols(NNZ_ROW_MAX)
            real(dp), intent(OUT) :: vals(NNZ_ROW_MAX)
            real(dp), intent(OUT) :: bval

            integer :: k

            nb   = 0
            bval = 0.0_dp

            call add_member(r,r,nb,cols,vals,bval)
            do k = mem_ptr(r), mem_ptr(r+1)-1
                call add_member(mem_list(k),r,nb,cols,vals,bval)
            end do

        end subroutine assemble_row

        subroutine add_member(mm,r,nb,cols,vals,bval)
            ! Terms of unknown mm (r itself or tied to r) added to the row of r
            integer,  intent(IN)    :: mm, r
            integer,  intent(INOUT) :: nb
            integer,  intent(INOUT) :: cols(NNZ_ROW_MAX)
            real(dp), intent(INOUT) :: vals(NNZ_ROW_MAX)
            real(dp), intent(INOUT) :: bval

            integer  :: i, j, ip1, jp1, qq, mask
            integer  :: ci, cj
            integer  :: d(4)
            real(dp) :: H(4,4)
            real(dp) :: beta_now

            qq = (mm+1)/2
            i  = lgs%n2i(qq)
            j  = lgs%n2j(qq)

            ! Neighbour index for the front orientation (periodic wrap, as before)
            ip1 = i+1; if (ip1 .eq. nx+1) ip1 = 1
            jp1 = j+1; if (jp1 .eq. ny+1) jp1 = 1

            if (mod(mm,2) .eq. 1) then
                ! ---- ux(i,j) ----

                ! Face: basal drag and driving stress, or boundary work at a front
                mask     = ssa_mask_acx(i,j)
                beta_now = beta_acx(i,j)
                if (mask .eq. 1 .and. beta_acx(i,j) .eq. 0.0_wp) beta_now = beta_min
                if (mask .eq. 3) then
                    ! Calving front: only the ice half of the face's control area has drag
                    call add_entry(r,0.5_dp*beta_now*dxdy,nb,cols,vals)
                    if (f_ice(i,j) .eq. 1.0_wp .and. f_ice(ip1,j) .lt. 1.0_wp) then
                        bval = bval + taul_int_acx(i,j)*real(dy,dp)
                    else
                        bval = bval - taul_int_acx(i,j)*real(dy,dp)
                    end if
                else
                    call add_entry(r,beta_now*dxdy,nb,cols,vals)
                    bval = bval - taud_acx(i,j)*dxdy
                end if

                ! Cells (i,j) and (i+1,j); corners (i,j) and (i,j-1)
                ci = i;   cj = j;   if (cell_dofs(ci,cj,d,H))   call add_element(mm,d,H,nb,cols,vals,bval)
                ci = i+1; cj = j;   if (cell_dofs(ci,cj,d,H))   call add_element(mm,d,H,nb,cols,vals,bval)
                ci = i;   cj = j;   if (corner_dofs(ci,cj,d,H)) call add_element(mm,d,H,nb,cols,vals,bval)
                ci = i;   cj = j-1; if (corner_dofs(ci,cj,d,H)) call add_element(mm,d,H,nb,cols,vals,bval)

            else
                ! ---- uy(i,j) ----

                mask     = ssa_mask_acy(i,j)
                beta_now = beta_acy(i,j)
                if (mask .eq. 1 .and. beta_acy(i,j) .eq. 0.0_wp) beta_now = beta_min
                if (mask .eq. 3) then
                    ! Calving front: only the ice half of the face's control area has drag
                    call add_entry(r,0.5_dp*beta_now*dxdy,nb,cols,vals)
                    if (f_ice(i,j) .eq. 1.0_wp .and. f_ice(i,jp1) .lt. 1.0_wp) then
                        bval = bval + taul_int_acy(i,j)*real(dx,dp)
                    else
                        bval = bval - taul_int_acy(i,j)*real(dx,dp)
                    end if
                else
                    call add_entry(r,beta_now*dxdy,nb,cols,vals)
                    bval = bval - taud_acy(i,j)*dxdy
                end if

                ! Cells (i,j) and (i,j+1); corners (i,j) and (i-1,j)
                ci = i;   cj = j;   if (cell_dofs(ci,cj,d,H))   call add_element(mm,d,H,nb,cols,vals,bval)
                ci = i;   cj = j+1; if (cell_dofs(ci,cj,d,H))   call add_element(mm,d,H,nb,cols,vals,bval)
                ci = i;   cj = j;   if (corner_dofs(ci,cj,d,H)) call add_element(mm,d,H,nb,cols,vals,bval)
                ci = i-1; cj = j;   if (corner_dofs(ci,cj,d,H)) call add_element(mm,d,H,nb,cols,vals,bval)

            end if

        end subroutine add_member

        logical function cell_dofs(ci,cj,d,H) result(ok)
            ! Membrane term of cell (ci,cj): N (2 u_x^2 + 2 v_y^2 + 2 u_x v_y) dx dy,
            ! u_x = (a2-a1)/dx, v_y = (b2-b1)/dy with
            ! a1 = ux(ci-1,cj), a2 = ux(ci,cj), b1 = uy(ci,cj-1), b2 = uy(ci,cj).
            integer,  intent(INOUT) :: ci, cj
            integer,  intent(OUT)   :: d(4)
            real(dp), intent(OUT)   :: H(4,4)
            real(dp) :: Nc, hx, hy, hc

            ok = element_index(ci,nx,per_x)
            if (ok) ok = element_index(cj,ny,per_y)
            if (.not. ok) return

            d(1) = dof_ux(ci-1,cj)
            d(2) = dof_ux(ci,  cj)
            d(3) = dof_uy(ci,  cj-1)
            d(4) = dof_uy(ci,  cj)

            Nc = real(N_aa(ci,cj),dp)
            hx = 4.0_dp*Nc*dyodx
            hy = 4.0_dp*Nc*dxody
            hc = 2.0_dp*Nc

            H(1,:) = [  hx, -hx,  hc, -hc ]
            H(2,:) = [ -hx,  hx, -hc,  hc ]
            H(3,:) = [  hc, -hc,  hy, -hy ]
            H(4,:) = [ -hc,  hc, -hy,  hy ]

        end function cell_dofs

        logical function corner_dofs(ci,cj,d,H) result(ok)
            ! Shear term of corner (ci,cj): 1/2 N_ab (u_y + v_x)^2 dx dy,
            ! u_y = (c2-c1)/dy, v_x = (e2-e1)/dx with
            ! c1 = ux(ci,cj), c2 = ux(ci,cj+1), e1 = uy(ci,cj), e2 = uy(ci+1,cj).
            integer,  intent(INOUT) :: ci, cj
            integer,  intent(OUT)   :: d(4)
            real(dp), intent(OUT)   :: H(4,4)
            real(dp) :: Nc, hx, hy

            ok = element_index(ci,nx,per_x)
            if (ok) ok = element_index(cj,ny,per_y)
            if (.not. ok) return

            d(1) = dof_ux(ci,  cj)
            d(2) = dof_ux(ci,  cj+1)
            d(3) = dof_uy(ci,  cj)
            d(4) = dof_uy(ci+1,cj)

            Nc = real(N_ab(ci,cj),dp)
            hx = Nc*dxody
            hy = Nc*dyodx

            H(1,:) = [  hx, -hx,  Nc, -Nc ]
            H(2,:) = [ -hx,  hx, -Nc,  Nc ]
            H(3,:) = [  Nc, -Nc,  hy, -hy ]
            H(4,:) = [ -Nc,  Nc, -hy,  hy ]

        end function corner_dofs

        subroutine add_element(mm,d,H,nb,cols,vals,bval)
            ! For every slot of the element holding unknown mm, add its Hessian row:
            ! free columns to their root, Dirichlet columns lifted to the RHS.
            integer,  intent(IN)    :: mm
            integer,  intent(IN)    :: d(4)
            real(dp), intent(IN)    :: H(4,4)
            integer,  intent(INOUT) :: nb
            integer,  intent(INOUT) :: cols(NNZ_ROW_MAX)
            real(dp), intent(INOUT) :: vals(NNZ_ROW_MAX)
            real(dp), intent(INOUT) :: bval
            integer :: l, k, c

            do l = 1, 4
                if (d(l) .ne. mm) cycle
                do k = 1, 4
                    c = d(k)
                    if (c .eq. 0) cycle                         ! zero ghost
                    select case(dkind(c))
                        case(DOF_FREE)
                            call add_entry(c,H(l,k),nb,cols,vals)
                        case(DOF_TIED)
                            call add_entry(droot(c),H(l,k),nb,cols,vals)
                        case(DOF_DIR)
                            bval = bval - H(l,k)*dval(c)
                    end select
                end do
            end do

        end subroutine add_element

        subroutine add_entry(c,v,nb,cols,vals)
            integer,  intent(IN)    :: c
            real(dp), intent(IN)    :: v
            integer,  intent(INOUT) :: nb
            integer,  intent(INOUT) :: cols(NNZ_ROW_MAX)
            real(dp), intent(INOUT) :: vals(NNZ_ROW_MAX)
            integer :: k
            do k = 1, nb
                if (cols(k) .eq. c) then
                    vals(k) = vals(k) + v
                    return
                end if
            end do
            nb = nb + 1
            if (nb .gt. NNZ_ROW_MAX) then
                write(io_unit_err,*) "linear_solver_matrix_ssa_ac_csr_2D_energy:: Error: &
                    &more than NNZ_ROW_MAX entries in a row."
                error stop 1
            end if
            cols(nb) = c
            vals(nb) = v
        end subroutine add_entry

        subroutine sort_row(nb,cols,vals)
            ! Insertion sort by column (rows are short)
            integer,  intent(IN)    :: nb
            integer,  intent(INOUT) :: cols(NNZ_ROW_MAX)
            real(dp), intent(INOUT) :: vals(NNZ_ROW_MAX)
            integer  :: k, l, c
            real(dp) :: v
            do k = 2, nb
                c = cols(k)
                v = vals(k)
                l = k-1
                do while (l .ge. 1)
                    if (cols(l) .le. c) exit
                    cols(l+1) = cols(l)
                    vals(l+1) = vals(l)
                    l = l-1
                end do
                cols(l+1) = c
                vals(l+1) = v
            end do
        end subroutine sort_row

    end subroutine linear_solver_matrix_ssa_ac_csr_2D_energy

end module solver_ssa_ac_energy
