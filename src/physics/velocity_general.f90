module velocity_general 
    ! This module contains general routines that are used by several solvers. 
    
    use yelmo_defs ,only  : sp, dp, wp, tol_underflow, io_unit_err, jacobian_3D_class, MASK_FRNT_ICE_FREE_LAND, &
                            A_FRONT_MIN
    use yelmo_tools, only : boundary_code, get_neighbor_indices_bc_codes, get_periodic_directions, &
                            integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax
    use gaussian_quadrature, only : gq2D_class, gq2D_init, gq2D_to_nodes_aa, &
                                    gq2D_to_nodes_acx, gq2D_to_nodes_acy, &
                                    gq3D_class, gq3D_init, gq3D_to_nodes_aa, &
                                    gq3D_to_nodes_acx, gq3D_to_nodes_acy

    use solver_linear, only : linear_solver_status_str
    use deformation, only : calc_strain_rate_horizontal_2D, calc_active_faces

    use solver_ssa_ac, only : ssa_diagnostics_write_init, ssa_diagnostics_write_step

    implicit none 

    private 
    public :: calc_uz_3D_jac
    public :: calc_uz_3D
    public :: calc_uz_3D_aa
    public :: calc_driving_stress
    public :: calc_driving_stress_gl
    public :: calc_lateral_bc_stress_2D
    public :: set_inactive_margins
    public :: calc_ice_flux
    public :: calc_grounding_line_flux
    public :: calc_vel_ratio

    public :: picard_calc_error 
    public :: picard_calc_error_angle 
    public :: picard_calc_convergence_l1rel_matrix
    public :: picard_calc_convergence_l2
    public :: picard_relax_vel
    public :: picard_relax_visc
    
contains 
    
    subroutine calc_uz_3D_jac(uz,uz_star,ux,uy,jvel,H_ice,f_ice,f_grnd,smb,bmb,dzbdt,dzsdt, &
                                    dzsdx,dzsdy,dzbdx,dzbdy,zeta_aa,zeta_ac,dx,dy,use_bmb,boundaries)
        ! Following algorithm outlined by the Glimmer ice sheet model:
        ! https://www.geos.ed.ac.uk/~mhagdorn/glide/glide-doc/glimmer_htmlse9.html#x17-660003.1.5

        ! Note: dzbdt and dzsdt are the kinematic rates of the base and surface of the
        ! ice column (tpo%now%dzbdt_kin, dzsdt_kin): vertical column change plus bedrock
        ! and sea-level rates, without lateral changes (calving, front advance, removals).

        implicit none 

        real(wp), intent(OUT) :: uz(:,:,:)          ! nx,ny,nz_ac
        real(wp), intent(OUT) :: uz_star(:,:,:)     ! nx,ny,nz_ac
        real(wp), intent(IN)  :: ux(:,:,:)          ! nx,ny,nz_aa
        real(wp), intent(IN)  :: uy(:,:,:)          ! nx,ny,nz_aa
        type(jacobian_3D_class), intent(IN) :: jvel 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: smb(:,:) 
        real(wp), intent(IN)  :: bmb(:,:) 
        real(wp), intent(IN)  :: dzbdt(:,:)         ! [m/a] Kinematic rate of the column base
        real(wp), intent(IN)  :: dzsdt(:,:) 
        real(wp), intent(IN)  :: dzsdx(:,:) 
        real(wp), intent(IN)  :: dzsdy(:,:) 
        real(wp), intent(IN)  :: dzbdx(:,:) 
        real(wp), intent(IN)  :: dzbdy(:,:) 
        real(wp), intent(IN)  :: zeta_aa(:)    ! z-coordinate, aa-nodes 
        real(wp), intent(IN)  :: zeta_ac(:)    ! z-coordinate, ac-nodes  
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        logical,  intent(IN)  :: use_bmb
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac
        integer :: im1, ip1, jm1, jp1  
        integer  :: im1m, ip1m, jm1m, jp1m 
        integer :: kup, kdn, kmid
        real(wp) :: f_bmb
        real(wp) :: H_now
        real(wp) :: H_inv
        real(wp) :: dzbdx_aa
        real(wp) :: dzbdy_aa
        real(wp) :: dzsdx_aa
        real(wp) :: dzsdy_aa 
        real(wp) :: dudx_aa
        real(wp) :: dvdy_aa
        real(wp) :: dudz_aa 
        real(wp) :: dvdz_aa 
        real(wp) :: ux_aa 
        real(wp) :: uy_aa 
        real(wp) :: uz_grid 
        real(wp) :: uz_srf 
        real(wp) :: zeta_now 
        real(wp) :: c_x 
        real(wp) :: c_y 
        real(wp) :: c_t 
        real(wp) :: c_z 

        real(wp) :: dzsdtn(4)
        real(wp) :: dzbdxn(4)
        real(wp) :: dzbdyn(4)
        real(wp) :: dzsdxn(4)
        real(wp) :: dzsdyn(4)
        real(wp) :: dudxn(4) 
        real(wp) :: dvdyn(4) 
        real(wp) :: uxn_up(4) 
        real(wp) :: uxn_dn(4) 
        real(wp) :: uxn(4) 
        real(wp) :: uyn_up(4) 
        real(wp) :: uyn_dn(4)
        real(wp) :: uyn(4)  
        real(wp) :: dudzn(4) 
        real(wp) :: dvdzn(4) 
        
        real(wp) :: dudxn8(8) 
        real(wp) :: dvdyn8(8) 

        real(wp) :: dzsdt_now
        real(wp) :: dzbdt_now 
        
        logical, allocatable :: is_ice(:,:)
        

        logical, allocatable :: act_acx(:,:), act_acy(:,:)

        type(gq2D_class) :: gq2D, gq2D_global
        type(gq3D_class) :: gq3D, gq3D_global
        real(wp) :: dz0, dz1
        integer  :: km1, kp1
        logical, parameter :: use_gq3D = .TRUE.

        integer  :: BC

        ! Initialize gaussian quadrature calculations
        call gq2D_init(gq2D_global)
        if (use_gq3D) call gq3D_init(gq3D_global)

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1) 
        
        allocate(is_ice(nx,ny))
        is_ice = (f_ice .eq. 1.0)
        
        ! Initialize vertical velocity to zero 
        uz = 0.0 

        ! Define switch for bmb
        if (use_bmb) then 
            f_bmb = 1.0 
        else 
            f_bmb = 0.0 
        end if 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Velocity faces with a solution (next to an ice-covered cell); the
        ! others hold zero velocity and do not enter the quadrature means
        allocate(act_acx(nx,ny),act_acy(nx,ny))
        call calc_active_faces(act_acx,act_acy,f_ice,BC)

        ! Next, calculate vertical velocity at each point through the column

        !$omp parallel &
        !$omp& private(i,j,im1,ip1,jm1,jp1,gq2D,gq3D) &
        !$omp& private(dzsdt_now,dzbdt_now,H_now,H_inv,dzbdxn,dzbdx_aa,dzbdyn,dzbdy_aa) &
        !$omp& private(dzsdxn,dzsdx_aa,dzsdyn,dzsdy_aa,uxn,ux_aa,uyn,uy_aa,uz_grid,k,kmid) &
        !$omp& private(dudxn,dudx_aa,dvdyn,dvdy_aa,km1,kp1,dz0,dudxn8,dvdyn8) &
        !$omp& private(kup,kdn,uxn_up,uxn_dn,uyn_up,uyn_dn,zeta_now,c_x,c_y,c_t,c_z) &
        !$omp& shared(gq2D_global,gq3D_global)
        gq2D = gq2D_global
        gq3D = gq3D_global

        !$omp do collapse(2) schedule(dynamic,64)
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! Kinematic rates of the column surface and base
            dzsdt_now = dzsdt(i,j)
            dzbdt_now = dzbdt(i,j)

            if (f_ice(i,j) .eq. 1.0) then

                H_now  = H_ice(i,j) 
                H_inv = 1.0/H_now 

                ! Get the centered ice-base gradient
                call gq2D_to_nodes_acx(gq2d,dzbdxn,dzbdx,dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                dzbdx_aa = sum(dzbdxn*gq2d%wt)/gq2d%wt_tot

                call gq2D_to_nodes_acy(gq2d,dzbdyn,dzbdy,dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                dzbdy_aa = sum(dzbdyn*gq2d%wt)/gq2d%wt_tot

                ! Get the centered surface gradient
                call gq2D_to_nodes_acx(gq2d,dzsdxn,dzsdx,dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                dzsdx_aa = sum(dzsdxn*gq2d%wt)/gq2d%wt_tot
                
                call gq2D_to_nodes_acy(gq2d,dzsdyn,dzsdy,dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                dzsdy_aa = sum(dzsdyn*gq2d%wt)/gq2d%wt_tot
                
                ! Get the aa-node centered horizontal velocity at the base
                call gq2D_to_nodes_acx(gq2d,uxn,ux(:,:,1),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                ux_aa = sum(uxn*gq2d%wt)/gq2d%wt_tot
                
                call gq2D_to_nodes_acy(gq2d,uyn,uy(:,:,1),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                uy_aa = sum(uyn*gq2d%wt)/gq2d%wt_tot

                ! Determine grid vertical velocity at the base due to sigma-coordinates 
                ! Glimmer, Eq. 3.35 
                ! ajr, 2020-01-27, untested:::
!                 uz_grid = dzsdt(i,j) + (ux_aa*dzsdx_aa + uy_aa*dzsdy_aa) &
!                             - ( (1.0_wp-zeta_ac(1))*dHdt(i,j) + ux_aa*dHdx_aa + uy_aa*dHdy_aa )
                uz_grid = 0.0_wp 

                ! ===================================================================
                ! Greve and Blatter (2009) style:

                ! Determine basal vertical velocity for this grid point 
                ! Following Eq. 5.31 of Greve and Blatter (2009)
                uz(i,j,1) = dzbdt_now + uz_grid + f_bmb*bmb(i,j) + ux_aa*dzbdx_aa + uy_aa*dzbdy_aa
                if (abs(uz(i,j,1)) .lt. TOL_UNDERFLOW) uz(i,j,1) = 0.0_wp 

                ! Determine surface vertical velocity following kinematic boundary condition 
                ! Glimmer, Eq. 3.10 [or Folwer, Chpt 10, Eq. 10.8]
                !uz_srf = dzsdt(i,j) + ux_aa*dzsdx_aa + uy_aa*dzsdy_aa - smb(i,j) 

                ! Integrate upward to each point above base until just below surface is reached 
                ! Integrate on vertical ac-nodes (ie, vertical cell borders between aa-node centers)
                do k = 2, nz_ac 

                    ! Note: center of cell below this one is zeta_aa(k-1)
                    kmid = k-1

                    ! Get dudz/dvdz values at vertical aa-nodes, in order
                    ! to vertically integrate each cell up to ac-node border.

if (.not. use_gq3D) then
    ! 2D QUADRATURE

                    call gq2D_to_nodes_acx(gq2d,dudxn,jvel%dxx(:,:,kmid),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                    dudx_aa = sum(dudxn*gq2d%wt)/gq2d%wt_tot

                    call gq2D_to_nodes_acy(gq2d,dvdyn,jvel%dyy(:,:,kmid),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                    dvdy_aa = sum(dvdyn*gq2d%wt)/gq2d%wt_tot

else 
    ! 3D QUADRATURE

                    km1 = kmid-1
                    kp1 = kmid+1
                    if (kmid .eq. 1)  km1 = 1
                    if (kmid .eq. nz_aa) kp1 = nz_aa

                    if (kmid .gt. 1) then
                        dz0 = H_ice(i,j)*(zeta_aa(kmid) - zeta_aa(km1))
                    else
                        dz0 = H_ice(i,j)*(zeta_aa(2) - zeta_aa(1))
                    end if

                    call gq3D_to_nodes_acx(gq3D,dudxn8,jvel%dxx,dx,dy,dz0,dz1,i,j,kmid,im1,ip1,jm1,jp1,km1,kp1,act=act_acx)
                    dudx_aa = sum(dudxn8*gq3D%wt)/gq3D%wt_tot

                    call gq3D_to_nodes_acy(gq3D,dvdyn8,jvel%dyy,dx,dy,dz0,dz1,i,j,kmid,im1,ip1,jm1,jp1,km1,kp1,act=act_acy)
                    dvdy_aa = sum(dvdyn8*gq3D%wt)/gq3D%wt_tot

end if 

                    ! Calculate vertical velocity of this layer
                    ! (Greve and Blatter, 2009, Eq. 5.95)
                    uz(i,j,k) = uz(i,j,k-1) - H_now*(zeta_ac(k)-zeta_ac(k-1))*(dudx_aa+dvdy_aa)

                    ! Apply correction to match kinematic boundary condition at surface 
                    !uz(i,j,k) = uz(i,j,k) - zeta_ac(k)*(uz(i,j,k)-uz_srf)

                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 

                end do 

                ! === Also calculate the sigma-coordinate advective vertical velocity uz_star ===
                ! uz_star = uz + (coordinate-transform corrections) = H*dzeta/dt following a parcel.
                ! It is the vertical velocity that MUST be paired with horizontal advection taken at
                ! constant sigma (constant zeta level). It is a general kinematic field, not a
                ! thermodynamics internal: thermodynamics is presently its sole consumer, but any
                ! constant-sigma scalar advection (e.g. age/tracer transport) requires uz_star in
                ! place of the true vertical velocity uz. Computed here, alongside uz and owned by
                ! dyn, to reuse the surface/base geometry derivatives already formed above for uz.
                ! Note: the 1/H factor is deferred to each consumer's advection step (see below).
                
                do k = 1, nz_ac 

                    ! Get the centered horizontal velocity of box on vertical ac-nodes at the level k 
                    ! ajr: Given that the correction is applied to uz, which is defined on 
                    ! ac-nodes, it seems the correction should also be calculated on ac-nodes.
                    ! Note: nz_ac = nz_aa + 1
    
                    if (k .eq. 1) then 
                        kup = k 
                        kdn = k 
                    else if (k .eq. nz_ac) then 
                        kup = k-1 
                        kdn = k-1 
                    else
                        kup = k 
                        kdn = k-1 
                    end if
                    
                    call gq2D_to_nodes_acx(gq2d,uxn_up,ux(:,:,kup),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                    call gq2D_to_nodes_acx(gq2d,uxn_dn,ux(:,:,kdn),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                    uxn = 0.5_wp*(uxn_up+uxn_dn)
                    ux_aa = sum(uxn*gq2d%wt)/gq2d%wt_tot
                    
                    call gq2D_to_nodes_acy(gq2d,uyn_up,uy(:,:,kup),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                    call gq2D_to_nodes_acy(gq2d,uyn_dn,uy(:,:,kdn),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                    uyn = 0.5_wp*(uyn_up+uyn_dn)
                    uy_aa = sum(uyn*gq2d%wt)/gq2d%wt_tot
                    
                    ! Take zeta directly at vertical cell edge where uz is calculated
                    ! (this is also where ux_aa and uy_aa are calculated above)
                    zeta_now = zeta_ac(k)

                    ! Calculate sigma-coordinate derivative correction factors
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132, 
                    !  also shown in 1D with Eq. 5.145)

                    ! Note: not dividing by H here, since the thermodynamics advection step
                    ! is in terms of d/dz, not d/dzeta. Units of uz_star should still be m/m
                    ! c_x = - H_inv * ( (1.0-zeta_now)*dzbdx_aa  + zeta_now*dzsdx_aa )
                    ! c_y = - H_inv * ( (1.0-zeta_now)*dzbdy_aa  + zeta_now*dzsdy_aa )
                    ! c_t = - H_inv * ( (1.0-zeta_now)*dzbdt_now + zeta_now*dzsdt_now )
                    ! c_z = H_inv 
                    c_x = -( (1.0-zeta_now)*dzbdx_aa  + zeta_now*dzsdx_aa )
                    c_y = -( (1.0-zeta_now)*dzbdy_aa  + zeta_now*dzsdy_aa )
                    c_t = -( (1.0-zeta_now)*dzbdt_now + zeta_now*dzsdt_now )
                    c_z = 1.0 

                    ! Calculate adjusted vertical velocity for advection 
                    ! of this layer - units of m/m
                    ! (e.g., Greve and Blatter, 2009, Eq. 5.148)
                    uz_star(i,j,k) = ux_aa*c_x + uy_aa*c_y + uz(i,j,k)*c_z + c_t

                    if (abs(uz_star(i,j,k)) .lt. TOL_UNDERFLOW) uz_star(i,j,k) = 0.0_wp

                end do 
                
            else 
                ! No ice here, set vertical velocity equal to negative accum and bedrock change 

                do k = 1, nz_ac

                    uz(i,j,k) = dzbdt_now - max(smb(i,j),0.0)
                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 

                    uz_star(i,j,k) = uz(i,j,k)

               end do 

            end if 

        end do 
        end do 
        !$omp end do 
        !$omp end parallel

        return 

    end subroutine calc_uz_3D_jac

    subroutine calc_uz_3D(uz,uz_star,ux,uy,H_ice,f_ice,f_grnd,smb,bmb,dzbdt,dzsdt, &
                                    dzsdx,dzsdy,dzbdx,dzbdy,zeta_aa,zeta_ac,dx,dy,use_bmb,boundaries)
        ! Following algorithm outlined by the Glimmer ice sheet model:
        ! https://www.geos.ed.ac.uk/~mhagdorn/glide/glide-doc/glimmer_htmlse9.html#x17-660003.1.5

        ! Note: dzbdt and dzsdt are the kinematic rates of the base and surface of the
        ! ice column (tpo%now%dzbdt_kin, dzsdt_kin): vertical column change plus bedrock
        ! and sea-level rates, without lateral changes (calving, front advance, removals).

        implicit none 

        real(wp), intent(OUT) :: uz(:,:,:)          ! nx,ny,nz_ac
        real(wp), intent(OUT) :: uz_star(:,:,:)     ! nx,ny,nz_ac
        real(wp), intent(IN)  :: ux(:,:,:)          ! nx,ny,nz_aa
        real(wp), intent(IN)  :: uy(:,:,:)          ! nx,ny,nz_aa
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: smb(:,:) 
        real(wp), intent(IN)  :: bmb(:,:) 
        real(wp), intent(IN)  :: dzbdt(:,:)         ! [m/a] Kinematic rate of the column base
        real(wp), intent(IN)  :: dzsdt(:,:) 
        real(wp), intent(IN)  :: dzsdx(:,:) 
        real(wp), intent(IN)  :: dzsdy(:,:) 
        real(wp), intent(IN)  :: dzbdx(:,:) 
        real(wp), intent(IN)  :: dzbdy(:,:) 
        real(wp), intent(IN)  :: zeta_aa(:)    ! z-coordinate, aa-nodes 
        real(wp), intent(IN)  :: zeta_ac(:)    ! z-coordinate, ac-nodes  
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        logical,  intent(IN)  :: use_bmb
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac
        integer :: im1, ip1, jm1, jp1   
        integer :: kup, kdn
        real(wp) :: f_bmb
        real(wp) :: H_now
        real(wp) :: H_inv
        real(wp) :: dzbdx_aa
        real(wp) :: dzbdy_aa
        real(wp) :: dzsdx_aa
        real(wp) :: dzsdy_aa 
        real(wp) :: dudx_aa
        real(wp) :: dvdy_aa
        real(wp) :: dudz_aa 
        real(wp) :: dvdz_aa 
        real(wp) :: ux_aa 
        real(wp) :: uy_aa 
        real(wp) :: uz_grid 
        real(wp) :: uz_srf 
        real(wp) :: zeta_now 
        real(wp) :: c_x 
        real(wp) :: c_y 
        real(wp) :: c_t 

        real(wp) :: dzsdtn(4)
        real(wp) :: dzbdxn(4)
        real(wp) :: dzbdyn(4)
        real(wp) :: dzsdxn(4)
        real(wp) :: dzsdyn(4)
        real(wp) :: dudxn(4) 
        real(wp) :: dvdyn(4) 
        real(wp) :: uxn_up(4) 
        real(wp) :: uxn_dn(4) 
        real(wp) :: uxn(4) 
        real(wp) :: uyn_up(4) 
        real(wp) :: uyn_dn(4)
        real(wp) :: uyn(4)  
        real(wp) :: dudzn(4) 
        real(wp) :: dvdzn(4) 
        
        real(wp) :: dzsdt_now
        real(wp) :: dzbdt_now 

        real(wp), allocatable :: dudx(:,:,:)
        real(wp), allocatable :: dvdy(:,:,:)
        real(wp), allocatable :: dudy(:,:)
        real(wp), allocatable :: dvdx(:,:)
        
        
        logical, allocatable :: act_acx(:,:), act_acy(:,:)

        type(gq2D_class) :: gq2D, gq2D_global
        
        integer  :: BC

        ! Initialize gaussian quadrature calculations
        call gq2D_init(gq2D_global)

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1) 
        
        ! Allocate arrays 
        allocate(dudx(nx,ny,nz_aa))
        allocate(dvdy(nx,ny,nz_aa))
        allocate(dudy(nx,ny))
        allocate(dvdx(nx,ny))
        
        ! Initialize vertical velocity to zero 
        uz = 0.0 

        ! Define switch for bmb
        if (use_bmb) then 
            f_bmb = 1.0 
        else 
            f_bmb = 0.0 
        end if 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Velocity faces with a solution (next to an ice-covered cell); the
        ! others hold zero velocity and do not enter the quadrature means
        allocate(act_acx(nx,ny),act_acy(nx,ny))
        call calc_active_faces(act_acx,act_acy,f_ice,BC)

        ! First calculate horizontal strain rates at each layer for later use,
        ! with no correction factor for sigma-transformation.
        ! Note: we only need dudx and dvdy, but routine also calculate cross terms, which will not be used.

        !$omp parallel do private(k,dudy,dvdx)
        do k = 1, nz_aa
            call calc_strain_rate_horizontal_2D(dudx(:,:,k),dudy,dvdx,dvdy(:,:,k),ux(:,:,k),uy(:,:,k),f_ice,dx,dy,boundaries)
        end do
        !$omp end parallel do

        ! Next, calculate vertical velocity at each point through the column

        !$omp parallel private(i,j,k,im1,ip1,jm1,jp1,dzsdt_now,dzbdt_now,H_now,H_inv) &
        !$omp& private(dzbdxn,dzbdx_aa,dzbdyn,dzbdy_aa,dzsdxn,dzsdx_aa,dzsdyn,dzsdy_aa,uxn,ux_aa,uyn,uy_aa) &
        !$omp& private(uz_grid,dudxn,dudx_aa,dvdyn,dvdy_aa) &
        !$omp& private(kup,kdn,uxn_up,uxn_dn,dudzn,dudz_aa,uyn_up,uyn_dn,dvdzn,dvdz_aa,zeta_now,c_x,c_y,c_t,gq2d) &
        !$omp& shared(gq2D_global)
        gq2D = gq2D_global

        !$omp do collapse(2)
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! Kinematic rates of the column surface and base
            dzsdt_now = dzsdt(i,j)
            dzbdt_now = dzbdt(i,j)

            if (f_ice(i,j) .eq. 1.0) then

                H_now  = H_ice(i,j) 
                H_inv = 1.0/H_now 

                ! Get the centered ice-base gradient
                call gq2D_to_nodes_acx(gq2d,dzbdxn,dzbdx,dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                dzbdx_aa = sum(dzbdxn*gq2d%wt)/gq2d%wt_tot
                
                call gq2D_to_nodes_acy(gq2d,dzbdyn,dzbdy,dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                dzbdy_aa = sum(dzbdyn*gq2d%wt)/gq2d%wt_tot

                ! Get the centered surface gradient
                call gq2D_to_nodes_acx(gq2d,dzsdxn,dzsdx,dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                dzsdx_aa = sum(dzsdxn*gq2d%wt)/gq2d%wt_tot
                
                call gq2D_to_nodes_acy(gq2d,dzsdyn,dzsdy,dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                dzsdy_aa = sum(dzsdyn*gq2d%wt)/gq2d%wt_tot
                
                ! Get the aa-node centered horizontal velocity at the base
                call gq2D_to_nodes_acx(gq2d,uxn,ux(:,:,1),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                ux_aa = sum(uxn*gq2d%wt)/gq2d%wt_tot
                
                call gq2D_to_nodes_acy(gq2d,uyn,uy(:,:,1),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                uy_aa = sum(uyn*gq2d%wt)/gq2d%wt_tot
                
                ! Determine grid vertical velocity at the base due to sigma-coordinates 
                ! Glimmer, Eq. 3.35 
                ! ajr, 2020-01-27, untested:::
!                 uz_grid = dzsdt(i,j) + (ux_aa*dzsdx_aa + uy_aa*dzsdy_aa) &
!                             - ( (1.0_wp-zeta_ac(1))*dHdt(i,j) + ux_aa*dHdx_aa + uy_aa*dHdy_aa )
                uz_grid = 0.0_wp 

                ! ===================================================================
                ! Greve and Blatter (2009) style:

                ! Determine basal vertical velocity for this grid point 
                ! Following Eq. 5.31 of Greve and Blatter (2009)
                uz(i,j,1) = dzbdt_now + uz_grid + f_bmb*bmb(i,j) + ux_aa*dzbdx_aa + uy_aa*dzbdy_aa
                if (abs(uz(i,j,1)) .lt. TOL_UNDERFLOW) uz(i,j,1) = 0.0_wp 

                ! Determine surface vertical velocity following kinematic boundary condition 
                ! Glimmer, Eq. 3.10 [or Folwer, Chpt 10, Eq. 10.8]
                !uz_srf = dzsdt(i,j) + ux_aa*dzsdx_aa + uy_aa*dzsdy_aa - smb(i,j) 
                
                ! Integrate upward to each point above base until just below surface is reached 
                ! Integrate on vertical ac-nodes (ie, vertical cell borders between aa-node centers)
                do k = 2, nz_ac 

                    ! Calculate sigma-coordinate derivative correction factors for this layer
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132, 
                    !  also shown in 1D with Eq. 5.145)
                    ! See src/physics/deformation.f90::calc_jacobian_vel_3D()

                    ! Take zeta at the center of the cell below the current vertical ac boundary
                    ! (this is also where dudz_aa and dvz_aa are calculated above)
                    zeta_now = zeta_aa(k-1)

                    c_x = -H_inv * ( (1.0-zeta_now)*dzbdx_aa + zeta_now*dzsdx_aa )
                    c_y = -H_inv * ( (1.0-zeta_now)*dzbdy_aa + zeta_now*dzsdy_aa )
                    
                    ! Get dudz/dvdz values at vertical aa-nodes, in order
                    ! to vertically integrate each cell up to ac-node border.
                    ! Note: nz_ac = nz_aa + 1
                    if (k .eq. 2) then
                        kup = k 
                        kdn = k-1 
                    else if (k .eq. nz_ac) then 
                        kup = k-1 
                        kdn = k-2 
                    else 
                        ! Centered on k-1 
                        kup = k 
                        kdn = k-2
                    end if

                    call gq2D_to_nodes_acx(gq2d,uxn_up,ux(:,:,kup),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                    call gq2D_to_nodes_acx(gq2d,uxn_dn,ux(:,:,kdn),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                    dudzn = (uxn_up - uxn_dn) / (zeta_aa(kup)-zeta_aa(kdn))
                    dudz_aa = sum(dudzn*gq2d%wt)/gq2d%wt_tot

                    call gq2D_to_nodes_acy(gq2d,uyn_up,uy(:,:,kup),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                    call gq2D_to_nodes_acy(gq2d,uyn_dn,uy(:,:,kdn),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                    dvdzn = (uyn_up - uyn_dn) / (zeta_aa(kup)-zeta_aa(kdn))
                    dvdz_aa = sum(dvdzn*gq2d%wt)/gq2d%wt_tot

                    ! Calculate sigma-corrected derivatives
                    call gq2D_to_nodes_acx(gq2d,dudxn,dudx(:,:,k-1),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                    dudx_aa = sum(dudxn*gq2d%wt)/gq2d%wt_tot  +  c_x*dudz_aa 

                    call gq2D_to_nodes_acy(gq2d,dvdyn,dvdy(:,:,k-1),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                    dvdy_aa = sum(dvdyn*gq2d%wt)/gq2d%wt_tot  +  c_y*dvdz_aa 

                    ! Calculate vertical velocity of this layer
                    ! (Greve and Blatter, 2009, Eq. 5.95)
                    uz(i,j,k) = uz(i,j,k-1) - H_now*(zeta_ac(k)-zeta_ac(k-1))*(dudx_aa+dvdy_aa)

                    ! Apply correction to match kinematic boundary condition at surface 
                    !uz(i,j,k) = uz(i,j,k) - zeta_ac(k)*(uz(i,j,k)-uz_srf)

                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 
                    
                end do 


                ! === Also calculate the sigma-coordinate advective vertical velocity uz_star ===
                ! uz_star = uz + (coordinate-transform corrections) = H*dzeta/dt following a parcel.
                ! It is the vertical velocity that MUST be paired with horizontal advection taken at
                ! constant sigma (constant zeta level). It is a general kinematic field, not a
                ! thermodynamics internal: thermodynamics is presently its sole consumer, but any
                ! constant-sigma scalar advection (e.g. age/tracer transport) requires uz_star in
                ! place of the true vertical velocity uz. Computed here, alongside uz and owned by
                ! dyn, to reuse the surface/base geometry derivatives already formed above for uz.
                ! Note: the 1/H factor is deferred to each consumer's advection step (see below).
                
                do k = 1, nz_ac 

                    ! Get the centered horizontal velocity of box on vertical ac-nodes at the level k 
                    ! ajr: Given that the correction is applied to uz, which is defined on 
                    ! ac-nodes, it seems the correction should also be calculated on ac-nodes.
                    ! Note: nz_ac = nz_aa + 1
    
                    if (k .eq. 1) then 
                        kup = k 
                        kdn = k 
                    else if (k .eq. nz_ac) then 
                        kup = k-1 
                        kdn = k-1 
                    else
                        kup = k 
                        kdn = k-1 
                    end if
                    
                    call gq2D_to_nodes_acx(gq2d,uxn_up,ux(:,:,kup),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                    call gq2D_to_nodes_acx(gq2d,uxn_dn,ux(:,:,kdn),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acx)
                    uxn = 0.5_wp*(uxn_up+uxn_dn)
                    ux_aa = sum(uxn*gq2d%wt)/gq2d%wt_tot
                    
                    call gq2D_to_nodes_acy(gq2d,uyn_up,uy(:,:,kup),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                    call gq2D_to_nodes_acy(gq2d,uyn_dn,uy(:,:,kdn),dx,dy,i,j,im1,ip1,jm1,jp1,act=act_acy)
                    uyn = 0.5_wp*(uyn_up+uyn_dn)
                    uy_aa = sum(uyn*gq2d%wt)/gq2d%wt_tot
                    
                    ! Take zeta directly at vertical cell edge where uz is calculated
                    ! (this is also where ux_aa and uy_aa are calculated above)
                    zeta_now = zeta_ac(k)

                    ! Calculate sigma-coordinate derivative correction factors
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132, 
                    !  also shown in 1D with Eq. 5.145)

                    ! Note: not dividing by H here, since this is done in the thermodynamics advection step
                    c_x = -( (1.0-zeta_now)*dzbdx_aa  + zeta_now*dzsdx_aa )
                    c_y = -( (1.0-zeta_now)*dzbdy_aa  + zeta_now*dzsdy_aa )
                    c_t = -( (1.0-zeta_now)*dzbdt_now + zeta_now*dzsdt_now )
                    
                    ! Calculate adjusted vertical velocity for advection 
                    ! of this layer
                    ! (e.g., Greve and Blatter, 2009, Eq. 5.148)
                    uz_star(i,j,k) = uz(i,j,k) + ux_aa*c_x + uy_aa*c_y + c_t 

                    if (abs(uz_star(i,j,k)) .lt. TOL_UNDERFLOW) uz_star(i,j,k) = 0.0_wp
                    
                end do 
                
            else 
                ! No ice here, set vertical velocity equal to negative accum and bedrock change 

                do k = 1, nz_ac

                    uz(i,j,k) = dzbdt_now - max(smb(i,j),0.0)
                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 

                    uz_star(i,j,k) = uz(i,j,k)

               end do 

            end if 

        end do 
        end do 
        !$omp end do 
        !$omp end parallel

        return 

    end subroutine calc_uz_3D

    subroutine calc_uz_3D_aa(uz,uz_star,ux,uy,H_ice,f_ice,f_grnd,z_bed,z_srf,smb,bmb,dzbdt,dzsdt, &
                                            dzsdx,dzsdy,dzbdx,dzbdy,zeta_aa,zeta_ac,dx,dy,use_bmb,boundaries)
        ! Following algorithm outlined by the Glimmer ice sheet model:
        ! https://www.geos.ed.ac.uk/~mhagdorn/glide/glide-doc/glimmer_htmlse9.html#x17-660003.1.5

        ! Note: dzbdt and dzsdt are the kinematic rates of the base and surface of the
        ! ice column (tpo%now%dzbdt_kin, dzsdt_kin): vertical column change plus bedrock
        ! and sea-level rates, without lateral changes (calving, front advance, removals).

        implicit none 

        real(wp), intent(OUT) :: uz(:,:,:)          ! nx,ny,nz_ac
        real(wp), intent(OUT) :: uz_star(:,:,:)     ! nx,ny,nz_ac
        real(wp), intent(IN)  :: ux(:,:,:)          ! nx,ny,nz_aa
        real(wp), intent(IN)  :: uy(:,:,:)          ! nx,ny,nz_aa
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: z_bed(:,:) 
        real(wp), intent(IN)  :: z_srf(:,:) 
        real(wp), intent(IN)  :: smb(:,:) 
        real(wp), intent(IN)  :: bmb(:,:) 
        real(wp), intent(IN)  :: dzbdt(:,:)         ! [m/a] Kinematic rate of the column base
        real(wp), intent(IN)  :: dzsdt(:,:)
        real(wp), intent(IN)  :: dzsdx(:,:) 
        real(wp), intent(IN)  :: dzsdy(:,:) 
        real(wp), intent(IN)  :: dzbdx(:,:) 
        real(wp), intent(IN)  :: dzbdy(:,:) 
        real(wp), intent(IN)  :: zeta_aa(:)    ! z-coordinate, aa-nodes 
        real(wp), intent(IN)  :: zeta_ac(:)    ! z-coordinate, ac-nodes  
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        logical,  intent(IN)  :: use_bmb
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, k, nx, ny, nz_aa, nz_ac
        integer  :: im1, ip1, jm1, jp1   
        real(wp) :: f_bmb 
        real(wp) :: H_now
        real(wp) :: H_inv
        real(wp) :: dzbdx_aa
        real(wp) :: dzbdy_aa
        real(wp) :: dzsdx_aa
        real(wp) :: dzsdy_aa
        real(wp) :: duxdx_aa
        real(wp) :: duydy_aa
        real(wp) :: duxdz_aa 
        real(wp) :: duydz_aa
        real(wp) :: duxdx_now 
        real(wp) :: duydy_now 
        real(wp) :: ux_aa 
        real(wp) :: uy_aa 
        real(wp) :: uz_grid 
        real(wp) :: uz_srf 
        real(wp) :: zeta_now 
        real(wp) :: c_x 
        real(wp) :: c_y 
        real(wp) :: c_t 

        real(wp) :: dzsdt_now
        real(wp) :: dzbdt_now 

        
        integer  :: BC

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1) 

        ! Initialize vertical velocity to zero 
        uz = 0.0 

        ! Define switch for bmb
        if (use_bmb) then 
            f_bmb = 1.0 
        else 
            f_bmb = 0.0 
        end if 
        
        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Next, calculate velocity 

        !$omp parallel do collapse(2) private(i,j,k,im1,ip1,jm1,jp1,dzsdt_now,dzbdt_now,H_now,H_inv) &
        !$omp& private(dzbdx_aa,dzbdy_aa,dzsdx_aa,dzsdy_aa,ux_aa,uy_aa) &
        !$omp& private(uz_grid,duxdx_aa,duydy_aa,duxdz_aa,duydz_aa,duxdx_now,duydy_now) &
        !$omp& private(zeta_now,c_x,c_y,c_t)
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! Kinematic rates of the column surface and base
            dzsdt_now = dzsdt(i,j)
            dzbdt_now = dzbdt(i,j)

            if (f_ice(i,j) .eq. 1.0) then

                ! Get weighted ice thickness for stability
!                 H_now = (4.0*H_ice(i,j) + 2.0*(H_ice(im1,j)+H_ice(ip1,j)+H_ice(i,jm1)+H_ice(i,jp1))) / 16.0 &
!                       + (H_ice(im1,jm1)+H_ice(ip1,jm1)+H_ice(ip1,jp1)+H_ice(im1,jp1)) / 16.0 

                H_now  = H_ice(i,j) 
                H_inv = 1.0/H_now 

                ! Get the centered bedrock gradient 
                dzbdx_aa = 0.5*(dzbdx(i,j)+dzbdx(im1,j))
                dzbdy_aa = 0.5*(dzbdy(i,j)+dzbdy(i,jm1))

                ! Get the centered surface gradient 
                dzsdx_aa = 0.5*(dzsdx(i,j)+dzsdx(im1,j))
                dzsdy_aa = 0.5*(dzsdy(i,j)+dzsdy(i,jm1))

                ! Get the centered horizontal velocity at the base
                ux_aa = 0.5_wp* (ux(im1,j,1) + ux(i,j,1))
                uy_aa = 0.5_wp* (uy(i,jm1,1) + uy(i,j,1))
                
                ! Determine grid vertical velocity at the base due to sigma-coordinates 
                ! Glimmer, Eq. 3.35 
                ! ajr, 2020-01-27, untested:::
!                 uz_grid = dzsdt(i,j) + (ux_aa*dzsdx_aa + uy_aa*dzsdy_aa) &
!                             - ( (1.0_wp-zeta_ac(1))*dHdt(i,j) + ux_aa*dHdx_aa + uy_aa*dHdy_aa )
                uz_grid = 0.0_wp 

                ! ===================================================================
                ! Greve and Blatter (2009) style:

                ! Determine basal vertical velocity for this grid point 
                ! Following Eq. 5.31 of Greve and Blatter (2009)
                uz(i,j,1) = dzbdt_now + uz_grid + f_bmb*bmb(i,j) + ux_aa*dzbdx_aa + uy_aa*dzbdy_aa
                if (abs(uz(i,j,1)) .lt. TOL_UNDERFLOW) uz(i,j,1) = 0.0_wp 

                ! Determine surface vertical velocity following kinematic boundary condition 
                ! Glimmer, Eq. 3.10 [or Folwer, Chpt 10, Eq. 10.8]
                !uz_srf = dzsdt(i,j) + ux_aa*dzsdx_aa + uy_aa*dzsdy_aa - smb(i,j) 
                
                ! Integrate upward to each point above base until just below surface is reached 
                do k = 2, nz_ac 

                    ! Greve and Blatter (2009), Eq. 5.72
                    ! Bueler and Brown  (2009), Eq. 4
                    
                    duxdx_aa  = (ux(i,j,k-1)   - ux(im1,j,k-1)  )/dx
                    duydy_aa  = (uy(i,j,k-1)   - uy(i,jm1,k-1)  )/dy

                    ! Note: nz_ac = nz_aa + 1
                    if (k .eq. 2) then
                        duxdz_aa  = ( 0.5*(ux(i,j,k)+ux(im1,j,k)) &
                            - 0.5*(ux(i,j,k-1)+ux(im1,j,k-1)) ) / (zeta_aa(k)-zeta_aa(k-1))
                        duydz_aa  = ( 0.5*(uy(i,j,k)+uy(i,jm1,k)) &
                            - 0.5*(uy(i,j,k-1)+uy(i,jm1,k-1)) ) / (zeta_aa(k)-zeta_aa(k-1))
                    else if (k .eq. nz_ac) then
                        duxdz_aa  = ( 0.5*(ux(i,j,k-1)+ux(im1,j,k-1)) &
                            - 0.5*(ux(i,j,k-2)+ux(im1,j,k-2)) ) / (zeta_aa(k-1)-zeta_aa(k-2))
                        duydz_aa  = ( 0.5*(uy(i,j,k-1)+uy(i,jm1,k-1)) &
                            - 0.5*(uy(i,j,k-2)+uy(i,jm1,k-2)) ) / (zeta_aa(k-1)-zeta_aa(k-2))
                    else
                        ! Centered difference vertically, centered on k-1
                        duxdz_aa  = ( 0.5*(ux(i,j,k)+ux(im1,j,k)) &
                            - 0.5*(ux(i,j,k-2)+ux(im1,j,k-2)) ) / (zeta_aa(k)-zeta_aa(k-2))
                        duydz_aa  = ( 0.5*(uy(i,j,k)+uy(i,jm1,k)) &
                            - 0.5*(uy(i,j,k-2)+uy(i,jm1,k-2)) ) / (zeta_aa(k)-zeta_aa(k-2))
                    end if 

                    ! Calculate sigma-coordinate derivative correction factors
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132)
                    c_x = -H_inv * ( (1.0-zeta_ac(k))*dzbdx_aa + zeta_ac(k)*dzsdx_aa )
                    c_y = -H_inv * ( (1.0-zeta_ac(k))*dzbdy_aa + zeta_ac(k)*dzsdy_aa )

                    ! Calculate derivatives
                    duxdx_now = duxdx_aa + c_x*duxdz_aa 
                    duydy_now = duydy_aa + c_y*duydz_aa 
                    
                    ! Calculate vertical velocity of this layer
                    ! (Greve and Blatter, 2009, Eq. 5.95)
                    uz(i,j,k) = uz(i,j,k-1) & 
                        - H_now*(zeta_ac(k)-zeta_ac(k-1))*(duxdx_now+duydy_now)

                    ! Apply correction to match kinematic boundary condition at surface 
                    !uz(i,j,k) = uz(i,j,k) - zeta_ac(k)*(uz(i,j,k)-uz_srf)

                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 
                    
                end do 
                
                ! === Also calculate the sigma-coordinate advective vertical velocity uz_star ===
                ! uz_star = uz + (coordinate-transform corrections) = H*dzeta/dt following a parcel.
                ! It is the vertical velocity that MUST be paired with horizontal advection taken at
                ! constant sigma (constant zeta level). It is a general kinematic field, not a
                ! thermodynamics internal: thermodynamics is presently its sole consumer, but any
                ! constant-sigma scalar advection (e.g. age/tracer transport) requires uz_star in
                ! place of the true vertical velocity uz. Computed here, alongside uz and owned by
                ! dyn, to reuse the surface/base geometry derivatives already formed above for uz.
                ! Note: the 1/H factor is deferred to each consumer's advection step (see below).
                
                do k = 1, nz_ac 

                    ! Get the centered horizontal velocity of box
                    ! on vertical ac-nodes at the level k 
                    ! ajr: Given that the correction is 
                    ! applied to uz, which is defined on 
                    ! ac-nodes, it seems the correction
                    ! should also be calculated on ac-nodes.
                    ! Note: nz_ac = nz_aa + 1
    
                    if (k .eq. 1) then 
                        ux_aa = 0.5_wp* (ux(im1,j,k) + ux(i,j,k))
                        uy_aa = 0.5_wp* (uy(i,jm1,k) + uy(i,j,k))
                    else if (k .eq. nz_ac) then  
                        ux_aa = 0.5_wp* (ux(im1,j,k-1) + ux(i,j,k-1))
                        uy_aa = 0.5_wp* (uy(i,jm1,k-1) + uy(i,j,k-1))
                    else
                        ux_aa = 0.25_wp* (ux(im1,j,k) + ux(i,j,k) + ux(im1,j,k-1) + ux(i,j,k-1))
                        uy_aa = 0.25_wp* (uy(i,jm1,k) + uy(i,j,k) + uy(i,jm1,k-1) + uy(i,j,k-1))
                    end if 

                    ! Take zeta directly at vertical cell edge where uz is calculated
                    ! (this is also where ux_aa and uy_aa are calculated above)
                    zeta_now = zeta_ac(k)

                    ! Calculate sigma-coordinate derivative correction factors
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132, 
                    !  also shown in 1D with Eq. 5.145)

                    ! Note: not dividing by H here, since this is done in the thermodynamics advection step
                    c_x = -ux_aa * ( (1.0_wp-zeta_now)*dzbdx_aa  + zeta_now*dzsdx_aa )
                    c_y = -uy_aa * ( (1.0_wp-zeta_now)*dzbdy_aa  + zeta_now*dzsdy_aa )
                    c_t =         -( (1.0_wp-zeta_now)*dzbdt_now + zeta_now*dzsdt_now )

                    ! Calculate adjusted vertical velocity for advection 
                    ! of this layer
                    ! (e.g., Greve and Blatter, 2009, Eq. 5.148)
                    uz_star(i,j,k) = uz(i,j,k) + c_x + c_y + c_t 

                    if (abs(uz_star(i,j,k)) .lt. TOL_UNDERFLOW) uz_star(i,j,k) = 0.0_wp 
                    
                end do 
                
            else 
                ! No ice here, set vertical velocity equal to negative accum and bedrock change 

                do k = 1, nz_ac 

                    uz(i,j,k) = dzbdt_now - max(smb(i,j),0.0)
                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 

                    uz_star(i,j,k) = uz(i,j,k)

               end do 

            end if 

        end do 
        end do 
        !$omp end parallel do 

        return 

    end subroutine calc_uz_3D_aa

    subroutine calc_driving_stress(taud_acx,taud_acy,H_ice,f_ice,dzsdx,dzsdy,dx,taud_lim,rho_ice,g,boundaries)
        ! Calculate driving stress on staggered grid points
        ! Units: taud [Pa] == [kg m-1 s-2]
        ! taud = rho_ice*g*H_ice*dzs/dx

        ! Note: interpolation to Ab nodes no longer used here.

        ! Note: use parameter taud_lim to limit maximum allowed driving stress magnitude applied in the model.
        ! Should be an extreme value. eg, if dzdx = taud_lim / (rho*g*H), 
        ! then for taud_lim=5e5 and H=3000m, dzdx = 2e5 / (910*9.81*3000) = 0.02, 
        ! which is a rather steep slope for a shallow-ice model.

        implicit none 

        real(wp), intent(OUT) :: taud_acx(:,:)      ! [Pa]
        real(wp), intent(OUT) :: taud_acy(:,:)      ! [Pa]
        real(wp), intent(IN)  :: H_ice(:,:)         ! [m]
        real(wp), intent(IN)  :: f_ice(:,:)         ! [--]
        real(wp), intent(IN)  :: dzsdx(:,:)         ! [--]
        real(wp), intent(IN)  :: dzsdy(:,:)         ! [--]
        real(wp), intent(IN)  :: dx                 ! [m] 
        real(wp), intent(IN)  :: taud_lim           ! [Pa]
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: g 
        character(len=*), intent(IN) :: boundaries  ! Boundary conditions to apply 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer    :: im1, ip1, jm1, jp1 
        real(wp) :: dy, rhog 
        real(wp) :: H_mid
        real(wp) :: taud_mean 
        real(wp) :: taud_eps

        real(wp) :: taud_eps_lim = 0.4

        real(wp), allocatable :: taud_acx_0(:,:)
        real(wp), allocatable :: taud_acy_0(:,:)

        integer  :: BC

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        allocate(taud_acx_0(nx,ny))
        allocate(taud_acy_0(nx,ny))
        
        ! Define shortcut parameter 
        rhog = rho_ice * g 

        ! Assume grid resolution is symmetrical 
        dy = dx 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,H_mid)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! x-direction
            if (f_ice(i,j) .eq. 1.0_wp .and. f_ice(ip1,j) .lt. 1.0_wp) then 
                H_mid = 0.5_wp*(H_ice(i,j)+0.0_wp)
            else if (f_ice(i,j) .lt. 1.0_wp .and. f_ice(ip1,j) .eq. 1.0_wp) then 
                H_mid = 0.5_wp*(0.0_wp+H_ice(ip1,j))
            else  
                H_mid = 0.5_wp*(H_ice(i,j)+H_ice(ip1,j)) 
            end if
            taud_acx(i,j) = rhog * H_mid * dzsdx(i,j) 
            
            ! Apply limit
            if (taud_acx(i,j) .gt.  taud_lim) taud_acx(i,j) =  taud_lim 
            if (taud_acx(i,j) .lt. -taud_lim) taud_acx(i,j) = -taud_lim 
            
            ! y-direction 
            if (f_ice(i,j) .eq. 1.0_wp .and. f_ice(i,jp1) .lt. 1.0_wp) then 
                H_mid = 0.5_wp*(H_ice(i,j)+0.0_wp)
            else if (f_ice(i,j) .lt. 1.0_wp .and. f_ice(i,jp1) .eq. 1.0_wp) then 
                H_mid = 0.5_wp*(0.0_wp+H_ice(i,jp1))
            else  
                H_mid = 0.5_wp*(H_ice(i,j)+H_ice(i,jp1)) 
            end if
            taud_acy(i,j) = rhog * H_mid * dzsdy(i,j) 

            ! Apply limit
            if (taud_acy(i,j) .gt.  taud_lim) taud_acy(i,j) =  taud_lim 
            if (taud_acy(i,j) .lt. -taud_lim) taud_acy(i,j) = -taud_lim 
            
        end do
        end do 
        !$omp end parallel do

        ! SPECIAL CASES: edge cases for stability...

if (.FALSE.) then
        taud_acx_0 = taud_acx 
        taud_acy_0 = taud_acy 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! === x-direction ===

            if (f_ice(i,j) .eq. 1.0_wp .and. f_ice(ip1,j) .eq. 1.0_wp) then 

                taud_mean = 0.5*(taud_acx_0(im1,j)+taud_acx_0(ip1,j))
                taud_eps  = (taud_acx_0(i,j) - taud_mean) / taud_lim

                if (abs(taud_acx_0(i,j)) .eq. taud_lim .and. abs(taud_eps) .gt. taud_eps_lim) then 
                    ! Set driving stress equal to weighted average of neighbors
                    taud_acx(i,j) = 0.20*taud_acx_0(i,j) + 0.40*taud_acx_0(im1,j)+ 0.40*taud_acx_0(ip1,j)
                end if 

            end if

            ! === y-direction ===

            if (f_ice(i,j) .eq. 1.0_wp .and. f_ice(i,jp1) .eq. 1.0_wp) then 

                taud_mean = 0.5*(taud_acy_0(i,jm1)+taud_acy_0(i,jp1))
                taud_eps  = (taud_acy_0(i,j) - taud_mean) / taud_lim

                if (abs(taud_acy_0(i,j)) .eq. taud_lim .and. abs(taud_eps) .gt. taud_eps_lim) then 
                    ! Set driving stress equal to weighted average of neighbors
                    taud_acy(i,j) = 0.20*taud_acy_0(i,j) + 0.40*taud_acy_0(i,jm1)+ 0.40*taud_acy_0(i,jp1)
                end if 

            end if
            
        end do
        end do 

end if 


        return 

    end subroutine calc_driving_stress

    subroutine calc_driving_stress_gl(taud_acx,taud_acy,H_ice,z_srf,z_bed,z_sl,H_grnd, &
                                      f_grnd,f_grnd_acx,f_grnd_acy,dx,rho_ice,rho_sw,g,method,beta_gl_stag)
        ! taud = rho_ice*g*H_ice
        ! Calculate driving stress on staggered grid points, with 
        ! special treatment of the grounding line 
        ! Units: taud [Pa] == [kg m-1 s-2]
        
        ! Note: interpolation to Ab nodes no longer used here.

        implicit none 

        real(wp), intent(OUT) :: taud_acx(:,:)
        real(wp), intent(OUT) :: taud_acy(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: z_srf(:,:)
        real(wp), intent(IN)  :: z_bed(:,:)
        real(wp), intent(IN)  :: z_sl(:,:)
        real(wp), intent(IN)  :: H_grnd(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: f_grnd_acx(:,:)
        real(wp), intent(IN)  :: f_grnd_acy(:,:)
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw 
        real(wp), intent(IN)  :: g 
        integer,    intent(IN)  :: method        ! Which driving stress calculation to use
        integer,    intent(IN)  :: beta_gl_stag  ! Method of grounding line staggering of beta 

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1  
        real(wp) :: dy, rhog 
        real(wp) :: taud_grnd, taud_flt, taud_now 
        real(wp) :: H_mid, H_gl, z_gl, H_grnd_mid 
        real(wp) :: dzsdx, dzsdy
        real(wp) :: dzsdx_1, dzsdx_2
        real(wp) :: H_1, H_2  
        real(wp) :: taud_old, fac_gl   

        real(wp), parameter :: slope_max = 0.05   ! Very high limit == 0.05, low limit < 0.01 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Define shortcut parameter 
        rhog = rho_ice * g 

        ! Assume grid resolution is symmetrical 
        dy = dx 

        ! === Subgrid treatment === 
        
        ! Next, refine the driving stress at the grounding line
        ! as desired 

        select case(method)

            case(-1)
                ! Unsupported option: hard failure

                write(*,*) "calc_driving_stress_gl:: Error: taud_gl_method=-1 not supported."
                error stop "taud_gl_method=-1 not supported"

            case(1)
                ! Weighted average using the grounded fraction (ac-nodes)
                ! or one-sided choice
                ! between surface slope and virtual slope of 
                ! floating ice (using ice thickness)

                ! x-direction 
                do j = 1, ny 
                do i = 1, nx-1 

                    if ( f_grnd_acx(i,j) .gt. 0.0 .and. f_grnd_acx(i,j) .lt. 1.0) then 
                        ! Grounding line point (ac-node)

                        ! Get the ice thickness at the ac-node as the average of two neighbors
                        H_gl    = 0.5_wp*(H_ice(i,j)+H_ice(i+1,j))

                        ! Get slope of grounded point and of floating point (taken as zero),
                        ! then assume slope is the weighted average of the two 
                        dzsdx_1 = (z_srf(i+1,j)-z_srf(i,j)) / dx 
                        dzsdx_2 = 0.0 
                        dzsdx   = f_grnd_acx(i,j)*dzsdx_1 + (1.0-f_grnd_acx(i,j))*dzsdx_2  
                        
                        ! Limit the slope
                        call minmax(dzsdx,slope_max)  
                                                    
                        ! Get the driving stress
                        taud_old = taud_acx(i,j) 
                        taud_acx(i,j) = rhog * H_gl * dzsdx

                    end if

                end do 
                end do 

                ! y-direction 
                do j = 1, ny-1 
                do i = 1, nx 

                    if ( f_grnd_acy(i,j) .gt. 0.0 .and. f_grnd_acy(i,j) .lt. 1.0) then 
                        ! Grounding line point (ac-node)

                        ! Get the ice thickness at the ac-node as the average of two neighbors
                        H_gl    = 0.5_wp*(H_ice(i,j)+H_ice(i,j+1))

                        ! Get slope of grounded point and of floating point (taken as zero),
                        ! then assume slope is the weighted average of the two 
                        dzsdx_1 = (z_srf(i,j+1)-z_srf(i,j)) / dx 
                        dzsdx_2 = 0.0 
                        dzsdy   = f_grnd_acy(i,j)*dzsdx_1 + (1.0-f_grnd_acy(i,j))*dzsdx_2  
                        
                        call minmax(dzsdy,slope_max)  

                        ! Get the driving stress
                        taud_acy(i,j) = rhog * H_gl * dzsdy
                        
                    end if 

                end do 
                end do 

            case(2)
                ! One-sided differences upstream and downstream of the grounding line
                ! analgous to Feldmann et al. (2014, JG)

                ! x-direction 
                do j = 1, ny 
                do i = 1, nx-1 

                    if ( f_grnd_acx(i,j) .gt. 0.0 .and. f_grnd_acx(i,j) .lt. 1.0) then 
                        ! Grounding line point (ac-node)

                        H_grnd_mid = 0.5_wp*(H_grnd(i,j) + H_grnd(i+1,j))

                        if (H_grnd_mid .gt. 0.0) then 
                            ! Consider grounded 
                            dzsdx = (z_srf(i+1,j)-z_srf(i,j)) / dx 
                        else 
                            ! Consider floating: slope of the floating surface z_sl+(1-rho_ice/rho_sw)*H_ice
                            dzsdx = ( (z_sl(i+1,j)-z_sl(i,j)) + (1.0_wp-rho_ice/rho_sw)*(H_ice(i+1,j)-H_ice(i,j)) ) / dx
                        end if 
                        call minmax(dzsdx,slope_max)  

                        ! Get the ice thickness at the ac-node
                        H_gl    = 0.5_wp*(H_ice(i,j)+H_ice(i+1,j))

                        ! Get the driving stress
                        taud_old = taud_acx(i,j) 
                        taud_acx(i,j) = rhog * H_gl * dzsdx

                    end if

                end do 
                end do 

                ! y-direction 
                do j = 1, ny-1 
                do i = 1, nx 

                    if ( f_grnd_acy(i,j) .gt. 0.0 .and. f_grnd_acy(i,j) .lt. 1.0) then 
                        ! Grounding line point (ac-node)

                        H_grnd_mid = 0.5_wp*(H_grnd(i,j) + H_grnd(i,j+1))

                        if (H_grnd_mid .gt. 0.0) then 
                            ! Consider grounded 
                            dzsdx = (z_srf(i,j+1)-z_srf(i,j)) / dx 
                        else 
                            ! Consider floating: slope of the floating surface z_sl+(1-rho_ice/rho_sw)*H_ice
                            dzsdx = ( (z_sl(i,j+1)-z_sl(i,j)) + (1.0_wp-rho_ice/rho_sw)*(H_ice(i,j+1)-H_ice(i,j)) ) / dx
                        end if 
                        call minmax(dzsdx,slope_max)  

                        ! Get the ice thickness at the ac-node
                        H_gl    = 0.5_wp*(H_ice(i,j)+H_ice(i,j+1))

                        ! Get the driving stress
                        taud_acy(i,j) = rhog * H_gl * dzsdx
                        
                    end if 

                end do 
                end do 

            case(3) 
                ! Linear interpolation following Gladstone et al. (2010, TC) Eq. 27

                 
                do j = 1, ny 
                do i = 1, nx

                    im1 = max(1, i-1)
                    ip1 = min(nx,i+1)
                    
                    jm1 = max(1, j-1)
                    jp1 = min(ny,j+1)

                    ! x-direction
                    if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(ip1,j) .le. 0.0) then 
                        ! Grounding line point 

                        ! (i,j) grounded; (ip1,j) floating
                        call integrate_gl_driving_stress_linear(taud_acx(i,j),H_ice(i,j),H_ice(ip1,j), &
                                            z_bed(i,j),z_bed(ip1,j),z_sl(i,j), z_sl(ip1,j),dx,rho_ice,rho_sw,g)
                    
                    else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(ip1,j) .gt. 0.0) then 
                        ! (i,j) floating; (ip1,j) grounded 

                        call integrate_gl_driving_stress_linear(taud_acx(i,j),H_ice(ip1,j),H_ice(i,j), &
                                            z_bed(ip1,j),z_bed(i,j),z_sl(ip1,j),z_sl(i,j),dx,rho_ice,rho_sw,g)
                        
                        ! Set negative for direction 
                        taud_acx(i,j) = -taud_acx(i,j)

                    end if 

                    ! y-direction 
                    if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i,jp1) .le. 0.0) then
                        ! Grounding line point
 
                        ! (i,j) grounded; (i,jp1) floating 
                        call integrate_gl_driving_stress_linear(taud_acy(i,j),H_ice(i,j),H_ice(i,jp1), &
                                            z_bed(i,j),z_bed(i,jp1),z_sl(i,j),z_sl(i,jp1),dx,rho_ice,rho_sw,g)

                    else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i,jp1) .gt. 0.0) then
                        ! (i,j) floating; (i,jp1) grounded

                        call integrate_gl_driving_stress_linear(taud_acy(i,j),H_ice(i,jp1),H_ice(i,j), &
                                            z_bed(i,jp1),z_bed(i,j),z_sl(i,jp1),z_sl(i,j),dx,rho_ice,rho_sw,g)
                        
                        ! Set negative for direction 
                        taud_acy(i,j) = -taud_acy(i,j) 
                    end if

                end do 
                end do 

            case DEFAULT  
                ! Do nothing, use the standard no-subgrid treatment 
        end select 

        return 

    end subroutine calc_driving_stress_gl

    subroutine calc_lateral_bc_stress_2D(tau_bc_int_acx,tau_bc_int_acy,mask_frnt,H_ice,f_ice,z_srf,z_sl, &
                                                                                rho_ice,rho_sw,g,boundaries)
            ! Calculate the vertically integrated lateral stress [Pa m] boundary condition
            ! at the ice front. 

        implicit none 

        real(wp), intent(OUT) :: tau_bc_int_acx(:,:) 
        real(wp), intent(OUT) :: tau_bc_int_acy(:,:) 
        integer,  intent(IN)  :: mask_frnt(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: z_srf(:,:) 
        real(wp), intent(IN)  :: z_sl(:,:) 
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw 
        real(wp), intent(IN)  :: g 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        integer  :: i1, j1  
        real(wp) :: H_ice_now 
        real(wp) :: z_srf_now 
        real(wp) :: z_sl_now 
        
        integer  :: BC

        nx = size(tau_bc_int_acx,1) 
        ny = size(tau_bc_int_acx,2) 

        ! Intialize boundary fields to zero 
        tau_bc_int_acx = 0.0 
        tau_bc_int_acy = 0.0 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,i1,j1,H_ice_now,z_srf_now,z_sl_now)
        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! == acx nodes == 

            ! if ( (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) .or. &
            !      (f_ice(i,j) .lt. 1.0 .and. f_ice(ip1,j) .eq. 1.0) ) then 
            !     ! Ice margin detected, proceed to calculating the boundary stress

            if ( (mask_frnt(i,j) .gt. 0 .and. mask_frnt(ip1,j) .lt. 0) .or. & 
                 (mask_frnt(i,j) .lt. 0 .and. mask_frnt(ip1,j) .gt. 0) ) then 
                ! Ice margin detected, proceed to calculating the boundary stress

                ! Determine index of the ice covered cell.
                if (mask_frnt(i,j) .lt. 0.0) then 
                    i1 = ip1 
                else
                    i1 = i 
                end if 

                ! Get current ice thickness, surface elevation 
                ! and sea level at ice front from aa-node values.
                ! (No f_ice scaling since f_ice=1)
                H_ice_now = H_ice(i1,j)     
                z_srf_now = z_srf(i1,j) 
                z_sl_now  = z_sl(i1,j) 

                ! No water back-pressure across a face to ice-free land
                if (mask_frnt(i,j) .eq. MASK_FRNT_ICE_FREE_LAND .or. &
                    mask_frnt(ip1,j) .eq. MASK_FRNT_ICE_FREE_LAND) z_sl_now = z_srf_now - H_ice_now

                ! Calculate the lateral stress bc for this point
                call calc_lateral_bc_stress(tau_bc_int_acx(i,j),H_ice_now, &
                                                z_srf_now,z_sl_now,rho_ice,rho_sw,g)

            end if 


            ! == acy nodes == 

            ! if ( (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) .or. &
            !      (f_ice(i,j) .lt. 1.0 .and. f_ice(i,jp1) .eq. 1.0) ) then 
            !     ! Ice margin detected, proceed to calculating the boundary stress

            if ( (mask_frnt(i,j) .gt. 0 .and. mask_frnt(i,jp1) .lt. 0) .or. & 
                 (mask_frnt(i,j) .lt. 0 .and. mask_frnt(i,jp1) .gt. 0) ) then 
                ! Ice margin detected, proceed to calculating the boundary stress
                
                ! Determine index of the ice covered cell.
                if (mask_frnt(i,j) .lt. 0.0) then 
                    j1 = jp1 
                else
                    j1 = j
                end if 

                ! Get current ice thickness, surface elevation 
                ! and sea level at ice front from aa-node values.
                ! (No f_ice scaling since f_ice=1)
                H_ice_now = H_ice(i,j1)     
                z_srf_now = z_srf(i,j1) 
                z_sl_now  = z_sl(i,j1) 

                ! No water back-pressure across a face to ice-free land
                if (mask_frnt(i,j) .eq. MASK_FRNT_ICE_FREE_LAND .or. &
                    mask_frnt(i,jp1) .eq. MASK_FRNT_ICE_FREE_LAND) z_sl_now = z_srf_now - H_ice_now

                ! Calculate the lateral stress bc for this point
                call calc_lateral_bc_stress(tau_bc_int_acy(i,j),H_ice_now, &
                                                z_srf_now,z_sl_now,rho_ice,rho_sw,g)

            end if 

        end do 
        end do
        !$omp end parallel do

        return

    end subroutine calc_lateral_bc_stress_2D
    
    subroutine calc_lateral_bc_stress(tau_bc_int,H_ice,z_srf,z_sl,rho_ice,rho_sw,g)
            ! Calculate the vertically integrated lateral stress [Pa m] boundary condition
            ! at the ice front for given conditions.

            ! =========================================================
            ! Generalized solution for all ice fronts (floating and grounded)
            ! See Lipscomb et al. (2019), Eqs. 11 & 12, and 
            ! Winkelmann et al. (2011), Eq. 27 

        implicit none 

        real(wp), intent(OUT) :: tau_bc_int
        real(wp), intent(IN)  :: H_ice 
        real(wp), intent(IN)  :: z_srf
        real(wp), intent(IN)  :: z_sl
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw
        real(wp), intent(IN)  :: g 

        ! Local variables 
        real(wp) :: f_submerged 
        real(wp) :: H_ocn 

        ! Get current ocean thickness bordering ice sheet
        ! (for bedrock above sea level, this will give zero)
        f_submerged = 1.d0 - min((z_srf-z_sl)/H_ice,1.d0)
        H_ocn       = H_ice*f_submerged
        
        tau_bc_int = 0.5d0*rho_ice*g*H_ice**2 &         ! tau_out_int                                                ! p_out
                   - 0.5d0*rho_sw *g*H_ocn**2           ! tau_in_int

        return

    end subroutine calc_lateral_bc_stress

    subroutine integrate_gl_driving_stress_linear(taud,H_a,H_b,zb_a,zb_b,z_sl_a,z_sl_b,dx,rho_ice,rho_sw,g)
        ! Compute the driving stress for the grounding line more precisely (subgrid)
        ! following Gladstone et al. (2010, TC), Eq. 27 
        ! Note: here cell i is grounded and cell i+1 is floating 
        ! Units: taud [Pa] 

        implicit none 

        real(wp), intent(OUT) :: taud
        real(wp), intent(IN) :: H_a, H_b          ! Ice thickness cell i and i+1, resp.
        real(wp), intent(IN) :: zb_a, zb_b        ! Bedrock elevation cell i and i+1, resp.
        real(wp), intent(IN) :: z_sl_a,  z_sl_b   ! Sea level cell i and i+1, resp.
        real(wp), intent(IN) :: dx  
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw
        real(wp), intent(IN)  :: g 

        ! Local variables 
        real(wp) :: Ha, Hb, Sa, Sb, Ba, Bb, sla, slb 
        real(wp) :: dl, dx_ab 
        real(wp) :: H_mid, dzsdx 
        real(wp) :: rho_sw_ice, rho_ice_sw 
        integer  :: n 
        integer, parameter :: ntot = 100 

        ! Parameters 
        rho_sw_ice = rho_sw / rho_ice 
        rho_ice_sw = rho_ice / rho_sw 

        ! Get step size (dimensionless) and step resolution,
        ! such that ntot segments span exactly one cell (from a to b)
        dl    = 1.0_wp / real(ntot,wp)
        dx_ab  = dx*dl 

        ! Initialize driving stress to zero 
        taud = 0.0_wp 

        ! Step through the grid cell and calculate
        ! each piecewise value of driving stress
        do n = 1, ntot 

            Ha  = H_a + (H_b-H_a)*dl*(n-1)
            Hb  = H_a + (H_b-H_a)*dl*(n)

            Ba  = zb_a + (zb_b-zb_a)*dl*(n-1)
            Bb  = zb_a + (zb_b-zb_a)*dl*(n)
            
            sla = z_sl_a + (z_sl_b-z_sl_a)*dl*(n-1)
            slb = z_sl_a + (z_sl_b-z_sl_a)*dl*(n)
            
            if (Ha < rho_sw_ice*(sla-Ba)) then 
                Sa = sla + (1.0-rho_ice_sw)*Ha
            else 
                Sa = Ba + Ha 
            end if 

            if (Hb < rho_sw_ice*(slb-Bb)) then 
                Sb = slb + (1.0-rho_ice_sw)*Hb
            else 
                Sb = Bb + Hb 
            end if 
            
            H_mid = 0.5_wp * (Ha+Hb)
            dzsdx = (Sb-Sa) / dx_ab 

            taud  = taud + (H_mid * dzsdx)*dl 

        end do 

        ! Finally multiply with rho_ice*g 
        taud = rho_ice*g *taud

        return 

    end subroutine integrate_gl_driving_stress_linear
    
    subroutine set_inactive_margins(ux,uy,f_ice,boundaries,a_front)
        ! Zero the velocity on faces between an ice-free cell and a cell that
        ! is not fully ice covered: a partial cell fills before ice flows
        ! beyond it. With a_front, the area fraction of each cell behind a
        ! prescribed front (the level set), a face into an ice-free cell stays
        ! open where the front covers at least A_FRONT_MIN of that cell, so
        ! cells fill as the front advances over them. Without a prescribed
        ! front the only front is the ice itself (a_front = f_ice).

        implicit none

        real(wp), intent(INOUT) :: ux(:,:) 
        real(wp), intent(INOUT) :: uy(:,:) 
        real(wp), intent(IN)    :: f_ice(:,:) 
        character(len=*), intent(IN) :: boundaries 
        real(wp), intent(IN), optional :: a_front(:,:)  ! [--] Area fraction behind a prescribed front

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1
        integer :: BC
        logical, allocatable :: fill(:,:)               ! Ice-free cells that may be filled

        nx = size(f_ice,1) 
        ny = size(f_ice,2) 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        allocate(fill(nx,ny))
        if (present(a_front)) then
            fill = a_front .ge. A_FRONT_MIN
        else
            fill = .FALSE.
        end if

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            if (face_closed(f_ice(i,j),f_ice(ip1,j),fill(i,j),fill(ip1,j))) ux(i,j) = 0.0_wp
            if (face_closed(f_ice(i,j),f_ice(i,jp1),fill(i,j),fill(i,jp1))) uy(i,j) = 0.0_wp

        end do
        end do
        !$omp end parallel do

        return

    contains

        pure logical function face_closed(f_a,f_b,fill_a,fill_b)
            ! Face between an ice-free cell that may not be filled and a
            ! cell that is not fully ice covered
            real(wp), intent(IN) :: f_a, f_b
            logical,  intent(IN) :: fill_a, fill_b
            face_closed = (f_b .eq. 0.0_wp .and. .not. fill_b .and. f_a .lt. 1.0_wp) .or. &
                          (f_a .eq. 0.0_wp .and. .not. fill_a .and. f_b .lt. 1.0_wp)
        end function face_closed

    end subroutine set_inactive_margins

    subroutine calc_ice_flux(qq_acx,qq_acy,ux_bar,uy_bar,H_ice,dx,dy,boundaries)
        ! Calculate the ice flux through each cell face (ac-nodes), using the
        ! upwind ice thickness, as in the advection solvers (impl-lis, expl-upwind).
        ! qq      [m3 a-1] 
        ! ux,uy   [m a-1]
        ! H_ice   [m] 

        implicit none 

        real(wp), intent(OUT) :: qq_acx(:,:)     ! [m3 a-1] Ice flux (acx nodes)
        real(wp), intent(OUT) :: qq_acy(:,:)     ! [m3 a-1] Ice flux (acy nodes)
        real(wp), intent(IN)  :: ux_bar(:,:)     ! [m a-1]  Vertically averaged velocity (acx nodes)
        real(wp), intent(IN)  :: uy_bar(:,:)     ! [m a-1]  Vertically averaged velocity (acy nodes)
        real(wp), intent(IN)  :: H_ice(:,:)      ! [m]      Ice thickness, aa-nodes
        real(wp), intent(IN)  :: dx              ! [m]      Horizontal resolution, x-dir
        real(wp), intent(IN)  :: dy              ! [m]      Horizontal resolution, y-dir 
        character(len=*), intent(IN) :: boundaries

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1
        integer :: i2, j2
        integer :: BC
        logical :: per_x, per_y
        real(wp) :: area_ac 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Set boundary condition code
        BC = boundary_code(boundaries)
        call get_periodic_directions(per_x,per_y,BC)

        ! The last ac-node is an interior node only in a periodic direction
        i2 = nx-1
        if (per_x) i2 = nx
        j2 = ny-1
        if (per_y) j2 = ny

        ! Reset fluxes to zero 
        qq_acx = 0.0 
        qq_acy = 0.0 

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, i2 
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            if (ux_bar(i,j) .ge. 0.0_wp) then
                area_ac = H_ice(i,j)   * dy
            else
                area_ac = H_ice(ip1,j) * dy
            end if
            qq_acx(i,j) = area_ac*ux_bar(i,j)
        end do 
        end do 

        ! acy-nodes 
        do j = 1, j2 
        do i = 1, nx 
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            if (uy_bar(i,j) .ge. 0.0_wp) then
                area_ac = H_ice(i,j)   * dx
            else
                area_ac = H_ice(i,jp1) * dx
            end if
            qq_acy(i,j) = area_ac*uy_bar(i,j)
        end do 
        end do 
        
        return 

    end subroutine calc_ice_flux

    subroutine calc_grounding_line_flux(qq_gl_acx,qq_gl_acy,qq_acx,qq_acy,f_grnd,f_ice,boundaries)
        ! Ice flux across the grounding line [m3 a-1], on ac-nodes: the ice flux
        ! (qq_acx/qq_acy, see calc_ice_flux) through faces between a (partially)
        ! grounded cell (f_grnd > 0) and a floating ice cell (f_grnd = 0, f_ice > 0),
        ! and zero elsewhere. As for qq, the sign gives the direction (+x/+y).
        ! Grounding-line cells are defined as in calc_distance_to_grounding_line.

        implicit none

        real(wp), intent(OUT) :: qq_gl_acx(:,:)  ! [m3 a-1] Grounding-line flux (acx nodes)
        real(wp), intent(OUT) :: qq_gl_acy(:,:)  ! [m3 a-1] Grounding-line flux (acy nodes)
        real(wp), intent(IN)  :: qq_acx(:,:)     ! [m3 a-1] Ice flux (acx nodes)
        real(wp), intent(IN)  :: qq_acy(:,:)     ! [m3 a-1] Ice flux (acy nodes)
        real(wp), intent(IN)  :: f_grnd(:,:)     ! [--]     Grounded fraction, aa-nodes
        real(wp), intent(IN)  :: f_ice(:,:)      ! [--]     Ice-covered fraction, aa-nodes
        character(len=*), intent(IN) :: boundaries

        ! Local variables
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1
        integer :: BC

        nx = size(f_grnd,1)
        ny = size(f_grnd,2)

        BC = boundary_code(boundaries)

        qq_gl_acx = 0.0_wp
        qq_gl_acy = 0.0_wp

        do j = 1, ny
        do i = 1, nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            if (is_gl_face(f_grnd(i,j),f_grnd(ip1,j),f_ice(i,j),f_ice(ip1,j))) then
                qq_gl_acx(i,j) = qq_acx(i,j)
            end if

            if (is_gl_face(f_grnd(i,j),f_grnd(i,jp1),f_ice(i,j),f_ice(i,jp1))) then
                qq_gl_acy(i,j) = qq_acy(i,j)
            end if

        end do
        end do

        return

    contains

        pure logical function is_gl_face(fg0,fg1,fi0,fi1)
            ! Face between a (partially) grounded cell and a floating ice cell

            implicit none

            real(wp), intent(IN) :: fg0, fg1, fi0, fi1

            is_gl_face = (fg0 .gt. 0.0_wp .and. fg1 .eq. 0.0_wp .and. fi1 .gt. 0.0_wp) .or. &
                         (fg1 .gt. 0.0_wp .and. fg0 .eq. 0.0_wp .and. fi0 .gt. 0.0_wp)

        end function is_gl_face

    end subroutine calc_grounding_line_flux

    elemental function calc_vel_ratio(uxy_base,uxy_srf) result(f_vbvs)

        implicit none 

        real(wp), intent(IN) :: uxy_base, uxy_srf  
        real(wp) :: f_vbvs 

        ! Calculate the basal to surface velocity ratio, f_vbvs
        if ( uxy_srf .gt. 0.0) then
            f_vbvs = min(1.0, uxy_base / uxy_srf) 
        else 
            ! No ice (or no velocity)
            f_vbvs = 1.0 
        end if 

        return 

    end function calc_vel_ratio

    subroutine picard_calc_error(corr,ux,uy,ux_prev,uy_prev)
        ! Calculate the current error, ie, the 'correction vector'
        ! as defined by De Smedt et al. (2010):
        ! corr = U_now - U_prev.

        implicit none

        real(wp), intent(OUT) :: corr(:) 
        real(wp), intent(IN)  :: ux(:,:) 
        real(wp), intent(IN)  :: uy(:,:) 
        real(wp), intent(IN)  :: ux_prev(:,:)
        real(wp), intent(IN)  :: uy_prev(:,:) 
        
        ! Local variables 
        integer :: i, j, nx, ny, k  

        nx = size(ux,1)
        ny = size(ux,2) 

        ! Consistency check 
        if (size(corr,1) .ne. 2*nx*ny) then 
            write(*,*) "calc_convergence_angle:: Error: corr(N) must have N=2*nx*ny."
            stop 
        end if 

        k = 0

        do j = 1, ny 
        do i = 1, nx  

            k = k+1 
            corr(k) = ux(i,j) - ux_prev(i,j) 

            k = k+1 
            corr(k) = uy(i,j) - uy_prev(i,j) 

        end do 
        end do 

        return

    end subroutine picard_calc_error

    subroutine picard_calc_error_angle(theta,corr_nm1,corr_nm2)
        ! Calculate the current error, ie, the 'correction vector'
        ! as defined by De Smedt et al. (2010):
        ! corr = U_now - U_prev.

        implicit none

        real(wp), intent(OUT) :: theta
        real(wp), intent(IN)  :: corr_nm1(:) 
        real(wp), intent(IN)  :: corr_nm2(:) 
        
        
        ! Local variables   
        real(dp) :: val 

        real(dp), parameter :: tol = 1e-5 

        val = sum(corr_nm1*corr_nm2) / &
                (sqrt(sum(corr_nm1**2))*sqrt(sum(corr_nm2**2))+tol)

        theta = acos(val) 

        return

    end subroutine picard_calc_error_angle

    subroutine picard_calc_convergence_l1rel_matrix(err_x,err_y,ux,uy,ux_prev,uy_prev)

        implicit none 

        real(wp), intent(OUT) :: err_x(:,:)
        real(wp), intent(OUT) :: err_y(:,:)
        real(wp), intent(IN)  :: ux(:,:) 
        real(wp), intent(IN)  :: uy(:,:) 
        real(wp), intent(IN)  :: ux_prev(:,:) 
        real(wp), intent(IN)  :: uy_prev(:,:)  

        ! Local variables
        integer :: i, j

        real(wp), parameter :: ssa_vel_tolerance = 1e-2   ! [m/a] only consider points with velocity above this tolerance limit
        real(wp), parameter :: tol = 1e-5 

        !$omp parallel do collapse(2) private(i,j)
        do j = 1, size(ux,2)
        do i = 1, size(ux,1)

            ! Error in x-direction
            if (abs(ux(i,j)) .gt. ssa_vel_tolerance) then
                err_x(i,j) = 2.0_wp * abs(ux(i,j) - ux_prev(i,j)) / abs(ux(i,j) + ux_prev(i,j) + tol)
            else
                err_x(i,j) = 0.0_wp
            end if

            ! Error in y-direction 
            if (abs(uy(i,j)) .gt. ssa_vel_tolerance) then
                err_y(i,j) = 2.0_wp * abs(uy(i,j) - uy_prev(i,j)) / abs(uy(i,j) + uy_prev(i,j) + tol)
            else
                err_y(i,j) = 0.0_wp
            end if

        end do
        end do
        !$omp end parallel do

        return 

    end subroutine picard_calc_convergence_l1rel_matrix

    subroutine picard_calc_convergence_l2(is_converged,resid,ux,uy,ux_prev,uy_prev, &
                                                mask_acx,mask_acy,resid_tol,iter,iter_max,log,lin_iter,lin_status)

        ! Calculate the current level of convergence using the 
        ! L2 relative error norm between the current and previous
        ! velocity solutions.

        ! Note for parameter norm_method defined below: 

        ! norm_method=1: as defined by De Smedt et al. (2010):
        ! conv = ||U_now - U_prev||/||U_prev||.

        ! norm_method=2: as defined by Gagliardini et al., GMD, 2013, Eq. 65:
        ! conv = 2*||U_now - U_prev||/||U_now+U_prev||.
        
        implicit none 

        logical,  intent(OUT) :: is_converged
        real(wp), intent(OUT) :: resid 
        real(wp), intent(IN) :: ux(:,:) 
        real(wp), intent(IN) :: uy(:,:) 
        real(wp), intent(IN) :: ux_prev(:,:) 
        real(wp), intent(IN) :: uy_prev(:,:)  
        logical,  intent(IN) :: mask_acx(:,:) 
        logical,  intent(IN) :: mask_acy(:,:) 
        real(wp), intent(IN) :: resid_tol 
        integer,  intent(IN) :: iter 
        integer,  intent(IN) :: iter_max 
        logical,  intent(IN) :: log 
        integer,  intent(IN) :: lin_iter                ! Linear solver iterations of this Picard iteration
        integer,  intent(IN) :: lin_status              ! Linear solver return status

        ! Local variables
        integer :: i, j, nx, ny, k
        real(dp) :: tmpx, tmpy
        real(dp) :: res1, res2, res3
        real(dp) :: res1_j(size(ux,2))              ! Row sums of res1, res2, res3
        real(dp) :: res2_j(size(ux,2))
        real(dp) :: res3_j(size(ux,2))
        
        real(wp) :: ux_resid_max 
        real(wp) :: uy_resid_max 
        integer  :: nx_check, ny_check  
        character(len=1) :: converged_txt 

        real(dp), parameter :: du_reg  = 1e-10              ! [m/yr] Small regularization factor to avoid divide-by-zero
        real(wp), parameter :: vel_tol = 1e-5               ! [m/yr] only consider points with velocity above this tolerance limit
        integer,  parameter :: norm_method = 1              ! See note above.
        
        nx = size(ux,1)
        ny = size(ux,2)

        ! One pass over the domain: count the points to check, sum the
        ! squared differences and get the maximum difference per direction.
        ! The squared differences are summed per row here and the rows are
        ! added serially below, so the residual (and the convergence decision)
        ! does not depend on the number of threads.
        nx_check = 0
        ny_check = 0
        ux_resid_max = 0.0
        uy_resid_max = 0.0

        !$omp parallel do private(i,j,tmpx,tmpy,res1,res2,res3) &
        !$omp& reduction(+:nx_check,ny_check) reduction(max:ux_resid_max,uy_resid_max)
        do j = 1, ny

        res1 = 0.0
        res2 = 0.0
        res3 = 0.0

        do i = 1, nx

            ! x-direction contribution (acx-node, own mask)
            if (abs(ux(i,j)) .gt. vel_tol .and. mask_acx(i,j)) then
                nx_check = nx_check + 1
                ux_resid_max = max(ux_resid_max,abs(ux(i,j)-ux_prev(i,j)))

                tmpx = ux(i,j)-ux_prev(i,j)
                if (dabs(tmpx) .lt. TOL_UNDERFLOW) tmpx = 0.0
                res1 = res1 + tmpx*tmpx

                tmpx = ux_prev(i,j)
                if (dabs(tmpx) .lt. TOL_UNDERFLOW) tmpx = 0.0
                res2 = res2 + tmpx*tmpx

                tmpx = ux(i,j)+ux_prev(i,j)
                res3 = res3 + tmpx*tmpx
            end if

            ! y-direction contribution (acy-node, own mask)
            if (abs(uy(i,j)) .gt. vel_tol .and. mask_acy(i,j)) then
                ny_check = ny_check + 1
                uy_resid_max = max(uy_resid_max,abs(uy(i,j)-uy_prev(i,j)))

                tmpy = uy(i,j)-uy_prev(i,j)
                if (dabs(tmpy) .lt. TOL_UNDERFLOW) tmpy = 0.0
                res1 = res1 + tmpy*tmpy

                tmpy = uy_prev(i,j)
                if (dabs(tmpy) .lt. TOL_UNDERFLOW) tmpy = 0.0
                res2 = res2 + tmpy*tmpy

                tmpy = uy(i,j)+uy_prev(i,j)
                res3 = res3 + tmpy*tmpy
            end if

        end do

        res1_j(j) = res1
        res2_j(j) = res2
        res3_j(j) = res3

        end do
        !$omp end parallel do

        res1 = sum(res1_j)
        res2 = sum(res2_j)
        res3 = sum(res3_j)

        if ( (nx_check+ny_check) .gt. 0 ) then

            select case(norm_method)

                case(1)

                    resid = sqrt(res1)/(sqrt(res2)+du_reg)

                case(2)

                    resid = 2.0_wp*sqrt(res1)/(sqrt(res3)+du_reg)

            end select 

        else 
            ! No points available for comparison, set residual equal to zero 

            resid = 0.0_wp 

        end if

        ! Check for convergence
        if (resid .lt. resid_tol) then 
            is_converged = .TRUE. 
            converged_txt = "C"
        else if (iter .eq. iter_max) then 
            is_converged = .TRUE. 
            converged_txt = "X" 
        else 
            is_converged = .FALSE. 
            converged_txt = ""
        end if 

        !if (log .and. is_converged) then
        if (log) then
            ! Write summary to log if desired and iterations have completed

            ! Write summary to log
            write(*,"(a,a2,i4,g12.4,a3,2i8,2g12.4,a,i6,1x,a)") &
                "ssa: ", trim(converged_txt), iter, resid, " | ", nx_check, ny_check, ux_resid_max, uy_resid_max, &
                " | lin", lin_iter, trim(linear_solver_status_str(lin_status))

        end if 


if (.TRUE.) then
        if (ux_resid_max .ge. 9999.0_wp .or. uy_resid_max .ge. 9999.0_wp) then 
            ! Strange case is occurring. Poor convergence with high error, investigate

            write(io_unit_err,*) "ssa: Error: strange case occurring."
            write(io_unit_err,"(a,a2,i4,g12.4,a3,2i8,2g12.4)") &
            "ssa: ", trim(converged_txt), iter, resid, " | ", nx_check, ny_check, ux_resid_max, uy_resid_max 

            write(io_unit_err,*) "Writing diagnostic file: ssa_check.nc."

            ! Initialize output file 
            call ssa_diagnostics_write_init("ssa_check.nc",nx=size(ux,1),ny=size(ux,2),time_init=real(iter,wp))

            ! Write file with dummy zero values for variables we don't have
            call ssa_diagnostics_write_step("ssa_check.nc",ux,uy,resid,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp, &
                                    int(ux*0.0_wp),int(ux*0.0_wp),ux-ux_prev,uy-uy_prev,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp, &
                                    ux*0.0_wp,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp,ux_prev,uy_prev,time=real(iter,wp))

            stop 

        end if 
end if 
        
        return 

    end subroutine picard_calc_convergence_l2

    subroutine picard_relax_vel(ux,uy,ux_prev,uy_prev,rel)
        ! Relax velocity solution with previous iteration 

        implicit none 

        real(wp), intent(INOUT) :: ux(:,:)
        real(wp), intent(INOUT) :: uy(:,:)
        real(wp), intent(IN)    :: ux_prev(:,:)
        real(wp), intent(IN)    :: uy_prev(:,:)
        real(wp), intent(IN)    :: rel

        ! Local variables
        integer :: i, j

        ! Apply relaxation 
        !$omp parallel do collapse(2) private(i,j)
        do j = 1, size(ux,2)
        do i = 1, size(ux,1)
            ux(i,j) = ux_prev(i,j) + rel*(ux(i,j)-ux_prev(i,j))
            uy(i,j) = uy_prev(i,j) + rel*(uy(i,j)-uy_prev(i,j))
        end do
        end do
        !$omp end parallel do
        
        return 

    end subroutine picard_relax_vel

    subroutine picard_relax_visc(visc,visc_prev,rel)
        ! Relax viscosity solution with previous iteration, in log-space

        implicit none 

        real(wp), intent(INOUT) :: visc(:,:,:)
        real(wp), intent(IN)    :: visc_prev(:,:,:)
        real(wp), intent(IN)    :: rel

        ! Local variables
        integer :: i, j, k

        ! Apply relaxation 
        !$omp parallel do collapse(2) private(i,j,k)
        do k = 1, size(visc,3)
        do j = 1, size(visc,2)
        do i = 1, size(visc,1)
            visc(i,j,k) = exp( (1.0-rel)*log(visc_prev(i,j,k)) + rel*log(visc(i,j,k)) )
        end do
        end do
        end do
        !$omp end parallel do
        
        return 

    end subroutine picard_relax_visc



    

end module velocity_general
