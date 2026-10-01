module calving_ac
    ! Definitions for various calving laws on ac nodes. 
    ! This will be used for flux calving (lsf) 

    use yelmo_defs, only : sp, dp, wp, prec, TOL_UNDERFLOW
    use yelmo_tools, only : boundary_code, get_neighbor_indices_bc_codes
    use thermodynamics, only : calc_T_freeze_sw

    implicit none 
    private 
    
    ! === Floating calving routines === 
    public :: calc_calving_threshold_lsf
    public :: calc_calving_rate_vonmises_m16
    
    ! === Grounded calving routines === 
    public :: calc_fmb_ismip7

    ! === CalvMIP calving rates ===
    public :: calvmip_exp1
    public :: calvmip_exp2
    public :: calvmip_exp5_aa

contains 

    ! ===================================================================
    !
    ! Calving - floating ice 
    !
    ! ===================================================================
    
    subroutine calc_calving_threshold_lsf(cr_acx,cr_acy,u_acx,v_acy,H_ice,H_ice_c,f_ice,boundaries)
        ! Threshold calving rate flux based on CalvingMIP experiment 5.
        ! Valid for floating and grounded ice.
            
        implicit none
            
        real(wp), intent(OUT) :: cr_acx(:,:), cr_acy(:,:)
        real(wp), intent(IN)  :: u_acx(:,:),  v_acy(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: H_ice_c
        real(wp), intent(IN)  :: f_ice(:,:)                ! Ocean mask. Extrapolate values into that mask.
        character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose
                
        ! local variables
        integer  :: i, j, ip1, im1, jp1, jm1, nx, ny
        real(wp) :: wv_acx, wv_acy, H_acx, H_acy
        integer  :: BC

        nx = size(u_acx,1)
        ny = size(u_acx,2)

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,H_acx,H_acy,wv_acx,wv_acy)
        do j = 1, ny
        do i = 1, nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! Stagger ice thickness into ac-nodes                        
            H_acx = 0.5*(H_ice(i,j)+H_ice(ip1,j))
            H_acy = 0.5*(H_ice(i,j)+H_ice(i,jp1))
                
            ! Special case for border
            ! x-axis  
            if ((f_ice(i,j) .gt. 0.0_wp) .and. (f_ice(ip1,j) .eq. 0.0_wp)) then
                H_acx = H_ice(i,j)
            else if ((f_ice(i,j) .eq. 0.0_wp) .and. (f_ice(ip1,j) .gt. 0.0_wp)) then
                H_acx = H_ice(ip1,j)
            end if

            ! y-axis  
            if ((f_ice(i,j) .gt. 0.0_wp) .and. (f_ice(i,jp1) .eq. 0.0_wp)) then
                H_acy = H_ice(i,j)
            else if ((f_ice(i,j) .eq. 0.0_wp) .and. (f_ice(i,jp1) .gt. 0.0_wp)) then
                H_acy = H_ice(i,jp1)
            end if

            ! Compute calving-rates on ac-nodes
            wv_acx      = MAX(0.0_wp,1.0_wp+(H_ice_c-H_acx)/H_ice_c)
            cr_acx(i,j) = -u_acx(i,j)*wv_acx
            wv_acy      = MAX(0.0_wp,1.0_wp+(H_ice_c-H_acy)/H_ice_c)
            cr_acy(i,j) = -v_acy(i,j)*wv_acy

        end do
        end do
        !$omp end parallel do
    
        return
        
    end subroutine calc_calving_threshold_lsf

    subroutine calc_calving_rate_vonmises_m16(cr_acx,cr_acy,u_acx,v_acy,tau_1,tau_ice_c,f_ice,boundaries)
        ! Calculate the calving rate [m/yr] based on the 
        ! von Mises stress approach, as outlined by Morlighem et al. (2016)
        ! DOI: 10.1002/2016gl067695
        ! Eq. 4: c = v*tau_1/tau_ice (tau_ice_flt or tau_ice_grnd)

        implicit none 

        real(wp), intent(INOUT) :: cr_acx(:,:), cr_acy(:,:) ! Simulated calving rate. ac-nodes.
        real(wp), intent(IN)    :: u_acx(:,:),  v_acy(:,:)  ! Velocity fields. ac-nodes.
        real(wp), intent(IN)    :: tau_1(:,:)               ! 1st principal stress [Pa]. aa-nodes.
        real(wp), intent(IN)    :: tau_ice_c                ! Ice fracture strength [Pa].
        real(wp), intent(IN)    :: f_ice(:,:)               ! Ocean mask. Extrapolate values into that mask.
        character(len=*), intent(IN) :: boundaries 

        ! local variables
        integer  :: i, j, ip1, im1, jp1, jm1, nx, ny
        real(wp) :: tau1_acx, tau1_acy, wv_acx, wv_acy
        integer  :: BC

        nx = size(u_acx,1)
        ny = size(u_acx,2)

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,tau1_acx,tau1_acy,wv_acx,wv_acy)
        do j = 1, ny
            do i = 1, nx
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
                    
                ! Stagger 1st ppal stress into ac-nodes                        
                tau1_acx = 0.5*(tau_1(i,j)+tau_1(ip1,j))
                tau1_acy = 0.5*(tau_1(i,j)+tau_1(i,jp1))
                        
                ! Special case for border
                ! x-axis  
                if ((f_ice(i,j) .gt. 0.0_wp) .and. (f_ice(ip1,j) .eq. 0.0_wp)) then
                    tau1_acx = tau_1(i,j)
                else if ((f_ice(i,j) .eq. 0.0_wp) .and. (f_ice(ip1,j) .gt. 0.0_wp)) then
                    tau1_acx = tau_1(ip1,j)
                end if

                ! y-axis  
                if ((f_ice(i,j) .gt. 0.0_wp) .and. (f_ice(i,jp1) .eq. 0.0_wp)) then
                    tau1_acy = tau_1(i,j)
                else if ((f_ice(i,j) .eq. 0.0_wp) .and. (f_ice(i,jp1) .gt. 0.0_wp)) then
                    tau1_acy = tau_1(i,jp1)
                end if

                ! Compute calving-rates on ac-nodes
                wv_acx      = MAX(0.0_wp,tau1_acx/tau_ice_c)
                cr_acx(i,j) = -u_acx(i,j)*wv_acx
                wv_acy      = MAX(0.0_wp,tau1_acy/tau_ice_c)
                cr_acy(i,j) = -v_acy(i,j)*wv_acy
        end do
        end do
        !$omp end parallel do

        return 

    end subroutine calc_calving_rate_vonmises_m16
       
    subroutine calc_calving_rate_eigen(mb_calv,H_ice,f_ice,f_grnd,eps_eff,dx,k2,boundaries)
        ! Calculate the 'horizontal' calving rate [m/yr] based on the 
        ! von Mises stress approach, as outlined by Lipscomb et al. (2019)
        ! Eqs. 73-75.
        ! L19: kt = 0.0025 m yr-1 Pa-1, w2=25

        implicit none 

        real(wp), intent(OUT) :: mb_calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)  
        real(wp), intent(IN)  :: eps_eff(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: k2
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j
        integer  :: im1, jm1, ip1, jp1
        integer  :: nx, ny
        integer  :: n_ocean 
        logical  :: is_margin 
        real(wp) :: dy   
        real(wp) :: calv_ref
        real(wp) :: calv_now
        real(wp) :: H_eff 

        real(wp) :: dxx, dyy, dxy 
        real(wp) :: eps_eig_1_now, eps_eig_2_now
        real(wp) :: eps_eff_neighb(4)
        real(wp) :: wt
        
        real(wp), parameter :: calv_lim = 2000.0_wp     ! To avoid really high calving values

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Assume square grid cells 
        dy = dx 

        mb_calv = 0.0_wp

        ! NOTE: Eigen calving (Levermann et al., 2012) is not fully implemented here.
        ! The lateral calving rate (calv_ref) is computed but never used, and calv_now
        ! was previously applied uninitialized. Fail loudly rather than return garbage.
        error stop "calc_calving_rate_eigen: Eigen calving is not implemented"

        return

    end subroutine calc_calving_rate_eigen
     
    ! ===================================================================
    !
    ! Calving - Grounded ice (Marine terminating glaciers)
    !
    ! ===================================================================

    subroutine calc_fmb_ismip7(cr_acx,cr_acy,lsf,z_bed,z_sl,Qd,T_ocn,T0,dx,f_ice,boundaries)
        ! Calculate the retreat rate of marine terminating glaciers based on ISMIP7 protocol
        ! 
        ! m = (a h_w q^alpha + b) TF^beta [m/d]
        ! q = 86400*Q/A [m/d]
        !
        ! a, alpha, b, beta: constants
        ! h_w: water depth at the terminus, z_sl - z_bed [m]
        ! Q: subglacial discharge [m3/s]
        ! A: submerged area of the terminus face, h_w*dx [m2]
        ! TF: thermal forcing, T_ocn - T_f(h_w) [K]
        !
        ! T_f is the seawater freezing point following Jenkins (1991),
        ! evaluated at the water depth h_w assuming a constant salinity 
        ! (see calc_T_freeze_sw). The retreat rate m is converted to [m/yr]
        ! and applied on ac-nodes along the inward front normal -n, with
        ! n = grad(lsf)/|grad(lsf)| the outward normal (lsf > 0 is ocean).
        ! The front then retreats at rate m independently of the ice flow
        ! (in the LSF velocity w = u + cr), so a stagnant front also retreats.
        ! Where |grad(lsf)| = 0 (away from the front) there is no normal and
        ! the retreat rate is zero; the level set is not moved there anyway.
        ! Note: with lsf_method="snap" the lsf is saturated to +-1 next to
        ! the front, so the normal is only resolved as axis-aligned or
        ! diagonal directions; "redist" gives a smoother normal.

        implicit none 

        real(wp), intent(INOUT) :: cr_acx(:,:), cr_acy(:,:) ! Simulated calving rate. ac-nodes.
        real(wp), intent(IN)    :: lsf(:,:)                 ! Level-set function (aa-nodes, lsf > 0 is ocean)
        real(wp), intent(IN)    :: z_bed(:,:)               ! Bedrock elevation [m]
        real(wp), intent(IN)    :: z_sl(:,:)                ! Sea level [m]
        real(wp), intent(IN)    :: Qd(:,:)                  ! subglacial discharge [m3/s]
        real(wp), intent(IN)    :: T_ocn(:,:)               ! Ocean temperature [K]
        real(wp), intent(IN)    :: T0                       ! Reference freezing temperature [K]
        real(wp), intent(IN)    :: dx                       ! Resolution [m]
        real(wp), intent(IN)    :: f_ice(:,:)               ! Ocean mask. Extrapolate values into that mask.
        character(len=*), intent(IN) :: boundaries 

        ! local variables
        integer  :: i, j, ip1, im1, jp1, jm1, nx, ny
        real(wp) :: a, b, alpha, beta, m_acx, m_acy
        real(wp) :: gx, gy, gxy
        real(wp), allocatable :: m_aa(:,:), h_w(:,:), TF(:,:)
        integer  :: BC

        nx = size(z_bed,1)
        ny = size(z_bed,2) 

        allocate(m_aa(nx,ny))
        allocate(h_w(nx,ny))
        allocate(TF(nx,ny))

        a     = 3.0e-4
        b     = 0.15
        alpha = 0.39
        beta  = 1.18
        m_aa  = 0.0_wp

        ! Water depth and thermal forcing relative to the local freezing point
        h_w   = MAX(0.0_wp, z_sl - z_bed)
        TF    = MAX(0.0_wp, T_ocn - calc_T_freeze_sw(h_w,T0))

        ! Discharge is a non-negative volume flux (q**alpha is NaN for q < 0)
        m_aa  = 365.25*(a*h_w*((86400.0*MAX(0.0_wp,Qd)/(h_w*dx+1e-8))**alpha)+b)*(TF**beta) ! is in m/yr
        where(f_ice .eq. 0.0) m_aa = 0.0_wp

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,m_acx,m_acy,gx,gy,gxy)
        do j = 1, ny
            do i = 1, nx
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
                    
                ! Stagger retreat rate into ac-nodes                        
                m_acx = 0.5*(m_aa(i,j)+m_aa(ip1,j))
                m_acy = 0.5*(m_aa(i,j)+m_aa(i,jp1))
                        
                ! Special case for border
                ! x-axis  
                if ((f_ice(i,j) .gt. 0.0_wp) .and. (f_ice(ip1,j) .eq. 0.0_wp)) then
                    m_acx = m_aa(i,j)
                else if ((f_ice(i,j) .eq. 0.0_wp) .and. (f_ice(ip1,j) .gt. 0.0_wp)) then
                    m_acx = m_aa(ip1,j)
                end if

                ! y-axis  
                if ((f_ice(i,j) .gt. 0.0_wp) .and. (f_ice(i,jp1) .eq. 0.0_wp)) then
                    m_acy = m_aa(i,j)
                else if ((f_ice(i,j) .eq. 0.0_wp) .and. (f_ice(i,jp1) .gt. 0.0_wp)) then
                    m_acy = m_aa(i,jp1)
                end if

                ! Compute calving-rates on ac-nodes (retreat along the inward
                ! front normal). grad(lsf) on each ac-node: normal component
                ! from the two adjacent aa-nodes, tangential component from
                ! the centred difference averaged over them. Only the
                ! direction is needed, so the common 1/dx factor is dropped
                ! (dx = dy).

                ! x-direction
                gx  = lsf(ip1,j) - lsf(i,j)
                gy  = 0.25_wp*(lsf(i,jp1)+lsf(ip1,jp1)-lsf(i,jm1)-lsf(ip1,jm1))
                gxy = sqrt(gx**2 + gy**2)
                if (gxy .gt. 0.0_wp) then
                    cr_acx(i,j) = -(gx/gxy)*MAX(0.0_wp,m_acx)
                else
                    cr_acx(i,j) = 0.0_wp
                end if

                ! y-direction
                gy  = lsf(i,jp1) - lsf(i,j)
                gx  = 0.25_wp*(lsf(ip1,j)+lsf(ip1,jp1)-lsf(im1,j)-lsf(im1,jp1))
                gxy = sqrt(gx**2 + gy**2)
                if (gxy .gt. 0.0_wp) then
                    cr_acy(i,j) = -(gy/gxy)*MAX(0.0_wp,m_acy)
                else
                    cr_acy(i,j) = 0.0_wp
                end if
            end do
        end do
        !$omp end parallel do

        deallocate(m_aa)
        deallocate(h_w)
        deallocate(TF)

        return 

    end subroutine calc_fmb_ismip7
    
    ! ===================================================================
    !
    !                      CalvMIP experiments
    !
    ! ===================================================================

    subroutine calvmip_exp1(cr_acx,cr_acy,u_acx,v_acy,lsf_aa,dx,boundaries)
        ! Experiment 1 & 3 of CalvMIP
        implicit none
    
        real(wp), intent(OUT) :: cr_acx(:,:),cr_acy(:,:)   ! Calving rates on ac-nodes
        real(wp), intent(IN)  :: u_acx(:,:),v_acy(:,:)     ! Velocities on ac-nodes
        real(wp), intent(IN)  :: lsf_aa(:,:)               ! LSF mask on aa-nodes
        real(wp), intent(IN)  :: dx                        ! Ice resolution
        character(len=*), intent(IN)  :: boundaries        ! Boundary conditions to impose
    
        ! Local variables
        integer  :: i, j, im1, ip1, jm1, jp1, nx, ny
        real(wp) :: r, rip1, rim1, rjp1, rjm1
        integer  :: BC

        nx = size(u_acx,1)
        ny = size(u_acx,2)

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Initialize calving rates to opposite as velocity
        cr_acx = -u_acx 
        cr_acy = -v_acy
    
        r    = 0.0_wp
        rip1 = 0.0_wp
        rim1 = 0.0_wp
        rjp1 = 0.0_wp
        rjm1 = 0.0_wp
        
        do j = 1, ny
        do i = 1, nx
    
            ! Below radius, no calving.
            ! aa-nodes indices
            r = sqrt((0.5*(nx+1)-i)*(0.5*(nx+1)-i) + (0.5*(ny+1)-j)*(0.5*(ny+1)-j))*dx

            ! Below radius
            if (r .lt. 750e3) then
                ! Now treat border points
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
                rip1 = sqrt((0.5*(nx+1)-ip1)*(0.5*(nx+1)-ip1) + (0.5*(ny+1)-j)*(0.5*(ny+1)-j))*dx
                rim1 = sqrt((0.5*(nx+1)-im1)*(0.5*(nx+1)-im1) + (0.5*(ny+1)-j)*(0.5*(ny+1)-j))*dx
                rjp1 = sqrt((0.5*(nx+1)-i)*(0.5*(nx+1)-i) + (0.5*(ny+1)-jp1)*(0.5*(ny+1)-jp1))*dx
                rjm1 = sqrt((0.5*(nx+1)-i)*(0.5*(nx+1)-i) + (0.5*(ny+1)-jm1)*(0.5*(ny+1)-jm1))*dx
        
                ! === Check direction ===
                ! x-direction
                if (rip1 .ge. 750e3) then
                    ! border point right
                    cr_acx(i,j)   = -u_acx(i,j)
                else
                    cr_acx(i,j)   = 0.0_wp
                end if

                if (rim1 .ge. 750e3) then
                    ! border point left
                    cr_acx(im1,j) = -u_acx(im1,j)
                else
                    cr_acx(im1,j) = 0.0_wp
                end if

                ! y-direction
                if (rjp1 .ge. 750e3) then
                    ! border point top
                    cr_acy(i,j)   = -v_acy(i,j)
                else
                    !border point top
                    cr_acy(i,j)   = 0.0_wp
                end if

                if (rjm1 .ge. 750e3) then
                    ! border point bottom
                    cr_acy(i,jm1) = -v_acy(i,jm1)
                else
                    ! border point bottom
                    cr_acy(i,jm1) = 0.0_wp
                end if
                
            end if
    
        end do
        end do
    
        return
    
    end subroutine calvmip_exp1
    
    subroutine calvmip_exp2(cr_acx,cr_acy,u_acx,v_acy,time,boundaries)
        ! Experiment 2 & 4 of CalvMIP
    
        implicit none
    
        real(wp), intent(OUT) :: cr_acx(:,:), cr_acy(:,:)
        real(wp), intent(IN)  :: u_acx(:,:),  v_acy(:,:)
        real(wp), intent(IN)  :: time
        character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose
    
        ! local variables
        integer  :: i, j, ip1, im1, jp1, jm1, nx, ny
        real(wp) :: wv, uxy_acx, uxy_acy, u_acy, v_acx
        integer  :: BC

        real(wp), parameter :: pi = acos(-1.0)  ! Calculate pi intrinsically

        nx = size(u_acx,1)
        ny = size(u_acx,2) 
        
        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Initialize    
        wv      = -300.0 * sin(2.0 * pi * time / 1000.0) 
        uxy_acx = 0.0_wp
        uxy_acy = 0.0_wp
        u_acy   = 0.0_wp
        v_acx   = 0.0_wp
        cr_acx  = 0.0_wp
        cr_acy  = 0.0_wp
    
        do j = 1, ny
        do i = 1, nx
            ! Stagger velocities x/y ac-velocities into y/x ac-nodes
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            u_acy = 0.25_wp*(u_acx(i,j)+u_acx(im1,j)+u_acx(im1,jp1)+u_acx(i,jp1))
            v_acx = 0.25_wp*(v_acy(i,j)+v_acy(i,jm1)+v_acy(ip1,jm1)+v_acy(ip1,j))
            ! x-direction
            uxy_acx     = MAX(1e-8,(u_acx(i,j)**2 + v_acx**2)**0.5)
            cr_acx(i,j) = -u_acx(i,j)+(u_acx(i,j)/uxy_acx)*wv
            ! y-direction
            uxy_acy     = MAX(1e-8,(v_acy(i,j)**2 + u_acy**2)**0.5)
            cr_acy(i,j) = -v_acy(i,j)+(v_acy(i,j)/uxy_acy)*wv
        end do
        end do

        return
    
    end subroutine calvmip_exp2

    subroutine calvmip_exp5_aa(cr_acx,cr_acy,u_acx,v_acy,H_ice,H_ice_c,f_ice,boundaries)
        ! Experiment 5 of CalvMIP
            
        implicit none
            
        real(wp), intent(OUT) :: cr_acx(:,:), cr_acy(:,:)
        real(wp), intent(IN)  :: u_acx(:,:),  v_acy(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: H_ice_c
        real(wp), intent(IN)  :: f_ice(:,:)
        character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose
                
        ! local variables
        integer  :: i, j, ip1, im1, jp1, jm1, nx, ny
        real(wp) :: uxy_aa, uxy_acx, uxy_acy, u_acy, v_acx
        real(wp), allocatable :: H_ice_fill(:,:), wv_aa(:,:) 
        integer  :: BC

        nx = size(u_acx,1)
        ny = size(u_acx,2) 
        allocate(H_ice_fill(nx,ny))
        allocate(wv_aa(nx,ny))

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Initialize    
        uxy_aa     = 0.0_wp
        H_ice_fill = H_ice
        wv_aa      = 0.0_wp
                
        do j = 1, ny
        do i = 1, nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            ! velocity on aa-node
            uxy_aa     = ((0.5*(u_acx(i,j)+u_acx(im1,j)))**2 + (0.5*(v_acy(i,j)+v_acy(i,jm1)))**2)**0.5
            wv_aa(i,j) = MAX(0.0_wp,1.0_wp+(H_ice_c-H_ice_fill(i,j))/H_ice_c)*uxy_aa                
        end do
        end do

        do j = 1, ny
        do i = 1, nx
            ! Stagger velocities x/y ac-velocities into y/x ac-nodes
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            u_acy = 0.25_wp*(u_acx(i,j)+u_acx(im1,j)+u_acx(im1,jp1)+u_acx(i,jp1))
            v_acx = 0.25_wp*(v_acy(i,j)+v_acy(i,jm1)+v_acy(ip1,jm1)+v_acy(ip1,j))
            ! x-direction
            uxy_acx     = MAX(1e-8,(u_acx(i,j)**2 + v_acx**2)**0.5)
            cr_acx(i,j) = -(u_acx(i,j)/uxy_acx)*0.5*(wv_aa(i,j)+wv_aa(ip1,j))
            ! y-direction
            uxy_acy     = MAX(1e-8,(v_acy(i,j)**2 + u_acy**2)**0.5)
            cr_acy(i,j) = -(v_acy(i,j)/uxy_acy)*0.5*(wv_aa(i,j)+wv_aa(i,jp1))
        end do
        end do

        deallocate(H_ice_fill)
        deallocate(wv_aa)
    
        return
        
    end subroutine calvmip_exp5_aa
    
end module calving_ac
