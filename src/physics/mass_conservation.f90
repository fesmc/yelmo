module mass_conservation

    use yelmo_defs, only : sp, dp, wp, TOL, TOL_UNDERFLOW, MISSING_VALUE, io_unit_err, &
                           MASK_ICE_NONE, MASK_ICE_FIXED, MASK_ICE_DYNAMIC, A_FRONT_MIN
    use yelmo_tools, only : boundary_code, get_neighbor_indices_bc_codes, &
                            fill_borders_2D, set_boundaries_2D_aa

    use solver_advection, only : calc_advec2D  
    use velocity_general, only : set_inactive_margins 
    use topography, only : calc_ice_fraction, calc_front_cells, calc_front_H_ref

    implicit none 

    private

    public :: check_mass_conservation
    public :: apply_tendency
    public :: calc_G_advec_simple
    public :: calc_G_advec
    public :: calc_G_mbal
    public :: calc_G_calv
    public :: calc_G_boundaries
    public :: set_tau_relax
    public :: calc_G_relaxation

    public :: calc_G_remove_fractional_ice
    public :: calc_G_front_advance
    public :: calc_G_calving_front
    public :: calc_G_lsf_front
    public :: remove_icebergs
    
contains 
    
    subroutine check_mass_conservation(H_ice,f_ice,f_grnd,dHidt,mb_net,cmb,dHidt_dyn,smb,bmb,fmb,dmb, &
                                                                mb_resid,dx,sec_year,time,dt,units,label)

        implicit none

        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: dHidt(:,:)
        real(wp), intent(IN) :: mb_net(:,:)
        real(wp), intent(IN) :: cmb(:,:)
        real(wp), intent(IN) :: dHidt_dyn(:,:)
        real(wp), intent(IN) :: smb(:,:)
        real(wp), intent(IN) :: bmb(:,:)
        real(wp), intent(IN) :: fmb(:,:)
        real(wp), intent(IN) :: dmb(:,:)
        real(wp), intent(IN) :: mb_resid(:,:)
        real(wp), intent(IN) :: dx 
        real(wp), intent(IN) :: sec_year
        real(wp), intent(IN) :: time
        real(wp), intent(IN) :: dt
        character(len=*), intent(IN) :: units 
        character(len=*), intent(IN) :: label 

        ! Local variables
        integer  :: npts
        ! Totals are accumulated in double precision, so that the check
        ! does not add summation round-off of its own
        real(dp) :: tot_dHidt
        real(dp) :: tot_components
        real(dp) :: tot_mb_net
        real(dp) :: tot_cmb
        real(dp) :: tot_dHidt_dyn
        real(dp) :: tot_gross
        real(dp) :: conv
        real(dp) :: resid
        real(dp) :: resid_rel
        character(len=4) :: flag

        ! Tolerance on the relative residual: with every tendency passed through
        ! apply_tendency, the residual is only working-precision round-off of the
        ! per-cell updates, found to be <~0.3*epsilon(wp) of the gross throughput
        ! (EISMINT, TROUGH). 100*epsilon leaves ample headroom while still
        ! flagging any real leak (1e-5 of the throughput for wp=sp).
        real(dp), parameter :: tol_rel = 100.0_dp*epsilon(1.0_wp)

        ! Determine conversion factor to units of interest from [m^3/yr]

        select case(trim(units))

            case("m^3/yr")

                conv = 1.0_dp

            case("km^3/yr")

                conv = 1e-9_dp

            case("Sv")

                conv = 1e-6_dp / real(sec_year,dp)

            case DEFAULT

                write(io_unit_err,*) "check_mass_conservation:: Error: units not recognized."
                write(io_unit_err,*) "units = ", trim(units)
                stop

        end select

        ! Calculate totals, initially [m^3/yr] => [units]
 
        tot_dHidt       = sum(real(dHidt,dp))     * real(dx,dp)**2 * conv
        tot_mb_net      = sum(real(mb_net,dp))    * real(dx,dp)**2 * conv
        tot_cmb         = sum(real(cmb,dp))       * real(dx,dp)**2 * conv
        tot_dHidt_dyn   = sum(real(dHidt_dyn,dp)) * real(dx,dp)**2 * conv

        ! Gross throughput: sum of the magnitudes of the same component fluxes
        ! per cell, so opposing fluxes do not cancel (non-zero at equilibrium)
        tot_gross       = sum(abs(real(dHidt_dyn,dp))+abs(real(mb_net,dp))+abs(real(cmb,dp))) &
                                                        * real(dx,dp)**2 * conv

        ! Get total of components and residual, absolute [units] and relative
        ! to the gross throughput
        ! (dHidt_dyn integrates to the net flux across the domain boundary,
        ! so it must be included for the budget to close)
        tot_components = tot_dHidt_dyn + tot_mb_net + tot_cmb
        resid          = tot_components - tot_dHidt
        resid_rel      = resid / max(tot_gross,tiny(tot_gross))

        flag = ""
        if (abs(resid_rel) .gt. tol_rel) flag = "FAIL"

        write(*,"(a8,a,2f9.3,a3,4g14.4,1x,a4,a3,3g13.4)") &
                    trim(label), " mbcheck ["//trim(units)//"]: ", time, dt, " | ", &
                    tot_dHidt, tot_components, resid, resid_rel, flag, " | ", &
                    tot_dHidt_dyn, tot_mb_net, tot_cmb

        return

    end subroutine check_mass_conservation

    subroutine apply_tendency(H_ice,mb_dot,dt,label,adjust_mb)

        implicit none

        real(wp), intent(INOUT) :: H_ice(:,:)
        real(wp), intent(INOUT) :: mb_dot(:,:)
        real(wp), intent(IN)    :: dt 
        character(len=*),  intent(IN) :: label 
        logical,  intent(IN), optional :: adjust_mb
        
        ! Local variables
        integer :: i, j, nx, ny
        real(wp) :: H_prev
        real(wp) :: dHdt 
        logical  :: allow_adjust_mb

        if (dt .gt. 0.0) then 
            ! Only apply this routine if dt > 0!

            allow_adjust_mb = .FALSE.
            if (present(adjust_mb)) allow_adjust_mb = adjust_mb 

            nx = size(H_ice,1)
            ny = size(H_ice,2)

            do j = 1, ny 
            do i = 1, nx 

                ! Store previous ice thickness
                H_prev = H_ice(i,j) 

                ! Now update ice thickness with tendency for this timestep 
                H_ice(i,j) = H_prev + dt*mb_dot(i,j)

                ! Limit ice thickness to zero 
                if (H_ice(i,j) .lt. 0.0) H_ice(i,j) = 0.0 

                ! Ensure tiny numeric ice thicknesses are removed
                if (abs(H_ice(i,j)) .lt. TOL) H_ice(i,j) = 0.0
                
                ! Calculate actual current rate of change
                dHdt = (H_ice(i,j) - H_prev) / dt 

                ! Update mb rate to match ice rate of change perfectly
                if (allow_adjust_mb) then
                    mb_dot(i,j) = dHdt
                end if 

            end do
            end do

        end if 

        return

    end subroutine apply_tendency

    subroutine calc_G_advec_simple(G_advec,H_ice,f_ice,ux,uy,mask_ice, &
                                                    solver,boundaries,dx,dt,F)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp),         intent(OUT)   :: G_advec(:,:)         ! [m/yr] Tendency due to advection
        real(wp),         intent(IN)    :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp),         intent(IN)    :: f_ice(:,:)           ! [--]  Ice area fraction 
        real(wp),         intent(IN)    :: ux(:,:)              ! [m/a] Depth-averaged velocity, x-direction (ac-nodes)
        real(wp),         intent(IN)    :: uy(:,:)              ! [m/a] Depth-averaged velocity, y-direction (ac-nodes)
        integer,          intent(IN)    :: mask_ice(:,:)        ! Advection mask
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries
        real(wp),         intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(wp),         intent(IN)    :: dt                   ! [a]   Timestep 
        real(wp),         intent(IN), optional :: F(:,:) 

        ! Local variables 
        integer :: i, j, nx, ny
        real(wp), allocatable :: F_now(:,:) 
        real(wp), allocatable :: ux_tmp(:,:) 
        real(wp), allocatable :: uy_tmp(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(F_now(nx,ny))
        allocate(ux_tmp(nx,ny))
        allocate(uy_tmp(nx,ny))

        ! Set local velocity fields with no margin treatment intially
        ux_tmp = ux
        uy_tmp = uy
        
        F_now = 0.0_wp 
        if (present(F)) F_now = F 

        ! Ensure that no velocity is defined for outer boundaries of partially-filled margin points
        call set_inactive_margins(ux_tmp,uy_tmp,f_ice,boundaries)

        ! Determine current advective rate of change (time=n)
        call calc_advec2D(G_advec,H_ice,f_ice,ux_tmp,uy_tmp,F_now,mask_ice,dx,dx,dt,solver,boundaries)

        return 

    end subroutine calc_G_advec_simple

    subroutine calc_G_advec(G_adv,dHdt_n,H_ice_n,H_ice_pred,H_ice,f_ice,ux,uy, &
                        mask_pred_new,mask_corr_new,solver,mask_ice,boundaries, &
                        dx,dt,beta,pc_step,F)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp),         intent(OUT)   :: G_adv(:,:)           ! [m/yr] Tendency due to advection
        real(wp),         intent(INOUT) :: dHdt_n(:,:)          ! [m/a] Advective rate of ice thickness change from previous=>current timestep 
        real(wp),         intent(INOUT) :: H_ice_n(:,:)         ! [m]   Ice thickness from previous=>current timestep 
        real(wp),         intent(IN)    :: H_ice_pred(:,:)      ! [m]   Ice thickness from predicted timestep 
        real(wp),         intent(IN)    :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp),         intent(IN)    :: f_ice(:,:)           ! [--]  Ice area fraction 
        real(wp),         intent(IN)    :: ux(:,:)              ! [m/a] Depth-averaged velocity, x-direction (ac-nodes)
        real(wp),         intent(IN)    :: uy(:,:)              ! [m/a] Depth-averaged velocity, y-direction (ac-nodes)
        integer,          intent(IN)    :: mask_pred_new(:,:)   
        integer,          intent(IN)    :: mask_corr_new(:,:)  
        integer,          intent(IN)    :: mask_ice(:,:)        ! Advection mask  
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries
        real(wp),         intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(wp),         intent(IN)    :: dt                   ! [a]   Timestep 
        real(wp),         intent(IN)    :: beta(4)              ! Timestep weighting parameters
        character(len=*), intent(IN)    :: pc_step              ! Current predictor-corrector step ('predictor' or 'corrector')
        real(wp),         intent(IN), optional :: F(:,:) 
        
        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1  
        real(wp), allocatable :: F_now(:,:) 
        real(wp), allocatable :: dHdt_advec(:,:) 
        real(wp), allocatable :: ux_tmp(:,:) 
        real(wp), allocatable :: uy_tmp(:,:) 

        real(wp), parameter :: dHdt_advec_lim = 10.0_wp     ! [m/a] Hard limit on advection rate

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(F_now(nx,ny))
        allocate(ux_tmp(nx,ny))
        allocate(uy_tmp(nx,ny))
        allocate(dHdt_advec(nx,ny))

        dHdt_advec = 0.0_wp 

        F_now = 0.0_wp 
        if (present(F)) F_now = F 

        ! Set local velocity fields with no margin treatment intially
        ux_tmp = ux
        uy_tmp = uy
        
        ! ===================================================================================
        ! Resolve the dynamic part (ice advection) using multistep method

        select case(trim(pc_step))
        
            case("predictor") 
                
                ! Fill velocity field for new cells 
                call fill_vel_new_cells(ux_tmp,uy_tmp,mask_corr_new,boundaries)

                ! Ensure that no velocity is defined for outer boundaries of partially-filled margin points
                call set_inactive_margins(ux_tmp,uy_tmp,f_ice,boundaries)

                ! Store ice thickness from time=n
                H_ice_n   = H_ice 

                ! Store advective rate of change from saved from previous timestep (now represents time=n-1)
                dHdt_advec = dHdt_n 

                ! Determine current advective rate of change (time=n)
                call calc_advec2D(dHdt_n,H_ice,f_ice,ux_tmp,uy_tmp,F_now,mask_ice,dx,dx,dt,solver,boundaries)

                ! Calculate rate of change using weighted advective rates of change 
                dHdt_advec = beta(1)*dHdt_n + beta(2)*dHdt_advec 
                
                ! Calculate predicted ice thickness (time=n+1,pred)
                !H_ice = H_ice_n + dt*dHdt_advec 

            case("corrector") ! corrector 

                ! Fill velocity field for new cells 
                call fill_vel_new_cells(ux_tmp,uy_tmp,mask_pred_new,boundaries)

                ! Ensure that no velocity is defined for outer boundaries of partially-filled margin points
                call set_inactive_margins(ux_tmp,uy_tmp,f_ice,boundaries)

                ! Determine advective rate of change based on predicted H,ux/y fields (time=n+1,pred)
                call calc_advec2D(dHdt_advec,H_ice_pred,f_ice,ux_tmp,uy_tmp,F_now,mask_ice,dx,dx,dt,solver,boundaries)

                ! Calculate rate of change using weighted advective rates of change 
                dHdt_advec = beta(3)*dHdt_advec + beta(4)*dHdt_n 
                
                ! Calculate corrected ice thickness (time=n+1)
                !H_ice = H_ice_n + dt*dHdt_advec 

                ! Finally, update dHdt_n with correct term to use as n-1 on next iteration
                dHdt_n = dHdt_advec 

        end select
        
        ! Store advective tendency 
        G_adv = dHdt_advec 

        return 

    end subroutine calc_G_advec

    subroutine fill_vel_new_cells(ux,uy,mask,boundaries)

        implicit none

        real(wp), intent(INOUT) :: ux(:,:) 
        real(wp), intent(INOUT) :: uy(:,:) 
        integer,  intent(IN)    :: mask(:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        integer :: BC

        nx = size(mask,1)
        ny = size(mask,2) 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            
            if (mask(i,j) .eq. 2) then 
                ! This site just filled with ice, so 
                ! velocity may not be defined on borders
                ! Check borders and fill in empty velocities

                ! x-direction
                if (ux(i,j) .eq. 0.0_wp .and. ux(im1,j) .ne. 0.0_wp) then 
                    ux(i,j) = ux(im1,j) 
                else if (ux(i,j) .ne. 0.0_wp .and. ux(im1,j) .eq. 0.0_wp) then
                    ux(im1,j) = ux(i,j)
                end if 

                ! y-direction
                if (uy(i,j) .eq. 0.0_wp .and. uy(i,jm1) .ne. 0.0_wp) then 
                    uy(i,j) = uy(i,jm1) 
                else if (uy(i,j) .ne. 0.0_wp .and. uy(i,jm1) .eq. 0.0_wp) then
                    uy(i,jm1) = uy(i,j)
                end if 
                
            end if 

        end do 
        end do


        return

    end subroutine fill_vel_new_cells

    subroutine calc_G_mbal(G_mb,H_ice,f_grnd,mbal,dt,f_ice)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp), intent(OUT)   :: G_mb(:,:)            ! [m/yr] Actual tendency due to mass balance
        real(wp), intent(IN)    :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp), intent(IN)    :: f_grnd(:,:)          ! [--]  Grounded fraction 
        real(wp), intent(IN)    :: mbal(:,:)            ! [m/yr] Net mass balance; mbal = smb+bmb+fmb+calv
        real(wp), intent(IN)    :: dt                   ! [a]   Timestep  
        real(wp), intent(IN), optional :: f_ice(:,:)    ! [--]  Ice area fraction: the rate acts on the covered area of partial cells

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! ==== MASS BALANCE =====

        ! Initialize G_mb object with diagnosed mass balance everywhere
        G_mb = mbal

        ! Partial cells: the per-area rate applies to the covered area only
        if (present(f_ice)) then
            where (f_ice .gt. 0.0_wp .and. f_ice .lt. 1.0_wp) G_mb = f_ice*G_mb
        end if

        ! Ensure melting is only counted where ice exists 
        where(G_mb .lt. 0.0 .and. H_ice .eq. 0.0) G_mb = 0.0 

        ! Additionally ensure ice cannot form in open ocean 
        where(f_grnd .eq. 0.0 .and. H_ice .eq. 0.0)  G_mb = 0.0  

        ! Ensure melt is limited to amount of available ice to melt  
        where((H_ice+dt*G_mb) .lt. 0.0) G_mb = -H_ice/dt

        return 

    end subroutine calc_G_mbal

    subroutine calc_G_calv(G_calv,H_ice,calv_flt,calv_grnd,dt,calv_flt_method,boundaries)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp), intent(OUT)   :: G_calv(:,:)          ! [m/yr] Actual calving rate applied to real ice points
        real(wp), intent(IN)    :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp), intent(IN)    :: calv_flt(:,:)        ! [m/a] Potential calving rate (floating)
        real(wp), intent(IN)    :: calv_grnd(:,:)       ! [m/a] Potential calving rate (grounded)
        real(wp), intent(IN)    :: dt                   ! [a]   Timestep   
        character(len=*), intent(IN) :: calv_flt_method
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        logical :: is_margin
        logical :: kill_floating
        real(wp) :: calv_flt_now 
        real(wp) :: calv_grnd_now 
        integer :: BC

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Determine whether a kill method is being applied
        kill_floating = .FALSE. 
        if (trim(calv_flt_method) .eq. "kill" .or. &
            trim(calv_flt_method) .eq. "kill-pos") then 

            kill_floating = .TRUE. 

        end if 

        ! ===== CALVING ======

        ! Combine grounded and floating calving into one field for output.
        ! It has already been scaled by area of ice in cell (f_ice).

        ! Note 1: Only allow calving at the current margin 
        ! If ice has retreated before applying calving, then H_ice is 
        ! zero and so G_calv will also be zero. But if ice has advanced,
        ! then calving should also go to zero. 

        ! Note 2: for floating ice, allow calving everywhere if kill_floating is active.

        G_calv = 0.0_wp 
        
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            
            is_margin = H_ice(i,j) .gt. 0.0 .and. &
                count([H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)].eq.0.0) .gt. 0

            if (is_margin .or. kill_floating) then
                calv_flt_now = calv_flt(i,j)
            else
                calv_flt_now = 0.0_wp
            end if 

            if (is_margin) then
                calv_grnd_now = calv_grnd(i,j) 
            else
                calv_grnd_now = 0.0_wp
            end if

            ! Calculate calving rate tendency (negative == mass loss)
            G_calv(i,j) = -(calv_flt_now + calv_grnd_now)
            
            ! Limit calving rate to available ice
            if (H_ice(i,j)+dt*G_calv(i,j) .lt. 0.0) G_calv(i,j) = -H_ice(i,j)/dt
            
        end do 
        end do

        return 

    end subroutine calc_G_calv

    subroutine calc_G_boundaries(mb_resid,H_ice,H_eff,mask_cf,f_grnd,uxy_b,mask_ice,boundaries, &
                                                            H_ice_ref,H_min_flt,H_min_grnd,tau,dt)

        implicit none

        real(wp),           intent(INOUT)   :: mb_resid(:,:)            ! [m/yr] Residual mass balance
        real(wp),           intent(IN)      :: H_ice(:,:)               ! [m] Ice thickness
        real(wp),           intent(IN)      :: H_eff(:,:)               ! [m] Effective ice thickness (tpo%now%H_eff)
        logical,            intent(IN)      :: mask_cf(:,:)             ! Subgrid front cells (ytopo.front_subgrid)
        real(wp),           intent(IN)      :: f_grnd(:,:)              ! [--] Grounded ice fraction
        real(wp),           intent(IN)      :: uxy_b(:,:)               ! [m/a] Basal sliding speed, aa-nodes
        integer,            intent(IN)      :: mask_ice(:,:)            ! Mask: MASK_ICE_NONE forced zero, MASK_ICE_FIXED imposed (=H_ice_ref), MASK_ICE_DYNAMIC active
        character(len=*),   intent(IN)      :: boundaries               ! Boundary condition choice
        real(wp),           intent(IN)      :: H_ice_ref(:,:)           ! [m]  Reference ice thickness to fill with for boundaries=="fixed"
        real(wp),           intent(IN)      :: H_min_flt                ! [m] Minimum allowed floating ice thickness 
        real(wp),           intent(IN)      :: H_min_grnd               ! [m] Minimum allowed grounded ice thickness 
        real(wp),           intent(IN)      :: tau                      ! [yr] Timescale for removing margin ice thinner than H_min_flt/H_min_grnd
        real(wp),           intent(IN)      :: dt                       ! [yr] Timestep

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        real(wp), allocatable :: H_ice_new(:,:)
        real(wp), allocatable :: H_tmp(:,:)
        real(wp) :: H_max 
        logical  :: is_margin 
        logical  :: is_island 
        logical  :: is_isthmus_x 
        logical  :: is_isthmus_y 
        real(wp) :: f_rm 
        integer  :: BC

        real(wp), parameter :: H_min_tol = 1e-6

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        allocate(H_tmp(nx,ny)) 
        allocate(H_ice_new(nx,ny)) 
        H_ice_new = H_ice 

        ! Apply special case for symmetric EISMINT domain when basal sliding is active
        ! (ensure summit thickness does not grow disproportionately)
        if (trim(boundaries) .eq. "EISMINT" .and. maxval(uxy_b) .gt. 0.0) then 
            i = (nx-1)/2 
            j = (ny-1)/2
            H_ice_new(i,j) = (H_ice(i-1,j)+H_ice(i+1,j) &
                                    +H_ice(i,j-1)+H_ice(i,j+1)) / 4.0 
        end if  
        
        ! Remove margin points that are too thin, or points that are below tolerance ====

        ! Too-thin margin ice is removed at the rate H/tau, so that the removal
        ! per year does not depend on the timestep (all of it when dt >= tau)
        f_rm = 1.0_wp
        if (tau .gt. dt) f_rm = dt/tau

        H_tmp = H_ice_new 

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,is_margin)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            
            is_margin = H_tmp(i,j) .gt. 0.0 .and. &
                count([H_tmp(im1,j),H_tmp(ip1,j),H_tmp(i,jm1),H_tmp(i,jp1)].eq.0.0) .gt. 0

            if (is_margin) then
                ! Ice covered point at the margin

                ! Remove ice that is too thin (effective thickness: the
                ! thickness of a partial front cell over its covered area)
                if ( (f_grnd(i,j) .eq. 0.0_wp .and. H_eff(i,j) .lt. H_min_flt) .or. &
                     (f_grnd(i,j) .gt. 0.0_wp .and. H_eff(i,j) .lt. H_min_grnd) ) then
                    H_ice_new(i,j) = (1.0_wp-f_rm)*H_ice_new(i,j)
                end if
 
            end if 

            ! Also remove very thin ice points (eg 1e-6 m thick - thicker than machine tolerance, but thinner than relevant)
            ! E.g., for a very small timestep of dt=1e-3 and an accumulation rate of 0.1 m/yr, 
            ! after one timestep, H = 1e-3*0.1 = 1e-4 m. 
            if (H_tmp(i,j) .lt. H_min_tol) H_ice_new(i,j) = 0.0_wp 

        end do 
        end do
        !$omp end parallel do

        ! Remove ice islands =====

        H_tmp = H_ice_new 

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,is_island)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            
            ! Check for ice islands
            is_island = H_tmp(i,j) .gt. 0.0 .and. &
                count([H_tmp(im1,j),H_tmp(ip1,j),H_tmp(i,jm1),H_tmp(i,jp1)].gt.0.0) .eq. 0

            if (is_island) then 
                ! Ice-covered island
                ! Remove ice completely. 

                H_ice_new(i,j)   = 0.0_wp 

            end if 

        end do 
        end do
        !$omp end parallel do

        ! Reduce ice thickness for margin points that are thicker 
        ! than inland neighbors. Not for subgrid front cells: there
        ! the excess above H_eff is moved on by the front advance ====

        H_tmp = H_ice_new

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,is_margin,H_max)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            
            is_margin = H_tmp(i,j) .gt. 0.0 .and. .not. mask_cf(i,j) .and. &
                count([H_tmp(im1,j),H_tmp(ip1,j),H_tmp(i,jm1),H_tmp(i,jp1)].eq.0.0) .gt. 0

            if (is_margin) then
                ! Ice covered point at the margin

                ! Calculate maximum thickness of neighbors 
                H_max = maxval([H_tmp(im1,j),H_tmp(ip1,j), &
                                H_tmp(i,jm1),H_tmp(i,jp1)])

                if (H_ice_new(i,j) .gt. H_max) H_ice_new(i,j) = H_max
                
            end if
            
        end do 
        end do
        !$omp end parallel do
        
        select case(trim(boundaries))

            case("MISMIP3D","TROUGH")

                ! Ensure that x-boundary ice-thickness matches "infinite" and "zero" case.
                ! This is not fully handled by ytopo.solver="impl-lis", since
                ! subsequently the mb forcing is applied.

                H_ice_new(1,:)    = H_ice_new(2,:)          ! x=0, Symmetry 
                H_ice_new(nx,:)   = 0.0                     ! x=max, no ice

            case("periodic","periodic-xy")

                ! Do nothing - this should be handled by the ice advection routine
                ! if the default choice ytopo.solver="impl-lis" is used.

            case("periodic-x")
                ! Periodic x: nothing to do (handled by the ice advection routine).
                ! Infinite y: set border points equal to inner neighbors, as
                ! for "infinite", since subsequently the mb forcing is applied.

                H_ice_new(:,1)  = H_ice_new(:,2)
                H_ice_new(:,ny) = H_ice_new(:,ny-1)

            case("infinite")
                ! Set border points equal to inner neighbors 
                ! This is not fully handled by ytopo.solver="impl-lis", since
                ! subsequently the mb forcing is applied.
                
                H_ice_new(1,:)    = H_ice_new(2,:)          ! x=0, copy inner neighbor
                H_ice_new(nx,:)   = H_ice_new(nx-1,:)       ! x=max, copy inner neighbor
                H_ice_new(:,1)  = H_ice_new(:,2)
                H_ice_new(:,ny) = H_ice_new(:,ny-1)

            case("fixed") 
                ! Set border points equal to prescribed values from array

                ! Do nothing - this should be handled by the ice advection routine
                ! if the default choice ytopo.solver="impl-lis" is used.

                !call fill_borders_2D(H_ice_new,nfill=1,fill=H_ice_ref)

            case("zeros","EISMINT")
                ! Force zero ice thickness on all four grid borders
                ! (overrides bnd%mask_ice on border cells)

                H_ice_new(1,:)  = 0.0
                H_ice_new(nx,:) = 0.0

                H_ice_new(:,1)  = 0.0
                H_ice_new(:,ny) = 0.0

            case("mask")
                ! Defer border treatment entirely to bnd%mask_ice
                ! (handled by the universal mask enforcement below).

            case DEFAULT    ! e.g., None/none
                ! No special border treatment; rely on bnd%mask_ice
                ! (handled by the universal mask enforcement below).

        end select

        ! Impose the per-cell ice-domain mask (bnd%mask_ice) universally,
        ! independent of the boundary BC choice above: MASK_ICE_NONE forces zero
        ! ice thickness and MASK_ICE_FIXED imposes the reference thickness.
        ! mask_ice is a per-cell constraint, not a border treatment, so it must
        ! be honored for every 'boundaries' choice - matching how the advection
        ! solver already treats mask_ice. Without this, a boundary choice like
        ! "zeros" only zeroes the outer domain borders, and positive smb re-grows
        ! ice in interior MASK_ICE_NONE cells that dynamics had zeroed each step.
        where (mask_ice .eq. MASK_ICE_NONE)  H_ice_new = 0.0_wp
        where (mask_ice .eq. MASK_ICE_FIXED) H_ice_new = H_ice_ref

        ! Determine rate of mass balance related to changes applied here.
        ! Where ice is removed completely (H_ice_new == 0), overshoot by 10%
        ! so that apply_tendency's clip-to-zero removes it robustly (no
        ! round-off remnant). Everywhere else (MASK_ICE_FIXED imposed values,
        ! margin reduction, boundary copies) use the exact rate, so that
        ! apply_tendency lands on H_ice_new (to round-off) without over- or
        ! undershooting it.
        if (dt .ne. 0.0) then
            where (H_ice_new .eq. 0.0_wp .and. mask_ice .ne. MASK_ICE_FIXED)
                mb_resid = 1.1_wp * (H_ice_new - H_ice) / dt
            elsewhere
                mb_resid = (H_ice_new - H_ice) / dt
            end where
        else
            mb_resid = 0.0
        end if

        return

    end subroutine calc_G_boundaries

    subroutine set_tau_relax(tau_relax,H_ice,f_grnd,mask_grz,H_ref,topo_rel,tau,boundaries)
        ! Build a per-cell relaxation timescale field from a topo_rel mask choice and
        ! a uniform tau. Only tau > 0 produces actual relaxation; to impose an ice
        ! thickness directly, set bnd%mask_ice == MASK_ICE_FIXED instead.

        implicit none 

        real(wp), intent(OUT)   :: tau_relax(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:)  
        integer,  intent(IN)    :: mask_grz(:,:) 
        real(wp), intent(IN)    :: H_ref(:,:) 
        integer,  intent(IN)    :: topo_rel 
        real(wp), intent(IN)    :: tau
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, nx, ny 
        integer  :: im1, ip1, jm1, jp1
        integer  :: BC

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            
            select case(topo_rel)

                case(1) 
                    ! Relax the shelf (floating) ice and ice-free points
                
                    if (f_grnd(i,j) .eq. 0.0 .or. H_ref(i,j) .eq. 0.0) then
                        tau_relax(i,j) = tau
                    else
                        tau_relax(i,j) = -1.0
                    end if
            
                case(2) 
                    ! Relax the shelf (floating) ice and ice-free points
                    ! and the grounding-line ice too
                    
                    if (f_grnd(i,j) .eq. 0.0 .or. H_ref(i,j) .eq. 0.0) then
                        tau_relax(i,j) = tau
                    
                    else if (f_grnd(i,j) .gt. 0.0 .and. &
                            (f_grnd(im1,j) .eq. 0.0 .or. f_grnd(ip1,j) .eq. 0.0 &
                            .or. f_grnd(i,jm1) .eq. 0.0 .or. f_grnd(i,jp1) .eq. 0.0)) then

                        tau_relax(i,j) = tau
                    
                    else

                        tau_relax(i,j) = -1.0

                    end if
            
                case(3)
                    ! Relax all points
                    
                    tau_relax(i,j) = tau
                
                case(4) 
                    ! Relax all grounded grounding-zone points 

                    if (mask_grz(i,j) .eq. 0 .or. mask_grz(i,j) .eq. 1) then 

                        tau_relax(i,j) = tau
                    
                    else

                        tau_relax(i,j) = -1.0

                    end if 

                case DEFAULT ! topo_rel == 0

                    ! No relaxation

                    tau_relax(i,j) = -1.0

            end select
            
        end do 
        end do 


        return 

    end subroutine set_tau_relax
    
    subroutine calc_G_relaxation(dHdt,H_ice,H_ref,tau_relax,dt)
        ! Relax ice toward a reference state with a finite timescale tau_relax > 0.
        ! Note: imposing the ice thickness directly is handled separately via
        ! mask_ice == MASK_ICE_FIXED in calc_G_boundaries (no longer via tau_relax == 0).

        implicit none

        real(wp), intent(OUT)   :: dHdt(:,:)
        real(wp), intent(IN)    :: H_ice(:,:)
        real(wp), intent(IN)    :: H_ref(:,:)
        real(wp), intent(IN)    :: tau_relax(:,:)
        real(wp), intent(IN)    :: dt

        ! Local variables
        integer  :: i, j, nx, ny

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        dHdt = 0.0

        !$omp parallel do collapse(2) private(i,j)
        do j = 1, ny
        do i = 1, nx

            if (tau_relax(i,j) .gt. 0.0) then
                ! Apply relaxation to reference state
                dHdt(i,j) = (H_ref(i,j) - H_ice(i,j)) / tau_relax(i,j)
            end if

        end do
        end do
        !$omp end parallel do

        return

    end subroutine calc_G_relaxation

    subroutine calc_G_front_advance(mb_adv,H_ice,H_eff,mask_cf,mask_ocn,ux,uy,dt,boundaries)
        ! Advance the ice front (CISM advance_calving_front): in a front cell
        ! that holds more ice than its effective thickness, move the excess
        ! (plus a small amount, so that the cell stays partial) to its
        ! ice-free ocean edge neighbours, split by the outward velocity across
        ! each face. Returns the rate [m/yr]; it sums to zero (transport).

        implicit none

        real(wp), intent(OUT) :: mb_adv(:,:)            ! [m/yr] Rate of thickness change
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m]    Ice thickness
        real(wp), intent(IN)  :: H_eff(:,:)             ! [m]    Effective thickness
        logical,  intent(IN)  :: mask_cf(:,:)           ! Front cells
        logical,  intent(IN)  :: mask_ocn(:,:)          ! Ice-free ocean cells
        real(wp), intent(IN)  :: ux(:,:)                ! [m/yr] Depth-averaged velocity (acx-nodes)
        real(wp), intent(IN)  :: uy(:,:)                ! [m/yr] Depth-averaged velocity (acy-nodes)
        real(wp), intent(IN)  :: dt                     ! [yr]
        character(len=*), intent(IN) :: boundaries

        ! Local variables
        integer  :: i, j, k, nx, ny, im1, ip1, jm1, jp1, BC
        integer  :: in(4), jn(4)
        real(wp) :: u_out(4), u_tot, dH_tot

        real(wp), parameter :: dH_small = 0.1_wp       ! [m] Leaves the cell slightly below H_eff

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        BC = boundary_code(boundaries)

        mb_adv = 0.0_wp

        if (dt .le. 0.0_wp) return

        ! Serial loop: a receiving cell can border several front cells
        do j = 1, ny
        do i = 1, nx

            if (.not. mask_cf(i,j) .or. H_ice(i,j) .le. H_eff(i,j)) cycle

            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            in    = [ip1,im1,i,i]
            jn    = [j,j,jp1,jm1]
            u_out = [ux(i,j), -ux(im1,j), uy(i,j), -uy(i,jm1)]

            u_tot = 0.0_wp
            do k = 1, 4
                if (mask_ocn(in(k),jn(k)) .and. u_out(k) .gt. 0.0_wp) u_tot = u_tot + u_out(k)
            end do

            if (u_tot .le. 0.0_wp) cycle

            dH_tot = min(H_ice(i,j) - H_eff(i,j) + dH_small, H_ice(i,j))

            do k = 1, 4
                if (mask_ocn(in(k),jn(k)) .and. u_out(k) .gt. 0.0_wp) then
                    mb_adv(in(k),jn(k)) = mb_adv(in(k),jn(k)) + dH_tot*(u_out(k)/u_tot)/dt
                end if
            end do
            mb_adv(i,j) = mb_adv(i,j) - dH_tot/dt

        end do
        end do

        return

    end subroutine calc_G_front_advance

    subroutine calc_G_calving_front(cmb,H_ice,cmb_dem,mask_cf,mask_elig,mask_ocn,ux,uy,dt,boundaries)
        ! Apply a calving demand at the ice front (CISM apply_calving_dthck).
        ! cmb_dem is the thinning rate for one exposed face (H_eff*c/dx); it
        ! is scaled by the front length, L/dx = 1, sqrt(2), 2 for 1, 2, >=3
        ! ocean faces. Demand beyond a front cell's ice is taken from its
        ! upstream eligible edge neighbours, split by inflow (one level; any
        ! rest is not calved). Returns the applied calving rate [m/yr, <= 0].

        implicit none

        real(wp), intent(OUT) :: cmb(:,:)               ! [m/yr] Applied calving rate
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m]    Ice thickness
        real(wp), intent(IN)  :: cmb_dem(:,:)           ! [m/yr] Calving demand (<= 0), one face
        logical,  intent(IN)  :: mask_cf(:,:)           ! Front cells
        logical,  intent(IN)  :: mask_elig(:,:)         ! Eligible (marine) ice cells
        logical,  intent(IN)  :: mask_ocn(:,:)          ! Ice-free ocean cells
        real(wp), intent(IN)  :: ux(:,:)                ! [m/yr] Depth-averaged velocity (acx-nodes)
        real(wp), intent(IN)  :: uy(:,:)                ! [m/yr] Depth-averaged velocity (acy-nodes)
        real(wp), intent(IN)  :: dt                     ! [yr]
        character(len=*), intent(IN) :: boundaries

        ! Local variables
        integer  :: i, j, k, nx, ny, im1, ip1, jm1, jp1, BC, n_ocn
        integer  :: in(4), jn(4)
        real(wp) :: q_in(4), q_tot, dH
        real(wp), allocatable :: H_now(:,:), rest(:,:)

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        BC = boundary_code(boundaries)

        cmb = 0.0_wp
        if (dt .le. 0.0_wp) return

        allocate(H_now(nx,ny), rest(nx,ny))
        H_now = H_ice
        rest  = 0.0_wp

        ! Front cells: calve up to the full column, keep the rest
        do j = 1, ny
        do i = 1, nx
            if (.not. mask_cf(i,j) .or. cmb_dem(i,j) .ge. 0.0_wp) cycle
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            n_ocn = count([mask_ocn(im1,j),mask_ocn(ip1,j),mask_ocn(i,jm1),mask_ocn(i,jp1)])
            dH = -cmb_dem(i,j)*dt
            if (n_ocn .eq. 2) dH = dH*sqrt(2.0_wp)
            if (n_ocn .ge. 3) dH = dH*2.0_wp
            if (dH .gt. H_now(i,j)) then
                rest(i,j)  = dH - H_now(i,j)
                H_now(i,j) = 0.0_wp
            else
                H_now(i,j) = H_now(i,j) - dH
            end if
        end do
        end do

        ! Rest: calve upstream eligible (non-front) edge neighbours, split by inflow
        do j = 1, ny
        do i = 1, nx
            if (rest(i,j) .le. 0.0_wp) cycle
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            in = [im1,ip1,i,i]
            jn = [j,j,jm1,jp1]
            ! Volume flux into (i,j) across each face [m2/yr per unit width]
            q_in = [ max( ux(im1,j),0.0_wp)*H_ice(im1,j), max(-ux(i,j),0.0_wp)*H_ice(ip1,j), &
                     max( uy(i,jm1),0.0_wp)*H_ice(i,jm1), max(-uy(i,j),0.0_wp)*H_ice(i,jp1) ]
            do k = 1, 4
                if (.not. mask_elig(in(k),jn(k)) .or. mask_cf(in(k),jn(k))) q_in(k) = 0.0_wp
            end do
            q_tot = sum(q_in)
            if (q_tot .le. 0.0_wp) cycle
            do k = 1, 4
                if (q_in(k) .le. 0.0_wp) cycle
                dH = min(rest(i,j)*q_in(k)/q_tot, H_now(in(k),jn(k)))
                H_now(in(k),jn(k)) = H_now(in(k),jn(k)) - dH
            end do
        end do
        end do

        cmb = (H_now - H_ice)/dt

        return

    end subroutine calc_G_calving_front

    subroutine calc_G_lsf_front(cmb,H_ice,a_lsf,z_bed,z_sl,rho_ice,rho_sw, &
                                    front_subgrid,H_eff_min,dHdx,dx,dt,boundaries)
        ! Thickness of the subgrid front cells follows the level set (CISM
        ! subgrid calving mask, H/H_eff = 1 - mask): eligible cells with less
        ! than A_FRONT_MIN of their area behind the front are emptied, and
        ! front cells (also eligible cells touching the ocean at a corner)
        ! hold at most a_lsf*H_ref. The reference H_ref comes from the
        ! remaining interior cells (calc_front_H_ref) and does not depend on
        ! the trimmed cell's own thickness, so repeated trimming does not
        ! compound. Cells without an interior neighbour are not trimmed.
        ! Returns the applied calving rate [m/yr, <= 0].

        implicit none

        real(wp), intent(OUT) :: cmb(:,:)               ! [m/yr] Applied calving rate
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m]    Ice thickness
        real(wp), intent(IN)  :: a_lsf(:,:)             ! [--]   Area fraction behind the level-set front
        real(wp), intent(IN)  :: z_bed(:,:)
        real(wp), intent(IN)  :: z_sl(:,:)
        real(wp), intent(IN)  :: rho_ice
        real(wp), intent(IN)  :: rho_sw
        character(len=*), intent(IN) :: front_subgrid   ! "floating" or "marine"
        real(wp), intent(IN)  :: H_eff_min              ! [m]    Minimum H_eff of eligible cells
        real(wp), intent(IN)  :: dHdx                   ! [m/m]  Thickness gradient assumed at a full front
        real(wp), intent(IN)  :: dx                     ! [m]    Grid resolution
        real(wp), intent(IN)  :: dt                     ! [yr]   Timestep
        character(len=*), intent(IN) :: boundaries

        ! Local variables
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1, BC
        logical  :: is_flt
        real(wp) :: H_max
        real(wp), allocatable :: H_now(:,:), H_ref(:,:)
        logical,  allocatable :: mask_cf(:,:), mask_elig(:,:), mask_ocn(:,:), mask_part(:,:), has_ref(:,:)

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        BC = boundary_code(boundaries)

        is_flt = trim(front_subgrid) .eq. "floating"

        allocate(H_now(nx,ny),H_ref(nx,ny))
        allocate(mask_cf(nx,ny),mask_elig(nx,ny),mask_ocn(nx,ny),mask_part(nx,ny),has_ref(nx,ny))

        H_now = H_ice

        ! Empty eligible cells (almost) entirely beyond the front
        call calc_front_cells(mask_cf,mask_elig,mask_ocn,H_now,z_bed,z_sl,rho_ice,rho_sw,front_subgrid,boundaries)
        where (mask_elig .and. a_lsf .lt. A_FRONT_MIN) H_now = 0.0_wp

        ! Partial cells: eligible cells not entirely behind the front that
        ! touch the ocean at an edge (front cells) or a corner
        call calc_front_cells(mask_cf,mask_elig,mask_ocn,H_now,z_bed,z_sl,rho_ice,rho_sw,front_subgrid,boundaries)

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
        do j = 1, ny
        do i = 1, nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            mask_part(i,j) = mask_elig(i,j) .and. a_lsf(i,j) .lt. 1.0_wp .and. &
                ( mask_cf(i,j) .or. mask_ocn(im1,jm1) .or. mask_ocn(ip1,jm1) .or. &
                                    mask_ocn(im1,jp1) .or. mask_ocn(ip1,jp1) )
        end do
        end do
        !$omp end parallel do

        ! Reference thickness from interior cells (eligible, not front, not partial)
        call calc_front_H_ref(H_ref,has_ref,H_now,z_bed,z_sl,rho_ice,rho_sw,mask_part, &
                              mask_elig .and. .not. (mask_cf .or. mask_part),front_subgrid,dHdx,dx,boundaries)

        ! Trim partial cells to a_lsf*H_ref (same bounds as H_eff: at least
        ! H_eff_min, at most flotation for "floating")
        !$omp parallel do collapse(2) private(i,j,H_max)
        do j = 1, ny
        do i = 1, nx
            if (.not. has_ref(i,j)) cycle
            H_max = max(H_ref(i,j), H_eff_min)
            if (is_flt) H_max = min(H_max, max( (z_sl(i,j)-z_bed(i,j))*rho_sw/rho_ice, 0.0_wp ))
            H_now(i,j) = min(H_now(i,j), a_lsf(i,j)*H_max)
        end do
        end do
        !$omp end parallel do

        if (dt .gt. 0.0_wp) then
            cmb = (H_now - H_ice)/dt
        else
            cmb = 0.0_wp
        end if

        return

    end subroutine calc_G_lsf_front

    subroutine calc_G_remove_fractional_ice(mb_diff,H_ice,f_ice,tau,dt,boundaries)
        ! Eliminate fractional ice covered points (icebergs) that have
        ! no fully ice-covered edge or diagonal neighbor, at the rate
        ! H/tau (all of it when dt >= tau). 

        implicit none 

        real(wp), intent(OUT) :: mb_diff(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: tau                ! [yr] Removal timescale
        real(wp), intent(IN)  :: dt 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        real(wp), allocatable :: H_new(:,:) 
        real(wp) :: f_rm 
        integer :: BC

        nx = size(H_ice,1) 
        ny = size(H_ice,2) 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        allocate(H_new(nx,ny)) 

        H_new = H_ice 

        f_rm = 1.0_wp
        if (tau .gt. dt) f_rm = dt/tau

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            if (f_ice(i,j) .gt. 0.0 .and. f_ice(i,j) .lt. 1.0) then 
                ! Fractional ice-covered point 

                if ( count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1), &
                            f_ice(im1,jm1),f_ice(ip1,jm1),f_ice(im1,jp1),f_ice(ip1,jp1)] &
                            .eq. 1.0) .eq. 0) then 
                    ! No fully ice-covered neighbors available.
                    ! Point should be removed. 

                    H_new(i,j) = (1.0_wp-f_rm)*H_ice(i,j) 

                end if

            end if

        end do 
        end do
        !$omp end parallel do

        ! Determine rate of mass balance related to changes applied here

        if (dt .ne. 0.0) then 
            mb_diff = (H_new - H_ice) / dt 
        else 
            mb_diff = 0.0
        end if

        return

    end subroutine calc_G_remove_fractional_ice

    subroutine remove_icebergs(H_ice)

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:) 

        return

    end subroutine remove_icebergs


end module mass_conservation
