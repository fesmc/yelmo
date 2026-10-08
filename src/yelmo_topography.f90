
module yelmo_topography

    use nml  
    use ncio
    
    use yelmo_defs
    use yelmo_tools 
    
    use mass_conservation
    use calving_aa
    use calving_ac
    use lsf_module
    use solver_linear, only : LGS_SUCCESS
    use topography 
    use discharge
    use velocity_general, only : set_inactive_margins

    use runge_kutta 
    use distances

    implicit none
    
    private
    public :: calc_ytopo_pc 
    public :: calc_ytopo_diagnostic 
    public :: calc_ytopo_rates
    public :: calc_transport_velocity
    public :: ytopo_par_load
    public :: ytopo_alloc
    public :: ytopo_dealloc
    
    ! Integers
    public :: mask_bed_ocean  
    public :: mask_bed_land  
    public :: mask_bed_frozen
    public :: mask_bed_stream
    public :: mask_bed_grline
    public :: mask_bed_float 
    public :: mask_bed_island
    
contains
    
    subroutine calc_ytopo_pc(tpo,dyn,mat,thrm,bnd,dta,time,topo_fixed,pc_step,use_H_pred,filter_vel)

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd
        type(ydata_class),  intent(IN)    :: dta 
        real(wp),           intent(IN)    :: time
        logical,            intent(IN)    :: topo_fixed  
        character(len=*),   intent(IN)    :: pc_step 
        logical, optional,  intent(IN)    :: use_H_pred
        logical, optional,  intent(IN)    :: filter_vel     ! Advect with mean of current and previous velocity solutions

        ! Local variables 
        integer  :: i, j, nx, ny
        real(wp) :: dt  
        real(wp), allocatable :: dHidt_now(:,:) 
        real(wp), allocatable :: H_prev(:,:)
        real(wp), allocatable :: ux_adv(:,:)            ! Transport velocity (calc_transport_velocity)
        real(wp), allocatable :: uy_adv(:,:)
        logical,  allocatable :: mask_cf(:,:), mask_elig(:,:), mask_ocn(:,:)
        integer  :: lin_iter, lin_status                ! Linear solver iterations and status of the advection
        logical  :: filt

        logical, parameter :: use_rk4 = .FALSE. 

        real(8)  :: cpu_time0, cpu_time1

        ! Store initial cpu time for the speed metric
        call yelmo_cpu_time(cpu_time0)

        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)

        allocate(dHidt_now(nx,ny))
        allocate(H_prev(nx,ny))
        allocate(ux_adv(nx,ny))
        allocate(uy_adv(nx,ny))
        allocate(mask_cf(nx,ny),mask_elig(nx,ny),mask_ocn(nx,ny))

        filt = .FALSE.
        if (present(filter_vel)) filt = filter_vel

        ! Initialize time if necessary 
        if (tpo%par%time .gt. dble(time)) then 
            tpo%par%time = dble(time) 
        end if 
        
        ! Initialize rk4 object if necessary 
        if (.not. allocated(tpo%rk4%tau)) then 
            call rk4_2D_init(tpo%rk4,nx,ny,dt_min=1e-3_wp)
        end if 

        ! Get time step
        dt = dble(time) - tpo%par%time 

        ! Step 0: Get some diagnostic quantities for mass balance calculation --------

        ! Get ice thickness entering routine
        H_prev = tpo%now%H_ice 

        ! Step 1: Go through predictor-corrector-advance steps

        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then 

            ! Ice thickness evolution from dynamics alone

            select case(trim(pc_step))

                case("predictor") 
                    ! Determine predicted ice thickness 

                    ! Store dynamic rate of change from previous timestep,
                    ! along with ice thickness and surface elevation
                    ! (the latter only for calculating rate of change later)
                    tpo%now%H_ice_n     = tpo%now%H_ice
                    tpo%now%H_ice_dyn_n = tpo%now%H_ice_dyn
                    tpo%now%z_srf_n     = tpo%now%z_srf
                    tpo%now%lsf_n       = tpo%now%lsf

                    ! Get ice-fraction mask for current ice thickness  
                    call update_ice_fraction(tpo,bnd)

                    call calc_transport_velocity(ux_adv,uy_adv,tpo,dyn,bnd,filt)

if (use_rk4) then
                    call rk4_2D_step(tpo%rk4,tpo%now%H_ice,tpo%now%f_ice,dHidt_now,ux_adv,uy_adv, &
                                                bnd%mask_ice,tpo%par%dx,dt,tpo%par%solver,tpo%par%boundaries)

else
                    call calc_G_advec_simple(dHidt_now,tpo%now%H_ice,tpo%now%f_ice,ux_adv,uy_adv, &
                                                 bnd%mask_ice,tpo%par%solver,tpo%par%boundaries,tpo%par%dx,dt, &
                                                 lin_iter=lin_iter,lin_status=lin_status)
                    tpo%par%adv_lin_iter = lin_iter
                    tpo%par%adv_lin_fail = merge(1,0,lin_status .ne. LGS_SUCCESS)
                 
end if

                    ! Store raw advective rate at current state (f_n) for the
                    ! corrector and for shifting at advance
                    tpo%now%dHidt_dyn_raw = dHidt_now

                    ! Calculate rate of change using weighted advective rates of change
                    ! depending on timestepping method chosen
                    tpo%now%dHidt_dyn = tpo%par%dt_beta(1)*dHidt_now + tpo%par%dt_beta(2)*tpo%now%dHidt_dyn_raw_n

                    ! Apply rate and update ice thickness (predicted).
                    ! The predictor-corrector mixing of dHdt with previous timesteps can
                    ! give small negative ice thicknesses at the margin. apply_tendency
                    ! clips them to zero (a mass source). The clip is booked in mb_clip,
                    ! so that dHidt_dyn stays pure transport and the budget shows it.
                    ! dHidt_vert holds the applied rate (dHidt_dyn + mb_clip), see below.
                    ! dHidt_dyn is then the applied transport rate (with the rounding of
                    ! the update), so the budget closes; the history is dHidt_dyn_raw.
                    tpo%now%H_ice      = tpo%now%H_ice_n
                    tpo%now%dHidt_vert = tpo%now%dHidt_dyn
                    call apply_tendency(tpo%now%H_ice,tpo%now%dHidt_vert,dt,"dyn_pred",adjust_mb=.TRUE., &
                                        mb_clip=tpo%now%mb_clip)
                    tpo%now%dHidt_dyn  = tpo%now%dHidt_vert - tpo%now%mb_clip

                case("corrector") 

                    ! Set current thickness to predicted thickness
                    tpo%now%H_ice = tpo%now%pred%H_ice 
                    tpo%now%lsf   = tpo%now%pred%lsf

                    ! Get ice-fraction mask for predicted ice thickness  
                    call update_ice_fraction(tpo,bnd)

                    call calc_transport_velocity(ux_adv,uy_adv,tpo,dyn,bnd,filt)

if (use_rk4) then
                    call rk4_2D_step(tpo%rk4,tpo%now%H_ice,tpo%now%f_ice,dHidt_now,ux_adv,uy_adv, &
                                                bnd%mask_ice,tpo%par%dx,dt,tpo%par%solver,tpo%par%boundaries)
else
                    call calc_G_advec_simple(dHidt_now,tpo%now%H_ice,tpo%now%f_ice,ux_adv,uy_adv, &
                                                bnd%mask_ice,tpo%par%solver,tpo%par%boundaries,tpo%par%dx,dt, &
                                                lin_iter=lin_iter,lin_status=lin_status)
                    tpo%par%adv_lin_iter = tpo%par%adv_lin_iter + lin_iter
                    tpo%par%adv_lin_fail = tpo%par%adv_lin_fail + merge(1,0,lin_status .ne. LGS_SUCCESS)
                 
end if

                    ! Calculate rate of change using weighted advective rates of change 
                    ! depending on timestepping method chosen 
                    tpo%now%dHidt_dyn = tpo%par%dt_beta(3)*dHidt_now + tpo%par%dt_beta(4)*tpo%now%dHidt_dyn_raw
                    
                    ! Apply rate and update ice thickness (corrected), clip booked in mb_clip,
                    ! dHidt_dyn the applied transport rate (as in the predictor)
                    tpo%now%H_ice      = tpo%now%H_ice_n
                    tpo%now%lsf        = tpo%now%lsf_n
                    tpo%now%dHidt_vert = tpo%now%dHidt_dyn
                    call apply_tendency(tpo%now%H_ice,tpo%now%dHidt_vert,dt,"dyn_corr",adjust_mb=.TRUE., &
                                        mb_clip=tpo%now%mb_clip)
                    tpo%now%dHidt_dyn  = tpo%now%dHidt_vert - tpo%now%mb_clip

            end select

            ! The mass balance below acts on the ice present after transport
            ! (area fraction of cells that received or lost ice)
            call update_ice_fraction(tpo,bnd)

            ! Note: at this point, mass has only been advected (moved around), plus the
            ! clip of negative thicknesses (mb_clip, above).

            select case(trim(pc_step))

                case("predictor","corrector")
                    ! For either predictor or corrector step, also calculate all mass balance changes

                    ! === smb =====
                    call calc_G_mbal(tpo%now%smb,tpo%now%H_ice,tpo%now%f_grnd,bnd%smb,dt,tpo%now%f_ice)

                    ! Apply rate and update ice thickness
                    call apply_tendency(tpo%now%H_ice,tpo%now%smb,dt,"smb",adjust_mb=.TRUE.)
                    
                    ! === bmb =====

                    ! Calculate grounded fraction on aa-nodes
                    ! (only to be used with basal mass balance, later all
                    !  f_grnd arrays will be calculated according to use choices)
                    call determine_grounded_fractions(tpo%now%f_grnd_bmb,H_grnd=tpo%now%H_grnd, &
                                                                        boundaries=tpo%par%boundaries)
                    
                    ! Combine basal mass balance into one field accounting for 
                    ! grounded/floating fraction of grid cells 
                    call calc_bmb_total(tpo%now%bmb_ref,thrm%now%bmb_grnd,bnd%bmb_shlf,tpo%now%H_ice, &
                                        tpo%now%H_grnd,tpo%now%f_grnd_bmb,tpo%par%gz_Hg0,tpo%par%gz_Hg1, &
                                        tpo%par%gz_nx,tpo%par%bmb_gl_method,tpo%par%boundaries)

                    if (tpo%par%use_bmb) then
                        call calc_G_mbal(tpo%now%bmb,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%bmb_ref,dt,tpo%now%f_ice)
                    else
                        ! Mainly for when running EISMINT1
                        tpo%now%bmb = 0.0
                    end if

                    ! Apply rate and update ice thickness
                    call apply_tendency(tpo%now%H_ice,tpo%now%bmb,dt,"bmb",adjust_mb=.TRUE.)

                    ! === fmb =====

                    ! Calculate frontal mass balance
                    call calc_fmb_total(tpo%now%fmb_ref,bnd%fmb_shlf,bnd%bmb_shlf,tpo%now%H_ice, &
                                    tpo%now%H_grnd,tpo%now%f_ice,tpo%par%fmb_method,tpo%par%fmb_scale, &
                                    tpo%par%fmb_lambda, bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%dx,tpo%par%boundaries, &
                                    bnd%Qd,bnd%tf_shlf)

                    if (tpo%par%use_bmb) then
                        call calc_G_mbal(tpo%now%fmb,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%fmb_ref,dt)
                    else
                        ! Mainly for when running EISMINT1
                        tpo%now%fmb = 0.0
                    end if

                    ! Apply rate and update ice thickness
                    call apply_tendency(tpo%now%H_ice,tpo%now%fmb,dt,"fmb",adjust_mb=.TRUE.)
                    
                    ! === dmb =====

                    call calc_mb_discharge(tpo%now%dmb_ref,tpo%now%H_ice,tpo%now%z_srf,bnd%z_bed_sd,tpo%now%dist_grline, &
                                tpo%now%dist_margin,tpo%par%dmb_method,tpo%par%dx,tpo%par%dmb_alpha_max, &
                                tpo%par%dmb_tau,tpo%par%dmb_sigma_ref,tpo%par%dmb_m_d,tpo%par%dmb_m_r)
                    
                    call calc_G_mbal(tpo%now%dmb,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%dmb_ref,dt)

                    ! Apply rate and update ice thickness
                    call apply_tendency(tpo%now%H_ice,tpo%now%dmb,dt,"dmb",adjust_mb=.TRUE.)
                    
                    ! === mb_net =====
                    tpo%now%mb_net = tpo%now%smb + tpo%now%bmb + tpo%now%fmb + tpo%now%dmb 

                    ! Vertical thickness change of the ice column (advection, smb, bmb;
                    ! relaxation added below). Lateral changes (fmb, dmb, calving, front
                    ! advance, removals) are not vertical motion of the column surface or base.
                    ! (dHidt_vert holds the applied dynamic rate dHidt_dyn + mb_clip)
                    tpo%now%dHidt_vert = tpo%now%dHidt_vert + tpo%now%smb + tpo%now%bmb

                    ! === calving ===
                    ! Calculate and apply calving
                    if (tpo%par%use_lsf) then
                        ! Level-set function as calving
                        call calc_ytopo_calving_lsf(tpo,dyn,mat,thrm,bnd,dt,H_prev,time)
                    else
                        ! Mass balance calving
                        call calc_ytopo_calving(tpo,dyn,mat,thrm,bnd,dt)
                    end if

                    ! Applied calving (all removal at the front) split by the state
                    ! of the cell: floating (f_grnd = 0) or grounded
                    where (tpo%now%f_grnd .eq. 0.0_wp)
                        tpo%now%cmb_flt  = tpo%now%cmb
                        tpo%now%cmb_grnd = 0.0_wp
                    elsewhere
                        tpo%now%cmb_flt  = 0.0_wp
                        tpo%now%cmb_grnd = tpo%now%cmb
                    end where

                    ! Get ice-fraction mask for ice thickness  
                    call update_ice_fraction(tpo,bnd)

                    ! If desired, finally relax solution to reference state
                    if (tpo%par%topo_rel .ne. 0) then 

                        if (tpo%par%topo_rel .eq. -1) then
                            ! Use externally defined tau_relax field
                            tpo%now%tau_relax = bnd%tau_relax
                        else
                            ! Define the relaxation timescale field, if needed
                            call set_tau_relax(tpo%now%tau_relax,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%mask_grz,bnd%H_ice_ref, &
                                               tpo%par%topo_rel,tpo%par%topo_rel_tau,tpo%par%boundaries)
                        end if

                        ! Relaxation is only meaningful where ice is actively solved
                        where (bnd%mask_ice .ne. MASK_ICE_DYNAMIC) tpo%now%tau_relax = -1.0_wp

                        select case(trim(tpo%par%topo_rel_field))

                            case("H_ref")
                                ! Relax towards reference ice thickness field H_ref

                                call calc_G_relaxation(tpo%now%mb_relax,tpo%now%H_ice,bnd%H_ice_ref,tpo%now%tau_relax,dt)

                            case("H_ice_n")
                                ! Relax towards previous iteration ice thickness 
                                ! (ie slow down changes)
                                ! ajr: needs testing, not sure if this works well or helps anything.

                                call calc_G_relaxation(tpo%now%mb_relax,tpo%now%H_ice,tpo%now%H_ice_n,tpo%now%tau_relax,dt)

                            case DEFAULT 

                                write(*,*) "calc_ytopo:: Error: topo_rel_field not recognized."
                                write(*,*) "topo_rel_field = ", trim(tpo%par%topo_rel_field)
                                error stop 1

                        end select

                        ! Apply rate and update ice thickness
                        call apply_tendency(tpo%now%H_ice,tpo%now%mb_relax,dt,"relax",adjust_mb=.TRUE.)

                        ! Add relaxation tendency to mb_net for proper accounting of mass change
                        tpo%now%mb_net = tpo%now%mb_net + tpo%now%mb_relax
                        tpo%now%dHidt_vert = tpo%now%dHidt_vert + tpo%now%mb_relax

                        ! Get ice-fraction mask for ice thickness  
                        call update_ice_fraction(tpo,bnd)

                    end if

                    ! Finally, apply all additional (generally artificial) ice thickness adjustments 
                    ! and store changes in residual mass balance field. 
                    call calc_front_cells(mask_cf,mask_elig,mask_ocn,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                    bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%front_subgrid,tpo%par%boundaries)
                    call calc_G_boundaries(tpo%now%mb_resid,tpo%now%H_ice,tpo%now%H_eff,mask_cf,tpo%now%f_grnd, &
                                            dyn%now%uxy_b,bnd%mask_ice,tpo%par%boundaries,bnd%H_ice_ref, &
                                            tpo%par%H_min_flt,tpo%par%H_min_grnd,tpo%par%H_min_tau,dt)

                    ! Apply rate and update ice thickness
                    call apply_tendency(tpo%now%H_ice,tpo%now%mb_resid,dt,"resid",adjust_mb=.TRUE.)

                    ! Add residual tendency to mb_net for proper accounting of mass change
                    tpo%now%mb_net = tpo%now%mb_net + tpo%now%mb_resid

                    ! Get ice-fraction mask for ice thickness  
                    call update_ice_fraction(tpo,bnd)
                    

            end select 

            select case(trim(pc_step))

                case("predictor") 

                    ! Save current predictor fields, 
                    ! proceed with predictor fields for calculating dynamics.
                    tpo%now%pred%H_ice      = tpo%now%H_ice 
                    tpo%now%pred%dHidt_dyn  = tpo%now%dHidt_dyn
                    tpo%now%pred%dHidt_vert = tpo%now%dHidt_vert
                    tpo%now%pred%mb_net     = tpo%now%mb_net 
                    tpo%now%pred%mb_relax   = tpo%now%mb_relax 
                    tpo%now%pred%mb_resid   = tpo%now%mb_resid 
                    tpo%now%pred%mb_clip    = tpo%now%mb_clip
                    tpo%now%pred%smb        = tpo%now%smb
                    tpo%now%pred%bmb        = tpo%now%bmb
                    tpo%now%pred%fmb        = tpo%now%fmb
                    tpo%now%pred%dmb        = tpo%now%dmb
                    tpo%now%pred%cmb        = tpo%now%cmb 
                    tpo%now%pred%cmb_flt    = tpo%now%cmb_flt 
                    tpo%now%pred%cmb_grnd   = tpo%now%cmb_grnd
                    tpo%now%pred%lsf        = tpo%now%lsf 
                    tpo%now%pred%cmb_flt_x      = tpo%now%cmb_flt_x
                    tpo%now%pred%cmb_flt_y      = tpo%now%cmb_flt_y
                    tpo%now%pred%cmb_grnd_x     = tpo%now%cmb_grnd_x
                    tpo%now%pred%cmb_grnd_y     = tpo%now%cmb_grnd_y
                    tpo%now%pred%cr_acx         = tpo%now%cr_acx
                    tpo%now%pred%cr_acy         = tpo%now%cr_acy
                    tpo%now%pred%calv_rate_flt  = tpo%now%calv_rate_flt
                    tpo%now%pred%calv_rate_grnd = tpo%now%calv_rate_grnd
                    
                case("corrector")
                    ! Determine corrected ice thickness 

                    ! Save current corrector fields
                    tpo%now%corr%H_ice      = tpo%now%H_ice 
                    tpo%now%corr%dHidt_dyn  = tpo%now%dHidt_dyn
                    tpo%now%corr%dHidt_vert = tpo%now%dHidt_vert
                    tpo%now%corr%mb_net     = tpo%now%mb_net 
                    tpo%now%corr%mb_relax   = tpo%now%mb_relax 
                    tpo%now%corr%mb_resid   = tpo%now%mb_resid 
                    tpo%now%corr%mb_clip    = tpo%now%mb_clip
                    tpo%now%corr%smb        = tpo%now%smb
                    tpo%now%corr%bmb        = tpo%now%bmb
                    tpo%now%corr%fmb        = tpo%now%fmb
                    tpo%now%corr%dmb        = tpo%now%dmb
                    tpo%now%corr%cmb        = tpo%now%cmb 
                    tpo%now%corr%cmb_flt    = tpo%now%cmb_flt 
                    tpo%now%corr%cmb_grnd   = tpo%now%cmb_grnd
                    tpo%now%corr%lsf        = tpo%now%lsf
                    tpo%now%corr%cmb_flt_x      = tpo%now%cmb_flt_x
                    tpo%now%corr%cmb_flt_y      = tpo%now%cmb_flt_y
                    tpo%now%corr%cmb_grnd_x     = tpo%now%cmb_grnd_x
                    tpo%now%corr%cmb_grnd_y     = tpo%now%cmb_grnd_y
                    tpo%now%corr%cr_acx         = tpo%now%cr_acx
                    tpo%now%corr%cr_acy         = tpo%now%cr_acy
                    tpo%now%corr%calv_rate_flt  = tpo%now%calv_rate_flt
                    tpo%now%corr%calv_rate_grnd = tpo%now%calv_rate_grnd
                    
                    ! Restore main ice thickness field to original 
                    ! value at the beginning of the timestep for 
                    ! calculation of remaining quantities (thermo, material)
                    tpo%now%H_ice = tpo%now%H_ice_n 
                    tpo%now%lsf   = tpo%now%lsf_n 

                case("advance")
                    ! Now let's actually advance the ice thickness field

                    if (.not. present(use_H_pred)) then 
                        write(*,*) "calc_ytopo_pc:: Error: &
                        & For step='advance', the argument use_H_pred&
                        & must be provided."
                        error stop 1
                    end if 

                    ! Determine which ice thickness to use going forward
                    if (use_H_pred) then 

                        ! Load predictor fields in current state variables
                        tpo%now%H_ice       = tpo%now%pred%H_ice 
                        tpo%now%dHidt_dyn   = tpo%now%pred%dHidt_dyn
                        tpo%now%dHidt_vert  = tpo%now%pred%dHidt_vert
                        tpo%now%mb_net      = tpo%now%pred%mb_net 
                        tpo%now%mb_relax    = tpo%now%pred%mb_relax 
                        tpo%now%mb_resid    = tpo%now%pred%mb_resid 
                        tpo%now%mb_clip     = tpo%now%pred%mb_clip
                        tpo%now%smb         = tpo%now%pred%smb 
                        tpo%now%bmb         = tpo%now%pred%bmb 
                        tpo%now%fmb         = tpo%now%pred%fmb 
                        tpo%now%dmb         = tpo%now%pred%dmb 
                        tpo%now%cmb         = tpo%now%pred%cmb 
                        tpo%now%cmb_flt     = tpo%now%pred%cmb_flt
                        tpo%now%cmb_grnd    = tpo%now%pred%cmb_grnd
                        tpo%now%lsf         = tpo%now%pred%lsf 
                        tpo%now%cmb_flt_x      = tpo%now%pred%cmb_flt_x
                        tpo%now%cmb_flt_y      = tpo%now%pred%cmb_flt_y
                        tpo%now%cmb_grnd_x     = tpo%now%pred%cmb_grnd_x
                        tpo%now%cmb_grnd_y     = tpo%now%pred%cmb_grnd_y
                        tpo%now%cr_acx         = tpo%now%pred%cr_acx
                        tpo%now%cr_acy         = tpo%now%pred%cr_acy
                        tpo%now%calv_rate_flt  = tpo%now%pred%calv_rate_flt
                        tpo%now%calv_rate_grnd = tpo%now%pred%calv_rate_grnd
                        
                    else
                        ! Load corrector fields in current state variables
                        tpo%now%H_ice       = tpo%now%corr%H_ice 
                        tpo%now%dHidt_dyn   = tpo%now%corr%dHidt_dyn
                        tpo%now%dHidt_vert  = tpo%now%corr%dHidt_vert
                        tpo%now%mb_net      = tpo%now%corr%mb_net 
                        tpo%now%mb_relax    = tpo%now%corr%mb_relax 
                        tpo%now%mb_resid    = tpo%now%corr%mb_resid 
                        tpo%now%mb_clip     = tpo%now%corr%mb_clip
                        tpo%now%smb         = tpo%now%corr%smb 
                        tpo%now%bmb         = tpo%now%corr%bmb 
                        tpo%now%fmb         = tpo%now%corr%fmb 
                        tpo%now%dmb         = tpo%now%corr%dmb 
                        tpo%now%cmb         = tpo%now%corr%cmb 
                        tpo%now%cmb_flt     = tpo%now%corr%cmb_flt
                        tpo%now%cmb_grnd    = tpo%now%corr%cmb_grnd
                        tpo%now%lsf         = tpo%now%corr%lsf
                        tpo%now%cmb_flt_x      = tpo%now%corr%cmb_flt_x
                        tpo%now%cmb_flt_y      = tpo%now%corr%cmb_flt_y
                        tpo%now%cmb_grnd_x     = tpo%now%corr%cmb_grnd_x
                        tpo%now%cmb_grnd_y     = tpo%now%corr%cmb_grnd_y
                        tpo%now%cr_acx         = tpo%now%corr%cr_acx
                        tpo%now%cr_acy         = tpo%now%corr%cr_acy
                        tpo%now%calv_rate_flt  = tpo%now%corr%calv_rate_flt
                        tpo%now%calv_rate_grnd = tpo%now%corr%calv_rate_grnd

                    end if

                    ! Shift raw advective rate f_n -> f_{n-1} for the next step.
                    ! Done once here, regardless of the use_H_pred branch above.
                    tpo%now%dHidt_dyn_raw_n = tpo%now%dHidt_dyn_raw

            end select

            ! Determine rates of change
            ! Note: dzsdt is deferred until after calc_ytopo_diagnostic below,
            ! because z_srf is only refreshed from the updated H_ice/f_ice there
            ! (calc_z_srf_max).
            tpo%now%dHidt  = (tpo%now%H_ice - tpo%now%H_ice_n) / dt
            tpo%now%dlsfdt = (tpo%now%lsf   - tpo%now%lsf_n) / dt

            ! Determine mass balance error as the residual of dHidt with respect to
            ! all applied tendencies (dynamics, clip, mb_net incl. relax and resid, calving).
            ! Since every tendency passes through apply_tendency (adjust_mb=.TRUE.),
            ! this should vanish to round-off.
            tpo%now%mb_err = tpo%now%dHidt - (tpo%now%dHidt_dyn + tpo%now%mb_clip + tpo%now%mb_net + tpo%now%cmb)

        end if

        ! Update fields and masks (refreshes z_srf from the current H_ice/f_ice)
        call calc_ytopo_diagnostic(tpo,dyn,mat,thrm,bnd)

        ! Surface-elevation rate: now that calc_ytopo_diagnostic has refreshed
        ! z_srf, difference it against the start-of-step snapshot z_srf_n.
        ! Same guard as the other rates; z_srf_n is set at step start (predictor)
        ! and is not modified by calc_ytopo_diagnostic.
        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then
            tpo%now%dzsdt = (tpo%now%z_srf - tpo%now%z_srf_n) / dt
        end if

        ! Kinematic rates of the ice column's surface and base, the boundary
        ! conditions of the vertical velocity: vertical column change plus the
        ! bedrock and sea-level rates (no vertical change when the ice is not advanced)
        if (topo_fixed .or. dt .le. 0.0) tpo%now%dHidt_vert = 0.0_wp
        call calc_column_kinematic_rates(tpo%now%dzsdt_kin,tpo%now%dzbdt_kin,tpo%now%mask_kin,tpo%now%dHidt_vert, &
                    (.not. topo_fixed .and. dt .gt. 0.0),tpo%now%f_grnd,tpo%now%f_ice,tpo%now%H_ice,tpo%now%H_ice_dyn,tpo%now%H_ice_n, &
                    tpo%now%H_ice_dyn_n,bnd%dz_bed_dt,bnd%dz_sl_dt,bnd%c%rho_ice,bnd%c%rho_sw)

        ! When the ice is not advanced this step -- initialization (pc_step="none")
        ! or any topo_fixed step -- the mass-balance block above is skipped, so the
        ! applied-mb diagnostics (smb/bmb/fmb) would otherwise stay zero. Diagnose
        ! them from the current boundary forcing so the held/initial state carries a
        ! meaningful mass balance in the output. dt=0 => calc_G_mbal returns the
        ! forcing masked to the current ice, nothing is applied and H_ice is left
        ! untouched. mb_net stays zero (nothing applied) and bmb_grnd may still be
        ! zero at initialization (thermodynamics not yet computed); both expected.
        if ( topo_fixed .or. dt .le. 0.0 ) then
            call calc_ytopo_mb_diagnostic(tpo,thrm,bnd)
        end if


        ! Computational performance (model speed in kyr/hr) of the topography
        ! calls of this step: the predictor starts the count, the advance ends it
        call yelmo_cpu_time(cpu_time1)
        if (trim(pc_step) .eq. "predictor") tpo%par%cpu_step = 0.0d0
        tpo%par%cpu_step = tpo%par%cpu_step + (cpu_time1 - cpu_time0)

        if (trim(pc_step) .eq. "advance") then 
            ! Advance timestep here whether topo_fixed was true or not...

            call yelmo_calc_speed(tpo%par%speed,real(tpo%par%time,wp),time,0.0d0,tpo%par%cpu_step)
            
            ! Update ytopo time to current time 
            tpo%par%time = dble(time)
            
        end if

        return

    end subroutine calc_ytopo_pc

    subroutine calc_ytopo_mb_diagnostic(tpo,thrm,bnd)
        ! Diagnose the applied mass-balance fields (smb, combined bmb, fmb) from the
        ! current boundary forcing and geometry WITHOUT advancing the ice. Used by
        ! calc_ytopo_pc when the ice is not evolved this step (initialization via
        ! pc_step="none", or any topo_fixed step): the main mass-balance block there
        ! is gated on (.not. topo_fixed .and. dt>0), so without this the output
        ! smb/bmb/fmb would remain zero at the initial time. The calls mirror the
        ! diagnostic half of that block (calc_G_mbal etc.) but omit apply_tendency,
        ! and use dt=0 so calc_G_mbal returns the forcing masked to the current ice
        ! (no melt-limiting, H_ice untouched). mb_net is deliberately left unchanged
        ! (nothing is applied); bmb_grnd may still be zero at init (thermodynamics
        ! not yet computed), which only zeroes the grounded part of the combined bmb.

        implicit none

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd

        ! Diagnostic-only: no timestep is taken, so nothing is applied to H_ice.
        real(wp), parameter :: dt = 0.0_wp

        ! === smb ===
        call calc_G_mbal(tpo%now%smb,tpo%now%H_ice,tpo%now%f_grnd,bnd%smb,dt,tpo%now%f_ice)

        ! === bmb (combined grounded + shelf) ===
        call determine_grounded_fractions(tpo%now%f_grnd_bmb,H_grnd=tpo%now%H_grnd, &
                                                            boundaries=tpo%par%boundaries)
        call calc_bmb_total(tpo%now%bmb_ref,thrm%now%bmb_grnd,bnd%bmb_shlf,tpo%now%H_ice, &
                            tpo%now%H_grnd,tpo%now%f_grnd_bmb,tpo%par%gz_Hg0,tpo%par%gz_Hg1, &
                            tpo%par%gz_nx,tpo%par%bmb_gl_method,tpo%par%boundaries)
        if (tpo%par%use_bmb) then
            call calc_G_mbal(tpo%now%bmb,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%bmb_ref,dt,tpo%now%f_ice)
        else
            tpo%now%bmb = 0.0
        end if

        ! === fmb ===
        call calc_fmb_total(tpo%now%fmb_ref,bnd%fmb_shlf,bnd%bmb_shlf,tpo%now%H_ice, &
                        tpo%now%H_grnd,tpo%now%f_ice,tpo%par%fmb_method,tpo%par%fmb_scale, &
                        tpo%par%fmb_lambda,bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%dx,tpo%par%boundaries, &
                        bnd%Qd,bnd%tf_shlf)
        if (tpo%par%use_bmb) then
            call calc_G_mbal(tpo%now%fmb,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%fmb_ref,dt)
        else
            tpo%now%fmb = 0.0
        end if

        return

    end subroutine calc_ytopo_mb_diagnostic

    subroutine calc_ytopo_calving(tpo,dyn,mat,thrm,bnd,dt)

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 
        real(wp),           intent(IN)    :: dt

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp), allocatable :: mbal_now(:,:) 
        real(wp), allocatable :: cmb_sd(:,:)
        real(wp), allocatable :: tau_eig_1(:,:), tau_eig_2(:,:)
        logical,  allocatable :: mask_cf(:,:), mask_elig(:,:), mask_ocn(:,:)
        logical,  allocatable :: mask_kill(:,:)

        nx = size(tpo%now%H_ice,1) 
        ny = size(tpo%now%H_ice,2) 

        allocate(mbal_now(nx,ny)) 
        allocate(mask_cf(nx,ny),mask_elig(nx,ny),mask_ocn(nx,ny))
        allocate(mask_kill(nx,ny))
        allocate(cmb_sd(nx,ny)) 
        allocate(tau_eig_1(nx,ny),tau_eig_2(nx,ny))


        ! Make sure current ice mask is correct
        call update_ice_fraction(tpo,bnd)

        ! === CALVING ===

        ! == Diagnose strains and stresses relevant to calving ==

        ! eps_eff = effective strain = eigencalving e+*e- following Levermann et al. (2012)
        call calc_eps_eff(tpo%now%eps_eff,dyn%now%strn2D%eps_eig_1,dyn%now%strn2D%eps_eig_2,tpo%now%f_ice)

        ! Principal stresses on the current ice geometry (cells that received ice
        ! since the last velocity solution take their neighbours' stresses)
        call fill_stress_new_ice(tau_eig_1,tau_eig_2,mat%now%strs2D%tau_eig_1,mat%now%strs2D%tau_eig_2, &
                                    tpo%now%f_ice,dyn%now%f_ice_solv,tpo%par%boundaries)

        ! tau_eff = effective stress ~ von Mises stress following Lipscomb et al. (2019)
        call calc_tau_eff(tpo%now%tau_eff,tau_eig_1,tau_eig_2,tpo%now%f_ice,tpo%par%w2)

        ! == Determine thickness threshold for calving spatially ==

        call define_calving_thickness_threshold(tpo%now%H_calv,tpo%now%z_bed_filt,tpo%par%Hc_ref_flt, &
                                            tpo%par%Hc_deep,tpo%par%zb_deep_0,tpo%par%zb_deep_1)

        ! Define factor for calving stress spatially
        call define_calving_stress_factor(tpo%now%kt,tpo%now%z_bed_filt,tpo%par%kt_ref, &
                                            tpo%par%kt_deep,tpo%par%zb_deep_0,tpo%par%zb_deep_1)

        ! == Calculate potential floating calving rate ==

        select case(trim(tpo%par%calv_flt_method))

            case("zero","none")

                tpo%now%cmb_flt = 0.0 

            case("threshold") 
                ! Use threshold method

                call calc_calving_rate_threshold(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                 tpo%now%H_calv,tpo%par%calv_tau,tpo%par%boundaries)
                
            case("vm-l19")
                ! Use von Mises calving as defined by Lipscomb et al. (2019)

                ! Next, diagnose calving
                call calc_calving_rate_vonmises_l19(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                                        tpo%now%tau_eff,tpo%par%dx,tpo%now%kt,tpo%par%boundaries)

                ! Scale calving with 'thin' calving rate to ensure 
                ! small ice thicknesses are removed (not needed with the subgrid front)
                if (trim(tpo%par%front_subgrid) .eq. "none") then
                    call apply_calving_rate_thin(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,tpo%par%calv_thin,tpo%par%Hc_ref_thin,tpo%par%boundaries)
                end if

            case("eigen")
                ! Use Eigen calving as defined by Levermann et al. (2012)

                ! Next, diagnose calving
                call calc_calving_rate_eigen(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                                        tpo%now%eps_eff,tpo%par%dx,tpo%par%k2,tpo%par%boundaries)

                ! Scale calving with 'thin' calving rate to ensure 
                ! small ice thicknesses are removed (not needed with the subgrid front)
                if (trim(tpo%par%front_subgrid) .eq. "none") then
                    call apply_calving_rate_thin(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,tpo%par%calv_thin,tpo%par%Hc_ref_thin,tpo%par%boundaries)
                end if

            case("kill","kill-pos")
                ! Floating ice is removed at the end of the calving step (below)

                tpo%now%cmb_flt = 0.0 

            case DEFAULT 

                write(*,*) "calc_ytopo:: Error: floating calving method not recognized."
                write(*,*) "calv_flt_method = ", trim(tpo%par%calv_flt_method)
                error stop 1

        end select
        
        ! Additionally ensure higher calving rate for floating tongues of
        ! one grid-point width, or lower calving rate for embayed points.

        select case(trim(tpo%par%calv_flt_method))

            case("zero","none","kill","kill-pos")
                
                ! Do nothing for these methods

            case DEFAULT

                if (trim(tpo%par%front_subgrid) .eq. "none") then
                    call calc_calving_rate_tongues(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice, &
                                                tpo%now%f_grnd,tpo%par%calv_tau,tpo%par%boundaries)
                else
                    ! Subgrid front: demand scaled by the front length, and
                    ! demand beyond a front cell's ice taken from upstream
                    call calc_front_cells(mask_cf,mask_elig,mask_ocn,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                    bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%front_subgrid,tpo%par%boundaries)
                    mbal_now = tpo%now%cmb_flt
                    call calc_G_calving_front(tpo%now%cmb_flt,tpo%now%H_ice,mbal_now,mask_cf,mask_elig,mask_ocn, &
                                    dyn%now%ux_bar,dyn%now%uy_bar,dt,tpo%par%boundaries)
                end if

        end select

        
        ! Apply rate and update ice thickness
        call apply_tendency(tpo%now%H_ice,tpo%now%cmb_flt,dt,"calv_flt",adjust_mb=.TRUE.)


        ! Diagnose potential grounded-ice calving rate [m/yr]

        select case(trim(tpo%par%calv_grnd_method))

            case("zero","none")

                tpo%now%cmb_grnd = 0.0 

            case("stress-b12") 
                ! Use simple threshold method

                call calc_calving_ground_rate_stress_b12(tpo%now%cmb_grnd,tpo%now%H_ice,tpo%now%f_ice, &
                                                tpo%now%f_grnd,bnd%z_bed,bnd%z_sl-bnd%z_bed,tpo%par%calv_tau, &
                                                bnd%c%rho_ice,bnd%c%rho_sw,bnd%c%g,tpo%par%boundaries)

            case DEFAULT 

                write(*,*) "calc_ytopo:: Error: grounded calving method not recognized."
                write(*,*) "calv_grnd_method = ", trim(tpo%par%calv_grnd_method)
                error stop 1

        end select
        
        ! Additionally include parameterized grounded calving 
        ! to account for grid resolution 
        select case(trim(tpo%par%calv_grnd_method))

            case("zero","none")
                
                ! Do nothing for these methods

            case DEFAULT 

                call calc_calving_ground_rate_stdev(cmb_sd,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                bnd%z_bed_sd,tpo%par%sd_min,tpo%par%sd_max,tpo%par%calv_grnd_max,tpo%par%calv_tau,tpo%par%boundaries)
                tpo%now%cmb_grnd = tpo%now%cmb_grnd + cmb_sd 

        end select
        
        ! Apply rate and update ice thickness
        call apply_tendency(tpo%now%H_ice,tpo%now%cmb_grnd,dt,"calv_grnd",adjust_mb=.TRUE.)

        ! Finally, get the total combined calving mass balance
        tpo%now%cmb = tpo%now%cmb_flt + tpo%now%cmb_grnd 

        ! Update ice fraction mask 
        call update_ice_fraction(tpo,bnd)

        if (trim(tpo%par%front_subgrid) .ne. "none") then
            ! Advance the front: excess ice above H_eff in front cells moves
            ! to the ocean neighbours (transport, booked with dHidt_dyn)
            call calc_front_cells(mask_cf,mask_elig,mask_ocn,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                    bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%front_subgrid,tpo%par%boundaries)
            call calc_G_front_advance(mbal_now,tpo%now%H_ice,tpo%now%H_eff,mask_cf,mask_ocn, &
                                    dyn%now%ux_bar,dyn%now%uy_bar,dt,tpo%par%boundaries)
            call apply_tendency(tpo%now%H_ice,mbal_now,dt,"advance",adjust_mb=.TRUE.)
            tpo%now%dHidt_dyn = tpo%now%dHidt_dyn + mbal_now
            call update_ice_fraction(tpo,bnd)

            ! Treat fractional points that are not connected to full ice-covered points
            ! (a donor of the advance, left just below H_eff, counts as full)
            call calc_G_remove_fractional_ice(mbal_now,tpo%now%H_ice,tpo%now%f_ice,tpo%par%H_min_tau,dt, &
                                                tpo%par%boundaries,H_eff=tpo%now%H_eff)
        else
            ! Treat fractional points that are not connected to full ice-covered points
            call calc_G_remove_fractional_ice(mbal_now,tpo%now%H_ice,tpo%now%f_ice,tpo%par%H_min_tau,dt,tpo%par%boundaries)
        end if

        ! Apply rate and update ice thickness
        call apply_tendency(tpo%now%H_ice,mbal_now,dt,"frac",adjust_mb=.TRUE.)

        ! Add this rate to calving tendency
        tpo%now%cmb = tpo%now%cmb + mbal_now

        call update_ice_fraction(tpo,bnd)

        ! Kill methods: no floating ice may remain in the kill region ("kill": all
        ! floating ice, "kill-pos": floating ice where bnd%calv_mask). Applied last,
        ! after the front advance: otherwise ice the advance moves into the region
        ! survives the step and is removed in the next one, so the front cell fills
        ! and empties on alternate steps.
        select case(trim(tpo%par%calv_flt_method))

            case("kill","kill-pos")

                mask_kill = tpo%now%f_grnd .eq. 0.0_wp .and. tpo%now%H_ice .gt. 0.0_wp
                if (trim(tpo%par%calv_flt_method) .eq. "kill-pos") mask_kill = mask_kill .and. bnd%calv_mask

                call calc_calving_rate_kill(mbal_now,tpo%now%H_ice,mask_kill,tau=0.0_wp,dt=dt)
                call apply_tendency(tpo%now%H_ice,mbal_now,dt,"calv_kill",adjust_mb=.TRUE.)

                tpo%now%cmb_flt = tpo%now%cmb_flt + mbal_now
                tpo%now%cmb     = tpo%now%cmb     + mbal_now

                call update_ice_fraction(tpo,bnd)

        end select

        return

    end subroutine calc_ytopo_calving

    subroutine calc_ytopo_calving_lsf(tpo,dyn,mat,thrm,bnd,dt,H_prev,time_now)
        ! Calving computed as a flux. LSF mask is updated with velocity - calving front velocity.
        ! Ice in cells whose centre is in the ocean domain (LSF > 0) is deleted with a melt
        ! equal to the ice thickness. With a subgrid front, cells eligible for it follow the
        ! level-set area instead (calc_G_lsf_front).

        implicit none
    
        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd
        real(wp),           intent(IN)    :: dt
        real(wp),           intent(IN)    :: H_prev(:,:)
        real(wp), optional, intent(IN)    :: time_now
    
        ! Local variables
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: dt_kill
        real(wp), allocatable :: mbal_now(:,:)
        real(wp), allocatable :: a_lsf(:,:)
        real(wp), allocatable :: tau_eig_1(:,:), tau_eig_2(:,:)
        logical,  allocatable :: mask_cf(:,:), mask_elig(:,:), mask_ocn(:,:)
        !real(wp), allocatable :: u_acx_fill(:,:), v_acy_fill(:,:)
        integer  :: BC
        logical  :: is_front
        real(wp) :: cr_x, cr_y
        character(len=256) :: bnd_lsf

        ! Make sure dt is not zero
        dt_kill = dt 
        if (dt_kill .eq. 0.0) dt_kill = 1.0_wp
    
        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)

        ! Set boundary condition code
        BC = boundary_code(tpo%par%boundaries)

        allocate(mbal_now(nx,ny))
        allocate(a_lsf(nx,ny))
        allocate(mask_cf(nx,ny),mask_elig(nx,ny),mask_ocn(nx,ny))
        allocate(tau_eig_1(nx,ny),tau_eig_2(nx,ny))

        ! Principal stresses on the current ice geometry (cells that received ice
        ! since the last velocity solution take their neighbours' stresses)
        call fill_stress_new_ice(tau_eig_1,tau_eig_2,mat%now%strs2D%tau_eig_1,mat%now%strs2D%tau_eig_2, &
                                    tpo%now%f_ice,dyn%now%f_ice_solv,tpo%par%boundaries)

        ! === Floating calving laws ===
        
        ! Initialize the calving rates
        tpo%now%cmb_flt_x = 0.0_wp
        tpo%now%cmb_flt_y = 0.0_wp
        tpo%now%cmb_flt   = 0.0_wp

        select case(trim(tpo%par%calv_flt_method))
    
            case("zero","none")
                ! Do nothing. No calving.
                
            case("equil")
                ! For an equilibrated ice sheet calving rates should be opposite to ice velocity
                tpo%now%cmb_flt_x = -1*dyn%now%ux_bar
                tpo%now%cmb_flt_y = -1*dyn%now%uy_bar
    
            case("threshold")
                call calc_calving_threshold_lsf(tpo%now%cmb_flt_x,tpo%now%cmb_flt_y,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice,tpo%par%Hc_ref_flt,tpo%now%f_ice,tpo%par%boundaries)
        
            case("vm-m16")
                call calc_calving_rate_vonmises_m16(tpo%now%cmb_flt_x,tpo%now%cmb_flt_y,dyn%now%ux_bar,dyn%now%uy_bar,tau_eig_1,tpo%par%tau_ice_flt,tpo%now%f_ice,tpo%par%boundaries)
                
            ! TO DO: Add new laws
    
            ! === CalvMIP laws ===
            case("exp1","exp3")
                call calvmip_exp1(tpo%now%cmb_flt_x,tpo%now%cmb_flt_y,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%lsf,tpo%par%dx,tpo%par%boundaries) 
    
            case("exp2","exp4")
                call calvmip_exp2(tpo%now%cmb_flt_x,tpo%now%cmb_flt_y,dyn%now%ux_bar,dyn%now%uy_bar,time_now,tpo%par%boundaries)
            
            case("exp5")
                call calvmip_exp5_aa(tpo%now%cmb_flt_x,tpo%now%cmb_flt_y,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice,tpo%par%Hc_ref_flt,tpo%now%f_ice,tpo%par%boundaries)

            case DEFAULT
    
                write(*,*) "calc_ytopo:: Error: floating calving method not recognized."
                write(*,*) "calv_flt_method = ", trim(tpo%par%calv_flt_method)
                error stop 1
    
        end select
    
        ! === Marine terminating calving laws ===

        ! Initialize the calving rates
        tpo%now%cmb_grnd_x = 0.0_wp
        tpo%now%cmb_grnd_y = 0.0_wp
        tpo%now%cmb_grnd   = 0.0_wp

        select case(trim(tpo%par%calv_grnd_method))

            case("zero","none")
                ! Do nothing. No calving.

            case("equil")
                ! For an equilibrated ice sheet calving rates should be opposite to ice velocity
                tpo%now%cmb_grnd_x = -1*dyn%now%ux_bar
                tpo%now%cmb_grnd_y = -1*dyn%now%uy_bar

            case("threshold")
                ! Ice thickness threshold.
                call calc_calving_threshold_lsf(tpo%now%cmb_grnd_x,tpo%now%cmb_grnd_y,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice,tpo%par%Hc_ref_grnd,tpo%now%f_ice,tpo%par%boundaries)
        
            case("vm-m16")
                call calc_calving_rate_vonmises_m16(tpo%now%cmb_grnd_x,tpo%now%cmb_grnd_y,dyn%now%ux_bar,dyn%now%uy_bar,tau_eig_1,tpo%par%tau_ice_grnd,tpo%now%f_ice,tpo%par%boundaries)    

            case("ismip7")
                ! Retreat of marine-terminating glaciers following ISMIP7 protocol
                call calc_fmb_ismip7(tpo%now%cmb_grnd_x,tpo%now%cmb_grnd_y,tpo%now%lsf, &
                                     bnd%z_bed,bnd%z_sl,bnd%Qd,bnd%T_shlf,bnd%c%T0,tpo%par%dx,tpo%now%f_ice,tpo%par%boundaries)

            case DEFAULT
                ! To do: Add new laws
                ! MICI should be a marine terminating calving law (only for grounding-line points?)
                write(*,*) "calc_ytopo:: Error: grounded calving method not recognized."
                write(*,*) "calv_grnd_method = ", trim(tpo%par%calv_grnd_method)
                error stop 1
    
        end select
        
        ! === Land terminating calving laws ===
        ! For the moment we will assume no calving laws for land-terminating ice points.
        ! Only deformation.
    
        ! === Merge all calving law ===
        ! Merge all calving-rates into a single velocity field.
        ! Using ac-nodes for indices now.
        tpo%now%cr_acx = 0.0_wp
        tpo%now%cr_acy = 0.0_wp
        
        do j=1,ny
        do i=1,nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            ! x-ac node
            if (tpo%now%f_grnd_acx(i,j) .eq. 0.0) then
                ! Floating point
                tpo%now%cr_acx(i,j) = tpo%now%cmb_flt_x(i,j)
            else
                if (bnd%z_bed(i,j) .ge. bnd%z_sl(i,j) .and. bnd%z_bed(ip1,j) .ge. bnd%z_sl(ip1,j)) then
                    ! Face between two land cells (lsf pinned to -1 in both): the
                    ! level set does not move here. A face between a land and a
                    ! marine cell stays open, so the level set can follow the ice
                    ! across the coast.
                    tpo%now%cr_acx(i,j) = -1*dyn%now%ux_bar(i,j)
                else
                    ! Marine-terminating point.
                    tpo%now%cr_acx(i,j) = tpo%now%cmb_grnd_x(i,j)
                end if
            end if
                
            ! y-ac node
            if (tpo%now%f_grnd_acy(i,j) .eq. 0.0) then
                ! Floating point
                tpo%now%cr_acy(i,j) = tpo%now%cmb_flt_y(i,j)
            else
                if (bnd%z_bed(i,j) .ge. bnd%z_sl(i,j) .and. bnd%z_bed(i,jp1) .ge. bnd%z_sl(i,jp1)) then
                    ! Face between two land cells (see x-direction)
                    tpo%now%cr_acy(i,j) = -1*dyn%now%uy_bar(i,j)
                else
                    ! Marine-terminating point.
                    tpo%now%cr_acy(i,j) = tpo%now%cmb_grnd_y(i,j)
                end if
            end if
    
        end do
        end do

        ! === LSF advection ===
        ! Boundaries for the LSF advection: the model-wide tpo%par%boundaries
        ! (periodic directions wrap, as in LSFsnap and the area fraction), but
        ! "infinite" (Neumann-zero) instead of zero (Dirichlet) borders: the LSF is
        ! a signed-distance field that must continue smoothly outside the
        ! domain. A Dirichlet-zero boundary would create a spurious LSF=0
        ! contour one cell from the boundary that the Sussman/Osher
        ! redistance fights every step (see issue #34 follow-up). Matches
        ! Yelmo.jl, whose Oceananigans `:bounded` BC zeros only the halo,
        ! leaving edge cells free.
        bnd_lsf = tpo%par%boundaries
        if (trim(bnd_lsf) .eq. "zeros") bnd_lsf = "infinite"

        call LSFupdate(tpo%now%dlsfdt,tpo%now%lsf,tpo%now%cr_acx,tpo%now%cr_acy,dyn%now%ux_bar,dyn%now%uy_bar, &
                       tpo%par%dx,tpo%par%dy,dt,bnd_lsf)

        ! Marine points where ice is not allowed (bnd%mask_ice = MASK_ICE_NONE, 
        ! where H_ice is held at zero) are ocean by definition: keep the LSF
        ! at its ocean value there, so that the front cannot advance into them.
        where(bnd%mask_ice .eq. MASK_ICE_NONE .and. bnd%z_bed .lt. bnd%z_sl) tpo%now%lsf = 1.0_wp

        ! LSF should not affect grounded land points, i.e. points whose bed
        ! is at or above sea level. The comparison is inclusive (.ge.) so that
        ! a flat bed sitting exactly at sea level (e.g. the EISMINT/HALFAR
        ! benchmarks, z_bed = z_sl = 0) is treated as land and not calved;
        ! a strict .gt. left such domains flagged as ocean (lsf = 1) and the
        ! cmb loop below deleted all ice every step. Pin BEFORE the LSF
        ! discipline so that phi0 carries the right (land = -1) value into the
        ! Sussman/Osher iteration — Yelmo.jl uses this order. The snap pass
        ! after the cmb loop also benefits from a pre-pinned lsf.
        where(bnd%z_bed .ge. bnd%z_sl) tpo%now%lsf = -1.0_wp

        ! Choose LSF discipline: "redist" runs here (palma-ice #34 placement);
        ! "snap" runs after the cmb loop below, matching the pre-#34 ordering
        ! where the neighbour-snap was interleaved with cmb.
        !
        ! For "redist": lsf is in normalized ±1 units (LSFupdate saturates
        ! it to that range), so we redistance in grid-cell units (dx=dy=1)
        ! so that the PDE drives |grad lsf| -> 1 per cell near the zero
        ! level set, producing lsf ≈ ±1 at adjacent cells. Passing physical
        ! dx (e.g. 25000 m) would make the smoothed sign function
        ! ≈ ±lsf/dx ≈ 0 several cells out from the front, freezing the
        ! front in place (see issue #34). Matches Yelmo.jl. Boundaries as
        ! for the LSF advection call above (periodic wraps, otherwise
        ! Neumann-zero, the only sensible non-periodic BC for the SO
        ! redistance of a signed-distance field).
        select case(trim(tpo%par%lsf_method))
            case("redist")
                if (tpo%par%lsf_redist_n_iter .le. 0) then
                    write(io_unit_err,*) "calc_ytopo_calving_lsf:: Error: &
                        &lsf_method = 'redist' requires lsf_redist_n_iter > 0; &
                        &got lsf_redist_n_iter = ", tpo%par%lsf_redist_n_iter
                    error stop 1
                end if
                call LSFredistance(tpo%now%lsf,1.0_wp,1.0_wp, &
                                   tpo%par%lsf_redist_n_iter,bnd_lsf)
            case("snap")
                ! Handled after cmb loop below.
            case default
                write(io_unit_err,*) "calc_ytopo_calving_lsf:: Error: &
                    &unknown lsf_method = '"//trim(tpo%par%lsf_method)//"'. &
                    &Expected 'snap' or 'redist'."
                error stop 1
        end select

        ! === Calving ===
        ! Apply calving as a melt rate equal to ice thickness where lsf is positive,
        ! except in cells eligible for the subgrid front (none with front_subgrid="none"),
        ! which are emptied by area below (calc_G_lsf_front)
        call calc_front_cells(mask_cf,mask_elig,mask_ocn,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                              bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%front_subgrid,tpo%par%boundaries)

        tpo%now%cmb = 0.0_wp
        do j=1,ny
        do i=1,nx
            if (tpo%now%lsf(i,j) .gt. 0.0_wp .and. .not. mask_elig(i,j)) then
                ! Calve ice outside LSF mask (cmb = H_ice / dt_kill)
                tpo%now%cmb(i,j) = -(tpo%now%H_ice(i,j) / dt_kill)
            end if
        end do
        end do

        ! Legacy snap-mode LSF discipline runs after the cmb loop so that
        ! cmb sees the same pre-snap lsf(i,j) as in the pre-#34 interleaved
        ! loop. cmb does not read neighbour lsf, so splitting the passes is
        ! bitwise-equivalent to the inline version.
        if (trim(tpo%par%lsf_method) .eq. "snap") then
            call LSFsnap(tpo%now%lsf,time_now,tpo%par%dt_lsf,tpo%par%boundaries)
        end if

        ! Apply rate and update ice thickness
        call apply_tendency(tpo%now%H_ice,tpo%now%cmb,dt,"calving_lsf",adjust_mb=.TRUE.)

        if (trim(tpo%par%front_subgrid) .ne. "none") then
            ! Subgrid front: front-cell thickness follows the level set
            ! (at most a_lsf*H_eff), so f_ice = H_ice/H_eff is the area
            ! behind the front
            call calc_lsf_area_fraction(a_lsf,tpo%now%lsf,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,tpo%par%boundaries)
            call calc_G_lsf_front(mbal_now,tpo%now%H_ice,a_lsf,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice,bnd%c%rho_sw, &
                                  tpo%par%front_subgrid,tpo%par%front_H_eff_min,tpo%par%front_dHdx, &
                                  tpo%par%dx,dt_kill,tpo%par%boundaries)
            call apply_tendency(tpo%now%H_ice,mbal_now,dt,"calving_lsf_front",adjust_mb=.TRUE.)
            tpo%now%cmb = tpo%now%cmb + mbal_now
        end if

        ! Update ice fraction mask 
        call update_ice_fraction(tpo,bnd)

        ! Treat fractional points that are not connected to full ice-covered points
        call calc_G_remove_fractional_ice(mbal_now,tpo%now%H_ice,tpo%now%f_ice,tpo%par%H_min_tau,dt,tpo%par%boundaries)

        ! Apply rate and update ice thickness
        !mbal_now = 0.0_wp  ! ajr, commented this out, as it is zeroed out above. Otherwise
        ! the above fractional ice is not calculated and nothing should happen.
        call apply_tendency(tpo%now%H_ice,mbal_now,dt,"frac",adjust_mb=.TRUE.)

        ! Add this rate to calving tendency
        tpo%now%cmb = tpo%now%cmb + mbal_now

        call update_ice_fraction(tpo,bnd)

        ! Ice-free marine cells behind the front (for example emptied by oceanic melt) are
        ! returned to the ocean side of the lsf mask, unless an edge neighbour holds ice
        ! that can fill them as the front advances
        select case(trim(tpo%par%calv_flt_method))
            case("equil")
                    ! Do nothing here
            case DEFAULT
                !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
                do j = 1, ny
                do i = 1, nx
                    call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
                    if (tpo%now%H_ice(i,j) .le. 0.0_wp .and. tpo%now%lsf(i,j) .lt. 0.0_wp .and. &
                        bnd%z_bed(i,j) .lt. bnd%z_sl(i,j) .and. &
                        .not. (tpo%now%H_ice(im1,j) .gt. 0.0_wp .or. tpo%now%H_ice(ip1,j) .gt. 0.0_wp .or. &
                               tpo%now%H_ice(i,jm1) .gt. 0.0_wp .or. tpo%now%H_ice(i,jp1) .gt. 0.0_wp)) then
                        tpo%now%lsf(i,j) = 1.0_wp
                    end if
                end do
                end do
                !$omp end parallel do
        end select 

        ! Diagnostic: calving speed of the front [m/yr], the magnitude of the face
        ! rates the level set used (cr_acx/cr_acy, law chosen by the face f_grnd_acx/acy)
        ! averaged to the cell centre, at front cells (ice cells with an ice-free
        ! ocean edge neighbour), floating (f_grnd = 0) or grounded
        tpo%now%calv_rate_flt  = 0.0_wp
        tpo%now%calv_rate_grnd = 0.0_wp

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,is_front,cr_x,cr_y)
        do j = 1, ny
        do i = 1, nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            is_front = tpo%now%H_ice(i,j) .gt. 0.0_wp .and. bnd%z_bed(i,j) .lt. bnd%z_sl(i,j) .and. &
                ( is_ocean(im1,j) .or. is_ocean(ip1,j) .or. is_ocean(i,jm1) .or. is_ocean(i,jp1) )

            if (is_front) then
                cr_x = 0.5_wp*(tpo%now%cr_acx(im1,j)+tpo%now%cr_acx(i,j))
                cr_y = 0.5_wp*(tpo%now%cr_acy(i,jm1)+tpo%now%cr_acy(i,j))
                if (tpo%now%f_grnd(i,j) .eq. 0.0_wp) then
                    tpo%now%calv_rate_flt(i,j)  = sqrt(cr_x**2 + cr_y**2)
                else
                    tpo%now%calv_rate_grnd(i,j) = sqrt(cr_x**2 + cr_y**2)
                end if
            end if
        end do
        end do
        !$omp end parallel do

        return

    contains

        logical function is_ocean(ii,jj)
            ! Ice-free cell with the bed below sea level
            integer, intent(IN) :: ii, jj
            is_ocean = tpo%now%H_ice(ii,jj) .eq. 0.0_wp .and. bnd%z_bed(ii,jj) .lt. bnd%z_sl(ii,jj)
        end function is_ocean
    
    end subroutine calc_ytopo_calving_lsf

    subroutine calc_ytopo_diagnostic(tpo,dyn,mat,thrm,bnd)
        ! Calculate adjustments to surface elevation, bedrock elevation
        ! and ice thickness 

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 

        ! Local variables
        integer  :: gz_nx, gz_ny
        logical  :: gz_perx, gz_pery
        integer,  allocatable :: mask_src(:,:)
        real(wp), allocatable :: dist_cells(:,:)

        ! Final update of ice fraction mask (or define it now for fixed topography)
        call update_ice_fraction(tpo,bnd)

        ! Geometry of the active ice column, used for the surface and the
        ! dynamics: cells with at least A_FRONT_MIN of their area ice covered
        ! are active, and partial front cells (ytopo.front_subgrid) take part
        ! as full cells with thickness H_eff. A full front cell holding more
        ! ice than its reference H_eff keeps its own column (H_eff remains the
        ! reference of the front advance and trim). Cells below A_FRONT_MIN
        ! keep their ice (they fill by transport) but are ice free for the
        ! dynamics: otherwise a film of a few millimetres would be a full
        ! H_eff column (front_H_eff_min at least) and turn the ocean face of
        ! the neighbouring front cell into an interior face (GRL-8KM: such
        ! cells moved at ~8 km/yr). With front_subgrid="none", f_ice is
        ! binary, H_ice_dyn == H_ice and f_ice_dyn == f_ice.
        where (tpo%now%f_ice .ge. A_FRONT_MIN)
            tpo%now%H_ice_dyn = max(tpo%now%H_eff,tpo%now%H_ice)
            tpo%now%f_ice_dyn = 1.0_wp
        elsewhere
            tpo%now%H_ice_dyn = 0.0_wp
            tpo%now%f_ice_dyn = 0.0_wp
        end where

        ! Calculate grounding overburden ice thickness (from the actual thickness)
        call calc_H_grnd(tpo%now%H_grnd,tpo%now%H_ice,tpo%now%f_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice,bnd%c%rho_sw, &
                                                                                            use_f_ice=.FALSE.)

        ! Calculate the surface elevation of the ice column (H_eff in partial front cells)
        call calc_z_srf_max(tpo%now%z_srf,tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice,bnd%c%rho_sw)
        
        ! Calculate the ice base elevation too
        ! Define z_base as the elevation at the base of the ice sheet 
        ! This is used for the basal derivative instead of bedrock so
        ! that it is valid for both grounded and floating ice. Note, 
        ! for grounded ice, z_base==z_bed.  
        tpo%now%z_base = tpo%now%z_srf - tpo%now%H_ice_dyn 

        ! 2. Calculate additional topographic properties ------------------

        ! Calculate the surface slope and the ice thickness gradient (on staggered acx/y nodes)
        call calc_gradient_acx(tpo%now%dzsdx,tpo%now%z_srf,tpo%now%f_ice_dyn,tpo%par%dx,tpo%par%grad_lim,zero_outside=.FALSE.,boundaries=tpo%par%boundaries,slope_bg=tpo%par%slope_bg_x)
        call calc_gradient_acy(tpo%now%dzsdy,tpo%now%z_srf,tpo%now%f_ice_dyn,tpo%par%dy,tpo%par%grad_lim,zero_outside=.FALSE.,boundaries=tpo%par%boundaries,slope_bg=tpo%par%slope_bg_y)
        
        call calc_gradient_acx(tpo%now%dHidx,tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%par%dx,tpo%par%grad_lim,zero_outside=.TRUE.,boundaries=tpo%par%boundaries)
        call calc_gradient_acy(tpo%now%dHidy,tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%par%dy,tpo%par%grad_lim,zero_outside=.TRUE.,boundaries=tpo%par%boundaries)
        
        call calc_gradient_acx(tpo%now%dzbdx,tpo%now%z_base,tpo%now%f_ice_dyn,tpo%par%dx,tpo%par%grad_lim,zero_outside=.FALSE.,boundaries=tpo%par%boundaries,slope_bg=tpo%par%slope_bg_x)
        call calc_gradient_acy(tpo%now%dzbdy,tpo%now%z_base,tpo%now%f_ice_dyn,tpo%par%dy,tpo%par%grad_lim,zero_outside=.FALSE.,boundaries=tpo%par%boundaries,slope_bg=tpo%par%slope_bg_y)

        ! 3. Calculate new masks ------------------------------

        ! Calculate the grounded fraction and grounding line mask of each grid cell
        select case(tpo%par%gl_sep)

            case(1) 
                ! Binary f_grnd, linear f_grnd_acx/acy based on H_grnd

                call calc_f_grnd_subgrid_linear(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%now%H_grnd, &
                                                                tpo%par%boundaries)

            case(3) 
                ! Grounded area of H_grnd interpolated bilinearly between cell centres,
                ! analytical solutions of Leguy et al. (2021)

                call determine_grounded_fractions(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy, &
                                                                tpo%now%f_grnd_ab,tpo%now%H_grnd,tpo%par%boundaries)

        end select
        
        ! Calculate grounded fraction due to pinning points 
        call calc_f_grnd_pinning_points(tpo%now%f_grnd_pin,tpo%now%H_ice,tpo%now%f_ice, &
                                                bnd%z_bed,bnd%z_bed_sd,bnd%z_sl,bnd%c%rho_ice,bnd%c%rho_sw)

        ! === Grounding-zone and ice-margin distances (signed, in km) ===========
        ! Strategy: use calc_distance_to_*(calc_distances=.FALSE.) only to locate
        ! the grounding line / ice margin (0 at location, signed sentinel
        ! elsewhere: negative floating/ice-free, positive grounded/ice-covered),
        ! then compute the actual signed distance field with an (optionally
        ! periodic) chamfer transform via compute_distance_to_mask.

        gz_nx = size(tpo%now%dist_grline,1)
        gz_ny = size(tpo%now%dist_grline,2)

        allocate(mask_src(gz_nx,gz_ny))
        allocate(dist_cells(gz_nx,gz_ny))

        ! Derive periodic flags from the boundary condition string.
        ! (bcs mapping per solver_advection.f90: MISMIP3D/TROUGH => y-periodic.)
        select case(trim(tpo%par%boundaries))
            case("periodic","periodic-xy")
                gz_perx = .TRUE.  ; gz_pery = .TRUE.
            case("periodic-x")
                gz_perx = .TRUE.  ; gz_pery = .FALSE.
            case("periodic-y","MISMIP3D","TROUGH")
                gz_perx = .FALSE. ; gz_pery = .TRUE.
            case DEFAULT
                gz_perx = .FALSE. ; gz_pery = .FALSE.
        end select

        ! --- Grounding-line distance ---

        ! Locate grounding line: 0 at GL, -sentinel floating, +sentinel grounded
        call calc_distance_to_grounding_line(tpo%now%dist_grline,tpo%now%f_grnd,tpo%par%dx, &
                                                    tpo%par%boundaries,calc_distances=.FALSE.)

        ! Build source mask for chamfer: -1 floating (inside), 0 at GL, +1 grounded (outside)
        mask_src = 0
        where(tpo%now%dist_grline < 0.0_wp) mask_src = -1
        where(tpo%now%dist_grline > 0.0_wp) mask_src =  1

        ! Signed distance in grid cells, then convert to km
        call compute_distance_to_mask(dist_cells,mask_src,periodic_x=gz_perx,periodic_y=gz_pery)
        tpo%now%dist_grline = dist_cells * (tpo%par%dx*1.0e-3_wp)

        ! Define the grounding-zone mask from the signed km distance
        call calc_grounding_line_zone(tpo%now%mask_grz,tpo%now%dist_grline,tpo%par%dist_grz)

        ! --- Ice-margin distance ---

        ! Locate ice margin: 0 at margin, -sentinel ice-free, +sentinel ice-covered
        call calc_distance_to_ice_margin(tpo%now%dist_margin,tpo%now%f_ice,tpo%par%dx, &
                                                    tpo%par%boundaries,calc_distances=.FALSE.)

        ! Build source mask: -1 ice-free (inside), 0 at margin, +1 ice-covered (outside)
        mask_src = 0
        where(tpo%now%dist_margin < 0.0_wp) mask_src = -1
        where(tpo%now%dist_margin > 0.0_wp) mask_src =  1

        ! Signed distance in grid cells, then convert to km
        call compute_distance_to_mask(dist_cells,mask_src,periodic_x=gz_perx,periodic_y=gz_pery)
        tpo%now%dist_margin = dist_cells * (tpo%par%dx*1.0e-3_wp)

        deallocate(mask_src)
        deallocate(dist_cells)

        ! Calculate the general bed mask
        call gen_mask_bed(tpo%now%mask_bed,tpo%now%f_ice,thrm%now%f_pmp, &
                                            tpo%now%f_grnd,tpo%now%mask_grz.eq.0)


        ! Calculate the ice-front mask (mainly for use in dynamics)
        ! Partial front cells are front cells, the ice-front boundary is on their ocean faces
        call calc_ice_front(tpo%now%mask_frnt,tpo%now%f_ice_dyn,tpo%now%f_grnd,bnd%z_bed,bnd%z_sl,tpo%par%boundaries)

        
        return 

    end subroutine calc_ytopo_diagnostic

    subroutine calc_ytopo_rates(tpo,bnd,time,dt,step,check_mb,overwrite)
        ! Calculate average rates over outer timestep

        implicit none

        type(ytopo_class),  intent(INOUT) :: tpo 
        type(ybound_class), intent(IN)    :: bnd
        real(wp),           intent(IN)    :: time
        real(wp),           intent(IN)    :: dt 
        character(len=*),   intent(IN)    :: step 
        logical,            intent(IN)    :: check_mb 
        logical, optional,  intent(IN)    :: overwrite 
        
        real(wp), parameter :: tol_dt = 1e-3

        select case(trim(step))

            case("init")
                ! Initialization of averaging fields - set to zero

                tpo%now%rates%dzsdt         = 0.0
                tpo%now%rates%dHidt         = 0.0
                tpo%now%rates%dHidt_dyn     = 0.0
                tpo%now%rates%mb_net        = 0.0
                tpo%now%rates%mb_relax      = 0.0
                tpo%now%rates%mb_resid      = 0.0
                tpo%now%rates%mb_clip       = 0.0
                tpo%now%rates%mb_err        = 0.0
                tpo%now%rates%smb           = 0.0
                tpo%now%rates%bmb           = 0.0
                tpo%now%rates%fmb           = 0.0
                tpo%now%rates%dmb           = 0.0
                tpo%now%rates%cmb           = 0.0
                tpo%now%rates%cmb_flt       = 0.0
                tpo%now%rates%cmb_grnd      = 0.0
                tpo%now%rates%dlsfdt        = 0.0

                tpo%now%rates%dt_tot = 0.0 

            case("step")
                ! Add current step to total

                tpo%now%rates%dzsdt         = tpo%now%rates%dzsdt       + tpo%now%dzsdt*dt
                tpo%now%rates%dHidt         = tpo%now%rates%dHidt       + tpo%now%dHidt*dt
                tpo%now%rates%dHidt_dyn     = tpo%now%rates%dHidt_dyn   + tpo%now%dHidt_dyn*dt
                tpo%now%rates%mb_net        = tpo%now%rates%mb_net      + tpo%now%mb_net*dt
                tpo%now%rates%mb_relax      = tpo%now%rates%mb_relax    + tpo%now%mb_relax*dt
                tpo%now%rates%mb_resid      = tpo%now%rates%mb_resid    + tpo%now%mb_resid*dt
                tpo%now%rates%mb_clip       = tpo%now%rates%mb_clip     + tpo%now%mb_clip*dt
                tpo%now%rates%mb_err        = tpo%now%rates%mb_err      + tpo%now%mb_err*dt
                tpo%now%rates%smb           = tpo%now%rates%smb         + tpo%now%smb*dt
                tpo%now%rates%bmb           = tpo%now%rates%bmb         + tpo%now%bmb*dt
                tpo%now%rates%fmb           = tpo%now%rates%fmb         + tpo%now%fmb*dt
                tpo%now%rates%dmb           = tpo%now%rates%dmb         + tpo%now%dmb*dt
                tpo%now%rates%cmb           = tpo%now%rates%cmb         + tpo%now%cmb*dt
                tpo%now%rates%cmb_flt       = tpo%now%rates%cmb_flt     + tpo%now%cmb_flt*dt
                tpo%now%rates%cmb_grnd      = tpo%now%rates%cmb_grnd    + tpo%now%cmb_grnd*dt
                tpo%now%rates%dlsfdt        = tpo%now%rates%dlsfdt      + tpo%now%dlsfdt*dt


                tpo%now%rates%dt_tot = tpo%now%rates%dt_tot + dt  
                
            case("final")
                ! Divide by total time to get average rate

                if (tpo%now%rates%dt_tot .gt. 0.0) then 

                    tpo%now%rates%dzsdt         = tpo%now%rates%dzsdt / tpo%now%rates%dt_tot
                    tpo%now%rates%dHidt         = tpo%now%rates%dHidt / tpo%now%rates%dt_tot
                    tpo%now%rates%dHidt_dyn     = tpo%now%rates%dHidt_dyn / tpo%now%rates%dt_tot
                    tpo%now%rates%mb_net        = tpo%now%rates%mb_net / tpo%now%rates%dt_tot
                    tpo%now%rates%mb_relax      = tpo%now%rates%mb_relax / tpo%now%rates%dt_tot
                    tpo%now%rates%mb_resid      = tpo%now%rates%mb_resid / tpo%now%rates%dt_tot
                    tpo%now%rates%mb_clip       = tpo%now%rates%mb_clip / tpo%now%rates%dt_tot
                    tpo%now%rates%mb_err        = tpo%now%rates%mb_err / tpo%now%rates%dt_tot
                    tpo%now%rates%smb           = tpo%now%rates%smb / tpo%now%rates%dt_tot
                    tpo%now%rates%bmb           = tpo%now%rates%bmb / tpo%now%rates%dt_tot
                    tpo%now%rates%fmb           = tpo%now%rates%fmb / tpo%now%rates%dt_tot
                    tpo%now%rates%dmb           = tpo%now%rates%dmb / tpo%now%rates%dt_tot
                    tpo%now%rates%cmb           = tpo%now%rates%cmb / tpo%now%rates%dt_tot
                    tpo%now%rates%cmb_flt       = tpo%now%rates%cmb_flt / tpo%now%rates%dt_tot
                    tpo%now%rates%cmb_grnd      = tpo%now%rates%cmb_grnd / tpo%now%rates%dt_tot
                    tpo%now%rates%dlsfdt        = tpo%now%rates%dlsfdt / tpo%now%rates%dt_tot
                    
                    ! Check that dt_tot matches outer dt value
                    if ( abs(dt - tpo%now%rates%dt_tot) .gt. tol_dt) then
                        write(*,*) "calc_ytopo_rates: dt, dt_tot : ", dt, tpo%now%rates%dt_tot
                    end if 

                end if

                if (present(overwrite)) then
                    if (overwrite) then
                        ! Overwrite the instantaneous rates with averaged rates for output

                        tpo%now%dzsdt       = tpo%now%rates%dzsdt
                        tpo%now%dHidt       = tpo%now%rates%dHidt
                        tpo%now%dHidt_dyn   = tpo%now%rates%dHidt_dyn
                        tpo%now%mb_net      = tpo%now%rates%mb_net
                        tpo%now%mb_relax    = tpo%now%rates%mb_relax
                        tpo%now%mb_resid    = tpo%now%rates%mb_resid
                        tpo%now%mb_clip     = tpo%now%rates%mb_clip
                        tpo%now%mb_err      = tpo%now%rates%mb_err
                        tpo%now%smb         = tpo%now%rates%smb
                        tpo%now%bmb         = tpo%now%rates%bmb
                        tpo%now%fmb         = tpo%now%rates%fmb
                        tpo%now%dmb         = tpo%now%rates%dmb
                        tpo%now%cmb         = tpo%now%rates%cmb
                        tpo%now%cmb_flt     = tpo%now%rates%cmb_flt
                        tpo%now%cmb_grnd    = tpo%now%rates%cmb_grnd
                        tpo%now%dlsfdt      = tpo%now%rates%dlsfdt

                    end if
                end if 

            case DEFAULT 

                write(io_unit_err,*) "calc_ytopo_rates:: Error: step name not recognized."
                write(io_unit_err,*) "step = ", trim(step)
                error stop 1

        end select

        if (check_mb) then 
            ! Perform mass balance check to make sure that mass is conserved

            call check_mass_conservation(tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,tpo%now%dHidt, &
                        tpo%now%mb_net,tpo%now%cmb,tpo%now%dHidt_dyn,tpo%now%smb,tpo%now%bmb, &
                        tpo%now%fmb,tpo%now%dmb,tpo%now%mb_resid,tpo%now%mb_clip,tpo%par%dx,bnd%c%sec_year,time,dt, &
                        units="km^3/yr",label=step)
                        
        end if 

        return

    end subroutine calc_ytopo_rates
    
    subroutine ytopo_par_load(par,filename,group_ytopo,group_ycalv,nx,ny,dx,init)

        type(ytopo_param_class), intent(OUT) :: par
        character(len=*),        intent(IN)  :: filename
        character(len=*),        intent(IN)  :: group_ytopo ! Usually "ytopo"
        character(len=*),        intent(IN)  :: group_ycalv ! calving group
        integer,                 intent(IN)  :: nx, ny
        real(wp),                intent(IN)  :: dx
        logical, optional,       intent(IN)  :: init

        ! Local variables
        logical :: init_pars

        ! Canonical defaults file (schema) and canonical group names.
        ! User par files only need to list overrides; missing parameters
        ! are taken from this file and unknowns are caught by nml_validate.
        character(len=*), parameter :: def_file  = "input/yelmo_defaults.nml"
        character(len=*), parameter :: def_ytopo = "ytopo"
        character(len=*), parameter :: def_ycalv = "ycalv"

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE.

        ! Reject typos in the user file under either group
        call nml_validate(filename,def_file,group_ytopo,defaults_group=def_ytopo)
        call nml_validate(filename,def_file,group_ycalv,defaults_group=def_ycalv)

        ! Store parameter values in output object
        call nml_read(filename,group_ytopo,"solver",            par%solver,           init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"grad_lim",          par%grad_lim,         init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"grad_lim_zb",       par%grad_lim_zb,      init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"slope_bg_x",        par%slope_bg_x,       init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"slope_bg_y",        par%slope_bg_y,       init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"front_subgrid",     par%front_subgrid,    init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"front_H_eff_min",   par%front_H_eff_min,  init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"front_dHdx",        par%front_dHdx,       init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"use_bmb",           par%use_bmb,          init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"topo_fixed",        par%topo_fixed,       init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"topo_rel",          par%topo_rel,         init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"topo_rel_tau",      par%topo_rel_tau,     init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"topo_rel_field",    par%topo_rel_field,   init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        ! Grounding line
        call nml_read(filename,group_ytopo,"bmb_gl_method",     par%bmb_gl_method,    init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"gl_sep",            par%gl_sep,           init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"gz_nx",             par%gz_nx,            init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        ! pmpt method
        call nml_read(filename,group_ytopo,"dist_grz",          par%dist_grz,         init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"gz_Hg0",            par%gz_Hg0,           init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"gz_Hg1",            par%gz_Hg1,           init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        ! dmb
        call nml_read(filename,group_ytopo,"dmb_method",        par%dmb_method,       init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"dmb_alpha_max",     par%dmb_alpha_max,    init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"dmb_tau",           par%dmb_tau,          init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"dmb_sigma_ref",     par%dmb_sigma_ref,    init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"dmb_m_d",           par%dmb_m_d,          init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"dmb_m_r",           par%dmb_m_r,          init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        ! fmb
        call nml_read(filename,group_ytopo,"fmb_method",        par%fmb_method,       init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"fmb_scale",         par%fmb_scale,        init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)
        call nml_read(filename,group_ytopo,"fmb_lambda",        par%fmb_lambda,       init=init_pars,defaults_file=def_file,defaults_group=def_ytopo)

        ! === read calving routine ===
        call nml_read(filename,group_ycalv,"use_lsf",           par%use_lsf,            init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"lsf_method",        par%lsf_method,         init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"dt_lsf",            par%dt_lsf,             init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"lsf_redist_n_iter", par%lsf_redist_n_iter,  init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"calv_flt_method",   par%calv_flt_method,    init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"calv_grnd_method",  par%calv_grnd_method,   init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        ! ?
        call nml_read(filename,group_ycalv,"H_min_grnd",        par%H_min_grnd,         init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"H_min_flt",         par%H_min_flt,          init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"H_min_tau",         par%H_min_tau,          init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"sd_min",            par%sd_min,             init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"sd_max",            par%sd_max,             init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"calv_grnd_max",     par%calv_grnd_max,      init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        !
        call nml_read(filename,group_ycalv,"calv_tau",          par%calv_tau,           init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"calv_thin",         par%calv_thin,          init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"k2",                par%k2,                 init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"w2",                par%w2,                 init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"kt_ref",            par%kt_ref,             init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"kt_deep",           par%kt_deep,            init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"tau_ice_flt",       par%tau_ice_flt,        init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"tau_ice_grnd",      par%tau_ice_grnd,       init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        ! Threshold method
        call nml_read(filename,group_ycalv,"Hc_ref_flt",        par%Hc_ref_flt,       init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"Hc_ref_grnd",       par%Hc_ref_grnd,      init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"Hc_ref_thin",       par%Hc_ref_thin,      init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"Hc_deep",           par%Hc_deep,          init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"zb_deep_0",         par%zb_deep_0,        init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"zb_deep_1",         par%zb_deep_1,        init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)
        call nml_read(filename,group_ycalv,"zb_sigma",          par%zb_sigma,         init=init_pars,defaults_file=def_file,defaults_group=def_ycalv)

        ! === Validate parameter values ====
        call yelmo_check_enum(group_ytopo,"solver",         par%solver,         "none|expl|expl-new|expl-upwind|impl-upwind|expl-sico|impl-sico|impl-sico-lis|impl-lis")
        call yelmo_check_enum(group_ytopo,"front_subgrid",  par%front_subgrid,  "none|floating|marine")
        call yelmo_check_enum(group_ytopo,"bmb_gl_method",  par%bmb_gl_method,  "fcmp|fmp|pmp|pmpt|nmp")
        call yelmo_check_enum(group_ytopo,"topo_rel_field", par%topo_rel_field, "H_ref|H_ice_n")
        call yelmo_check_enum(group_ycalv,"lsf_method",     par%lsf_method,     "snap|redist")
        ! Allowed calv_flt_method values differ between the level-set and mass-balance calving paths
        ! (see calc_ytopo_calving_lsf and calc_ytopo_calving, respectively).
        if (par%use_lsf) then
            call yelmo_check_enum(group_ycalv,"calv_flt_method",par%calv_flt_method, &
                                  "zero|none|equil|threshold|vm-m16|exp1|exp2|exp3|exp4|exp5")
        else
            call yelmo_check_enum(group_ycalv,"calv_flt_method",par%calv_flt_method, &
                                  "zero|none|threshold|vm-l19|eigen|kill|kill-pos")
        end if
        ! As with calv_flt_method, the allowed calv_grnd_method values differ between the
        ! level-set and mass-balance calving paths.
        if (par%use_lsf) then
            call yelmo_check_enum(group_ycalv,"calv_grnd_method",par%calv_grnd_method, &
                                  "zero|none|equil|threshold|vm-m16|ismip7")
        else
            call yelmo_check_enum(group_ycalv,"calv_grnd_method",par%calv_grnd_method, &
                                  "zero|none|stress-b12")
        end if

        if (par%gl_sep .eq. 2) then
            ! The grounded area of gl_sep = 2 interpolated between the corner means of
            ! H_grnd only, so a cell grounded at its centre next to deep ocean had f_grnd = 0
            write(io_unit_err,*) "ytopo_par_load:: error: ytopo.gl_sep = 2 is deprecated; use gl_sep = 3."
            error stop 1
        else if (par%gl_sep .ne. 1 .and. par%gl_sep .ne. 3) then
            write(io_unit_err,*) "ytopo_par_load:: error: ytopo.gl_sep must be 1 or 3; got ", par%gl_sep
            error stop 1
        end if
        if (par%dt_lsf .gt. 0.0_wp .and. par%dt_lsf .lt. 0.01_wp) then
            ! LSFsnap checks the reflag time on a 0.01 yr resolution (nint(time*100)),
            ! so smaller positive intervals are not representable (and nint(dt_lsf*100)=0).
            write(io_unit_err,*) "ytopo_par_load:: error: ycalv.dt_lsf must be <= 0 (disabled) &
                                 &or >= 0.01 yr; got ", par%dt_lsf
            error stop 1
        end if
        if (par%front_H_eff_min .lt. 0.0_wp .or. par%front_dHdx .lt. 0.0_wp) then
            write(io_unit_err,*) "ytopo_par_load:: error: front_H_eff_min and front_dHdx must be >= 0; got ", &
                                 par%front_H_eff_min, par%front_dHdx
            error stop 1
        end if
        if (par%grad_lim .le. 0.0_wp) then
            write(io_unit_err,*) "ytopo_par_load:: error: grad_lim must be > 0; got ", par%grad_lim
            error stop 1
        end if
        if (par%grad_lim_zb .le. 0.0_wp) then
            write(io_unit_err,*) "ytopo_par_load:: error: grad_lim_zb must be > 0; got ", par%grad_lim_zb
            error stop 1
        end if
        if (par%H_min_tau .lt. 0.0_wp) then
            write(io_unit_err,*) "ytopo_par_load:: error: ycalv.H_min_tau must be >= 0; got ", par%H_min_tau
            error stop 1
        end if
        if (par%tau_ice_flt .le. 0.0_wp .or. par%tau_ice_grnd .le. 0.0_wp) then
            write(io_unit_err,*) "ytopo_par_load:: error: ycalv.tau_ice_flt and tau_ice_grnd must be > 0; got ", &
                                 par%tau_ice_flt, par%tau_ice_grnd
            error stop 1
        end if
        if (par%sd_min .ge. par%sd_max) then
            write(io_unit_err,*) "ytopo_par_load:: error: ycalv.sd_min must be < ycalv.sd_max; got ", &
                                 par%sd_min, par%sd_max
            error stop 1
        end if
        if (par%zb_deep_0 .lt. par%zb_deep_1) then
            write(io_unit_err,*) "ytopo_par_load:: error: ycalv.zb_deep_0 must be >= ycalv.zb_deep_1 &
                                 &(both negative; transition starts at zb_deep_0); got ", &
                                 par%zb_deep_0, par%zb_deep_1
            error stop 1
        end if

        ! === Set internal parameters ====
        par%nx  = nx 
        par%ny  = ny 
        par%dx  = dx 
        par%dy  = dx 

        ! Define how boundaries of grid should be treated 
        ! This should only be modified by the dom%par%experiment variable
        ! in yelmo_init. By default set boundaries to zero 
        par%boundaries = "zeros" 
        
        ! Define current time as unrealistic value
        par%time      = 1000000000   ! [a] 1 billion years in the future 

        ! Intialize timestepping parameters to Forward Euler (beta2=beta4=0: no contribution from previous timestep)
        par%dt_zeta     = 1.0 
        par%dt_beta(1)  = 1.0 
        par%dt_beta(2)  = 0.0 
        par%dt_beta(3)  = 1.0 
        par%dt_beta(4)  = 0.0 

        ! Set some additional values to start out right
        par%speed      = 0.0_wp 
        par%cpu_step   = 0.0d0

        par%adv_lin_iter = 0
        par%adv_lin_fail = 0 

        return

    end subroutine ytopo_par_load
    
    subroutine ytopo_alloc(now,nx,ny)

        implicit none 

        type(ytopo_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny  

        call ytopo_dealloc(now)

        call ytopo_pc_alloc(now%pred,nx,ny)
        call ytopo_pc_alloc(now%corr,nx,ny)

        ! Rates (for timestep averages)
        allocate(now%rates%dzsdt(nx,ny))
        allocate(now%rates%dHidt(nx,ny))
        allocate(now%rates%dHidt_dyn(nx,ny))
        allocate(now%rates%mb_net(nx,ny))
        allocate(now%rates%mb_relax(nx,ny))
        allocate(now%rates%mb_resid(nx,ny))
        allocate(now%rates%mb_clip(nx,ny))
        allocate(now%rates%mb_err(nx,ny))
        allocate(now%rates%smb(nx,ny))
        allocate(now%rates%bmb(nx,ny))
        allocate(now%rates%fmb(nx,ny))
        allocate(now%rates%dmb(nx,ny))
        allocate(now%rates%cmb(nx,ny))
        allocate(now%rates%cmb_flt(nx,ny))
        allocate(now%rates%cmb_grnd(nx,ny))
        allocate(now%rates%dlsfdt(nx,ny))
        
        ! Remaining ytopo fields...

        allocate(now%H_ice(nx,ny))
        allocate(now%z_srf(nx,ny))
        allocate(now%z_base(nx,ny))

        allocate(now%dzsdt(nx,ny))
        allocate(now%dHidt(nx,ny))
        allocate(now%dHidt_dyn(nx,ny))
        allocate(now%mb_net(nx,ny))
        allocate(now%mb_relax(nx,ny))
        allocate(now%mb_resid(nx,ny))
        allocate(now%mb_clip(nx,ny))
        allocate(now%mb_err(nx,ny))
        allocate(now%smb(nx,ny))
        allocate(now%bmb(nx,ny))
        allocate(now%fmb(nx,ny))
        allocate(now%dmb(nx,ny))
        allocate(now%cmb(nx,ny))
        allocate(now%cmb_flt(nx,ny))
        allocate(now%cmb_flt_x(nx,ny))
        allocate(now%cmb_flt_y(nx,ny))
        allocate(now%cmb_grnd(nx,ny))
        allocate(now%calv_rate_flt(nx,ny))
        allocate(now%calv_rate_grnd(nx,ny))
        allocate(now%cmb_grnd_x(nx,ny))
        allocate(now%cmb_grnd_y(nx,ny))
        allocate(now%cr_acx(nx,ny))
        allocate(now%cr_acy(nx,ny))

        allocate(now%lsf(nx,ny))       
        allocate(now%dlsfdt(nx,ny))
        
        allocate(now%bmb_ref(nx,ny))
        allocate(now%fmb_ref(nx,ny))
        allocate(now%dmb_ref(nx,ny))

        allocate(now%eps_eff(nx,ny))
        allocate(now%tau_eff(nx,ny))
        
        allocate(now%dzsdx(nx,ny))
        allocate(now%dzsdy(nx,ny))
        allocate(now%dHidx(nx,ny))
        allocate(now%dHidy(nx,ny))
        allocate(now%dzbdx(nx,ny))
        allocate(now%dzbdy(nx,ny))

        allocate(now%dzsdx_aa(nx,ny))
        allocate(now%dzsdy_aa(nx,ny))
        allocate(now%dHidx_aa(nx,ny))
        allocate(now%dHidy_aa(nx,ny))
        allocate(now%dzbdx_aa(nx,ny))
        allocate(now%dzbdy_aa(nx,ny))

        allocate(now%H_eff(nx,ny))
        allocate(now%H_grnd(nx,ny))
        allocate(now%H_calv(nx,ny))
        allocate(now%kt(nx,ny))
        allocate(now%z_bed_filt(nx,ny))

        ! Masks 
        allocate(now%f_grnd(nx,ny))
        allocate(now%f_grnd_acx(nx,ny))
        allocate(now%f_grnd_acy(nx,ny))
        allocate(now%f_grnd_ab(nx,ny))
        allocate(now%f_ice(nx,ny))

        allocate(now%f_grnd_bmb(nx,ny))
        allocate(now%f_grnd_pin(nx,ny))

        allocate(now%dist_margin(nx,ny))
        allocate(now%dist_grline(nx,ny))

        allocate(now%mask_bed(nx,ny))
        allocate(now%mask_grz(nx,ny))
        allocate(now%mask_frnt(nx,ny))
        
        allocate(now%dHidt_dyn_raw(nx,ny))
        allocate(now%dHidt_dyn_raw_n(nx,ny))
        allocate(now%dHidt_vert(nx,ny))
        allocate(now%dzsdt_kin(nx,ny))
        allocate(now%dzbdt_kin(nx,ny))
        allocate(now%mask_kin(nx,ny))
        allocate(now%H_ice_n(nx,ny))
        allocate(now%H_ice_dyn_n(nx,ny))
        allocate(now%z_srf_n(nx,ny))
        allocate(now%lsf_n(nx,ny))
        
        allocate(now%H_ice_dyn(nx,ny))
        allocate(now%f_ice_dyn(nx,ny))

        allocate(now%tau_relax(nx,ny))

        now%rates%dzsdt         = 0.0
        now%rates%dHidt         = 0.0
        now%rates%dHidt_dyn     = 0.0
        now%rates%mb_net        = 0.0
        now%rates%mb_relax      = 0.0
        now%rates%mb_resid      = 0.0
        now%rates%mb_clip       = 0.0
        now%rates%mb_err        = 0.0
        now%rates%smb           = 0.0
        now%rates%bmb           = 0.0
        now%rates%fmb           = 0.0
        now%rates%dmb           = 0.0
        now%rates%cmb           = 0.0
        now%rates%cmb_flt       = 0.0
        now%rates%cmb_grnd      = 0.0
        now%rates%dlsfdt        = 0.0

        now%H_ice       = 0.0 
        now%z_srf       = 0.0
        now%z_base      = 0.0  
        now%dzsdt       = 0.0 
        now%dHidt       = 0.0
        now%dHidt_dyn   = 0.0
        now%mb_net      = 0.0 
        now%mb_relax    = 0.0
        now%mb_resid    = 0.0
        now%mb_clip     = 0.0
        now%mb_err      = 0.0
        now%smb         = 0.0 
        now%bmb         = 0.0  
        now%fmb         = 0.0
        now%dmb         = 0.0
        now%cmb         = 0.0
        now%cmb_flt     = 0.0
        now%cmb_flt_x   = 0.0
        now%cmb_flt_y   = 0.0
        now%cmb_grnd    = 0.0
        now%calv_rate_flt  = 0.0
        now%calv_rate_grnd = 0.0
        now%lsf         = 1.0 ! init to 0.0?       
        now%dlsfdt      = 0.0
        
        now%bmb_ref     = 0.0  
        now%fmb_ref     = 0.0
        now%dmb_ref     = 0.0

        now%eps_eff     = 0.0
        now%tau_eff     = 0.0
        
        now%dzsdx       = 0.0 
        now%dzsdy       = 0.0 
        now%dHidx       = 0.0 
        now%dHidy       = 0.0
        now%dzbdx       = 0.0 
        now%dzbdy       = 0.0

        now%dzsdx_aa    = 0.0 
        now%dzsdy_aa    = 0.0 
        now%dHidx_aa    = 0.0 
        now%dHidy_aa    = 0.0
        now%dzbdx_aa    = 0.0 
        now%dzbdy_aa    = 0.0

        now%H_eff       = 0.0 
        now%H_grnd      = 0.0  
        now%H_calv      = 0.0  
        now%kt          = 0.0  
        now%z_bed_filt  = 0.0  

        now%f_grnd      = 0.0  
        now%f_grnd_acx  = 0.0  
        now%f_grnd_acy  = 0.0  
        now%f_grnd_ab   = 0.0
        now%f_grnd_bmb  = 0.0
        now%f_grnd_pin  = 0.0
        now%f_ice       = 0.0  
        now%dist_margin = 0.0
        now%dist_grline = 0.0 

        now%mask_bed    = 0 
        now%mask_grz    = 0 
        now%mask_frnt   = 0

        now%dHidt_dyn_raw   = 0.0
        now%dHidt_dyn_raw_n = 0.0
        now%dHidt_vert      = 0.0
        now%dzsdt_kin       = 0.0
        now%dzbdt_kin       = 0.0
        now%mask_kin        = 0
        now%H_ice_n     = 0.0
        now%H_ice_dyn_n = 0.0
        now%z_srf_n     = 0.0
        now%lsf_n     = 0.0 

        now%H_ice_dyn   = 0.0 
        now%f_ice_dyn   = 0.0 
        
        now%tau_relax   = 0.0

        return 

    end subroutine ytopo_alloc

    subroutine ytopo_dealloc(now)

        implicit none 

        type(ytopo_state_class), intent(INOUT) :: now

        call ytopo_pc_dealloc(now%pred)
        call ytopo_pc_dealloc(now%corr)

        if (allocated(now%rates%dzsdt))         deallocate(now%rates%dzsdt)
        if (allocated(now%rates%dHidt))         deallocate(now%rates%dHidt)
        if (allocated(now%rates%dHidt_dyn))     deallocate(now%rates%dHidt_dyn)
        if (allocated(now%rates%mb_net))        deallocate(now%rates%mb_net)
        if (allocated(now%rates%mb_relax))      deallocate(now%rates%mb_relax)
        if (allocated(now%rates%mb_resid))      deallocate(now%rates%mb_resid)
        if (allocated(now%rates%mb_clip))       deallocate(now%rates%mb_clip)
        if (allocated(now%rates%mb_err))        deallocate(now%rates%mb_err)
        if (allocated(now%rates%smb))           deallocate(now%rates%smb)
        if (allocated(now%rates%bmb))           deallocate(now%rates%bmb)
        if (allocated(now%rates%fmb))           deallocate(now%rates%fmb)
        if (allocated(now%rates%dmb))           deallocate(now%rates%dmb)
        if (allocated(now%rates%cmb))           deallocate(now%rates%cmb)
        if (allocated(now%rates%cmb_flt))       deallocate(now%rates%cmb_flt)
        if (allocated(now%rates%cmb_grnd))      deallocate(now%rates%cmb_grnd)
        if (allocated(now%rates%dlsfdt))        deallocate(now%rates%dlsfdt)
        
        if (allocated(now%H_ice))       deallocate(now%H_ice)
        if (allocated(now%z_srf))       deallocate(now%z_srf)
        if (allocated(now%z_base))      deallocate(now%z_base)

        if (allocated(now%dzsdt))       deallocate(now%dzsdt)
        if (allocated(now%dHidt))       deallocate(now%dHidt)
        if (allocated(now%dHidt_dyn))   deallocate(now%dHidt_dyn)
        if (allocated(now%mb_net))      deallocate(now%mb_net)
        if (allocated(now%mb_relax))    deallocate(now%mb_relax)
        if (allocated(now%mb_resid))    deallocate(now%mb_resid)
        if (allocated(now%mb_clip))     deallocate(now%mb_clip)
        if (allocated(now%mb_err))      deallocate(now%mb_err)
        if (allocated(now%smb))         deallocate(now%smb)
        if (allocated(now%bmb))         deallocate(now%bmb)
        if (allocated(now%fmb))         deallocate(now%fmb)
        if (allocated(now%dmb))         deallocate(now%dmb)
        if (allocated(now%cmb))         deallocate(now%cmb)
        if (allocated(now%cmb_flt))     deallocate(now%cmb_flt)
        if (allocated(now%cmb_flt_x))   deallocate(now%cmb_flt_x)
        if (allocated(now%cmb_flt_y))   deallocate(now%cmb_flt_y)
        if (allocated(now%cmb_grnd))    deallocate(now%cmb_grnd)
        if (allocated(now%calv_rate_flt))  deallocate(now%calv_rate_flt)
        if (allocated(now%calv_rate_grnd)) deallocate(now%calv_rate_grnd)
        if (allocated(now%cmb_grnd_x))  deallocate(now%cmb_grnd_x)
        if (allocated(now%cmb_grnd_y))  deallocate(now%cmb_grnd_y)
        if (allocated(now%cr_acx))      deallocate(now%cr_acx)
        if (allocated(now%cr_acy))      deallocate(now%cr_acy)
        if (allocated(now%lsf))         deallocate(now%lsf)       
        if (allocated(now%dlsfdt))      deallocate(now%dlsfdt)
        
        if (allocated(now%bmb_ref))     deallocate(now%bmb_ref)
        if (allocated(now%fmb_ref))     deallocate(now%fmb_ref)
        if (allocated(now%dmb_ref))     deallocate(now%dmb_ref)

        if (allocated(now%eps_eff))     deallocate(now%eps_eff)
        if (allocated(now%tau_eff))     deallocate(now%tau_eff)
        
        if (allocated(now%dzsdx))       deallocate(now%dzsdx)
        if (allocated(now%dzsdy))       deallocate(now%dzsdy)
        if (allocated(now%dHidx))       deallocate(now%dHidx)
        if (allocated(now%dHidy))       deallocate(now%dHidy)
        if (allocated(now%dzbdx))       deallocate(now%dzbdx)
        if (allocated(now%dzbdy))       deallocate(now%dzbdy)
        
        if (allocated(now%dzsdx_aa))       deallocate(now%dzsdx_aa)
        if (allocated(now%dzsdy_aa))       deallocate(now%dzsdy_aa)
        if (allocated(now%dHidx_aa))       deallocate(now%dHidx_aa)
        if (allocated(now%dHidy_aa))       deallocate(now%dHidy_aa)
        if (allocated(now%dzbdx_aa))       deallocate(now%dzbdx_aa)
        if (allocated(now%dzbdy_aa))       deallocate(now%dzbdy_aa)
        
        if (allocated(now%H_eff))       deallocate(now%H_eff)
        if (allocated(now%H_grnd))      deallocate(now%H_grnd)
        if (allocated(now%H_calv))      deallocate(now%H_calv)
        if (allocated(now%kt))          deallocate(now%kt)
        if (allocated(now%z_bed_filt))  deallocate(now%z_bed_filt)

        if (allocated(now%f_grnd))      deallocate(now%f_grnd)
        if (allocated(now%f_grnd_acx))  deallocate(now%f_grnd_acx)
        if (allocated(now%f_grnd_acy))  deallocate(now%f_grnd_acy)
        if (allocated(now%f_grnd_ab))   deallocate(now%f_grnd_ab)
        if (allocated(now%f_grnd_bmb))  deallocate(now%f_grnd_bmb)
        if (allocated(now%f_grnd_pin))  deallocate(now%f_grnd_pin)

        if (allocated(now%f_ice))       deallocate(now%f_ice)

        if (allocated(now%dist_margin)) deallocate(now%dist_margin)
        if (allocated(now%dist_grline)) deallocate(now%dist_grline)

        if (allocated(now%mask_bed))    deallocate(now%mask_bed)
        if (allocated(now%mask_grz))    deallocate(now%mask_grz)
        if (allocated(now%mask_frnt))   deallocate(now%mask_frnt)
        
        if (allocated(now%dHidt_dyn_raw))   deallocate(now%dHidt_dyn_raw)
        if (allocated(now%dHidt_dyn_raw_n)) deallocate(now%dHidt_dyn_raw_n)
        if (allocated(now%dHidt_vert))      deallocate(now%dHidt_vert)
        if (allocated(now%dzsdt_kin))       deallocate(now%dzsdt_kin)
        if (allocated(now%dzbdt_kin))       deallocate(now%dzbdt_kin)
        if (allocated(now%mask_kin))        deallocate(now%mask_kin)
        if (allocated(now%H_ice_n))     deallocate(now%H_ice_n)
        if (allocated(now%H_ice_dyn_n)) deallocate(now%H_ice_dyn_n)
        if (allocated(now%z_srf_n))     deallocate(now%z_srf_n)
        if (allocated(now%lsf_n))       deallocate(now%lsf_n)
        
        if (allocated(now%H_ice_dyn))   deallocate(now%H_ice_dyn)
        if (allocated(now%f_ice_dyn))   deallocate(now%f_ice_dyn)
        
        if (allocated(now%tau_relax))   deallocate(now%tau_relax)
        
        return 

    end subroutine ytopo_dealloc
    
    subroutine ytopo_pc_alloc(pc,nx,ny)

        implicit none

        type(ytopo_pc_class), intent(INOUT) :: pc 
        integer, intent(IN) :: nx, ny  

        ! First deallocate everything for safety
        call ytopo_pc_dealloc(pc)

        ! Allocate fields 
        allocate(pc%H_ice(nx,ny))
        allocate(pc%dHidt_dyn(nx,ny))
        allocate(pc%dHidt_vert(nx,ny))
        allocate(pc%mb_net(nx,ny))
        allocate(pc%mb_relax(nx,ny))
        allocate(pc%mb_resid(nx,ny))
        allocate(pc%mb_clip(nx,ny))
        allocate(pc%smb(nx,ny))
        allocate(pc%bmb(nx,ny))
        allocate(pc%fmb(nx,ny))
        allocate(pc%dmb(nx,ny))
        allocate(pc%cmb(nx,ny))      
        allocate(pc%cmb_flt(nx,ny))
        allocate(pc%cmb_grnd(nx,ny))
        allocate(pc%lsf(nx,ny))
        allocate(pc%cmb_flt_x(nx,ny))
        allocate(pc%cmb_flt_y(nx,ny))
        allocate(pc%cmb_grnd_x(nx,ny))
        allocate(pc%cmb_grnd_y(nx,ny))
        allocate(pc%cr_acx(nx,ny))
        allocate(pc%cr_acy(nx,ny))
        allocate(pc%calv_rate_flt(nx,ny))
        allocate(pc%calv_rate_grnd(nx,ny))
        
        ! Initialize to zero
        pc%H_ice        = 0.0
        pc%dHidt_dyn    = 0.0
        pc%dHidt_vert   = 0.0
        pc%mb_net       = 0.0
        pc%mb_relax     = 0.0
        pc%mb_resid     = 0.0
        pc%mb_clip      = 0.0
        pc%smb          = 0.0
        pc%bmb          = 0.0
        pc%fmb          = 0.0
        pc%dmb          = 0.0
        pc%cmb          = 0.0      
        pc%cmb_flt      = 0.0 
        pc%cmb_grnd     = 0.0
        pc%lsf          = 0.0            
        pc%cmb_flt_x      = 0.0
        pc%cmb_flt_y      = 0.0
        pc%cmb_grnd_x     = 0.0
        pc%cmb_grnd_y     = 0.0
        pc%cr_acx         = 0.0
        pc%cr_acy         = 0.0
        pc%calv_rate_flt  = 0.0
        pc%calv_rate_grnd = 0.0
        
        return

    end subroutine ytopo_pc_alloc

    subroutine ytopo_pc_dealloc(pc)

        implicit none

        type(ytopo_pc_class), intent(INOUT) :: pc 
        
        if (allocated(pc%H_ice))        deallocate(pc%H_ice)
        if (allocated(pc%dHidt_dyn))    deallocate(pc%dHidt_dyn)
        if (allocated(pc%dHidt_vert))   deallocate(pc%dHidt_vert)
        if (allocated(pc%mb_net))       deallocate(pc%mb_net)
        if (allocated(pc%mb_relax))     deallocate(pc%mb_relax)
        if (allocated(pc%mb_resid))     deallocate(pc%mb_resid)
        if (allocated(pc%mb_clip))      deallocate(pc%mb_clip)
        if (allocated(pc%smb))          deallocate(pc%smb)
        if (allocated(pc%bmb))          deallocate(pc%bmb)
        if (allocated(pc%fmb))          deallocate(pc%fmb)
        if (allocated(pc%dmb))          deallocate(pc%dmb)
        if (allocated(pc%cmb))          deallocate(pc%cmb)
        if (allocated(pc%cmb_flt))      deallocate(pc%cmb_flt)
        if (allocated(pc%cmb_grnd))     deallocate(pc%cmb_grnd)
        if (allocated(pc%lsf))          deallocate(pc%lsf)
        if (allocated(pc%cmb_flt_x)) deallocate(pc%cmb_flt_x)
        if (allocated(pc%cmb_flt_y)) deallocate(pc%cmb_flt_y)
        if (allocated(pc%cmb_grnd_x)) deallocate(pc%cmb_grnd_x)
        if (allocated(pc%cmb_grnd_y)) deallocate(pc%cmb_grnd_y)
        if (allocated(pc%cr_acx)) deallocate(pc%cr_acx)
        if (allocated(pc%cr_acy)) deallocate(pc%cr_acy)
        if (allocated(pc%calv_rate_flt)) deallocate(pc%calv_rate_flt)
        if (allocated(pc%calv_rate_grnd)) deallocate(pc%calv_rate_grnd)
        
        return

    end subroutine ytopo_pc_dealloc

    subroutine calc_transport_velocity(ux_t,uy_t,tpo,dyn,bnd,filter_vel)
        ! Depth-averaged velocity that transports ice thickness: the current
        ! solution, or with filter_vel the mean of the current and previous
        ! solutions (dyn%now fields remain the true solution). Faces that
        ! carry no ice are zero (set_inactive_margins, from the current
        ! tpo%now%f_ice): faces from partial cells into ice-free cells, except,
        ! with a subgrid front following the level set, into cells the front
        ! covers by at least A_FRONT_MIN.

        implicit none

        real(wp),           intent(OUT) :: ux_t(:,:)
        real(wp),           intent(OUT) :: uy_t(:,:)
        type(ytopo_class),  intent(IN)  :: tpo
        type(ydyn_class),   intent(IN)  :: dyn
        type(ybound_class), intent(IN)  :: bnd
        logical,            intent(IN)  :: filter_vel

        real(wp), allocatable :: a_front(:,:)

        if (filter_vel) then
            ux_t = 0.5_wp*(dyn%now%ux_bar + dyn%now%ux_bar_prev)
            uy_t = 0.5_wp*(dyn%now%uy_bar + dyn%now%uy_bar_prev)
        else
            ux_t = dyn%now%ux_bar
            uy_t = dyn%now%uy_bar
        end if

        if (tpo%par%use_lsf .and. trim(tpo%par%front_subgrid) .ne. "none") then
            allocate(a_front(size(tpo%now%H_ice,1),size(tpo%now%H_ice,2)))
            call calc_lsf_area_fraction(a_front,tpo%now%lsf,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,tpo%par%boundaries)
            call set_inactive_margins(ux_t,uy_t,tpo%now%f_ice,tpo%par%boundaries,a_front)
        else
            call set_inactive_margins(ux_t,uy_t,tpo%now%f_ice,tpo%par%boundaries)
        end if

        return

    end subroutine calc_transport_velocity

    subroutine update_ice_fraction(tpo,bnd)
        ! Update the ice area fraction tpo%now%f_ice and effective thickness
        ! tpo%now%H_eff from tpo%now%H_ice (CISM-style front scheme,
        ! ytopo.front_subgrid). With the level set, the front cells'
        ! thickness is trimmed to the level-set area first (calc_G_lsf_front),
        ! so the same f_ice = H_ice/H_eff holds in both calving paths. The
        ! front cells are classified as in the trim, from the current level
        ! set (a_lsf: cells cut by the front at an ocean corner are front cells).

        implicit none

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ybound_class), intent(IN)    :: bnd

        real(wp), allocatable :: a_lsf(:,:)

        if (tpo%par%use_lsf .and. trim(tpo%par%front_subgrid) .ne. "none") then
            allocate(a_lsf(size(tpo%now%H_ice,1),size(tpo%now%H_ice,2)))
            call calc_lsf_area_fraction(a_lsf,tpo%now%lsf,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,tpo%par%boundaries)
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_eff,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                   bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%front_subgrid, &
                                   tpo%par%front_H_eff_min,tpo%par%front_dHdx,tpo%par%dx,tpo%par%boundaries,a_lsf)
        else
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_eff,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                   bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%front_subgrid, &
                                   tpo%par%front_H_eff_min,tpo%par%front_dHdx,tpo%par%dx,tpo%par%boundaries)
        end if

        return

    end subroutine update_ice_fraction

end module yelmo_topography
