
module yelmo_tracers
    ! ytrc: passive-tracer subsystem controller.
    !
    ! Drives up to three interchangeable englacial-tracing backends, any subset
    ! of which may run simultaneously so their capabilities can be compared:
    !
    !   * euler  -- the in-tree Eulerian age tracer (ice_tracer::calc_tracer_3D),
    !               relocated here from ymat for v2.0.
    !   * trc    -- the Lagrangian particle model (tracer library).
    !   * elsa   -- the Lagrangian layer model (elsa library).
    !
    ! Each backend produces a gridded deposition-time field on the host sigma
    ! grid (t_dep_euler / t_dep_trc / t_dep_elsa); one is designated authoritative
    ! (t_dep_source) and copied into t_dep, from which the isochrone depths are
    ! diagnosed. elsa and tracer additionally own their native state objects and
    ! (later) their own NetCDF output/restart.
    !
    ! Note: the enh_bnd "*-tracer" advection stays in ymat -- it is a genuine
    ! material property and only needs the Eulerian solver.

    use nml

    use yelmo_defs
    use yelmo_tools, only : stagger_acx_aa, stagger_acy_aa
    use ice_tracer,  only : calc_tracer_3D, calc_isochrones

    use elsa,   only : elsa_init, elsa_update, elsa_end
    use tracer, only : tracer_init, tracer_update, tracer_end

    implicit none

    private
    public :: ytrc_par_load, ytrc_alloc, ytrc_dealloc
    public :: ytrc_init, calc_ytrc, ytrc_end

contains

    subroutine calc_ytrc(trc,tpo,dyn,thrm,bnd,grd,time)
        ! Advance all enabled passive-tracer backends by one host step and
        ! refresh the harmonized deposition-time diagnostics.

        implicit none

        type(ytrc_class),   intent(INOUT) :: trc
        type(ytopo_class),  intent(IN)    :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd
        type(grid_class),   intent(IN)    :: grd
        real(wp),           intent(IN)    :: time

        ! Local variables
        integer  :: nx, ny, nz_aa
        real(wp) :: dt
        real(wp), allocatable :: X_srf(:,:)
        logical,  allocatable :: mask_tracers(:,:)

        nx    = trc%par%nx
        ny    = trc%par%ny
        nz_aa = trc%par%nz_aa

        ! Initialize time if necessary
        if (trc%par%time .gt. time) trc%par%time = time

        ! Get time step and advance current time
        dt           = real(time,wp) - real(trc%par%time,wp)
        trc%par%time = real(time,dp)

        ! === 1. Eulerian backend ==========================================
        if (trc%par%use_euler .and. trc%par%calc_age .and. dt .gt. 0.0_wp) then

            allocate(X_srf(nx,ny))
            allocate(mask_tracers(nx,ny))

            ! Surface boundary condition is the current time
            X_srf = time

            ! Restrict tracing away from very fast-flowing ice (surface value imposed there)
            mask_tracers = .TRUE.
            where (dyn%now%uxy_bar .gt. 500.0_wp) mask_tracers = .FALSE.

            call calc_tracer_3D(trc%now%t_dep_euler,X_srf,dyn%now%ux,dyn%now%uy,dyn%now%uz_star, &
                tpo%now%H_ice,tpo%now%bmb,trc%par%zeta_aa,trc%par%zeta_ac,trc%par%tracer_method, &
                trc%par%tracer_impl_kappa,dt,thrm%par%dx,time,mask=mask_tracers)

            deallocate(X_srf,mask_tracers)

        end if

        ! === 2. Lagrangian layer backend (elsa) ===========================
        if (trc%par%use_elsa) then
            call ytrc_update_elsa(trc,tpo,dyn,time)
            ! TODO (M2): harmonize elsa layer stack -> trc%now%t_dep_elsa
            trc%now%t_dep_elsa = MV
        end if

        ! === 3. Lagrangian particle backend (tracer) ======================
        if (trc%par%use_tracer) then
            call ytrc_update_tracer(trc,tpo,dyn,grd,time)
            ! TODO (M3): harmonize particle cloud -> trc%now%t_dep_trc
            trc%now%t_dep_trc = MV
        end if

        ! === 4. Authoritative field + isochrones ==========================
        select case(trim(trc%par%t_dep_source))
            case("euler")
                trc%now%t_dep = trc%now%t_dep_euler
            case("trc")
                trc%now%t_dep = trc%now%t_dep_trc
            case("elsa")
                trc%now%t_dep = trc%now%t_dep_elsa
            case DEFAULT
                write(io_unit_err,*) "calc_ytrc:: Error: t_dep_source not recognized: "//trim(trc%par%t_dep_source)
                stop "Program stopped."
        end select

        ! Diagnose isochrone depths from the authoritative deposition-time field
        call calc_isochrones(trc%now%depth_iso,trc%now%t_dep,tpo%now%H_ice,trc%par%age_iso, &
                                                                        trc%par%zeta_aa,time)

        return

    end subroutine calc_ytrc

    subroutine ytrc_update_elsa(trc,tpo,dyn,time)
        ! Drive elsa with the host's native (acx/acy, sigma) velocities. elsa
        ! maps them onto its own grid internally and self-paces via dt_coupling.

        implicit none

        type(ytrc_class),  intent(INOUT) :: trc
        type(ytopo_class), intent(IN)    :: tpo
        type(ydyn_class),  intent(IN)    :: dyn
        real(wp),          intent(IN)    :: time

        ! elsa is double-precision internally; pass dp copies to the _dp interface.
        call elsa_update(trc%elsa,real(time,dp),real(tpo%now%H_ice,dp), &
                         real(dyn%now%ux,dp),real(dyn%now%uy,dp), &
                         real(tpo%now%smb,dp),real(tpo%now%bmb,dp))

        return

    end subroutine ytrc_update_elsa

    subroutine ytrc_update_tracer(trc,tpo,dyn,grd,time)
        ! Drive the Lagrangian particle model. tracer wants cell-centred (aa)
        ! velocities on the zeta_aa sigma axis, so destagger ux/uy from acx/acy
        ! and interpolate uz from zeta_ac to zeta_aa. Deposition/stats cadence
        ! is governed by the tracer namelist (dt_dep, dt_write_stats).

        implicit none

        type(ytrc_class),  intent(INOUT) :: trc
        type(ytopo_class), intent(IN)    :: tpo
        type(ydyn_class),  intent(IN)    :: dyn
        type(grid_class),  intent(IN)    :: grd
        real(wp),          intent(IN)    :: time

        ! Local variables
        integer :: k, nx, ny, nz_aa
        logical :: dep_now, stats_now
        real(wp), allocatable :: ux_aa(:,:,:), uy_aa(:,:,:), uz_aa(:,:,:)

        nx    = trc%par%nx
        ny    = trc%par%ny
        nz_aa = trc%par%nz_aa

        allocate(ux_aa(nx,ny,nz_aa))
        allocate(uy_aa(nx,ny,nz_aa))
        allocate(uz_aa(nx,ny,nz_aa))

        ! Destagger horizontal velocities acx/acy -> aa, per level
        do k = 1, nz_aa
            ux_aa(:,:,k) = stagger_acx_aa(dyn%now%ux(:,:,k))
            uy_aa(:,:,k) = stagger_acy_aa(dyn%now%uy(:,:,k))
        end do

        ! Interpolate vertical velocity zeta_ac -> zeta_aa (yelmo uz is positive up)
        call interp_zeta_ac_to_aa(uz_aa,dyn%now%uz,trc%par%zeta_ac,trc%par%zeta_aa)

        ! Decide whether to deposit this step (cadence dt_dep from tracer namelist).
        ! M1: gridded stats harmonization not yet wired, so stats output is off.
        dep_now   = .FALSE.
        stats_now = .FALSE.
        if (real(time,dp) .ge. trc%par%time_dep_next) then
            dep_now = .TRUE.
            trc%par%time_dep_next = trc%par%time_dep_next + real(trc%trc%par%dt_dep,dp)
        end if

        call tracer_update(trc%trc,real(time,wp),real(grd%G%x,wp),real(grd%G%y,wp),trc%par%zeta_aa, &
                           tpo%now%z_srf,tpo%now%H_ice,ux_aa,uy_aa,uz_aa, &
                           dep_now=dep_now,stats_now=stats_now,sigma_srf=1.0_wp)

        deallocate(ux_aa,uy_aa,uz_aa)

        return

    end subroutine ytrc_update_tracer

    subroutine interp_zeta_ac_to_aa(u_aa,u_ac,zeta_ac,zeta_aa)
        ! Linearly interpolate a 3D field from the ac (interface) vertical axis
        ! onto the aa (cell-centre) axis, column by column. Endpoints clamp.

        implicit none

        real(wp), intent(OUT) :: u_aa(:,:,:)
        real(wp), intent(IN)  :: u_ac(:,:,:)
        real(wp), intent(IN)  :: zeta_ac(:)
        real(wp), intent(IN)  :: zeta_aa(:)

        ! Local variables
        integer  :: k, kac, nz_aa, nz_ac
        real(wp) :: wt

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        do k = 1, nz_aa
            if (zeta_aa(k) .le. zeta_ac(1)) then
                u_aa(:,:,k) = u_ac(:,:,1)
            else if (zeta_aa(k) .ge. zeta_ac(nz_ac)) then
                u_aa(:,:,k) = u_ac(:,:,nz_ac)
            else
                ! Find the ac interval containing zeta_aa(k)
                do kac = 1, nz_ac-1
                    if (zeta_aa(k) .ge. zeta_ac(kac) .and. zeta_aa(k) .le. zeta_ac(kac+1)) exit
                end do
                wt = (zeta_aa(k) - zeta_ac(kac)) / (zeta_ac(kac+1) - zeta_ac(kac))
                u_aa(:,:,k) = (1.0_wp-wt)*u_ac(:,:,kac) + wt*u_ac(:,:,kac+1)
            end if
        end do

        return

    end subroutine interp_zeta_ac_to_aa

    subroutine ytrc_init(trc,grd,time,H_ice)
        ! Initialize the enabled backends. Called once from yelmo_init_state,
        ! after the topography (H_ice) has been set.

        implicit none

        type(ytrc_class), intent(INOUT) :: trc
        type(grid_class), intent(IN)    :: grd
        real(wp),         intent(IN)    :: time
        real(wp),         intent(IN)    :: H_ice(:,:)

        ! Eulerian backend: seed deposition time to the current time
        trc%now%t_dep_euler = real(time,wp)

        ! Layer backend (elsa): allocate the layer stack sized by par%time_end
        if (trc%par%use_elsa) then
            call elsa_init(trc%elsa,trim(trc%par%elsa_nml),trim(trc%par%elsa_group), &
                           real(time,dp),real(trc%par%time_end,dp),grd%G%x,grd%G%y, &
                           real(trc%par%zeta_aa,dp),real(H_ice,dp),"acx_acy")
        end if

        ! Particle backend (tracer): allocate the fixed particle pool
        if (trc%par%use_tracer) then
            call tracer_init(trc%trc,trim(trc%par%tracer_nml),real(time,wp), &
                             real(grd%G%x,wp),real(grd%G%y,wp),is_sigma=.TRUE.,grid=grd)

            ! Initialize deposition/stats cadence trackers from the tracer namelist
            trc%par%time_dep_next   = real(time,dp)
            trc%par%time_stats_next = real(time,dp)
        end if

        return

    end subroutine ytrc_init

    subroutine ytrc_par_load(par,filename,group,zeta_aa,zeta_ac,nx,ny,dx,init)

        implicit none

        type(ytrc_param_class), intent(OUT) :: par
        character(len=*),       intent(IN)  :: filename
        character(len=*),       intent(IN)  :: group        ! Usually "ytrc"
        real(wp),               intent(IN)  :: zeta_aa(:)
        real(wp),               intent(IN)  :: zeta_ac(:)
        integer,                intent(IN)  :: nx, ny
        real(wp),               intent(IN)  :: dx
        logical, optional,      intent(IN)  :: init

        ! Local variables
        logical  :: init_pars
        real(wp) :: age_iso(10)

        character(len=*), parameter :: def_file = "input/yelmo_defaults.nml"
        character(len=*), parameter :: def_ytrc = "ytrc"

        age_iso = 0.0

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE.

        call nml_validate(filename,def_file,group,defaults_group=def_ytrc)

        call nml_read(filename,group,"use_euler",         par%use_euler,         init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"use_tracer",        par%use_tracer,        init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"use_elsa",          par%use_elsa,          init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"t_dep_source",      par%t_dep_source,      init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"time_end",          par%time_end,          init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"calc_age",          par%calc_age,          init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"age_iso",           age_iso,               init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"tracer_method",     par%tracer_method,     init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"tracer_impl_kappa", par%tracer_impl_kappa, init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"elsa_nml",          par%elsa_nml,          init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"elsa_group",        par%elsa_group,        init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)
        call nml_read(filename,group,"tracer_nml",        par%tracer_nml,        init=init_pars,defaults_file=def_file,defaults_group=def_ytrc)

        ! Validate parameter values
        call yelmo_check_enum(group,"t_dep_source",  par%t_dep_source,  "euler|trc|elsa")
        call yelmo_check_enum(group,"tracer_method", par%tracer_method, "expl|impl")

        ! Set internal parameters
        par%nx    = nx
        par%ny    = ny
        par%dx    = dx
        par%dy    = dx
        par%nz_aa = size(zeta_aa,1)
        par%nz_ac = size(zeta_ac,1)

        if (allocated(par%zeta_aa)) deallocate(par%zeta_aa)
        allocate(par%zeta_aa(par%nz_aa))
        par%zeta_aa = zeta_aa

        if (allocated(par%zeta_ac)) deallocate(par%zeta_ac)
        allocate(par%zeta_ac(par%nz_ac))
        par%zeta_ac = zeta_ac

        ! Number of isochrones follows the target ages (as in the former ymat path)
        if ( (.not. par%calc_age) .or. count(age_iso .eq. 0.0) .eq. size(age_iso)) then
            par%n_iso = 1
        else
            par%n_iso = count(age_iso .ne. 0.0)
        end if

        if (allocated(par%age_iso)) deallocate(par%age_iso)
        allocate(par%age_iso(par%n_iso))
        par%age_iso = age_iso(1:par%n_iso)

        ! Define current time as unrealistic value
        par%time = 1000000000   ! [a] 1 billion years in the future

        return

    end subroutine ytrc_par_load

    subroutine ytrc_alloc(now,nx,ny,nz_aa,n_iso)

        implicit none

        type(ytrc_state_class), intent(INOUT) :: now
        integer :: nx, ny, nz_aa, n_iso

        call ytrc_dealloc(now)

        allocate(now%t_dep(nx,ny,nz_aa))
        allocate(now%t_dep_euler(nx,ny,nz_aa))
        allocate(now%t_dep_trc(nx,ny,nz_aa))
        allocate(now%t_dep_elsa(nx,ny,nz_aa))
        allocate(now%depth_iso(nx,ny,n_iso))

        now%t_dep       = 0.0
        now%t_dep_euler = 0.0
        now%t_dep_trc   = MV
        now%t_dep_elsa  = MV
        now%depth_iso   = 0.0

        return

    end subroutine ytrc_alloc

    subroutine ytrc_dealloc(now)

        implicit none

        type(ytrc_state_class), intent(INOUT) :: now

        if (allocated(now%t_dep))        deallocate(now%t_dep)
        if (allocated(now%t_dep_euler))  deallocate(now%t_dep_euler)
        if (allocated(now%t_dep_trc))    deallocate(now%t_dep_trc)
        if (allocated(now%t_dep_elsa))   deallocate(now%t_dep_elsa)
        if (allocated(now%depth_iso))    deallocate(now%depth_iso)

        return

    end subroutine ytrc_dealloc

    subroutine ytrc_end(trc)

        implicit none

        type(ytrc_class), intent(INOUT) :: trc

        if (trc%par%use_elsa)   call elsa_end(trc%elsa)
        if (trc%par%use_tracer) call tracer_end(trc%trc)

        call ytrc_dealloc(trc%now)

        return

    end subroutine ytrc_end

end module yelmo_tracers
