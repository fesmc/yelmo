
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
    use ice_tracer,  only : calc_tracer_3D, calc_isochrones

    use elsa,   only : elsa_init, elsa_update, elsa_end, elsa_restart_write
    use tracer, only : tracer_init, tracer_update, tracer_end, &
                       tracer_write_init, tracer_write, tracer_read

    implicit none

    private
    public :: ytrc_par_load, ytrc_alloc, ytrc_dealloc
    public :: ytrc_init, calc_ytrc, ytrc_end
    public :: ytrc_restart_write

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
            call ytrc_harmonize_elsa(trc)
        end if

        ! === 3. Lagrangian particle backend (tracer) ======================
        if (trc%par%use_tracer) then
            call ytrc_update_tracer(trc,tpo,dyn,grd,time)
            call ytrc_harmonize_tracer(trc,grd)
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

    subroutine ytrc_harmonize_elsa(trc)
        ! Map elsa's isochrone layer stack onto the host sigma grid as a gridded
        ! deposition-time field (t_dep_elsa), column by column. Requires
        ! grid_factor=1 (checked at init), so elsa's grid is the host grid.
        !
        ! elsa layers (bottom to top):
        !   1 .. n_layers_init         initial equal-thickness column (ice deposited
        !                              at or before time_init -> floored to time_init)
        !   n_layers_init+1 .. n_top-1 layer L closed by the isochrone at its top,
        !                              deposition time = time_add(L - n_layers_init)
        !   n_top                      current accumulating layer, top = surface = now
        !
        ! dsum_iso(:,:,L) is the height of layer L's top above the bed; the column
        ! total dsum_iso(:,:,n_top) equals H_ice. Deposition time at a given height
        ! is linearly interpolated between the bounding isochrone knots.

        implicit none

        type(ytrc_class), intent(INOUT) :: trc

        ! Local variables
        integer  :: i, j, k, L, nlt, ntop, nx, ny, nz_aa, nk
        real(dp) :: h, tt, H_col, t_init
        real(dp), allocatable :: hk(:), tk(:)

        real(wp), parameter :: H_MIN_ELSA = 1.0e-6_wp

        nx     = trc%par%nx
        ny     = trc%par%ny
        nz_aa  = trc%par%nz_aa
        nlt    = trc%elsa%par%n_layers_init
        ntop   = trc%elsa%now%n_top
        t_init = trc%par%elsa_time_init

        trc%now%t_dep_elsa = MV

        ! Knots: a bed knot (h=0) plus every boundary from n_layers_init to n_top
        nk = (ntop - nlt + 1) + 1
        allocate(hk(nk),tk(nk))

        do j = 1, ny
        do i = 1, nx

            H_col = trc%elsa%now%dsum_iso(i,j,ntop)
            if (H_col .le. real(H_MIN_ELSA,dp)) cycle   ! no ice -> leave MV

            ! Build monotonic (height, time) knots for this column
            hk(1) = 0.0_dp
            tk(1) = t_init
            do L = nlt, ntop
                hk(L-nlt+2) = trc%elsa%now%dsum_iso(i,j,L)
                if (L .lt. nlt+1) then
                    tk(L-nlt+2) = t_init                              ! base of accumulation
                else if (L .lt. ntop) then
                    tk(L-nlt+2) = trc%elsa%par%time_add(L-nlt)        ! closed isochrone
                else
                    tk(L-nlt+2) = trc%elsa%now%time                  ! surface (now)
                end if
            end do

            do k = 1, nz_aa
                h = real(trc%par%zeta_aa(k),dp) * H_col
                call interp_monotonic(tt,h,hk,tk,nk)
                trc%now%t_dep_elsa(i,j,k) = real(tt,wp)
            end do

        end do
        end do

        deallocate(hk,tk)

        return

    end subroutine ytrc_harmonize_elsa

    subroutine interp_monotonic(y,x,xk,yk,nk)
        ! Linear interpolation of y(x) from monotonic non-decreasing knots
        ! (xk,yk). Clamps outside the range; skips zero-width (tied-xk) intervals.

        implicit none

        real(dp), intent(OUT) :: y
        real(dp), intent(IN)  :: x
        real(dp), intent(IN)  :: xk(:), yk(:)
        integer,  intent(IN)  :: nk

        integer  :: m
        real(dp) :: wt

        if (x .le. xk(1)) then
            y = yk(1)
            return
        else if (x .ge. xk(nk)) then
            y = yk(nk)
            return
        end if

        ! Find the interval [xk(m), xk(m+1)] containing x
        do m = 1, nk-1
            if (x .ge. xk(m) .and. x .le. xk(m+1)) then
                if (xk(m+1) .gt. xk(m)) then
                    wt = (x - xk(m)) / (xk(m+1) - xk(m))
                    y  = (1.0_dp-wt)*yk(m) + wt*yk(m+1)
                else
                    y = yk(m+1)   ! zero-width interval -> take upper knot
                end if
                return
            end if
        end do

        y = yk(nk)   ! fallback (should not be reached)

        return

    end subroutine interp_monotonic

    subroutine ytrc_update_tracer(trc,tpo,dyn,grd,time)
        ! Drive the Lagrangian particle model. tracer samples the host velocities
        ! on their native Arakawa-C locations and interpolates directly to each
        ! particle point (fesmc/tracer#1), so pass ux/uy/uz staggered together
        ! with their axes: ux on acx (x_ux = x+dx/2), uy on acy (y_uy = y+dy/2),
        ! uz on the zeta_ac interfaces (z_uz = zeta_ac). Deposition cadence is
        ! governed by the tracer namelist (dt_dep).

        implicit none

        type(ytrc_class),  intent(INOUT) :: trc
        type(ytopo_class), intent(IN)    :: tpo
        type(ydyn_class),  intent(IN)    :: dyn
        type(grid_class),  intent(IN)    :: grd
        real(wp),          intent(IN)    :: time

        ! Local variables
        integer :: nx, ny
        logical :: dep_now, stats_now
        real(wp) :: dx, dy
        real(wp), allocatable :: x_ux(:), y_uy(:)

        nx = trc%par%nx
        ny = trc%par%ny
        dx = real(grd%G%dx,wp)
        dy = real(grd%G%dy,wp)

        ! Staggered velocity axes: acx is the right cell border (x+dx/2), acy the
        ! top border (y+dy/2); uz already sits on the zeta_ac interfaces.
        allocate(x_ux(nx),y_uy(ny))
        x_ux = real(grd%G%x,wp) + 0.5_wp*dx
        y_uy = real(grd%G%y,wp) + 0.5_wp*dy

        ! Decide whether to deposit this step (cadence dt_dep from tracer namelist).
        ! Gridded-stats output is off; the harmonized t_dep_trc is built afterwards.
        dep_now   = .FALSE.
        stats_now = .FALSE.
        if (real(time,dp) .ge. trc%par%time_dep_next) then
            dep_now = .TRUE.
            trc%par%time_dep_next = trc%par%time_dep_next + real(trc%trc%par%dt_dep,dp)
        end if

        call tracer_update(trc%trc,real(time,wp), &
                           real(grd%G%x,wp),real(grd%G%y,wp),trc%par%zeta_aa, &
                           tpo%now%z_srf,tpo%now%H_ice, &
                           dyn%now%ux,dyn%now%uy,dyn%now%uz, &
                           x_ux=x_ux,y_uy=y_uy,z_uz=trc%par%zeta_ac, &
                           dep_now=dep_now,stats_now=stats_now,sigma_srf=1.0_wp)

        deallocate(x_ux,y_uy)

        return

    end subroutine ytrc_update_tracer

    subroutine ytrc_harmonize_tracer(trc,grd)
        ! Grid the Lagrangian particle cloud onto the host sigma grid as a gridded
        ! deposition-time field (t_dep_trc). Each active particle is binned to the
        ! nearest host cell (i,j) and nearest sigma level k by its position, and
        ! the mean deposition time of the particles in a cell is stored. Cells
        ! with no particle stay MV -- the field is intentionally gappy, its sparsity
        ! reflecting how well the particle cloud samples the ice at that point.

        implicit none

        type(ytrc_class), intent(INOUT) :: trc
        type(grid_class), intent(IN)    :: grd

        ! Local variables
        integer  :: p, np, i, j, k, nx, ny, nz_aa
        real(wp) :: x0, y0, dx, dy, sig
        real(wp), allocatable :: sum_t(:,:,:)
        integer,  allocatable :: cnt(:,:,:)

        nx    = trc%par%nx
        ny    = trc%par%ny
        nz_aa = trc%par%nz_aa

        x0 = real(grd%G%x(1),wp)
        y0 = real(grd%G%y(1),wp)
        dx = real(grd%G%dx,wp)
        dy = real(grd%G%dy,wp)

        allocate(sum_t(nx,ny,nz_aa)); sum_t = 0.0_wp
        allocate(cnt(nx,ny,nz_aa));   cnt   = 0

        np = size(trc%trc%now%active,1)

        do p = 1, np

            if (trc%trc%now%active(p) .lt. 1) cycle   ! only activated/live particles

            ! Nearest host cell centre (aa) in the horizontal
            i = nint((trc%trc%now%x(p) - x0)/dx) + 1
            j = nint((trc%trc%now%y(p) - y0)/dy) + 1
            if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny) cycle

            ! Nearest sigma level (sigma=1 at surface, matching zeta_aa)
            sig = trc%trc%now%sigma(p)
            k   = minloc(abs(trc%par%zeta_aa - sig),1)

            sum_t(i,j,k) = sum_t(i,j,k) + trc%trc%dep%time(p)
            cnt(i,j,k)   = cnt(i,j,k) + 1

        end do

        where (cnt .gt. 0)
            trc%now%t_dep_trc = sum_t / real(cnt,wp)
        elsewhere
            trc%now%t_dep_trc = MV
        end where

        deallocate(sum_t,cnt)

        return

    end subroutine ytrc_harmonize_tracer

    subroutine ytrc_init(trc,grd,time,H_ice,restart)
        ! Initialize the enabled backends. Called once from yelmo_init_state,
        ! after the topography (H_ice) has been set. On a restart, `restart` is
        ! the yelmo restart path; each backend restores its own native state from
        ! a sidecar file derived from it.

        implicit none

        type(ytrc_class), intent(INOUT) :: trc
        type(grid_class), intent(IN)    :: grd
        real(wp),         intent(IN)    :: time
        real(wp),         intent(IN)    :: H_ice(:,:)
        character(len=*), intent(IN), optional :: restart

        ! Local variables
        logical            :: is_restart
        integer            :: n_add
        character(len=512) :: elsa_rst

        is_restart = .FALSE.
        if (present(restart)) then
            if (trim(restart) .ne. "None" .and. len_trim(restart) .gt. 0) is_restart = .TRUE.
        end if

        ! Eulerian backend: seed deposition time to the current time on a cold
        ! start. On a restart, t_dep_euler is read from the yelmo restart file by
        ! yelmo_restart_read (which runs before this), so do not overwrite it.
        if (.not. is_restart) trc%now%t_dep_euler = real(time,wp)

        ! Layer backend (elsa): allocate the layer stack (sized by par%time_end on
        ! a cold start, or read from the sidecar on a restart).
        if (trc%par%use_elsa) then

            if (is_restart) then
                elsa_rst = ytrc_restart_filename(restart,"elsa")
                call elsa_init(trc%elsa,trim(trc%par%elsa_nml),trim(trc%par%elsa_group), &
                               real(time,dp),real(trc%par%time_end,dp),grd%G%x,grd%G%y, &
                               real(trc%par%zeta_aa,dp),real(H_ice,dp),"acx_acy", &
                               restart=trim(elsa_rst))
            else
                call elsa_init(trc%elsa,trim(trc%par%elsa_nml),trim(trc%par%elsa_group), &
                               real(time,dp),real(trc%par%time_end,dp),grd%G%x,grd%G%y, &
                               real(trc%par%zeta_aa,dp),real(H_ice,dp),"acx_acy")
            end if

            ! The t_dep_elsa diagnostic maps elsa's layer stack onto the host
            ! sigma grid column-by-column; this only holds when elsa shares the
            ! host grid (grid_factor=1). Enforce it.
            if (abs(trc%elsa%par%grid_factor - 1.0_wp) .gt. 1e-6_wp) then
                write(io_unit_err,*) "ytrc_init:: Error: use_elsa requires grid_factor=1 in the elsa &
                                     &namelist (for the t_dep_elsa diagnostic); got ", trc%elsa%par%grid_factor
                stop "Program stopped."
            end if

            ! Floor time for elsa's pre-init initial layers (see ytrc_harmonize_elsa).
            ! On a restart the original init time is not stored, so recover it from
            ! the isochrone schedule (exact for a regular layer_resolution).
            if (is_restart) then
                n_add = size(trc%elsa%par%time_add)
                if (n_add .ge. 2) then
                    trc%par%elsa_time_init = trc%elsa%par%time_add(1) &
                                           - (trc%elsa%par%time_add(2) - trc%elsa%par%time_add(1))
                else if (n_add .eq. 1) then
                    trc%par%elsa_time_init = trc%elsa%par%time_add(1)
                else
                    trc%par%elsa_time_init = real(time,dp)
                end if
            else
                trc%par%elsa_time_init = real(time,dp)
            end if

        end if

        ! Particle backend (tracer): allocate the fixed particle pool (tracer_init),
        ! then on a restart overlay the saved particle cloud from the sidecar
        ! (tracer_read must run after tracer_init; it only overwrites state).
        if (trc%par%use_tracer) then
            call tracer_init(trc%trc,trim(trc%par%tracer_nml),real(time,wp), &
                             real(grd%G%x,wp),real(grd%G%y,wp),is_sigma=.TRUE.,grid=grd)

            if (is_restart) then
                call tracer_read(trc%trc,trim(ytrc_restart_filename(restart,"tracer")),real(time,wp))
            end if

            ! Initialize deposition/stats cadence trackers from the tracer namelist
            trc%par%time_dep_next   = real(time,dp)
            trc%par%time_stats_next = real(time,dp)
        end if

        return

    end subroutine ytrc_init

    subroutine ytrc_restart_write(trc,filename,time)
        ! Write the native restart state of each enabled backend to a sidecar
        ! file alongside the yelmo restart (e.g. yelmo_restart.nc ->
        ! yelmo_restart_elsa.nc / _tracer.nc). The Eulerian backend needs no
        ! sidecar: its t_dep_euler field is written into the yelmo restart via
        ! the io%trc table.

        implicit none

        type(ytrc_class), intent(IN) :: trc
        character(len=*), intent(IN) :: filename
        real(wp),         intent(IN) :: time

        ! Local variables
        integer            :: is
        character(len=512) :: path, fldr, fname
        type(tracer_class) :: trc_tmp

        ! Guard on allocation: a restart file may be written before the backends
        ! are initialized (e.g. the diagnostic write inside yelmo_restart_read,
        ! which runs before ytrc_init), and there is nothing to save yet.
        if (trc%par%use_elsa .and. allocated(trc%elsa%now%d_iso)) then
            call elsa_restart_write(trc%elsa,trim(ytrc_restart_filename(filename,"elsa")))
        end if

        ! tracer_write takes a folder + filename and appends to a fresh archive;
        ! split the sidecar path and create it (tracer_write_init) before writing.
        ! tracer_write is INTENT(INOUT) (it stamps time_write), so operate on a
        ! local copy to keep this routine and its callers read-only in `trc`.
        if (trc%par%use_tracer .and. allocated(trc%trc%now%x)) then
            path = ytrc_restart_filename(filename,"tracer")
            is   = index(trim(path),"/",back=.TRUE.)
            if (is .gt. 0) then
                fldr  = path(1:is-1)
                fname = path(is+1:)
            else
                fldr  = "."
                fname = trim(path)
            end if
            trc_tmp = trc%trc
            call tracer_write_init(trc_tmp,trim(fldr),trim(fname))
            call tracer_write(trc_tmp,real(time,wp),trim(fldr),trim(fname))
        end if

        return

    end subroutine ytrc_restart_write

    function ytrc_restart_filename(base,tag) result(path)
        ! Build a sidecar restart path by inserting "_<tag>" before the ".nc"
        ! extension of `base` (e.g. ("dir/yelmo_restart.nc","elsa") ->
        ! "dir/yelmo_restart_elsa.nc"). Falls back to appending if no ".nc".

        implicit none

        character(len=*), intent(IN) :: base, tag
        character(len=512) :: path

        integer :: ie

        ie = index(base,".nc",back=.TRUE.)
        if (ie .gt. 0) then
            path = base(1:ie-1)//"_"//trim(tag)//".nc"
        else
            path = trim(base)//"_"//trim(tag)//".nc"
        end if

        return

    end function ytrc_restart_filename

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
