
module yelmo_boundaries

    use nml 
    use ncio 
    use yelmo_defs 
    use coords,  only : grid_class
    use regions, only : regions_class, region_mask_class, regions_init_nml, regions_select, &
                        regions_basin_ids, regions_end
    use yelmo_tools, only : boundary_code, get_periodic_directions

    implicit none
    
    private
    public :: ybound_define_physical_constants
    public :: ybound_load_masks
    public :: ybound_define_mask_ice
    public :: ybound_alloc, ybound_dealloc
    public :: ybound_update_rates
    
contains

    subroutine ybound_define_physical_constants(c,phys_const,domain,grid_name,cnst,cnst_out)
        ! Fill Yelmo's physical-constants record.
        !
        ! The shared quantities come from a phys_const_class (fesm-utils), either
        ! supplied by the driver or, when absent, loaded here from Yelmo's own
        ! parameter file as before. A driver that owns the constants passes them
        ! in, so every component of a coupled program works from one set; the
        ! no-argument path keeps the standalone drivers and the C API working.
        !
        ! c is Yelmo's working-precision mirror of that record. It is kept rather
        ! than replaced so the several hundred bnd%c%<name> uses across the
        ! physics are untouched, and because Yelmo carries two things that are
        ! deliberately not shared physical constants: sec_year, which is a
        ! calendar choice, and conv_mmdwe_maie, which embeds a day count.

        implicit none

        type(ybound_const_class), intent(OUT) :: c
        character(len=*), intent(IN) :: phys_const
        character(len=*), intent(IN) :: domain
        character(len=*), intent(IN) :: grid_name
        type(phys_const_class), intent(IN),  optional :: cnst
        ! The record actually used, returned so that a caller can hand the same
        ! set to the other components of a coupled program without loading the
        ! parameter file again or knowing which group was selected.
        type(phys_const_class), intent(OUT), optional :: cnst_out

        ! Local variables
        logical :: init_pars
        character(len=56) :: group
        type(phys_const_class) :: cn
        character(len=512), parameter :: filename = "input/yelmo_phys_const.nml"

        ! Determine physical constants group name to use based on parameter choice
        select case(trim(phys_const))
            case("Earth")
                group = "Earth"
            case("EISMINT","EISMINT1","EISMINT2")
                group = "EISMINT"        
            case("MISMIP","MISMIP3D")
                group = "MISMIP3D"
            case("MISMIP+","MISMIPplus")
                group = "MISMIPplus"
            case("ISMIPHOM","ISMIP-HOM")
                group = "ISMIPHOM"
            case("CALVINGMIP","CalvingMIP")
                group = "CALVINGMIP"
            case("TROUGH")
                group = "TROUGH"
            case DEFAULT
                write(*,*) "ybound_define_physical_constants:: Error: yelmo.phys_const not recognized."
                write(*,*) "yelmo.phys_const = ", trim(phys_const)
                error stop 1
        end select

        ! Obtain the shared physical constants

        if (present(cnst)) then
            ! The driver owns them; use its set as-is.
            cn = cnst
        else
            ! Standalone: load them from Yelmo's own parameter file.
            call phys_const_load(cn,filename,group=group)
        end if

        call phys_const_require(cn,"ybound_define_physical_constants")

        ! Mirror the shared record into Yelmo's working precision. The derived
        ! factors are taken from the record rather than recomputed here, so that
        ! the ratios cannot drift from the densities they follow from.

        call phys_const_get(cn,"g",                  c%g)
        call phys_const_get(cn,"T0",                 c%T0)
        call phys_const_get(cn,"rho_ice",            c%rho_ice)
        call phys_const_get(cn,"rho_w",              c%rho_w)
        call phys_const_get(cn,"rho_sw",             c%rho_sw)
        call phys_const_get(cn,"rho_asth",           c%rho_a)
        call phys_const_get(cn,"L_ice",              c%L_ice)
        call phys_const_get(cn,"T_pmp_beta",         c%T_pmp_beta)
        call phys_const_get(cn,"area_seasurf",       c%area_seasurf)

        call phys_const_get(cn,"conv_we_ie",         c%conv_we_ie)
        call phys_const_get(cn,"conv_mmawe_maie",    c%conv_mmawe_maie)
        call phys_const_get(cn,"conv_m3_Gt",         c%conv_m3_Gt)
        call phys_const_get(cn,"conv_km3_Gt",        c%conv_km3_Gt)
        call phys_const_get(cn,"conv_millionkm3_Gt", c%conv_millionkm3_Gt)
        call phys_const_get(cn,"conv_km3_sle",       c%conv_km3_sle)

        ! Yelmo's own, not part of the shared record.

        ! Year length is a calendar choice, not a physical constant, so it stays
        ! with Yelmo. fesm-utils names the standard conventions (sec_year_365d,
        ! sec_year_tropical, ...) if this is ever to be tied to one of them.
        init_pars = .TRUE.
        call nml_read(filename,group,"sec_year",c%sec_year,init=init_pars)

        ! [mm d-1 w.e.] => [m a-1 i.e.]. Carries a day count, hence local: 365
        ! here, against smbpal's 360 and the CMIP forcing's 365.2422.
        c%conv_mmdwe_maie = 1e-3*365*c%conv_we_ie

        ! Hand the record back if asked.
        if (present(cnst_out)) cnst_out = cn

        if (yelmo_log) then
            write(*,*) ""
            write(*,*) "yelmo:: loaded physical constants for: ", trim(domain), " : ", trim(grid_name)
            call phys_const_log(cn)
            write(*,*) "  yelmo-specific:"
            write(*,*) "    sec_year           = ", c%sec_year
            write(*,*) "    conv_mmdwe_maie    = ", c%conv_mmdwe_maie
            write(*,*) ""
        end if

        return

    end subroutine ybound_define_physical_constants
    
    subroutine ybound_load_masks(bnd,nml_path,nml_group,domain,grid_name,grid,reg,mask_ice,masks)
        ! Masks of the boundary from the regions, zones and basins of FesmData
        ! (fesm-utils regions), made once here from selection expressions of
        ! the masks group (nml_group):
        !   regions_group     group of the regions ("None" = no regions)
        !   basins            basin ids of bnd%basins: "<set>" or "<set>.group" of the
        !                     regions, "None" (no basins) or "domain" (one basin, 1)
        !   mask_ice_dynamic  where ice is dynamic (MASK_ICE_DYNAMIC; else none)
        !   mask_ice_fixed    where ice is prescribed (MASK_ICE_FIXED; over dynamic)
        !   relax, relax_tau  where ice relaxes to H_ice_ref (bnd%tau_relax; topo_rel = -1)
        !   mask_rmse         where the error metrics are computed (bnd%mask_rmse)
        ! A driver may supply the regions (reg, on the Yelmo grid) instead of
        ! regions_group, and mask_ice instead of the mask_ice expressions.
        ! Without regions, the expressions can only be "all" or "none".
        ! masks: the named masks of the regions (for regional output).

        implicit none 

        type(ybound_class), intent(INOUT) :: bnd 
        character(len=*), intent(IN)      :: nml_path, nml_group
        character(len=*), intent(IN)      :: domain, grid_name 
        type(grid_class), intent(IN)      :: grid
        type(regions_class), intent(IN), optional :: reg
        integer,  intent(IN), optional    :: mask_ice(:,:)
        type(region_mask_class), allocatable, intent(OUT), optional :: masks(:)

        ! Local variables
        type(regions_class) :: reg_now
        logical             :: with_reg
        logical             :: extent(size(bnd%basins,1),size(bnd%basins,2))
        character(len=256)  :: regions_group
        character(len=56)   :: basins
        character(len=1000) :: expr_dynamic, expr_fixed, expr_relax, expr_rmse
        real(wp)            :: relax_tau

        character(len=*), parameter :: def_file  = "input/yelmo_defaults.nml"
        character(len=*), parameter :: def_masks = "yelmo_masks"

        call nml_validate(nml_path,def_file,nml_group,defaults_group=def_masks)

        call nml_read(nml_path,nml_group,"regions_group",   regions_group,   defaults_file=def_file,defaults_group=def_masks)
        call nml_read(nml_path,nml_group,"basins",          basins,          defaults_file=def_file,defaults_group=def_masks)
        call nml_read(nml_path,nml_group,"mask_ice_dynamic",expr_dynamic,defaults_file=def_file,defaults_group=def_masks)
        call nml_read(nml_path,nml_group,"mask_ice_fixed",  expr_fixed,  defaults_file=def_file,defaults_group=def_masks)
        call nml_read(nml_path,nml_group,"relax",           expr_relax,           defaults_file=def_file,defaults_group=def_masks)
        call nml_read(nml_path,nml_group,"relax_tau",       relax_tau,       defaults_file=def_file,defaults_group=def_masks)
        call nml_read(nml_path,nml_group,"mask_rmse",       expr_rmse,       defaults_file=def_file,defaults_group=def_masks)

        if (trim(expr_relax) .ne. "none" .and. relax_tau .le. 0.0_wp) then
            write(io_unit_err,*) "ybound_load_masks:: Error: relax needs relax_tau > 0."
            write(io_unit_err,*) "relax = ", trim(expr_relax), ", relax_tau = ", relax_tau
            error stop 1
        end if

        ! The regions: from the driver, from regions_group, or none
        with_reg = .TRUE.
        if (present(reg)) then
            reg_now = reg
        else if (trim(regions_group) .ne. "None") then
            call regions_init_nml(reg_now,nml_path,regions_group,domain=domain,grid_name=grid_name,grid=grid)
        else
            with_reg = .FALSE.
        end if

        ! Region codes (deepest level) and basins (ids of the set `basins`)
        bnd%regions    = 0.0_wp
        bnd%basins     = 0.0_wp
        bnd%basin_mask = 0.0_wp
        if (with_reg) then
            bnd%regions = real(reg_now%region_3,wp)

            bnd%basins = real(regions_basin_ids(reg_now,basins,extent=extent),wp)
            where (extent) bnd%basin_mask = 1.0_wp
        else
            select case(trim(basins))
                case("None")
                    ! No basins
                case("domain")
                    bnd%basins     = 1.0_wp
                    bnd%basin_mask = 1.0_wp
                case DEFAULT
                    write(io_unit_err,*) "ybound_load_masks:: Error: without regions, basins can only be None or domain."
                    write(io_unit_err,*) "basins = ", trim(basins)
                    error stop 1
            end select
        end if

        ! Where ice is allowed in the domain
        if (present(mask_ice)) then

            if (any(mask_ice .ne. MASK_ICE_NONE .and. mask_ice .ne. MASK_ICE_FIXED &
                                                .and. mask_ice .ne. MASK_ICE_DYNAMIC)) then
                write(io_unit_err,*) "ybound_load_masks:: Error: mask_ice values must be &
                                     &MASK_ICE_NONE, MASK_ICE_FIXED or MASK_ICE_DYNAMIC."
                write(io_unit_err,*) "range(mask_ice): ", minval(mask_ice), maxval(mask_ice)
                error stop 1
            end if

            bnd%mask_ice = mask_ice

        else

            bnd%mask_ice = MASK_ICE_NONE
            where (select_mask(expr_dynamic)) bnd%mask_ice = MASK_ICE_DYNAMIC
            where (select_mask(expr_fixed))   bnd%mask_ice = MASK_ICE_FIXED

        end if

        ! Relaxation timescale (used with ytopo.topo_rel = -1)
        bnd%tau_relax = -1.0_wp
        where (select_mask(expr_relax)) bnd%tau_relax = relax_tau

        ! Region of the error metrics
        bnd%mask_rmse = select_mask(expr_rmse)

        ! The named masks, for regional output
        if (present(masks)) then
            if (with_reg) then
                masks = reg_now%masks
            else
                allocate(masks(0))
            end if
        end if

        if (with_reg) call regions_end(reg_now)

        write(*,*) "ybound_load_masks:: range(basins):  ", minval(bnd%basins),  maxval(bnd%basins)
        write(*,*) "ybound_load_masks:: range(regions): ", minval(bnd%regions), maxval(bnd%regions)

        return 

    contains

        function select_mask(expr) result(mask)
            character(len=*), intent(IN) :: expr
            logical :: mask(size(bnd%mask_ice,1),size(bnd%mask_ice,2))

            if (with_reg) then
                mask = regions_select(reg_now,expr)
            else
                select case(trim(expr))
                    case("all")
                        mask = .TRUE.
                    case("none")
                        mask = .FALSE.
                    case DEFAULT
                        write(io_unit_err,*) "ybound_load_masks:: Error: without regions, &
                                             &a mask expression can only be all or none."
                        write(io_unit_err,*) "expression = ", trim(expr)
                        error stop 1
                end select
            end if

        end function select_mask

    end subroutine ybound_load_masks

    subroutine ybound_define_mask_ice(bnd,domain,boundaries,mask_border)
        ! Treatment of the domain border in mask_ice (ice dynamic, MASK_ICE_DYNAMIC,
        ! prescribed, MASK_ICE_FIXED, or forced to zero, MASK_ICE_NONE), after
        ! ybound_load_masks has set where ice is allowed in the domain.

        implicit none

        type(ybound_class), intent(INOUT) :: bnd
        character(len=*),   intent(IN)    :: domain
        character(len=*),   intent(IN)    :: boundaries     ! Topography boundary conditions
        character(len=*),   intent(IN)    :: mask_border    ! yelmo.mask_border

        ! Also set calv_mask false everywhere (no imposed calving front)
        bnd%calv_mask   = .FALSE.

        call define_mask_ice_border(bnd%mask_ice,domain,boundaries,mask_border)

        return

    end subroutine ybound_define_mask_ice

    subroutine define_mask_ice_border(mask_ice,domain,boundaries,mask_border)
        ! Treatment of the domain border (yelmo.mask_border):
        !   "auto":    by domain (below)
        !   "none":    no ice on the border
        !   "fixed":   ice thickness prescribed on the border (= H_ice_ref)
        !   "dynamic": the border is left as the domain mask defines it
        ! "none" and "fixed" skip periodic directions, where the border points
        ! are interior points.

        implicit none

        integer,          intent(INOUT) :: mask_ice(:,:)
        character(len=*), intent(IN)    :: domain
        character(len=*), intent(IN)    :: boundaries
        character(len=*), intent(IN)    :: mask_border

        ! Local variables
        integer :: nx, ny
        logical :: per_x, per_y

        nx = size(mask_ice,1)
        ny = size(mask_ice,2)

        call get_periodic_directions(per_x,per_y,boundary_code(boundaries))

        select case(trim(mask_border))

            case("auto")

                select case(trim(domain))

                    case ("North","Eurasia","Antarctica","EISMINT")
                        ! No ice on any border

                        call set_border(mask_ice,MASK_ICE_NONE,per_x=.FALSE.,per_y=.FALSE.)

                    case ("Greenland")
                        ! The border is outside of the allowed regions

                    case ("MISMIP","MISMIP3D","MISMIP+","TROUGH","TROUGH-F17")

                        ! Ice can grow everywhere, except farthest x-border.
                        !
                        ! "MISMIP3D" must be listed here explicitly. yelmo_init already
                        ! maps experiment="MISMIP3D" onto the MISMIP3D (y-periodic)
                        ! tpo/dyn/thrm boundaries, but this select case used to omit it,
                        ! so a domain named "MISMIP3D" fell through to case DEFAULT and
                        ! had all four borders marked MASK_ICE_FIXED. With bnd%H_ice_ref
                        ! left at its zero default, calc_G_boundaries then reset
                        ! H_ice = H_ice_ref = 0 on the whole perimeter every timestep --
                        ! silently draining ice at the flowband divide (i=1) and along
                        ! both lateral edges (j=1, j=ny), which no MISMIP-type setup wants.
                        mask_ice(nx,:) = MASK_ICE_NONE

                    case DEFAULT
                        ! Unknown domain: prescribed borders in non-periodic directions
                        ! (mask_ice can always be modified later)

                        call set_border(mask_ice,MASK_ICE_FIXED,per_x,per_y)

                end select

            case("none")

                call set_border(mask_ice,MASK_ICE_NONE,per_x,per_y)

            case("fixed")

                call set_border(mask_ice,MASK_ICE_FIXED,per_x,per_y)

            case("dynamic")

                ! Pass, border as defined by the domain mask

            case DEFAULT

                write(io_unit_err,*) "ybound_define_mask_ice:: Error: yelmo.mask_border not recognized."
                write(io_unit_err,*) "mask_border = ", trim(mask_border)
                write(io_unit_err,*) "Options: auto, none, fixed, dynamic."
                error stop 1

        end select

        return

    end subroutine define_mask_ice_border

    subroutine set_border(mask_ice,val,per_x,per_y)
        ! Set the border points of mask_ice to val, except in periodic directions.

        implicit none

        integer, intent(INOUT) :: mask_ice(:,:)
        integer, intent(IN)    :: val
        logical, intent(IN)    :: per_x, per_y

        ! Local variables
        integer :: nx, ny

        nx = size(mask_ice,1)
        ny = size(mask_ice,2)

        if (.not. per_x) then
            mask_ice(1,:)  = val
            mask_ice(nx,:) = val
        end if
        if (.not. per_y) then
            mask_ice(:,1)  = val
            mask_ice(:,ny) = val
        end if

        return

    end subroutine set_border

    subroutine ybound_update_rates(bnd,time)
        ! Rates of bedrock elevation and sea level since the previous call of
        ! yelmo_update (both are set by the driver between calls). Zero on the
        ! first call after initialisation. The previous-call state (z_bed_n,
        ! z_sl_n, time_n, rates_init) is in the restart, so a continued run (model
        ! time at initialisation = restart time) gets the same rates as a straight
        ! run; other restarts give zero rates on the first call.

        implicit none

        type(ybound_class), intent(INOUT) :: bnd
        real(dp),           intent(IN)    :: time

        if (bnd%rates_init .and. time .gt. bnd%time_n) then
            bnd%dz_bed_dt = (bnd%z_bed - bnd%z_bed_n) / real(time-bnd%time_n,wp)
            bnd%dz_sl_dt  = (bnd%z_sl  - bnd%z_sl_n)  / real(time-bnd%time_n,wp)
        else
            bnd%dz_bed_dt = 0.0_wp
            bnd%dz_sl_dt  = 0.0_wp
        end if

        bnd%z_bed_n    = bnd%z_bed
        bnd%z_sl_n     = bnd%z_sl
        bnd%time_n     = time
        bnd%rates_init = .TRUE.

        return

    end subroutine ybound_update_rates

    subroutine ybound_alloc(now,nx,ny)

        implicit none 

        type(ybound_class) :: now 
        integer :: nx, ny 

        call ybound_dealloc(now)

        allocate(now%z_bed(nx,ny))
        allocate(now%z_bed_sd(nx,ny))
        allocate(now%z_sl(nx,ny))
        allocate(now%H_sed(nx,ny))
        allocate(now%smb(nx,ny))
        allocate(now%T_srf(nx,ny))
        allocate(now%bmb_shlf(nx,ny))
        allocate(now%fmb_shlf(nx,ny))
        allocate(now%T_shlf(nx,ny))
        allocate(now%tf_shlf(nx,ny))
        allocate(now%Q_geo(nx,ny))
        allocate(now%Qd(nx,ny))

        allocate(now%enh_srf(nx,ny)) 

        allocate(now%basins(nx,ny))
        allocate(now%basin_mask(nx,ny))
        allocate(now%regions(nx,ny))
        allocate(now%mask_rmse(nx,ny))
        
        allocate(now%calv_mask(nx,ny))
        
        allocate(now%H_ice_ref(nx,ny))
        allocate(now%z_bed_ref(nx,ny))
        
        allocate(now%mask_ice(nx,ny))
        allocate(now%tau_relax(nx,ny))

        allocate(now%z_bed_corr(nx,ny))
        allocate(now%dzbdt_corr (nx,ny))

        allocate(now%z_bed_n(nx,ny))
        allocate(now%z_sl_n(nx,ny))
        allocate(now%dz_bed_dt(nx,ny))
        allocate(now%dz_sl_dt(nx,ny))
        
        now%z_bed       = 0.0_wp 
        now%z_bed_sd    = 0.0_wp
        now%z_sl        = 0.0_wp 
        now%H_sed       = 0.0_wp 
        now%smb         = 0.0_wp 
        now%T_srf       = 0.0_wp 
        now%bmb_shlf    = 0.0_wp 
        now%fmb_shlf    = 0.0_wp 
        now%T_shlf      = 0.0_wp
        now%tf_shlf     = 0.0_wp
        now%Q_geo       = 0.0_wp
        now%Qd          = 0.0_wp  

        now%enh_srf     = 1.0_wp 

        now%basins      = 0.0_wp 
        now%basin_mask  = 0.0_wp 
        now%regions     = 0.0_wp 
        now%mask_rmse   = .TRUE. 
        
        now%calv_mask   = .FALSE. ! By default no, no calving mask

        now%H_ice_ref   = 0.0_wp 
        now%z_bed_ref   = 0.0_wp

        now%mask_ice    = MASK_ICE_DYNAMIC   ! By default, ice is solved everywhere
        now%tau_relax   = -1.0_wp   ! By default, no relaxation anywhere

        now%z_bed_corr  = 0.0_wp
        now%dzbdt_corr  = 0.0_wp

        now%z_bed_n     = 0.0_wp
        now%z_sl_n      = 0.0_wp
        now%dz_bed_dt   = 0.0_wp
        now%dz_sl_dt    = 0.0_wp
        now%time_n      = 0.0_dp
        now%rates_init  = .FALSE.
        
        return 

    end subroutine ybound_alloc

    subroutine ybound_dealloc(now)

        implicit none 

        type(ybound_class) :: now

        if (allocated(now%z_bed))       deallocate(now%z_bed)
        if (allocated(now%z_bed_sd))    deallocate(now%z_bed_sd)
        if (allocated(now%z_sl))        deallocate(now%z_sl)
        if (allocated(now%H_sed))       deallocate(now%H_sed)
        if (allocated(now%smb))         deallocate(now%smb)
        if (allocated(now%T_srf))       deallocate(now%T_srf)
        if (allocated(now%bmb_shlf))    deallocate(now%bmb_shlf)
        if (allocated(now%fmb_shlf))    deallocate(now%fmb_shlf)
        if (allocated(now%T_shlf))      deallocate(now%T_shlf)
        if (allocated(now%tf_shlf))     deallocate(now%tf_shlf)
        if (allocated(now%Q_geo))       deallocate(now%Q_geo)
        if (allocated(now%Qd))          deallocate(now%Qd)

        if (allocated(now%enh_srf))     deallocate(now%enh_srf)
        
        if (allocated(now%basins))      deallocate(now%basins)
        if (allocated(now%basin_mask))  deallocate(now%basin_mask)
        if (allocated(now%regions))     deallocate(now%regions)
        if (allocated(now%mask_rmse))   deallocate(now%mask_rmse)
        
        if (allocated(now%calv_mask))   deallocate(now%calv_mask)
        
        if (allocated(now%H_ice_ref))   deallocate(now%H_ice_ref) 
        if (allocated(now%z_bed_ref))   deallocate(now%z_bed_ref)

        if (allocated(now%mask_ice))    deallocate(now%mask_ice)
        if (allocated(now%tau_relax))   deallocate(now%tau_relax)
        
        if (allocated(now%z_bed_corr))  deallocate(now%z_bed_corr)
        if (allocated(now%dzbdt_corr )) deallocate(now%dzbdt_corr )

        if (allocated(now%z_bed_n))     deallocate(now%z_bed_n)
        if (allocated(now%z_sl_n))      deallocate(now%z_sl_n)
        if (allocated(now%dz_bed_dt))   deallocate(now%dz_bed_dt)
        if (allocated(now%dz_sl_dt))    deallocate(now%dz_sl_dt)

        return 

    end subroutine ybound_dealloc

end module yelmo_boundaries
