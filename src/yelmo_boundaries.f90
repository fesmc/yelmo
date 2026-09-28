
module yelmo_boundaries

    use nml 
    use ncio 
    use yelmo_defs 
    use yelmo_tools, only : boundary_code, get_periodic_directions

    implicit none
    
    private
    public :: ybound_define_physical_constants
    public :: ybound_load_masks
    public :: ybound_define_mask_ice
    public :: ybound_alloc, ybound_dealloc
    
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
                stop
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
    
    subroutine ybound_load_masks(bnd,nml_path,nml_group,domain,grid_name)
        ! Load masks for managing regions and basins, etc. 

        implicit none 

        type(ybound_class), intent(INOUT) :: bnd 
        character(len=*), intent(IN)      :: nml_path, nml_group
        character(len=*), intent(IN)      :: domain, grid_name 

        ! Local variables
        logical            :: load_var
        character(len=512) :: filename
        character(len=56)  :: vnames(2)

        character(len=*), parameter :: def_file  = "input/yelmo_defaults.nml"
        character(len=*), parameter :: def_masks = "yelmo_masks"

        call nml_validate(nml_path,def_file,nml_group,defaults_group=def_masks)

        ! ====================================
        !
        ! basins
        !
        ! ====================================

        ! Specify default values
        bnd%basin_mask = 1.0
        bnd%basins     = 1.0

        call nml_read(nml_path,nml_group,"basins_load",load_var,defaults_file=def_file,defaults_group=def_masks)

        if (load_var) then

            call nml_read(nml_path,nml_group, "basins_path",filename,defaults_file=def_file,defaults_group=def_masks)
            call yelmo_parse_path(filename,domain,grid_name)
            call yelmo_check_file(nml_group,"basins_path",filename)

            call nml_read(nml_path,nml_group,"basins_nms",vnames,defaults_file=def_file,defaults_group=def_masks)
            ! Load basin information from a file
            call nc_read(filename,vnames(1),bnd%basins)

            if (trim(vnames(2)) .ne. "None") then
                ! If basins have been extrapolated, also load the original basin extent mask
                call nc_read(filename,vnames(2),bnd%basin_mask)
            end if

        end if

        ! ====================================
        !
        ! Read in the regions
        !
        ! ====================================

        ! First assign default region values (used if no regions file is loaded)
        bnd%region_mask = 1.0
        select case(trim(domain))
            case("Greenland")
                ! The Greenland domain is centered on the Greenland region
                bnd%regions = bnd%index_grl
            case DEFAULT
                ! No region information: mark all points as unclassified (0).
                ! Note: the hemisphere codes (North=1.0, Antarctica=2.0) cannot be
                ! used as defaults, since in the REGIONS files they denote points
                ! outside of any land subregion (open ocean), which is where
                ! ybound_define_mask_ice (and user code) forbids ice.
                bnd%regions = 0.0

        end select

        call nml_read(nml_path,nml_group,"regions_load",load_var,defaults_file=def_file,defaults_group=def_masks)

        if (load_var) then

            call nml_read(nml_path,nml_group, "regions_path",filename,defaults_file=def_file,defaults_group=def_masks)
            call yelmo_parse_path(filename,domain,grid_name)
            call yelmo_check_file(nml_group,"regions_path",filename)

            ! Load region information from a file
            call nml_read(nml_path,nml_group,"regions_nms",vnames,defaults_file=def_file,defaults_group=def_masks)
            call nc_read(filename,vnames(1),bnd%regions)

            if (trim(vnames(2)) .ne. "None") then
                ! If regions have been extrapolated, also load the original region extent mask
                call nc_read(filename,vnames(2),bnd%region_mask)
            end if

        end if 
        
        write(*,*) "ybound_load_masks:: range(basins):  ", minval(bnd%basins),  maxval(bnd%basins)
        write(*,*) "ybound_load_masks:: range(regions): ", minval(bnd%regions), maxval(bnd%regions)

        return 

    end subroutine ybound_load_masks

    subroutine ybound_define_mask_ice(bnd,domain,boundaries)
        ! Update mask defining where ice is dynamic (MASK_ICE_DYNAMIC),
        ! prescribed (MASK_ICE_FIXED), or forced to zero (MASK_ICE_NONE).

        implicit none

        type(ybound_class), intent(INOUT) :: bnd
        character(len=*),   intent(IN)    :: domain
        character(len=*),   intent(IN)    :: boundaries     ! Topography boundary conditions

        ! Local variables
        integer :: i, nx, ny
        logical :: per_x, per_y

        nx = size(bnd%mask_ice,1)
        ny = size(bnd%mask_ice,2)

        ! Initially mark all points as dynamic (ice is solved)
        bnd%mask_ice = MASK_ICE_DYNAMIC

        ! Also set calv_mask false everywhere (no imposed calving front)
        bnd%calv_mask   = .FALSE.

        ! Determine allowed regions based on domain
        select case(trim(domain))

            case ("North")
                ! Allow ice everywhere except the open ocean (region 1.0 in the
                ! REGIONS file; without a file, regions=0 and ice is allowed everywhere)

                where (bnd%regions .eq. 1.0) bnd%mask_ice = MASK_ICE_NONE
                bnd%mask_ice(1,:)  = MASK_ICE_NONE
                bnd%mask_ice(nx,:) = MASK_ICE_NONE
                bnd%mask_ice(:,1)  = MASK_ICE_NONE
                bnd%mask_ice(:,ny) = MASK_ICE_NONE

            case ("Eurasia")
                ! Allow ice only in the Eurasia domain (1.2*)

                if (count(bnd%regions .ge. 1.2 .and. bnd%regions .le. 1.29) .eq. 0) then
                    ! Without a regions file (regions=0), ice would be forbidden everywhere
                    write(io_unit_err,*) "ybound_define_mask_ice:: Error: domain='Eurasia' requires a regions &
                                         &field with Eurasia codes (1.2 <= regions <= 1.29), but none were found."
                    write(io_unit_err,*) "range(regions): ", minval(bnd%regions), maxval(bnd%regions)
                    stop
                end if

                where (bnd%regions .lt. 1.2 .or. bnd%regions .gt. 1.29) bnd%mask_ice = MASK_ICE_NONE
                bnd%mask_ice(1,:)  = MASK_ICE_NONE
                bnd%mask_ice(nx,:) = MASK_ICE_NONE
                bnd%mask_ice(:,1)  = MASK_ICE_NONE
                bnd%mask_ice(:,ny) = MASK_ICE_NONE

            case ("Greenland")

                bnd%mask_ice = MASK_ICE_NONE
                where (bnd%regions .eq. 1.3)  bnd%mask_ice = MASK_ICE_DYNAMIC   ! Main Greenland region
                where (bnd%regions .eq. 1.11) bnd%mask_ice = MASK_ICE_DYNAMIC   ! Ellesmere Island
                where (bnd%regions .eq. 1.0)  bnd%mask_ice = MASK_ICE_DYNAMIC   ! Open ocean (included some connections between 1.3 and 1.11)

            case ("Antarctica")
                ! Allow ice everywhere except the open ocean (region 2.0 in the
                ! REGIONS file; without a file, regions=0 and ice is allowed everywhere)

                where (bnd%regions .eq. 2.0) bnd%mask_ice = MASK_ICE_NONE
                bnd%mask_ice(1,:)  = MASK_ICE_NONE
                bnd%mask_ice(nx,:) = MASK_ICE_NONE
                bnd%mask_ice(:,1)  = MASK_ICE_NONE
                bnd%mask_ice(:,ny) = MASK_ICE_NONE


            case ("EISMINT")

                ! Ice can grow everywhere, except borders
                bnd%mask_ice       = MASK_ICE_DYNAMIC
                bnd%mask_ice(1,:)  = MASK_ICE_NONE
                bnd%mask_ice(nx,:) = MASK_ICE_NONE
                bnd%mask_ice(:,1)  = MASK_ICE_NONE
                bnd%mask_ice(:,ny) = MASK_ICE_NONE

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
                bnd%mask_ice       = MASK_ICE_DYNAMIC
                bnd%mask_ice(nx,:) = MASK_ICE_NONE

            case DEFAULT
                ! Unknown domain: dynamic interior, prescribed borders
                ! in non-periodic directions (in a periodic direction the
                ! border points are interior points)
                ! (mask_ice can always be modified later)

                call get_periodic_directions(per_x,per_y,boundary_code(boundaries))

                bnd%mask_ice       = MASK_ICE_DYNAMIC
                if (.not. per_x) then
                    bnd%mask_ice(1,:)  = MASK_ICE_FIXED
                    bnd%mask_ice(nx,:) = MASK_ICE_FIXED
                end if
                if (.not. per_y) then
                    bnd%mask_ice(:,1)  = MASK_ICE_FIXED
                    bnd%mask_ice(:,ny) = MASK_ICE_FIXED
                end if

        end select

        return

    end subroutine ybound_define_mask_ice

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
        allocate(now%region_mask(nx,ny))
        
        allocate(now%calv_mask(nx,ny))
        
        allocate(now%H_ice_ref(nx,ny))
        allocate(now%z_bed_ref(nx,ny))
        
        allocate(now%mask_ice(nx,ny))
        allocate(now%tau_relax(nx,ny))

        allocate(now%z_bed_corr(nx,ny))
        allocate(now%dzbdt_corr (nx,ny))
        
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
        now%region_mask = 0.0_wp 
        
        now%calv_mask   = .FALSE. ! By default no, no calving mask

        now%H_ice_ref   = 0.0_wp 
        now%z_bed_ref   = 0.0_wp

        now%mask_ice    = MASK_ICE_DYNAMIC   ! By default, ice is solved everywhere
        now%tau_relax   = -1.0_wp   ! By default, no relaxation anywhere

        now%z_bed_corr  = 0.0_wp
        now%dzbdt_corr  = 0.0_wp
        
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
        if (allocated(now%region_mask)) deallocate(now%region_mask)
        
        if (allocated(now%calv_mask))   deallocate(now%calv_mask)
        
        if (allocated(now%H_ice_ref))   deallocate(now%H_ice_ref) 
        if (allocated(now%z_bed_ref))   deallocate(now%z_bed_ref)

        if (allocated(now%mask_ice))    deallocate(now%mask_ice)
        if (allocated(now%tau_relax))   deallocate(now%tau_relax)
        
        if (allocated(now%z_bed_corr))  deallocate(now%z_bed_corr)
        if (allocated(now%dzbdt_corr )) deallocate(now%dzbdt_corr )

        return 

    end subroutine ybound_dealloc

end module yelmo_boundaries
