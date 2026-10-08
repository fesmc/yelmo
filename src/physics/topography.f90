module topography 

    use yelmo_defs, only : wp, dp, io_unit_err, pi, TOL, is_equal, &
                           MASK_FRNT_ICE_FREE, MASK_FRNT_ICE_FREE_LAND, MASK_FRNT_NONE, MASK_FRNT_FLOAT, &
                           MASK_FRNT_MARINE, MASK_FRNT_GRND
    use yelmo_tools, only : boundary_code, get_neighbor_indices_bc_codes, get_periodic_directions
    use subgrid, only : calc_subgrid_array_quad

    implicit none 

    ! Key for matching bed types given by mask_bed 
    integer, parameter :: mask_bed_ocean   = 0 
    integer, parameter :: mask_bed_land    = 1
    integer, parameter :: mask_bed_frozen  = 2
    integer, parameter :: mask_bed_stream  = 3
    integer, parameter :: mask_bed_grline  = 4
    integer, parameter :: mask_bed_float   = 5
    integer, parameter :: mask_bed_island  = 6
    integer, parameter :: mask_bed_partial = 7

    ! CISM limits of the effective surface at marine-grounded fronts
    ! ("AIS testing showed that values of 25 m and 0.001 prevent large
    ! ice speeds that can lead to instability")
    ! Ice thinner than H_ice_eps is ice free for the area fraction (f_ice = 0)
    ! and the front classification. A cell holding round-off amounts of ice
    ! would otherwise become a front cell with the full reference thickness
    ! H_eff >= front_H_eff_min, so that the front force would jump by one cell
    ! on round-off (CISM uses thck > eps11 for the same masks, in double
    ! precision).
    real(wp), parameter :: H_ice_eps     = 1e-3_wp  ! [m]
    real(wp), parameter :: dz_srf_max    = 25.0_wp  ! [m]   z_srf_eff - z_srf
    real(wp), parameter :: dz_srf_dx_max = 0.001_wp ! [m/m] upward surface slope at the front

    private  

    public :: gen_mask_bed
    public :: calc_column_kinematic_rates

    public :: calc_ice_fraction
    public :: calc_front_H_ref
    public :: calc_lsf_area_fraction
    public :: calc_front_cells
    public :: calc_ice_front

    public :: calc_z_srf_max
    public :: calc_H_eff
    public :: calc_H_grnd
    public :: calc_H_af
    public :: calc_f_grnd_subgrid_linear
    public :: calc_f_grnd_pinning_points
    public :: remove_englacial_lakes
    public :: calc_distance_to_ice_margin
    public :: calc_distance_to_grounding_line
    public :: calc_grounding_line_zone
    public :: calc_bmb_total
    public :: calc_fmb_total
    public :: calc_melt_rate_rignot16
    
    ! ajr: these routines are slow, do not use...
    !public :: distance_to_grline
    !public :: distance_to_margin
    
    public :: determine_grounded_fractions

    ! Integers
    public :: mask_bed_ocean  
    public :: mask_bed_land  
    public :: mask_bed_frozen
    public :: mask_bed_stream
    public :: mask_bed_grline
    public :: mask_bed_float 
    public :: mask_bed_island
    
contains 

    subroutine calc_column_kinematic_rates(dzsdt_kin,dzbdt_kin,mask_kin,dHidt_vert,advanced,f_grnd,f_ice, &
                                            H_ice,H_ice_dyn,H_ice_n,H_ice_dyn_n,dz_bed_dt,dz_sl_dt,rho_ice,rho_sw)
        ! Kinematic rates of the surface and base of the ice column, the
        ! boundary conditions of the vertical velocity. Grounded ice: the base
        ! follows the bedrock; floating ice: the column floats at sea level
        ! (z_b = z_sl - rho_ice/rho_sw*H). Blended by the grounded fraction.
        ! dHidt_vert holds only vertical changes of the column (not calving,
        ! front advance or removals) and follows the actual thickness, so it
        ! is only the column's rate where the column is the actual ice
        ! (H_ice_dyn == H_ice) at the start and the end of the step. Elsewhere
        ! (partial front cells, cells on the H_eff floor, newly ice-covered
        ! cells) the column is re-derived: zero rates.
        ! mask_kin marks the columns whose rate is given by the thickness step
        ! applied this step (valid column and ice advanced), where the vertical
        ! velocity can be closed against the applied thickness tendency.

        implicit none

        real(wp), intent(OUT) :: dzsdt_kin(:,:)     ! [m/a]
        real(wp), intent(OUT) :: dzbdt_kin(:,:)     ! [m/a]
        integer,  intent(OUT) :: mask_kin(:,:)
        real(wp), intent(IN)  :: dHidt_vert(:,:)    ! [m/a]
        logical,  intent(IN)  :: advanced           ! Ice thickness advanced this step
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)         ! [m] Thickness at the end of the step
        real(wp), intent(IN)  :: H_ice_dyn(:,:)     ! [m] Active column thickness at the end of the step
        real(wp), intent(IN)  :: H_ice_n(:,:)       ! [m] Thickness at the start of the step
        real(wp), intent(IN)  :: H_ice_dyn_n(:,:)   ! [m] Active column thickness at the start of the step
        real(wp), intent(IN)  :: dz_bed_dt(:,:)     ! [m/a]
        real(wp), intent(IN)  :: dz_sl_dt(:,:)      ! [m/a]
        real(wp), intent(IN)  :: rho_ice
        real(wp), intent(IN)  :: rho_sw

        integer  :: i, j, nx, ny
        real(wp) :: rho_frac, dzb_flt

        nx = size(dzsdt_kin,1)
        ny = size(dzsdt_kin,2)

        rho_frac = rho_ice/rho_sw

        !$omp parallel do collapse(2) private(i,j,dzb_flt)
        do j = 1, ny
        do i = 1, nx
            if (f_ice(i,j) .eq. 1.0_wp .and. H_ice_n(i,j) .gt. 0.0_wp .and. &
                H_ice_dyn(i,j) .eq. H_ice(i,j) .and. H_ice_dyn_n(i,j) .eq. H_ice_n(i,j)) then
                dzb_flt        = dz_sl_dt(i,j) - rho_frac*dHidt_vert(i,j)
                dzbdt_kin(i,j) = f_grnd(i,j)*dz_bed_dt(i,j) + (1.0_wp-f_grnd(i,j))*dzb_flt
                dzsdt_kin(i,j) = dzbdt_kin(i,j) + dHidt_vert(i,j)
                mask_kin(i,j)  = merge(1,0,advanced)
            else
                dzbdt_kin(i,j) = 0.0_wp
                dzsdt_kin(i,j) = 0.0_wp
                mask_kin(i,j)  = 0
            end if
        end do
        end do
        !$omp end parallel do

        return

    end subroutine calc_column_kinematic_rates


    elemental subroutine gen_mask_bed(mask,f_ice,f_pmp,f_grnd,mask_grline)
        ! Generate an output mask for model conditions at bed
        ! based on input masks 
        ! 0: ocean, 1: land, 2: sia, 3: streams, grline: 4, floating: 5, islands: 6
        ! 7: partially-covered ice cell.

        implicit none 

        integer,  intent(OUT) :: mask 
        real(wp), intent(IN)  :: f_ice, f_pmp, f_grnd
        logical,  intent(IN)  :: mask_grline

        if (mask_grline) then
            ! Grounding line

            mask = mask_bed_grline

        else if ( f_ice .eq. 0.0_wp ) then 
            ! Ice-free points 

            if (f_grnd .gt. 0.0) then
                ! Ice-free land

                mask = mask_bed_land

            else
                ! Ice-free ocean

                mask = mask_bed_ocean

            end if 

        else if (f_ice .gt. 0.0 .and. f_ice .lt. 1.0) then 
            ! Partially ice-covered points 

            mask = mask_bed_partial

        else
            ! Fully ice-covered points 

            if (f_grnd .gt. 0.0) then
                ! Grounded ice-covered points 

                if (f_pmp .gt. 0.5) then 
                    ! Temperate points

                    mask = mask_bed_stream 

                else
                    ! Frozen points 

                    mask = mask_bed_frozen 

                end if 

            else
                ! Floating ice-covered points 

                mask = mask_bed_float

            end if 

        end if 

        return 

    end subroutine gen_mask_bed

    subroutine calc_ice_fraction(f_ice,H_eff,H_ice,z_bed,z_sl,rho_ice,rho_sw, &
                                    front_subgrid,H_eff_min,dHdx,dx,boundaries,a_lsf)
        ! Ice area fraction f_ice and effective thickness H_eff of each cell,
        ! following the CISM subgrid calving-front scheme (which_ho_calving_front).
        !
        ! front_subgrid = "none":     f_ice binary, H_eff = H_ice.
        ! front_subgrid = "floating": floating cells can be partial front cells.
        ! front_subgrid = "marine":   floating and marine-grounded cells can be.
        !
        ! A front cell is an eligible ice cell with an ice-free ocean edge
        ! neighbour (with the level set a_lsf, also one cut by the front that
        ! touches the ocean at a corner; calc_front_cells). Its H_eff is the
        ! reference thickness from its interior (eligible, not front)
        ! neighbours (calc_front_H_ref); for "marine" the
        ! effective surface is also at most dz_srf_max above the actual
        ! surface. Front cells without an interior neighbour keep H_eff = H_ice.
        ! With the level set, front cells entirely behind the front (a_lsf = 1)
        ! take no reference either: the level set sets their area fraction,
        ! so they are full (H_eff = H_ice, above H_eff_min).
        ! H_eff >= H_eff_min in all eligible cells, and <= flotation in
        ! floating front cells ("floating"). f_ice = min(H_ice/H_eff,1) in front
        ! cells, 1 in other ice cells (H_ice > H_ice_eps), 0 elsewhere.

        implicit none 

        real(wp), intent(OUT) :: f_ice(:,:)             ! [--] Ice covered fraction (aa-nodes)
        real(wp), intent(OUT) :: H_eff(:,:)             ! [m]  Effective ice thickness (aa-nodes)
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m]  Ice thickness (aa-nodes)
        real(wp), intent(IN)  :: z_bed(:,:)             ! [m]  Bedrock elevation
        real(wp), intent(IN)  :: z_sl(:,:)              ! [m]  Sea-level elevation
        real(wp), intent(IN)  :: rho_ice
        real(wp), intent(IN)  :: rho_sw
        character(len=*), intent(IN) :: front_subgrid   ! "none", "floating" or "marine"
        real(wp), intent(IN)  :: H_eff_min              ! [m]  Minimum H_eff of eligible cells
        real(wp), intent(IN)  :: dHdx                   ! [m/m] Thickness gradient assumed at a full front
        real(wp), intent(IN)  :: dx                     ! [m]  Grid resolution
        character(len=*), intent(IN) :: boundaries
        real(wp), optional, intent(IN) :: a_lsf(:,:)    ! [--] Level-set area fraction (LSF calving)

        ! Local variables 
        integer  :: i, j, nx, ny
        real(wp) :: z_srf_eff, z_srf_max
        logical  :: is_none, is_flt, is_mar             ! front_subgrid choice
        logical, allocatable  :: mask_elig(:,:)         ! Eligible (marine) ice cells
        logical, allocatable  :: mask_cf(:,:)           ! Front cells (eligible, ocean edge neighbour)
        logical, allocatable  :: mask_ocn(:,:)          ! Ice-free ocean cells
        logical, allocatable  :: has_ref(:,:)           ! Front cells with an interior neighbour
        logical, allocatable  :: mask_ref(:,:)          ! Front cells that take a reference
        real(wp), allocatable :: H_ref(:,:)             ! [m] Reference thickness of front cells
        real(wp), allocatable :: H_flot(:,:)            ! [m] Flotation thickness

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        is_none = trim(front_subgrid) .eq. "none"
        is_flt  = trim(front_subgrid) .eq. "floating"
        is_mar  = trim(front_subgrid) .eq. "marine"

        if (.not. is_none) allocate(H_flot(nx,ny))

        ! Binary defaults (and flotation thickness)
        !$omp parallel do collapse(2) private(i,j)
        do j = 1, ny
        do i = 1, nx
            if (H_ice(i,j) .gt. H_ice_eps) then
                f_ice(i,j) = 1.0_wp
                H_eff(i,j) = H_ice(i,j)
            else
                f_ice(i,j) = 0.0_wp
                H_eff(i,j) = 0.0_wp
            end if
            if (.not. is_none) H_flot(i,j) = max( (z_sl(i,j)-z_bed(i,j))*rho_sw/rho_ice, 0.0_wp )
        end do
        end do
        !$omp end parallel do

        if (is_none) return

        allocate(mask_elig(nx,ny))
        allocate(mask_cf(nx,ny))
        allocate(mask_ocn(nx,ny))
        allocate(has_ref(nx,ny))
        allocate(H_ref(nx,ny))
        allocate(mask_ref(nx,ny))

        call calc_front_cells(mask_cf,mask_elig,mask_ocn,H_ice,z_bed,z_sl,rho_ice,rho_sw,front_subgrid,boundaries,a_lsf)

        ! Front cells entirely behind the level-set front take no reference
        mask_ref = mask_cf
        if (present(a_lsf)) mask_ref = mask_cf .and. a_lsf .lt. 1.0_wp

        call calc_front_H_ref(H_ref,has_ref,H_ice,z_bed,z_sl,rho_ice,rho_sw,mask_ref, &
                              mask_elig .and. .not. mask_cf,front_subgrid,dHdx,dx,boundaries)

        !$omp parallel do collapse(2) private(i,j,z_srf_eff,z_srf_max)
        do j = 1, ny
        do i = 1, nx

            if (has_ref(i,j)) then

                H_eff(i,j) = H_ref(i,j)

                if (is_mar) then
                    ! Limit the effective surface to dz_srf_max above the actual surface
                    z_srf_eff = srf_elev(H_eff(i,j),z_bed(i,j),z_sl(i,j),rho_ice,rho_sw)
                    z_srf_max = srf_elev(H_ice(i,j),z_bed(i,j),z_sl(i,j),rho_ice,rho_sw) + dz_srf_max
                    if (z_srf_eff .gt. z_srf_max) &
                        H_eff(i,j) = srf_thickness(z_srf_max,z_bed(i,j),z_sl(i,j),rho_ice,rho_sw)
                end if

            end if

            ! Lower limit in all eligible cells (most fronts are at least a few tens of metres thick)
            if (mask_elig(i,j)) H_eff(i,j) = max(H_eff(i,j), H_eff_min)

            ! Floating fronts do not exceed flotation (allows H_eff < H_eff_min in shallow water)
            if (is_flt .and. mask_cf(i,j)) H_eff(i,j) = min(H_eff(i,j), H_flot(i,j))

            ! Area fraction of front cells
            if (mask_cf(i,j) .and. H_eff(i,j) .gt. 0.0_wp) f_ice(i,j) = max( min(H_ice(i,j)/H_eff(i,j), 1.0_wp), TOL )

        end do
        end do
        !$omp end parallel do

        return 

    end subroutine calc_ice_fraction

    subroutine calc_front_H_ref(H_ref,has_ref,H_ice,z_bed,z_sl,rho_ice,rho_sw,mask_tgt,mask_int, &
                                    front_subgrid,dHdx,dx,boundaries)
        ! Reference (full-column) thickness of subgrid front cells, independent
        ! of the cell's own thickness: the thickest interior edge neighbour, or
        ! diagonal neighbour if there is none, minus dHdx*distance. Floating
        ! neighbours are capped at their flotation thickness ("floating"); for
        ! "marine" the effective surface rises at most dz_srf_dx_max*distance
        ! above the neighbour's surface. has_ref is false for target cells
        ! without an interior neighbour and outside mask_tgt (H_ref = 0).

        implicit none

        real(wp), intent(OUT) :: H_ref(:,:)             ! [m]  Reference thickness
        logical,  intent(OUT) :: has_ref(:,:)           ! Reference found
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m]  Ice thickness
        real(wp), intent(IN)  :: z_bed(:,:)             ! [m]  Bedrock elevation
        real(wp), intent(IN)  :: z_sl(:,:)              ! [m]  Sea-level elevation
        real(wp), intent(IN)  :: rho_ice
        real(wp), intent(IN)  :: rho_sw
        logical,  intent(IN)  :: mask_tgt(:,:)          ! Cells that need a reference
        logical,  intent(IN)  :: mask_int(:,:)          ! Interior cells (references)
        character(len=*), intent(IN) :: front_subgrid   ! "floating" or "marine"
        real(wp), intent(IN)  :: dHdx                   ! [m/m] Thickness gradient assumed at a full front
        real(wp), intent(IN)  :: dx                     ! [m]  Grid resolution
        character(len=*), intent(IN) :: boundaries

        ! Local variables
        integer  :: i, j, k, nx, ny, BC
        integer  :: im1, ip1, jm1, jp1
        integer  :: in(8), jn(8)
        integer  :: k_max
        real(wp) :: H_nb, H_max, dist
        real(wp) :: z_srf_eff, z_srf_max, z_srf_nb
        logical  :: is_flt, is_mar

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        BC = boundary_code(boundaries)

        is_flt = trim(front_subgrid) .eq. "floating"
        is_mar = trim(front_subgrid) .eq. "marine"

        !$omp parallel do collapse(2) private(i,j,k,im1,ip1,jm1,jp1,in,jn,k_max,H_nb,H_max,dist) &
        !$omp& private(z_srf_eff,z_srf_max,z_srf_nb)
        do j = 1, ny
        do i = 1, nx

            H_ref(i,j)   = 0.0_wp
            has_ref(i,j) = .FALSE.

            if (.not. mask_tgt(i,j)) cycle

            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! Edge neighbours (1-4), then diagonal neighbours (5-8)
            in = [im1,ip1,i,i,  im1,ip1,im1,ip1]
            jn = [j,j,jm1,jp1,  jm1,jm1,jp1,jp1]

            ! Thickest interior edge neighbour, else thickest interior diagonal
            k_max = 0
            H_max = 0.0_wp
            do k = 1, 8
                if (k .eq. 5 .and. k_max .gt. 0) exit
                if (.not. mask_int(in(k),jn(k))) cycle
                H_nb = H_ice(in(k),jn(k))
                if (is_flt) H_nb = min(H_nb, max( (z_sl(in(k),jn(k))-z_bed(in(k),jn(k)))*rho_sw/rho_ice, 0.0_wp ))
                if (H_nb .gt. H_max) then
                    H_max = H_nb
                    k_max = k
                end if
            end do

            if (k_max .eq. 0) cycle

            dist = dx
            if (k_max .gt. 4) dist = sqrt(2.0_wp)*dx

            H_ref(i,j)   = H_max - dHdx*dist
            has_ref(i,j) = .TRUE.

            if (is_mar) then
                ! Limit the upward slope of the effective surface from the neighbour
                z_srf_nb  = srf_elev(H_ice(in(k_max),jn(k_max)),z_bed(in(k_max),jn(k_max)),z_sl(in(k_max),jn(k_max)), &
                                     rho_ice,rho_sw)
                z_srf_eff = srf_elev(H_ref(i,j),z_bed(i,j),z_sl(i,j),rho_ice,rho_sw)
                z_srf_max = z_srf_nb + dz_srf_dx_max*dist
                if (z_srf_eff .gt. z_srf_max) H_ref(i,j) = srf_thickness(z_srf_max,z_bed(i,j),z_sl(i,j),rho_ice,rho_sw)
            end if

        end do
        end do
        !$omp end parallel do

        return

    end subroutine calc_front_H_ref

    pure function srf_elev(H,zb,zsl,rho_ice,rho_sw) result(zs)
        ! Surface elevation of a column: floating or grounded
        real(wp), intent(IN) :: H, zb, zsl, rho_ice, rho_sw
        real(wp) :: zs
        if (zb - zsl .lt. -rho_ice/rho_sw*H) then
            zs = zsl + (1.0_wp - rho_ice/rho_sw)*H
        else
            zs = zb + H
        end if
    end function srf_elev

    pure function srf_thickness(zs,zb,zsl,rho_ice,rho_sw) result(H)
        ! Thickness of a column with surface zs: floating, or grounded if thinner
        real(wp), intent(IN) :: zs, zb, zsl, rho_ice, rho_sw
        real(wp) :: H
        H = min( (zs - zsl)*rho_sw/(rho_sw-rho_ice), zs - zb )
    end function srf_thickness

    subroutine calc_front_cells(mask_cf,mask_elig,mask_ocn,H_ice,z_bed,z_sl,rho_ice,rho_sw,front_subgrid,boundaries,a_lsf)
        ! Front cells of the subgrid front scheme (ytopo.front_subgrid):
        ! eligible ice cells (floating, or floating and marine-grounded)
        ! with at least one ice-free ocean edge neighbour. With the level-set
        ! area fraction a_lsf (LSF calving), eligible cells not entirely behind
        ! the front (a_lsf < 1) that touch the ocean only at a corner are front
        ! cells too. This is the one front/interior classification of the
        ! level-set trim (calc_G_lsf_front) and of calc_ice_fraction. With
        ! "none" no cell is eligible.

        implicit none

        logical,  intent(OUT) :: mask_cf(:,:)           ! Front cells
        logical,  intent(OUT) :: mask_elig(:,:)         ! Eligible ice cells
        logical,  intent(OUT) :: mask_ocn(:,:)          ! Ice-free ocean cells
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: z_bed(:,:)
        real(wp), intent(IN)  :: z_sl(:,:)
        real(wp), intent(IN)  :: rho_ice
        real(wp), intent(IN)  :: rho_sw
        character(len=*), intent(IN) :: front_subgrid   ! "none", "floating" or "marine"
        character(len=*), intent(IN) :: boundaries
        real(wp), optional, intent(IN) :: a_lsf(:,:)    ! [--] Area fraction behind the level-set front

        ! Local variables
        integer :: i, j, nx, ny, im1, ip1, jm1, jp1, BC
        integer :: elig                                 ! 0: none, 1: floating, 2: marine

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        BC = boundary_code(boundaries)

        select case(trim(front_subgrid))
            case("none")
                elig = 0
            case("floating")
                elig = 1
            case("marine")
                elig = 2
            case DEFAULT
                write(io_unit_err,*) "calc_front_cells:: Error: front_subgrid not recognized: ", trim(front_subgrid)
                error stop 1
        end select

        !$omp parallel do collapse(2) private(i,j)
        do j = 1, ny
        do i = 1, nx
            mask_ocn(i,j) = H_ice(i,j) .le. H_ice_eps .and. z_bed(i,j) .lt. z_sl(i,j)
            select case(elig)
                case(1)
                    mask_elig(i,j) = H_ice(i,j) .gt. H_ice_eps .and. H_ice(i,j)*rho_ice/rho_sw .le. (z_sl(i,j)-z_bed(i,j))
                case(2)
                    mask_elig(i,j) = H_ice(i,j) .gt. H_ice_eps .and. z_bed(i,j) .lt. z_sl(i,j)
                case DEFAULT
                    mask_elig(i,j) = .FALSE.
            end select
        end do
        end do
        !$omp end parallel do

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
        do j = 1, ny
        do i = 1, nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            mask_cf(i,j) = mask_elig(i,j) .and. &
                ( mask_ocn(im1,j) .or. mask_ocn(ip1,j) .or. mask_ocn(i,jm1) .or. mask_ocn(i,jp1) )
        end do
        end do
        !$omp end parallel do

        if (present(a_lsf)) then
            ! Corner cells cut by the level-set front
            !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
            do j = 1, ny
            do i = 1, nx
                if (mask_cf(i,j) .or. .not. mask_elig(i,j) .or. a_lsf(i,j) .ge. 1.0_wp) cycle
                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
                mask_cf(i,j) = mask_ocn(im1,jm1) .or. mask_ocn(ip1,jm1) .or. &
                               mask_ocn(im1,jp1) .or. mask_ocn(ip1,jp1)
            end do
            end do
            !$omp end parallel do
        end if

        return

    end subroutine calc_front_cells

    subroutine calc_lsf_area_fraction(a_lsf,lsf,H_ice,z_bed,z_sl,boundaries)
        ! Area fraction of each cell behind the level-set front (lsf < 0).
        ! The cell is split into four quadrants with vertices at the cell
        ! centre, two edge midpoints and a corner (edge midpoints: mean of two
        ! aa-nodes, corners: mean of four), and the negative-LSF area of each
        ! quadrant is computed via marching squares. Including the centre value
        ! keeps the fraction continuous for 1-cell-wide features, where all
        ! four corner means are >= 0.
        ! Ice-free land neighbours take the value of the cell itself: their
        ! lsf is pinned to -1 to keep the level set off land, which says
        ! nothing about the front position in the cell. Counted as ice side,
        ! they gave an ocean cell (lsf > 0) between land cells an area fraction
        ! of ~0.4, so that it filled from a neighbouring front cell.

        implicit none

        real(wp), intent(OUT) :: a_lsf(:,:)             ! [--] Area fraction behind the front
        real(wp), intent(IN)  :: lsf(:,:)               ! [--] Level-set function (< 0: ice side)
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m]  Ice thickness
        real(wp), intent(IN)  :: z_bed(:,:)             ! [m]  Bedrock elevation
        real(wp), intent(IN)  :: z_sl(:,:)              ! [m]  Sea-level elevation
        character(len=*), intent(IN) :: boundaries

        integer  :: i, j, nx, ny, BC
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: l_c, l_W, l_E, l_S, l_N, l_SW, l_SE, l_NE, l_NW
        real(wp) :: phi_c, phi_W, phi_E, phi_S, phi_N
        real(wp) :: phi_BL, phi_BR, phi_TR, phi_TL
        logical, allocatable :: is_land(:,:)            ! Ice-free land cells

        nx = size(lsf,1)
        ny = size(lsf,2)
        BC = boundary_code(boundaries)

        allocate(is_land(nx,ny))
        is_land = z_bed .ge. z_sl .and. H_ice .le. H_ice_eps

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,phi_c,phi_W,phi_E,phi_S,phi_N) &
        !$omp& private(phi_BL,phi_BR,phi_TR,phi_TL,l_c,l_W,l_E,l_S,l_N,l_SW,l_SE,l_NE,l_NW)
        do j = 1, ny
        do i = 1, nx

            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! Neighbour values (ice-free land: the cell's own value)
            l_c  = lsf(i,j)
            l_W  = merge(l_c, lsf(im1,j),   is_land(im1,j))
            l_E  = merge(l_c, lsf(ip1,j),   is_land(ip1,j))
            l_S  = merge(l_c, lsf(i,jm1),   is_land(i,jm1))
            l_N  = merge(l_c, lsf(i,jp1),   is_land(i,jp1))
            l_SW = merge(l_c, lsf(im1,jm1), is_land(im1,jm1))
            l_SE = merge(l_c, lsf(ip1,jm1), is_land(ip1,jm1))
            l_NE = merge(l_c, lsf(ip1,jp1), is_land(ip1,jp1))
            l_NW = merge(l_c, lsf(im1,jp1), is_land(im1,jp1))

            ! Centre and edge midpoints
            phi_c = l_c
            phi_W = 0.5_wp*(l_W + l_c)
            phi_E = 0.5_wp*(l_c + l_E)
            phi_S = 0.5_wp*(l_S + l_c)
            phi_N = 0.5_wp*(l_c + l_N)

            ! Corners
            phi_BL = 0.25_wp*(l_SW + l_S + l_W + l_c)
            phi_BR = 0.25_wp*(l_S  + l_SE + l_c + l_E)
            phi_TR = 0.25_wp*(l_c  + l_E + l_N + l_NE)
            phi_TL = 0.25_wp*(l_W  + l_c + l_NW + l_N)

            ! Quadrants (BL, BR, TR, TL vertices of each)
            a_lsf(i,j) = 0.25_wp*( lsf_negative_area_fraction(phi_BL,phi_S, phi_c, phi_W ) &
                                 + lsf_negative_area_fraction(phi_S, phi_BR,phi_E, phi_c ) &
                                 + lsf_negative_area_fraction(phi_c, phi_E, phi_TR,phi_N ) &
                                 + lsf_negative_area_fraction(phi_W, phi_c, phi_N, phi_TL) )

        end do
        end do
        !$omp end parallel do

        return

    end subroutine calc_lsf_area_fraction

    function lsf_negative_area_fraction(phi_BL,phi_BR,phi_TR,phi_TL) result(f)
        ! Compute the area fraction of a unit square where a bilinear
        ! interpolant of the four corner phi values is negative. Uses the
        ! standard marching-squares topology with linear zero-crossing
        ! locations along each edge; the saddle case is split into two
        ! disjoint triangles (a centre-value test would give a sharper
        ! resolution but is unnecessary for the smooth LSF expected here).

        implicit none

        real(wp), intent(IN) :: phi_BL, phi_BR, phi_TR, phi_TL
        real(wp) :: f

        integer  :: mask
        real(wp) :: a, b, c, d

        ! Encode sign pattern: bit i set when corner phi > 0.
        mask = 0
        if (phi_BL .gt. 0.0_wp) mask = mask + 1
        if (phi_BR .gt. 0.0_wp) mask = mask + 2
        if (phi_TR .gt. 0.0_wp) mask = mask + 4
        if (phi_TL .gt. 0.0_wp) mask = mask + 8

        select case(mask)
        case(0)
            f = 1.0_wp
        case(15)
            f = 0.0_wp
        case(1)      ! BL positive (corner cut from ice)
            a = phi_BL/(phi_BL - phi_BR)
            b = phi_BL/(phi_BL - phi_TL)
            f = 1.0_wp - 0.5_wp*a*b
        case(2)      ! BR positive
            a = phi_BR/(phi_BR - phi_BL)
            b = phi_BR/(phi_BR - phi_TR)
            f = 1.0_wp - 0.5_wp*a*b
        case(4)      ! TR positive
            a = phi_TR/(phi_TR - phi_BR)
            b = phi_TR/(phi_TR - phi_TL)
            f = 1.0_wp - 0.5_wp*a*b
        case(8)      ! TL positive
            a = phi_TL/(phi_TL - phi_BL)
            b = phi_TL/(phi_TL - phi_TR)
            f = 1.0_wp - 0.5_wp*a*b
        case(14)     ! BL negative only (triangle of ice)
            a = phi_BL/(phi_BL - phi_BR)
            b = phi_BL/(phi_BL - phi_TL)
            f = 0.5_wp*a*b
        case(13)     ! BR negative only
            a = phi_BR/(phi_BR - phi_BL)
            b = phi_BR/(phi_BR - phi_TR)
            f = 0.5_wp*a*b
        case(11)     ! TR negative only
            a = phi_TR/(phi_TR - phi_BR)
            b = phi_TR/(phi_TR - phi_TL)
            f = 0.5_wp*a*b
        case(7)      ! TL negative only
            a = phi_TL/(phi_TL - phi_BL)
            b = phi_TL/(phi_TL - phi_TR)
            f = 0.5_wp*a*b
        case(3)      ! bottom positive, top negative (trapezoid on top)
            a = phi_TL/(phi_TL - phi_BL)
            b = phi_TR/(phi_TR - phi_BR)
            f = 0.5_wp*(a + b)
        case(12)     ! top positive, bottom negative
            a = phi_BL/(phi_BL - phi_TL)
            b = phi_BR/(phi_BR - phi_TR)
            f = 0.5_wp*(a + b)
        case(6)      ! right positive, left negative
            a = phi_BL/(phi_BL - phi_BR)
            b = phi_TL/(phi_TL - phi_TR)
            f = 0.5_wp*(a + b)
        case(9)      ! left positive, right negative
            a = phi_BR/(phi_BR - phi_BL)
            b = phi_TR/(phi_TR - phi_TL)
            f = 0.5_wp*(a + b)
        case(5)      ! saddle: BL+, BR-, TR+, TL-  (negative in BR and TL)
            a = phi_BR/(phi_BR - phi_BL)
            b = phi_BR/(phi_BR - phi_TR)
            c = phi_TL/(phi_TL - phi_TR)
            d = phi_TL/(phi_TL - phi_BL)
            f = 0.5_wp*a*b + 0.5_wp*c*d
        case(10)     ! saddle: BL-, BR+, TR-, TL+  (negative in BL and TR)
            a = phi_BL/(phi_BL - phi_BR)
            b = phi_BL/(phi_BL - phi_TL)
            c = phi_TR/(phi_TR - phi_TL)
            d = phi_TR/(phi_TR - phi_BR)
            f = 0.5_wp*a*b + 0.5_wp*c*d
        case default
            f = 0.0_wp
        end select

        if (f .lt. 0.0_wp) f = 0.0_wp
        if (f .gt. 1.0_wp) f = 1.0_wp

        return

    end function lsf_negative_area_fraction

    subroutine calc_ice_front(mask_frnt,f_ice,f_grnd,z_bed,z_sl,boundaries)
        ! Calculate a mask of ice front points that 
        ! demarcates both the ice front (last ice-covered point)
        ! and the ice-free front (first ice-free point), 
        ! and distinguishes between floating, marine and grounded fronts.

        implicit none 

        integer,  intent(OUT) :: mask_frnt(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: f_grnd(:,:) 
        real(wp), intent(IN)  :: z_bed(:,:)
        real(wp), intent(IN)  :: z_sl(:,:)
        character(len=*), intent(IN) :: boundaries 
        
        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        integer  :: n 
        real(wp) :: f_neighb(4) 
        integer  :: BC

        nx = size(mask_frnt,1) 
        ny = size(mask_frnt,2) 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Initialize mask to non-front everywhere to start 
        mask_frnt = MASK_FRNT_NONE

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,n,f_neighb)
        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            f_neighb = [f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)]
            n = count(f_neighb .lt. 1.0)

            if ( f_ice(i,j) .eq. 1.0_wp .and. n .gt. 0) then 
                ! This point is an ice front. 

                if (f_grnd(i,j) .gt. 0.0 .and. (z_sl(i,j) .le. z_bed(i,j)) ) then 
                    ! Ice front grounded above sea level 

                    mask_frnt(i,j) = MASK_FRNT_GRND 

                else if (f_grnd(i,j) .gt. 0.0) then 
                    ! Ice front grounded below sea level 

                    mask_frnt(i,j) = MASK_FRNT_MARINE 

                else
                    ! Floating ice front 

                    mask_frnt(i,j) = MASK_FRNT_FLOAT 

                end if 

                ! Ensure adjacent ice-free points are marked too, as ocean
                ! or land depending on their own bed, so that the front type
                ! can be decided per face (see set_ssa_masks)
                if (f_ice(im1,j) .lt. 1.0) mask_frnt(im1,j) = ice_free_code(z_bed(im1,j),z_sl(im1,j))
                if (f_ice(ip1,j) .lt. 1.0) mask_frnt(ip1,j) = ice_free_code(z_bed(ip1,j),z_sl(ip1,j))
                if (f_ice(i,jm1) .lt. 1.0) mask_frnt(i,jm1) = ice_free_code(z_bed(i,jm1),z_sl(i,jm1))
                if (f_ice(i,jp1) .lt. 1.0) mask_frnt(i,jp1) = ice_free_code(z_bed(i,jp1),z_sl(i,jp1))

            end if 

        end do
        end do 
        !$omp end parallel do

        return 

    contains

        pure integer function ice_free_code(z_bed_now,z_sl_now)
            ! Ice-free point next to a front: ocean if the bed is below sea level

            implicit none

            real(wp), intent(IN) :: z_bed_now, z_sl_now

            if (z_bed_now .lt. z_sl_now) then
                ice_free_code = MASK_FRNT_ICE_FREE
            else
                ice_free_code = MASK_FRNT_ICE_FREE_LAND
            end if

        end function ice_free_code

    end subroutine calc_ice_front

    elemental subroutine calc_z_srf_max(z_srf,H_ice,f_ice,z_bed,z_sl,rho_ice,rho_sw)
        ! Calculate surface elevation
        ! Adapted from Pattyn (2017), Eq. 1
        
        implicit none 

        real(wp), intent(INOUT) :: z_srf 
        real(wp), intent(IN)    :: H_ice
        real(wp), intent(IN)    :: f_ice
        real(wp), intent(IN)    :: z_bed
        real(wp), intent(IN)    :: z_sl
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw
        
        ! Local variables
        integer :: i, j, nx, ny 
        real(wp) :: rho_ice_sw
        real(wp) :: H_eff

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Get effective ice thickness
        call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)

        ! Initially calculate surface elevation everywhere 
        z_srf = max(z_bed + H_eff, z_sl + (1.0-rho_ice_sw)*H_eff)
        
        return 

    end subroutine calc_z_srf_max

    elemental subroutine calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero)
        ! Calculate ice-thickness, scaled at margins to actual thickness
        ! but as if it covered the whole grid cell.
        
        implicit none

        real(wp), intent(OUT) :: H_eff 
        real(wp), intent(IN)  :: H_ice 
        real(wp), intent(IN)  :: f_ice 
        logical,    intent(IN), optional :: set_frac_zero 

        if (f_ice .gt. 0.0) then 
            H_eff = H_ice / f_ice
        else 
            H_eff = H_ice 
        end if

        if (present(set_frac_zero)) then 
            if (set_frac_zero .and. f_ice .lt. 1.0) H_eff = 0.0_wp 
        end if 

        return

    end subroutine calc_H_eff

    elemental subroutine calc_H_grnd(H_grnd,H_ice,f_ice,z_bed,z_sl,rho_ice,rho_sw,use_f_ice)
        ! Calculate ice thickness overburden, H_grnd
        ! When H_grnd >= 0, grounded, when H_grnd < 0, floating 
        ! Also calculate rate of change for diagnostic related to grounding line 

        ! Note that this will not be identical to thickness above flotation,
        ! since it accounts for height above the water via z_bed too. 

        implicit none 

        real(wp), intent(INOUT) :: H_grnd
        real(wp), intent(IN)    :: H_ice
        real(wp), intent(IN)    :: f_ice
        real(wp), intent(IN)    :: z_bed
        real(wp), intent(IN)    :: z_sl 
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw
        logical,  intent(IN), optional :: use_f_ice

        ! Local variables   
        real(wp) :: rho_sw_ice 
        real(wp) :: H_eff 
        logical  :: use_f_ice_now 

        ! By default use f_ice to set points with f_ice < 1 to H_eff=0.0 
        use_f_ice_now = .TRUE. 
        if (present(use_f_ice)) use_f_ice_now = use_f_ice 

        rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
        ! Get effective ice thickness
        if (use_f_ice_now) then 
            call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)
        else 
            H_eff = H_ice 
        end if 

        ! Calculate new H_grnd (ice thickness overburden)
        !H_grnd = H_eff - rho_sw_ice*max(z_sl-z_bed,0.0_wp)

        ! ajr: testing. This ensures that ice-free ground above sea level
        ! also has H_grnd > 0.
        if (z_sl-z_bed .gt. 0.0) then 
            ! Grounded below sea level, diagnose overburden minus water thickness
            H_grnd = H_eff - rho_sw_ice*(z_sl-z_bed)
        else
            ! Grounded above sea level, simply sum elevation above sea level and ice thickness
            H_grnd = H_eff + (z_bed-z_sl)
        end if 

        ! ajr: to test somewhere eventually, more closely follows Gladstone et al (2010), Leguy et al (2021)
        !H_grnd = -( (z_sl-z_bed) - rho_ice_sw*H_eff )

        return 

    end subroutine calc_H_grnd

    elemental subroutine calc_H_af(H_af,H_ice,f_ice,z_bed,z_sl,rho_ice,rho_sw,use_f_ice)
        ! Calculate ice thickness above flotation, H_af

        implicit none 

        real(wp), intent(INOUT) :: H_af
        real(wp), intent(IN)    :: H_ice
        real(wp), intent(IN)    :: f_ice
        real(wp), intent(IN)    :: z_bed
        real(wp), intent(IN)    :: z_sl
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw
        
        logical,  intent(IN), optional :: use_f_ice

        ! Local variables   
        real(wp) :: rho_sw_ice 
        real(wp) :: H_eff 
        logical  :: use_f_ice_now 

        ! By default use f_ice to set points with f_ice < 1 to H_eff=0.0 
        use_f_ice_now = .TRUE. 
        if (present(use_f_ice)) use_f_ice_now = use_f_ice 

        rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
        ! Get effective ice thickness
        if (use_f_ice_now) then 
            call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)
        else 
            H_eff = H_ice 
        end if 

        ! Calculate ice thickness above flotation 
        ! ie, total ice column thickness minus the corresponding depth of seawater
        H_af = H_eff - rho_sw_ice*max(z_sl-z_bed,0.0_wp)

        ! Limit to positive values, since when it is negative, it is completely floating
        H_af = max(H_af,0.0_wp)

        return 

    end subroutine calc_H_af

    subroutine calc_f_grnd_subgrid_linear(f_grnd,f_grnd_x,f_grnd_y,H_grnd,boundaries)
        ! Calculate the grounded fraction of a cell in the x- and y-directions
        ! at the ac nodes
        !
        ! Given point 1 is grounded and point 2 is floating, we are looking for
        ! when H_grnd=0. The equation along an x-axis between point 1 (grounded point) 
        ! and point 2 (floating point) for H_grnd is:
        ! H_grnd = H_grnd_1 + f*(H_grnd_2-H_grnd_1)

        ! To find fraction f for when H_grnd == 0, rearrange equation:
        ! 0 = H_grnd_1 + f*(H_grnd_2-H_grnd_1)
        ! f = -H_grnd_1 / (H_grnd_2-H_grnd_1)
        ! where f is the distance along the 0:1 axis between the first and second points. 
        
        implicit none 

        real(wp), intent(OUT) :: f_grnd(:,:)
        real(wp), intent(OUT) :: f_grnd_x(:,:)
        real(wp), intent(OUT) :: f_grnd_y(:,:)
        real(wp), intent(IN)  :: H_grnd(:,:)
        character(len=*), intent(IN) :: boundaries

        ! Local variables  
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        real(wp) :: H_grnd_1, H_grnd_2
        integer  :: BC
        logical  :: per_x, per_y

        nx = size(f_grnd,1)
        ny = size(f_grnd,2)

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Central aa-node
        f_grnd = 1.0
        where (H_grnd < 0.0) f_grnd = 0.0
        
        ! x-direction, ac-node
        f_grnd_x = 1.0
        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,H_grnd_1,H_grnd_2)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(ip1,j) .le. 0.0) then 
                ! Point is grounded, neighbor is floating 

                H_grnd_1 = H_grnd(i,j) 
                H_grnd_2 = H_grnd(ip1,j) 

                ! Calculate fraction 
                f_grnd_x(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)

            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(ip1,j) .gt. 0.0) then 
                ! Point is floating, neighbor is grounded 

                H_grnd_1 = H_grnd(ip1,j) 
                H_grnd_2 = H_grnd(i,j) 

                ! Calculate fraction 
                f_grnd_x(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)

            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(ip1,j) .le. 0.0) then 
                ! Point is floating, neighbor is floating
                f_grnd_x(i,j) = 0.0 

            else 
                ! Point is grounded, neighbor is grounded
                f_grnd_x(i,j) = 1.0 

            end if 

        end do 
        end do 
        !$omp end parallel do

        ! y-direction, ac-node
        f_grnd_y = 1.0
        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,H_grnd_1,H_grnd_2)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i,jp1) .le. 0.0) then 
                ! Point is grounded, neighbor is floating 

                H_grnd_1 = H_grnd(i,j) 
                H_grnd_2 = H_grnd(i,jp1) 

                ! Calculate fraction 
                f_grnd_y(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)

            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i,jp1) .gt. 0.0) then 
                ! Point is floating, neighbor is grounded 

                H_grnd_1 = H_grnd(i,jp1) 
                H_grnd_2 = H_grnd(i,j) 

                ! Calculate fraction 
                f_grnd_y(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)
                
            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i,jp1) .le. 0.0) then 
                ! Point is floating, neighbor is floating
                f_grnd_y(i,j) = 0.0 

            else 
                ! Point is grounded, neighbor is grounded
                f_grnd_y(i,j) = 1.0 

            end if 

        end do 
        end do 
        !$omp end parallel do

        ! Set non-periodic boundary points equal to neighbor for aesthetics 
        call get_periodic_directions(per_x,per_y,BC)
        if (.not. per_x) f_grnd_x(nx,:) = f_grnd_x(nx-1,:) 
        if (.not. per_y) f_grnd_y(:,ny) = f_grnd_y(:,ny-1) 
        
        return 

    end subroutine calc_f_grnd_subgrid_linear
    
    subroutine calc_f_grnd_pinning_points(f_grnd,H_ice,f_ice,z_bed,z_bed_sd,z_sl,rho_ice,rho_sw)
        ! For floating points, determine how much of bed could be
        ! touching the base of the ice shelf due to subgrid pinning
        ! points following the distribution z = N(z_bed,z_bed_sd)

        implicit none

        real(wp), intent(OUT) :: f_grnd(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: z_bed(:,:) 
        real(wp), intent(IN)  :: z_bed_sd(:,:) 
        real(wp), intent(IN)  :: z_sl(:,:) 
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp) :: x, mu, sigma
        real(wp) :: rho_ice_sw
        real(wp) :: z_base_now
        real(wp) :: H_eff_now

        nx = size(f_grnd,1) 
        ny = size(f_grnd,2) 

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Initialize f_grnd to zero everywhere
        f_grnd = 0.0_wp 

        !$omp parallel do collapse(2) private(i,j,H_eff_now,z_base_now,x,mu,sigma)
        do j = 1, ny 
        do i = 1, nx 

            ! Get effective ice thickness
            call calc_H_eff(H_eff_now,H_ice(i,j),f_ice(i,j))

            ! Determine depth of base of ice shelf, assuming
            ! for now that it is floating. 
            z_base_now = z_sl(i,j) - rho_ice_sw*H_eff_now

            ! Check if it is actually floating: 

            if (z_base_now .gt. z_bed(i,j)) then 
                ! Floating point, calculate pinning fraction 

                ! Define parameters to calculate inverse CDF
                x     = z_base_now 
                mu    = z_bed(i,j) 
                sigma = z_bed_sd(i,j) 

                if ( is_equal(sigma,0.0_wp) ) then 
                    ! sigma not available, set f_grnd to zero 

                    f_grnd(i,j) = 0.0_wp 

                else

                    ! Calculate inverse cdf at z = ice_base, ie,
                    ! probability to have points at or higher than z,
                    ! which would be grounded. 

                    !f_grnd(i,j) = cdf(x,mu,sigma,inv=.TRUE.)

                    ! Alternatively, use approximation given by
                    ! Pollard and DeConto (2012), Eq. 13:

                    f_grnd(i,j) = 0.5_wp*max(0.0_wp,(1.0_wp-(x-mu)/sigma))

                end if 

            end if
            
        end do 
        end do 
        !$omp end parallel do

        return

    end subroutine calc_f_grnd_pinning_points

    subroutine remove_englacial_lakes(H_ice,z_bed,z_srf,z_sl,rho_ice,rho_sw)
        ! Diagnose where ice should be grounded, but a gap exists.
        ! In these locations, increase ice thickness in order
        ! to fill in the englacial lake. 
        ! Also check if ice thickness is below bedrock, in these
        ! places reduce ice thickness. 
        ! Used when loading datasets. 

        implicit none

        real(wp), intent(INOUT) :: H_ice(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: z_srf(:,:) 
        real(wp), intent(IN)    :: z_sl(:,:) 
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw 

        ! Local variables 
        integer :: i, j, nx , ny 
        real(wp) :: H_grnd_now 
        real(wp) :: rho_sw_ice
        
        rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Surface elevation is available, use to remove englacial lakes
        ! from ice thickness field. 
        do j = 1, ny 
        do i = 1, nx 

            ! Calculate H_grnd (ice thickness overburden)
            ! (ie, diagnose whether ice is grounded or floating)
            H_grnd_now = H_ice(i,j) - rho_sw_ice*max(z_sl(i,j)-z_bed(i,j),0.0_wp)

            if (H_grnd_now .gt. 0.0_wp .and. z_srf(i,j) - H_ice(i,j) .gt. z_bed(i,j)) then 
                ! Ice should be grounded, but a gap exists, 
                ! so a subglacial lake exists, increase thickness.

                H_ice(i,j) = z_srf(i,j) - z_bed(i,j) 

            else if (z_srf(i,j) - H_ice(i,j) .lt. z_bed(i,j)) then 
                ! Ice is erroneously thick, reduce it

                H_ice(i,j) = z_srf(i,j) - z_bed(i,j) 

            end if 

        end do 
        end do

        return

    end subroutine remove_englacial_lakes

    subroutine calc_distance_to_ice_margin(dist_mrgn,f_ice,dx,boundaries,calc_distances)
        ! Calculate distance to the ice margin
        
        ! Note: this subroutine is a wrapper that calls the
        ! grounding-line distance routine, since the algorithm
        ! works the same way. Simply substitute f_ice for f_grnd. 

        implicit none 

        real(wp), intent(OUT) :: dist_mrgn(:,:) ! [km] Distance to grounding line
        real(wp), intent(IN)  :: f_ice(:,:)     ! [1]  Fraction of grid-cell ice coverage 
        real(wp), intent(IN)  :: dx             ! [m]  Grid resolution (assume dy=dx)
        character(len=*), intent(IN) :: boundaries
        logical,  intent(IN)  :: calc_distances

        call calc_distance_to_grounding_line(dist_mrgn,f_ice,dx,boundaries,calc_distances)

        return 

    end subroutine calc_distance_to_ice_margin
    
    subroutine calc_distance_to_grounding_line(dist_gl,f_grnd,dx,boundaries,calc_distances)
        ! Calculate distance to the grounding line 
        
        implicit none 

        real(wp), intent(OUT) :: dist_gl(:,:)   ! [km] Distance to grounding line
        real(wp), intent(IN)  :: f_grnd(:,:)    ! [1]  Grounded grid-cell fraction 
        real(wp), intent(IN)  :: dx             ! [m]  Grid resolution (assume dy=dx)
        character(len=*), intent(IN) :: boundaries
        logical,  intent(IN)  :: calc_distances

        ! Local variables 
        integer  :: i, j, nx, ny, q
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: dist_direct_min 
        real(wp) :: dist_corners_min
        real(wp) :: dx_km 
        real(wp) :: dists(8) 
        integer  :: BC

        real(wp), parameter :: dist_max = 1e10          ! [km]
        real(wp), parameter :: sqrt_2   = sqrt(2.0_wp) 
        integer,  parameter :: iter_max = 1000

        real(wp), allocatable :: dist_gl_ref(:,:) 

        nx = size(dist_gl,1)
        ny = size(dist_gl,2)

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        allocate(dist_gl_ref(nx,ny)) 

        ! Units
        dx_km = dx*1e-3 

        ! 0: Assign a big distance to all points to start  ======================

        dist_gl  = dist_max

        ! 1. Next, determine grounding line =====================================

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! Grounded point or partially floating point with floating neighbors
            if (f_grnd(i,j) .gt. 0.0 .and. &
                ( is_equal(f_grnd(im1,j),0.0_wp) .or. is_equal(f_grnd(ip1,j),0.0_wp) .or. &
                 is_equal(f_grnd(i,jm1),0.0_wp) .or. is_equal(f_grnd(i,jp1),0.0_wp) ) ) then 
                
                dist_gl(i,j)  = 0.0_wp 

            end if 

        end do 
        end do 
        !$omp end parallel do

        ! 2. Next, determine distances to grounding line ======================
        
        if (calc_distances) then

            do q = 1, iter_max  
                ! Iterate distance of one neighbor at a time until grid is filled in

                dist_gl_ref = dist_gl

                !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,dists,dist_direct_min,dist_corners_min)
                do j = 1, ny 
                do i = 1, nx

                    if ( is_equal(dist_gl(i,j),dist_max) ) then 
                        ! Distance needs to be determined for this point 

                        ! Get neighbor indices
                        call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

                        ! Get distances to direct and corner neighbors
                        dists = [dist_gl_ref(im1,j),dist_gl_ref(ip1,j), &       ! Direct neighbors
                                dist_gl_ref(i,jm1),dist_gl_ref(i,jp1), &       ! Direct neighbors
                                dist_gl_ref(im1,jp1),dist_gl_ref(ip1,jp1), &   ! Corner neighbors
                                dist_gl_ref(im1,jm1),dist_gl_ref(ip1,jm1)]     ! Corner neighbors

                        if (count(dists .lt. dist_max) .ge. 2) then 
                            ! Some neighbors have been calculated already 
                            ! Note: check for at least 2 neighbors already
                            ! calculated to reduce errors in points with 
                            ! neighbors on both sides with similar distances. This 
                            ! Waiting for even more neighbors would reduce errors
                            ! further but greatly increases iterations. 

                            ! Determine minimum distance to grounding line for 
                            ! direct and diagonal neighbors, separately.
                            dist_direct_min  = minval(dists(1:4))
                            dist_corners_min = minval(dists(5:8))

                            if (dist_direct_min .le. dist_corners_min) then 
                                ! Assume nearest path to grounding line is along
                                ! direct neighbor path - add dx to dist for this point 

                                dist_gl(i,j) = dist_direct_min + dx_km 

                            else
                                ! Assume nearest path to grounding line is via 
                                ! a diagonal neighbor - add sqrt(2) to dist for this point

                                dist_gl(i,j) = dist_corners_min + sqrt_2*dx_km

                            end if 

                        end if 
                    end if

                end do 
                end do
                !$omp end parallel do
                
                if (count(dist_gl .eq. dist_max) .eq. 0) then 
                    ! No more points to check 
                    exit 
                end if 

            end do

        end if 

        ! Set all floating-point distances to negative values 
        where (f_grnd .eq. 0.0_wp) 
            dist_gl  = -dist_gl
        end where

        
        return 

    end subroutine calc_distance_to_grounding_line
    
    subroutine calc_grounding_line_zone(mask_grz,dist_gl,dist_grz)
        ! Define a mask of the grounding-zone region 
        ! mask_grz == -2: floating-point out of grounding zone
        ! mask_grz == -1: floating-point in grounding zone
        ! mask_grz ==  0: grounding line 
        ! mask_grz ==  1: grounded-point in grounding zone
        ! mask_grz ==  2: grounded-point out of grounding zone

        implicit none

        integer,  intent(OUT) :: mask_grz(:,:)  ! [1]  Grounding-zone mask 
        real(wp), intent(IN)  :: dist_gl(:,:)   ! [km] Distance to grounding line
        real(wp), intent(IN)  :: dist_grz       ! [km] Distance to define as grounding zone

        ! Finally define the grounding-zone mask too
        where(dist_gl .eq. 0.0_wp)
            mask_grz =  0 
        else where(dist_gl .lt. 0.0_wp .and. abs(dist_gl) .le. dist_grz)
            mask_grz = -1
        else where(dist_gl .gt. 0.0_wp .and.     dist_gl  .le. dist_grz)
            mask_grz =  1
        else where(dist_gl .lt. 0.0_wp)
            mask_grz = -2
        elsewhere
            mask_grz =  2
        end where

        return

    end subroutine calc_grounding_line_zone

    subroutine calc_bmb_total(bmb,bmb_grnd,bmb_shlf,H_ice,H_grnd,f_grnd,gz_Hg0,gz_Hg1, &
                                                            gz_nx,bmb_gl_method,boundaries)

        implicit none 

        real(wp),         intent(OUT) :: bmb(:,:) 
        real(wp),         intent(IN)  :: bmb_grnd(:,:) 
        real(wp),         intent(IN)  :: bmb_shlf(:,:) 
        real(wp),         intent(IN)  :: H_ice(:,:)
        real(wp),         intent(IN)  :: H_grnd(:,:)
        real(wp),         intent(IN)  :: f_grnd(:,:) 
        real(wp),         intent(IN)  :: gz_Hg0
        real(wp),         intent(IN)  :: gz_Hg1
        integer,          intent(IN)  :: gz_nx
        character(len=*), intent(IN)  :: bmb_gl_method 
        character(len=*), intent(IN)  :: boundaries 

        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: n_float  
        real(wp) :: bmb_shlf_now 

        nx = size(bmb,1)
        ny = size(bmb,2) 

        ! Combine floating and grounded parts into one field =========================
        ! Apply the floating basal mass balance according 
        ! to different subgridding options at the grounding line
        ! (following the notation of Leguy et al., 2021 - see Fig. 3)
        
        select case(bmb_gl_method)

            case("fcmp")
                ! Flotation criterion melt parameterization 
                ! Apply full bmb_shlf value where flotation criterion is met 

                where(H_grnd .le. 0.0_wp)
                    
                    bmb = bmb_shlf

                elsewhere

                    bmb = bmb_grnd 

                end where 

            case("fmp")
                ! Full melt parameterization
                ! Apply full bmb_shlf value to any cell that is at least
                ! partially floating. Perhaps unrealistic.

                where(f_grnd .lt. 1.0_wp)

                    bmb = bmb_shlf 

                elsewhere

                    bmb = bmb_grnd 

                end where 

            case("pmp")
                ! Partial melt parameterization
                ! Apply bmb_shlf to floating fraction of cell 

                where(f_grnd .lt. 1.0_wp)

                    bmb = f_grnd*bmb_grnd + (1.0_wp-f_grnd)*bmb_shlf 

                elsewhere

                    bmb = bmb_grnd 

                end where 

            case("pmpt")
                ! Partial melt parameterization with tidal grounding zone

                call calc_bmb_gl_pmpt(bmb,bmb_grnd,bmb_shlf,H_grnd,gz_Hg0,gz_Hg1,gz_nx,boundaries)

            case("nmp")
                ! No melt parameterization
                ! Apply bmb_shlf only where fully floating diagnosed

                where(f_grnd .eq. 0.0_wp)

                    bmb = bmb_shlf 

                elsewhere

                    bmb = bmb_grnd 

                end where 

        end select

        ! For aesthetics, also make sure that bmb is zero on ice-free land
        where (H_grnd .gt. 0.0_wp .and. H_ice .eq. 0.0_wp) bmb = 0.0_wp 

        return 

    end subroutine calc_bmb_total

    subroutine calc_fmb_total(fmb,fmb_shlf,bmb_shlf,H_ice,H_grnd,f_ice, &
                                fmb_method,fmb_scale,fmb_lambda,rho_ice,rho_sw,dx,boundaries, &
                                Q_sg, tf_shlf)

        implicit none 

        real(wp), intent(OUT) :: fmb(:,:) 
        real(wp), intent(IN)  :: fmb_shlf(:,:) 
        real(wp), intent(IN)  :: bmb_shlf(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: H_grnd(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        integer,  intent(IN)  :: fmb_method 
        real(wp), intent(IN)  :: fmb_scale
        real(wp), intent(IN)  :: fmb_lambda         
        real(wp), intent(IN)  :: rho_ice
        real(wp), intent(IN)  :: rho_sw
        real(wp), intent(IN)  :: dx
        character(len=*), intent(IN) :: boundaries
        real(wp), intent(IN), optional :: Q_sg(:,:)      ! Subglacial discharge [m3/s]
        real(wp), intent(IN), optional :: tf_shlf(:,:)  
    
        ! Local variables
        integer    :: i, j, nx, ny, n_margin
        integer    :: im1, ip1, jm1, jp1
        real(wp) :: H_eff, dz  
        real(wp) :: area_flt
        real(wp) :: area_tot
        real(wp) :: bmb_eff 
        logical  :: mask(4) 
        integer  :: BC


        real(wp) :: rho_ice_sw 
        
        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Total cell area 
        area_tot = dx*dx 

        nx = size(fmb,1)
        ny = size(fmb,2) 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        select case(fmb_method)

            case(0)
                ! fmb provided by boundary field fmb_shlf 

                fmb = fmb_shlf

            case(1,2) 
                ! Calculate fmb (1) as proportional to local bmb_shlf value, 
                ! or (2) from a reference field, but scaled to the area
                ! of the grid cell itself where it will be applied

                !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,mask,n_margin,H_eff,dz,area_flt,bmb_eff)
                do j = 1, ny 
                do i = 1, nx 

                    ! Get neighbor indices
                    call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

                    ! Get mask of neighbors that are ice free 
                    mask = ( [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)].eq.0.0 )

                    ! Count how many neighbors are ice free 
                    n_margin = count(mask) 

                    if (H_ice(i,j) .gt. 0.0_wp .and. &
                        H_grnd(i,j) .lt. H_ice(i,j) .and. &
                                            n_margin .gt. 0) then 
                        ! Cell is ice-covered, [grounded below sea level or floating] and at the ice margin

                        ! Get effective ice thickness
                        call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

                        ! Determine depth of adjacent water using centered cell information
                        if (H_grnd(i,j) .lt. 0.0_wp) then 
                            ! Cell is floating, calculate submerged ice thickness 
                            
                            dz = (H_eff*rho_ice_sw)
                            
                        else 
                            ! Cell is grounded, recover depth of seawater

                            dz = max( (H_eff - H_grnd(i,j)) * rho_ice_sw, 0.0_wp)

                        end if 

                        ! Get area of ice submerged and adjacent to seawater
                        area_flt = real(n_margin,wp)*dz*dx 

                        if (fmb_method .eq. 1) then
                            ! Also calculate the mean bmb_shlf value for the ice-free neighbors 
                            bmb_eff = sum([bmb_shlf(im1,j),bmb_shlf(ip1,j),bmb_shlf(i,jm1),bmb_shlf(i,jp1)], &
                                            mask=mask) / real(n_margin,wp)

                            ! Finally calculate the effective front mass balance rate 

                            fmb(i,j) = bmb_eff*(area_flt/area_tot)*fmb_scale

                        else if (fmb_method .eq. 2) then

                            ! Finally calculate the scaled front mass balance rate 

                            fmb(i,j) = fmb_shlf(i,j)*(area_flt/area_tot)
                        end if

                    else 
                        ! Set front mass balance equal to zero 

                        fmb(i,j) = 0.0_wp 

                    end if 

                end do 
                end do
                !$omp end parallel do

            case(3)
                ! Calculate fmb using the Rignot et al. (2016) parameterization
                ! (see calc_melt_rate_rignot16), scaled by fmb_lambda, with
                ! h the depth of the submerged face, A its area, Q_sg the
                ! subglacial discharge [m3/s] and TF = tf_shlf [K]
                
                !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,mask,n_margin,H_eff,dz,area_flt)
                do j = 1, ny 
                do i = 1, nx 
 
                    ! Get neighbor indices
                    call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
 
                    ! Get mask of neighbors that are ice free 
                    mask = ( [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)].eq.0.0 )
 
                    ! Count how many neighbors are ice free 
                    n_margin = count(mask) 
 
                    if (H_ice(i,j) .gt. 0.0_wp .and. &
                        H_grnd(i,j) .lt. H_ice(i,j) .and. &
                                            n_margin .gt. 0) then 
                        ! Cell is ice-covered, [grounded below sea level or floating] and at the ice margin
 
                        ! Get effective ice thickness
                        call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))
 
                        ! Determine depth of adjacent water (h in Rignot et al., 2016)
                        if (H_grnd(i,j) .lt. 0.0_wp) then 
                            ! Cell is floating, calculate submerged ice thickness 
                            dz = (H_eff*rho_ice_sw)
                            
                        else 
                            ! Cell is grounded, recover depth of seawater
                            dz = max( (H_eff - H_grnd(i,j)) * rho_ice_sw, 0.0_wp)
 
                        end if 
 
                        ! Get area of ice submerged and adjacent to seawater
                        area_flt = real(n_margin,wp)*dz*dx 

                        if (area_flt .gt. 0.0_wp) then
                            ! Melt rate [m/yr] of the submerged face, A = area_flt [m2],
                            ! scaled to the cell area
                            fmb(i,j) = - fmb_lambda * calc_melt_rate_rignot16(dz,Q_sg(i,j),area_flt,tf_shlf(i,j)) &
                                            * (area_flt/area_tot)
                        else
                            fmb(i,j) = 0.0_wp
                        end if 
 
                    else 
                        ! Set front mass balance equal to zero 
                        fmb(i,j) = 0.0_wp 
 
                    end if 
 
                end do 
                end do
                !$omp end parallel do

            case DEFAULT 

                write(*,*) "calc_fmb_total:: Error: fmb_method not recognized."
                write(*,*) "fmb_method = ", fmb_method 
                error stop 1

        end select

        return 

    end subroutine calc_fmb_total

    elemental function calc_melt_rate_rignot16(h,Q,area,TF) result(m)
        ! Submarine melt rate of a marine-terminating glacier front
        ! following Rignot et al. (2016) (ISMIP7 protocol):
        !
        !   m = (a*h*q^alpha + b)*TF^beta [m/d],  q = 86400*Q/area [m/d]
        !
        ! h:  water depth at the front [m]
        ! Q:  subglacial discharge [m3/s], a non-negative volume flux
        !     (q**alpha is NaN for q < 0)
        ! area: submerged area of the front face [m2]
        ! TF: thermal forcing [K]
        ! The rate is returned in [m/yr].

        implicit none

        real(wp), intent(IN) :: h
        real(wp), intent(IN) :: Q
        real(wp), intent(IN) :: area
        real(wp), intent(IN) :: TF
        real(wp) :: m

        ! Local variables
        real(wp) :: q_now

        real(wp), parameter :: a       = 3.0e-4_wp
        real(wp), parameter :: b       = 0.15_wp
        real(wp), parameter :: alpha   = 0.39_wp
        real(wp), parameter :: beta    = 1.18_wp
        real(wp), parameter :: days_yr = 365.25_wp

        if (area .gt. 0.0_wp) then
            q_now = 86400.0_wp*max(Q,0.0_wp)/area
        else
            q_now = 0.0_wp
        end if

        m = days_yr * (a*h*q_now**alpha + b) * max(TF,0.0_wp)**beta

        return

    end function calc_melt_rate_rignot16

    subroutine calc_bmb_gl_pmpt(bmb,bmb_grnd,bmb_shlf,H_grnd,gz_Hg0,gz_Hg1,nxi,boundaries)
        ! Calculate basal mass balance, with bmb at the grounding line
        ! determined via subgrid calculation of flotation and parameterization
        ! for tidal-induced grounded melt. 

        implicit none
        
        real(wp), intent(OUT) :: bmb(:,:)           ! aa-nodes 
        real(wp), intent(IN)  :: bmb_grnd(:,:)      ! aa-nodes 
        real(wp), intent(IN)  :: bmb_shlf(:,:)      ! aa-nodes 
        real(wp), intent(IN)  :: H_grnd(:,:)        ! aa-nodes
        real(wp), intent(IN)  :: gz_Hg0             ! Lower limit in H_grnd for grounding zone
        real(wp), intent(IN)  :: gz_Hg1             ! Upper limit in H_grnd for grounding zone
        integer,  intent(IN)  :: nxi                ! Number of interpolation points per side (nxi*nxi)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables
        integer  :: i, j, i1, j1, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: Hg_nb(9)
        real(wp) :: wt 
        integer  :: BC

        real(wp), allocatable :: Hg_int(:,:)
        real(wp), allocatable :: bmb_int(:,:)

        ! Consistency check
        if (gz_Hg0 .gt. 0.0) then 
            write(io_unit_err,*) "calc_bmb_gl_pmpt:: Error: lower limit on grounding zone must be <= 0.0."
            write(io_unit_err,*) "gz_Hg0 = ", gz_Hg0
            error stop 1
        end if 
        if (gz_Hg1 .lt. 0.0) then 
            write(io_unit_err,*) "calc_bmb_gl_pmpt:: Error: upper limit on grounding zone must be >= 0.0."
            write(io_unit_err,*) "gz_Hg1 = ", gz_Hg1
            error stop 1
        end if 

        nx = size(H_grnd,1)
        ny = size(H_grnd,2) 

        ! Set boundary condition code
        BC = boundary_code(boundaries)

        ! Allocate subgrid arrays
        allocate(Hg_int(nxi,nxi))
        allocate(bmb_int(nxi,nxi))

        !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,Hg_nb,Hg_int,bmb_int,i1,j1,wt)
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

            ! Gather the 3x3 neighborhood explicitly (a slice im1:ip1 is
            ! empty where the neighbor indices wrap around the domain)
            Hg_nb = [H_grnd(im1,jm1),H_grnd(i,jm1),H_grnd(ip1,jm1), &
                     H_grnd(im1,j),  H_grnd(i,j),  H_grnd(ip1,j),   &
                     H_grnd(im1,jp1),H_grnd(i,jp1),H_grnd(ip1,jp1)]

            if (minval(Hg_nb) .ge. gz_Hg1) then
                ! Entire cell is grounded

                bmb(i,j) = bmb_grnd(i,j)
            
            else if (maxval(Hg_nb) .lt. gz_Hg0) then
                ! Entire cell is floating

                bmb(i,j) = bmb_shlf(i,j) 

            else
                ! Point contains the grounding zone
                
                ! Calculate subgrid values of H_grnd (bilinear between cell
                ! centres, as the grounded fraction of gl_sep = 3)
                call calc_subgrid_array_quad(Hg_int,H_grnd,nxi,i,j,im1,ip1,jm1,jp1)

                ! Calculate individual bmb values for each subgrid point
                do j1 = 1, nxi
                do i1 = 1, nxi

                    if (Hg_int(i1,j1) .lt. gz_Hg0) then
                        ! Floating, outside of grounding zone 
                        wt = 0.0
                    else if (Hg_int(i1,j1) .ge. gz_Hg1) then
                        ! Grounded, outside of grounding zone
                        wt = 1.0 
                    else
                        ! Within grounding zone
                        wt = (Hg_int(i1,j1)-gz_Hg0) / (gz_Hg1 - gz_Hg0)
                    end if 

                    ! Get subgrid bmb weighted between floating and grounded contributions
                    bmb_int(i1,j1) = wt*bmb_grnd(i,j) + (1.0-wt)*bmb_shlf(i,j) 

                end do
                end do

                ! Get the mean bmb rate for the entire cell
                bmb(i,j) = sum(bmb_int) / real(nxi*nxi,wp)

            end if 

        end do 
        end do
        !$omp end parallel do

        return
        
    end subroutine calc_bmb_gl_pmpt
    
!! f_grnd calculations from IMAU-ICE / CISM 

! == Routines for determining the grounded fraction on all four grids
  
  subroutine determine_grounded_fractions(f_grnd,f_grnd_acx,f_grnd_acy,f_grnd_ab,H_grnd,boundaries)
    ! Determine the grounded fraction of centered and staggered grid points
    ! Uses the bilinear interpolation scheme (with analytical solutions) 
    ! from CISM (Leguy et al., 2021), as adapted from IMAU-ICE v2.0 code (rev. 4776833b)
    
    implicit none
    
    real(wp), intent(OUT) :: f_grnd(:,:) 
    real(wp), intent(OUT), optional :: f_grnd_acx(:,:) 
    real(wp), intent(OUT), optional :: f_grnd_acy(:,:) 
    real(wp), intent(OUT), optional :: f_grnd_ab(:,:) 
    real(wp), intent(IN)  :: H_grnd(:,:) 
    character(len=*), intent(IN) :: boundaries
    
    ! Local variables
    integer :: i, j, nx, ny 
    integer :: im1, ip1, jm1, jp1
    integer :: BC

    real(wp), allocatable :: f_grnd_NW(:,:)
    real(wp), allocatable :: f_grnd_NE(:,:)
    real(wp), allocatable :: f_grnd_SW(:,:)
    real(wp), allocatable :: f_grnd_SE(:,:)
    real(wp), allocatable :: f_flt(:,:) 

    nx = size(f_grnd,1)
    ny = size(f_grnd,2) 
    
    ! Set boundary condition code
    BC = boundary_code(boundaries)

    allocate(f_grnd_NW(nx,ny))
    allocate(f_grnd_NE(nx,ny))
    allocate(f_grnd_SW(nx,ny))
    allocate(f_grnd_SE(nx,ny))
    
    allocate(f_flt(nx,ny))

    ! Define aa-node variable f_flt as the flotation function
    ! following Leguy et al. (2021), Eq. 6. 
    ! Note: -H_grnd is not exactly the same, since it is in ice thickness,
    ! whereas the L21 equation is in water equivalent thickness, and it includes
    ! bedrock above sea level. But test this as it is first. 
    f_flt = -H_grnd
    
    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    call determine_grounded_fractions_CISM_quads(f_grnd_NW,f_grnd_NE,f_grnd_SW,f_grnd_SE,f_flt,boundaries)
    
    ! Get grounded fractions on all four grids by averaging over the quadrants
    !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
    do j = 1, ny
    do i = 1, nx 
        
        ! Get neighbor indices
        call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

      ! aa-nodes
      f_grnd(i,j)     = 0.25_wp * (f_grnd_NW(i,j) + f_grnd_NE(i,j) + f_grnd_SW(i,j) + f_grnd_SE(i,j))
      
      if (present(f_grnd_acx)) then
        ! acx-nodes
        f_grnd_acx(i,j) = 0.25_wp * (f_grnd_NE(i,j) + f_grnd_SE(i,j) + f_grnd_NW(ip1,j) + f_grnd_SW(ip1,j))
      end if 

      if (present(f_grnd_acy)) then
        ! acy-nodes
        f_grnd_acy(i,j) = 0.25_wp * (f_grnd_NE(i,j) + f_grnd_NW(i,j) + f_grnd_SE(i,jp1) + f_grnd_SW(i,jp1))
      end if 

      if (present(f_grnd_ab)) then
        ! ab-nodes
        f_grnd_ab(i,j)  = 0.25_wp * (f_grnd_NE(i,j) + f_grnd_NW(ip1,j) + f_grnd_SE(i,jp1) + f_grnd_SW(ip1,jp1))
      end if

    end do
    end do
    !$omp end parallel do
    
    return 

  end subroutine determine_grounded_fractions

  subroutine determine_grounded_fractions_CISM_quads(f_grnd_NW,f_grnd_NE,f_grnd_SW,f_grnd_SE,f_flt,boundaries)
    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    ! (using the approach from CISM, where grounded fractions are calculated
    !  based on analytical solutions to the bilinear interpolation)
    
    implicit none
    
    real(wp), intent(OUT) :: f_grnd_NW(:,:)
    real(wp), intent(OUT) :: f_grnd_NE(:,:)
    real(wp), intent(OUT) :: f_grnd_SW(:,:)
    real(wp), intent(OUT) :: f_grnd_SE(:,:)
    real(wp), intent(IN)  :: f_flt(:,:) 
    character(len=*), intent(IN) :: boundaries

    ! Local variables:
    integer  :: i, j, ii, jj, nx, ny
    integer  :: im1, ip1, jm1, jp1  
    real(wp) :: f_NW, f_N, f_NE, f_W, f_m, f_E, f_SW, f_S, f_SE
    real(wp) :: fq_NW, fq_NE, fq_SW, fq_SE
    integer  :: BC

    nx = size(f_flt,1)
    ny = size(f_flt,2)
    
    ! Set boundary condition code
    BC = boundary_code(boundaries)

    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    !$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1, f_NW,f_N,f_NE,f_W,f_m,f_E,f_SW,f_S,f_SE, fq_NW,fq_NE,fq_SW,fq_SE)
    do j = 1, ny
    do i = 1, nx
        
        ! Get neighbor indices
        call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

      f_NW = 0.25_wp * (f_flt(im1,jp1) + f_flt(i,jp1)   + f_flt(im1,j)   + f_flt(i,j))
      f_N  = 0.50_wp * (f_flt(i,jp1)   + f_flt(i,j))
      f_NE = 0.25_wp * (f_flt(i,jp1)   + f_flt(ip1,jp1) + f_flt(i,j)     + f_flt(ip1,j))
      f_W  = 0.50_wp * (f_flt(im1,j)   + f_flt(i,j))
      f_m  = f_flt(i,j) 
      f_E  = 0.50_wp * (f_flt(i,j)     + f_flt(ip1,j))
      f_SW = 0.25_wp * (f_flt(im1,j)   + f_flt(i,j)     + f_flt(im1,jm1) + f_flt(i,jm1))
      f_S  = 0.50_wp * (f_flt(i,j)     + f_flt(i,jm1))
      f_SE = 0.25_wp * (f_flt(i,j)     + f_flt(ip1,j)   + f_flt(i,jm1)   + f_flt(ip1,jm1))

      ! NW
      fq_NW = f_NW
      fq_NE = f_N
      fq_SW = f_W
      fq_SE = f_m
      call calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_NW(i,j) )
      
      ! NE
      fq_NW = f_N
      fq_NE = f_NE
      fq_SW = f_m
      fq_SE = f_E
      call calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_NE(i,j) )
      
      ! SW
      fq_NW = f_W
      fq_NE = f_m
      fq_SW = f_SW
      fq_SE = f_S
      call calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_SW(i,j) )
      
      ! SE
      fq_NW = f_m
      fq_NE = f_E
      fq_SW = f_S
      fq_SE = f_SE
      call calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_SE(i,j) )
      
    end do
    end do
    !$omp end parallel do

    return 

  end subroutine determine_grounded_fractions_CISM_quads

  subroutine calc_fraction_above_zero( f_NW, f_NE, f_SW, f_SE, phi)
    ! Given a square with function values at the four corners,
    ! calculate the fraction phi of the square where the function is larger than zero.
    
    ! Note: calculations below require double precision!! 
    
    implicit none
    
    ! In/output variables:
    real(wp), intent(IN)    :: f_NW, f_NE, f_SW, f_SE
    real(wp), intent(OUT)   :: phi
    
    ! Local variables:
    real(dp) :: f_NWp, f_NEp, f_SWp, f_SEp
    real(dp) :: aa,bb,cc,dd,x,f1,f2
    integer  :: scen

    real(wp), parameter :: ftol = 1e-4_dp
    
    ! The analytical solutions sometime give problems when one or more of the corner
    ! values is VERY close to zero; avoid this.
    if (f_NW == 0.0_dp) then
      f_NWp = ftol
    else if (f_NW > 0.0_dp) then
      f_NWp = MAX(  ftol, f_NW)
    else if (f_NW < 0.0_dp) then
      f_NWp = MIN( -ftol, f_NW)
    else
      f_NWp = f_NW
    end if
    if (f_NE == 0.0_dp) then
      f_NEp = ftol
    else if (f_NE > 0.0_dp) then
      f_NEp = MAX(  ftol, f_NE)
    else if (f_NE < 0.0_dp) then
      f_NEp = MIN( -ftol, f_NE)
    else
      f_NEp = f_NE
    end if
    if (f_SW == 0.0_dp) then
      f_SWp = ftol
    else if (f_SW > 0.0_dp) then
      f_SWp = MAX(  ftol, f_SW)
    else if (f_SW < 0.0_dp) then
      f_SWp = MIN( -ftol, f_SW)
    else
      f_SWp = f_SW
    end if
    if (f_SE == 0.0_dp) then
      f_SEp = ftol
    else if (f_SE > 0.0_dp) then
      f_SEp = MAX(  ftol, f_SE)
    else if (f_SE < 0.0_dp) then
      f_SEp = MIN( -ftol, f_SE)
    else
      f_SEp = f_SE
    end if
    
    if (f_NWp <= 0.0_dp .AND. f_NEp <= 0.0_dp .AND. f_SWp <= 0.0_dp .AND. f_SEp <= 0.0_dp) then
      ! All four corners are grounded.
      
      phi = 1.0_wp
      
    else if (f_NWp >= 0.0_dp .AND. f_NEp >= 0.0_dp .AND. f_SWp >= 0.0_dp .AND. f_SEp >= 0.0_dp) then
      ! All four corners are floating
      
      phi = 0.0_wp
      
    else
      ! At least one corner is grounded and at least one is floating;
      ! the grounding line must pass through this square!
      
      ! Only four "scenarios" exist (with rotational symmetries):
      ! 1: SW grounded, rest floating
      ! 2: SW floating, rest grounded
      ! 3: south grounded, north floating
      ! 4: SW & NE grounded, SE & NW floating
      ! Rotate the four-corner world until it matches one of these scenarios.
      call rotate_quad_until_match( f_NWp, f_NEp, f_SWp, f_SEp, scen)
    
      ! Calculate initial values of coefficients, and make correction
      ! for when d=0 (to avoid problems)

      aa  = f_SWp
      bb  = f_SEp - f_SWp
      cc  = f_NWp - f_SWp
      dd  = f_NEp + f_SWp - f_NWp - f_SEp

      ! Exception for when d=0
      if (ABS(dd) < ftol) then
        if (f_SWp > 0.0_dp) then
          f_SWp = f_SWp + 0.1_dp
        else
          f_SWp = f_SWp - 0.1_dp
        end if
        aa  = f_SWp
        bb  = f_SEp - f_SWp
        cc  = f_NWp - f_SWp
        dd  = f_NEp + f_SWp - f_NWp - f_SEp
      end if
        
      if (scen == 1) then
        ! 1: SW grounded, rest floating
        
        phi = ((bb*cc - aa*dd) * log(abs(1.0_dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
         
      else if (scen == 2) then
        ! 2: SW floating, rest grounded
        ! Assign negative coefficients to calculate floating fraction,
        ! then get complement to obtain grounded fraction. 

        aa  = -(f_SWp)
        bb  = -(f_SEp - f_SWp)
        cc  = -(f_NWp - f_SWp)
        dd  = -(f_NEp + f_SWp - f_NWp - f_SEp)

        ! Exception for when d=0
        if (ABS(dd) < 1e-4_dp) then
          if (f_SWp > 0.0_dp) then
            f_SWp = f_SWp + 0.1_dp
          else
            f_SWp = f_SWp - 0.1_dp
          end if
          aa  = -(f_SWp)
          bb  = -(f_SEp - f_SWp)
          cc  = -(f_NWp - f_SWp)
          dd  = -(f_NEp + f_SWp - f_NWp - f_SEp)
        end if
        
        phi = 1.0_dp - ((bb*cc - aa*dd) * log(abs(1.0_dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
        
      else if (scen == 3) then
        ! 3: south grounded, north floating
        
        ! Exception for when the GL runs parallel to the x-axis
        if (abs( 1.0_dp - f_NWp/f_NEp) < 1e-6_dp .and. abs( 1.0_dp - f_SWp/f_SEp) < 1e-6_dp) then
          
          phi = f_SWp / (f_SWp - f_NWp)
          
        else
            
          x   = 0.0_dp
          f1  = ((bb*cc - aa*dd) * log(abs(cc+dd*x)) - bb*dd*x) / (dd**2)
          x   = 1.0_dp
          f2  = ((bb*cc - aa*dd) * log(abs(cc+dd*x)) - bb*dd*x) / (dd**2)
          phi = f2-f1
                  
        end if
        
      else if (scen == 4) then
        ! 4: SW & NE grounded, SE & NW floating
        ! (recalculate coefficients here explicitly for two cases)

        ! SW corner
        aa  = f_SWp
        bb  = f_SEp - f_SWp
        cc  = f_NWp - f_SWp
        dd  = f_NEp + f_SWp - f_NWp - f_SEp
        phi = ((bb*cc - aa*dd) * log(abs(1.0_dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
        
        ! NE corner
        call rotate_quad( f_NWp, f_NEp, f_SWp, f_SEp)
        call rotate_quad( f_NWp, f_NEp, f_SWp, f_SEp)
        aa  = f_SWp
        bb  = f_SEp - f_SWp
        cc  = f_NWp - f_SWp
        dd  = f_NEp + f_SWp - f_NWp - f_SEp
        phi = phi + ((bb*cc - aa*dd) * log(abs(1.0_dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
        
      else
        write(io_unit_err,*) 'determine_grounded_fractions_CISM_quads - calc_fraction_above_zero - ERROR: unknown scenario [', scen, ']!'
        error stop 1
      end if
      
    end if
    
    if (phi < -0.01_wp .OR. phi > 1.01_wp .OR. phi /= phi) then
      write(io_unit_err,*) 'calc_fraction_above_zero - ERROR: phi = ', phi
      write(io_unit_err,*) 'scen = ', scen
      write(io_unit_err,*) 'f = [', f_NWp, ',', f_NEp, ',', f_SWp, ',', f_SEp, ']'
      write(io_unit_err,*) 'aa = ', aa, ', bb = ', bb, ', cc = ', cc, ', dd = ', dd, ', f1 = ', f1, ',f2 = ', f2
      error stop 1
    end if
    
    phi = MAX( 0.0_wp, MIN( 1.0_wp, phi))
    
    return

  end subroutine calc_fraction_above_zero

  subroutine rotate_quad_until_match( f_NW, f_NE, f_SW, f_SE, scen)
    ! Rotate the four corners until one of the four possible scenarios is found.
    ! 1: SW grounded, rest floating
    ! 2: SW floating, rest grounded
    ! 3: south grounded, north floating
    ! 4: SW & NE grounded, SE & NW floating
    
    implicit none
    
    ! In/output variables:
    real(dp), intent(INOUT) :: f_NW, f_NE, f_SW, f_SE
    integer,  intent(OUT)   :: scen
    
    ! Local variables:
    logical :: found_match
    integer :: nit
    
    found_match = .FALSE.
    scen        = 0
    nit         = 0
    
    do while (.not. found_match)
      
      nit = nit+1
      
      call rotate_quad( f_NW, f_NE, f_SW, f_SE)
      
      if     (f_SW < 0.0_wp .AND. f_SE > 0.0_wp .AND. f_NE > 0.0_wp .AND. f_NW > 0.0_wp) then
        ! 1: SW grounded, rest floating
        scen = 1
        found_match = .TRUE.
      else if (f_SW > 0.0_wp .AND. f_SE < 0.0_wp .AND. f_NE < 0.0_wp .AND. f_NW < 0.0_wp) then
        ! 2: SW floating, rest grounded
        scen = 2
        found_match = .TRUE.
      else if (f_SW < 0.0_wp .AND. f_SE < 0.0_wp .AND. f_NE > 0.0_wp .AND. f_NW > 0.0_wp) then
        ! 3: south grounded, north floating
        scen = 3
        found_match = .TRUE.
      else if (f_SW < 0.0_wp .AND. f_SE > 0.0_wp .AND. f_NE < 0.0_wp .AND. f_NW > 0.0_wp) then
        ! 4: SW & NE grounded, SE & NW floating
        scen = 4
        found_match = .TRUE.
      end if
      
      if (nit > 4) then
        write(io_unit_err,*) 
        write(io_unit_err,*) 'determine_grounded_fractions_CISM_quads - rotate_quad_until_match - ERROR: couldnt find matching scenario!'
        write(io_unit_err,*) 'f_SW, f_SE, f_NE, f_NW: ', f_SW, f_SE, f_NE, f_NW
        error stop 1
      end if
      
    end do
    
    return 

  end subroutine rotate_quad_until_match

  subroutine rotate_quad( f_NW, f_NE, f_SW, f_SE)
    ! Rotate the four corners anticlockwise by 90 degrees
    
    implicit none
    
    ! In/output variables:
    real(dp), intent(INOUT) :: f_NW, f_NE, f_SW, f_SE
    
    ! Local variables:
    real(dp) :: fvals(4)
    
    fvals = [f_NW,f_NE,f_SE,f_SW]
    f_NW = fvals( 2)
    f_NE = fvals( 3)
    f_SE = fvals( 4)
    f_SW = fvals( 1)
    
    return 

  end subroutine rotate_quad


    elemental function cdf(x,mu,sigma,inv) result(F)
        ! Solve for cumulative probability below (cdf)
        ! or above (cdf(inv=TRUE)) the value x
        ! given a normal distribution N(mu,sigma)

        ! See equation for F(x) and Q(x) in, e.g.:
        ! https://en.wikipedia.org/wiki/Normal_distribution

        implicit none

        real(wp), intent(IN)  :: x
        real(wp), intent(IN)  :: mu
        real(wp), intent(IN)  :: sigma
        logical,  intent(IN), optional :: inv
        real(wp) :: F
        
        real(wp), parameter :: sqrt2 = sqrt(2.0_wp)

        ! Calculate CDF at value x
        F = 0.5_wp*( 1.0_wp + error_function((x-mu)/(sqrt2*sigma)) )

        if (present(inv)) then 
            if (inv) then 
                ! Calculate inverse CDF at value x 
                F = 1.0_wp - F 
            end if 
        end if 

        return

    end function cdf

    elemental function error_function(X) result(ERR)
        ! Purpose: Compute error function erf(x)
        ! Input:   x   --- Argument of erf(x)
        ! Output:  ERR --- erf(x)
        
        ! Note: also separately defined in thermodynamics.f90

        implicit none 

        real(wp), intent(IN)  :: X
        real(wp) :: ERR
        
        ! Local variables:
        real(wp)              :: EPS
        real(wp)              :: X2
        real(wp)              :: ER
        real(wp)              :: R
        real(wp)              :: C0
        integer                 :: k
        
        EPS = 1.0e-15
        X2  = X * X
        if (abs(X) < 3.5) then
            ER = 1.0
            R  = 1.0
            do k = 1, 50
                R  = R * X2 / (real(k, wp) + 0.5)
                ER = ER+R
                if(abs(R) < abs(ER) * EPS) then
                    C0  = 2.0 / sqrt(pi) * X * exp(-X2)
                    ERR = C0 * ER
                    EXIT
                end if
            end do
        else
            ER = 1.0
            R  = 1.0
            do k = 1, 12
                R  = -R * (real(k, wp) - 0.5) / X2
                ER = ER + R
                C0  = EXP(-X2) / (abs(X) * sqrt(pi))
                ERR = 1.0 - C0 * ER
                if(X < 0.0) ERR = -ERR
            end do
        end if

        return

    end function error_function

end module topography


