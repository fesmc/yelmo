module yelmo_grid
    ! This module holds routines to manage the gridded coordinates of a given Yelmo domain.
    !
    ! The horizontal grid is represented by the fesm-utils/coords grid_class
    ! (dom%grd). This module provides yelmo-facing entry points that build a
    ! grid_class from a predefined name, a netcdf file, explicit options or
    ! explicit axes, all delegating the heavy lifting to coords grid_init.

    use yelmo_defs, only : sp, dp, wp, pi
    use yelmo_tools, only : get_region_indices

    use coords, only : grid_class, grid_init

    use ncio

    implicit none

    interface yelmo_init_grid
        module procedure yelmo_init_grid_fromname
        module procedure yelmo_init_grid_fromfile
        module procedure yelmo_init_grid_fromopt
        module procedure yelmo_init_grid_fromaxes
    end interface

    private

    public :: calc_zeta

    public :: yelmo_init_grid
    public :: yelmo_init_grid_fromfile
    public :: yelmo_init_grid_fromname
    public :: yelmo_init_grid_fromaxes
    public :: yelmo_init_grid_fromopt
    public :: yelmo_grid_write

contains

    subroutine calc_zeta(zeta_aa,zeta_ac,nz_ac,nz_aa,zeta_scale,zeta_exp)
        ! Calculate the vertical axis cell-centers first (aa-nodes),
        ! including a cell centered at the base and the surface.
        ! Then calculate the vertical axis cell-edges (ac-nodes).
        ! The base and surface ac-nodes coincide with the cell centers.
        ! There is one more cell-edge than cell-centers (nz_ac=nz_aa+1)

        implicit none

        real(wp), allocatable, intent(INOUT) :: zeta_aa(:)
        real(wp), allocatable, intent(INOUT) :: zeta_ac(:)
        integer,                 intent(OUT)   :: nz_ac
        integer,      intent(IN)   :: nz_aa
        character(*), intent(IN)   :: zeta_scale
        real(wp),   intent(IN)   :: zeta_exp

        ! Local variables
        integer :: k
        real(wp), allocatable :: tmp(:)

        ! Define size of zeta ac-node vector
        nz_ac = nz_aa + 1

        ! First allocate arrays
        if (allocated(zeta_aa)) deallocate(zeta_aa)
        if (allocated(zeta_ac)) deallocate(zeta_ac)
        allocate(zeta_aa(nz_aa))
        allocate(zeta_ac(nz_ac))

        ! Initially define a linear zeta scale
        ! Base = 0.0, Surface = 1.0
        do k = 1, nz_aa
            zeta_aa(k) = 0.0 + 1.0*(k-1)/real(nz_aa-1)
        end do

        ! Scale zeta to produce different resolution through column if desired
        ! zeta_scale = ["linear","exp","wave"]
        select case(trim(zeta_scale))

            case("exp")
                ! Increase resolution at the base

                zeta_aa = zeta_aa**(zeta_exp)

            case("exp-inv")
                ! Increase resolution at the surface

                zeta_aa = 1.0_wp - zeta_aa**(zeta_exp)

                ! Reverse order
                allocate(tmp(nz_aa))
                tmp = zeta_aa
                do k = 1, nz_aa
                    zeta_aa(k) = tmp(nz_aa-k+1)
                end do

            case("tanh")
                ! Increase resolution at base and surface

                zeta_aa = tanh(1.0*pi*(zeta_aa-0.5))
                zeta_aa = zeta_aa - minval(zeta_aa)
                zeta_aa = zeta_aa / maxval(zeta_aa)

            case DEFAULT
            ! Do nothing, scale should be linear as defined above

        end select

        ! Get zeta_ac (between zeta_aa values, as well as at the base and surface)
        zeta_ac(1) = 0.0
        do k = 2, nz_ac-1
            zeta_ac(k) = 0.5 * (zeta_aa(k-1)+zeta_aa(k))
        end do
        zeta_ac(nz_ac) = 1.0

        return

    end subroutine calc_zeta

    subroutine yelmo_init_grid_fromname(grd,grid_name)
        ! Build a predefined yelmo domain grid by name. The catalog stores each
        ! domain's geometry (corner/spacing/size, in km) and its map projection;
        ! the grid_class is built in meters by coords grid_init.

        implicit none

        type(grid_class), intent(INOUT) :: grd
        character(len=*), intent(IN)    :: grid_name

        ! Local variables
        real(dp) :: x0, y0, dx, dy
        real(dp) :: lambda, phi, alpha
        integer  :: nx, ny
        logical  :: center
        character(len=56) :: mtype

        ! All predefined domains use a WGS84 polar-stereographic projection.
        ! Defaults are the ESPG-3413 North projection (also used for the
        ! Eurasia and Greenland sub-domains); Antarctica overrides these below.
        mtype  = "polar_stereographic"
        lambda = -45.0_dp
        phi    =  70.0_dp
        alpha  =  20.0_dp
        x0     =   0.0_dp
        y0     =   0.0_dp
        center = .FALSE.

        select case(trim(grid_name))

            ! NORTH DOMAINS =======================
            ! ESPG-3413 polar stereographic (lambda=-45, phi=70). The Eurasia
            ! and Greenland sub-domains use the same projection for consistency.

            case("NH-40KM")
                x0=-4900.0_dp; y0=5400.0_dp; dx=40.0_dp; dy=40.0_dp; nx=221; ny=221
            case("NH-20KM")
                x0=-4900.0_dp; y0=5400.0_dp; dx=20.0_dp; dy=20.0_dp; nx=441; ny=441
            case("NH-10KM")
                x0=-4900.0_dp; y0=5400.0_dp; dx=10.0_dp; dy=10.0_dp; nx=881; ny=881
            case("NH-5KM")
                x0=-4900.0_dp; y0=5400.0_dp; dx=5.0_dp;  dy=5.0_dp;  nx=1761; ny=1761

            ! EURASIA DOMAINS =======================

            case("EIS-40KM")
                x0=380.0_dp; y0=-5000.0_dp; dx=40.0_dp; dy=40.0_dp; nx=89;  ny=161
            case("EIS-20KM")
                x0=380.0_dp; y0=-5000.0_dp; dx=20.0_dp; dy=20.0_dp; nx=177; ny=321
            case("EIS-10KM")
                x0=380.0_dp; y0=-5000.0_dp; dx=10.0_dp; dy=10.0_dp; nx=353; ny=641
            case("EIS-5KM")
                x0=380.0_dp; y0=-5000.0_dp; dx=5.0_dp;  dy=5.0_dp;  nx=705; ny=1281

            ! GREENLAND DOMAINS =======================

            case("GRL-40KM")
                x0=-720.0_dp; y0=-3450.0_dp; dx=40.0_dp; dy=40.0_dp; nx=43;   ny=73
            case("GRL-20KM")
                x0=-720.0_dp; y0=-3450.0_dp; dx=20.0_dp; dy=20.0_dp; nx=85;   ny=145
            case("GRL-10KM")
                x0=-720.0_dp; y0=-3450.0_dp; dx=10.0_dp; dy=10.0_dp; nx=169;  ny=289
            case("GRL-5KM")
                x0=-720.0_dp; y0=-3450.0_dp; dx=5.0_dp;  dy=5.0_dp;  nx=337;  ny=577
            case("GRL-2KM")
                x0=-720.0_dp; y0=-3450.0_dp; dx=2.0_dp;  dy=2.0_dp;  nx=841;  ny=1441
            case("GRL-1KM")
                x0=-720.0_dp; y0=-3450.0_dp; dx=1.0_dp;  dy=1.0_dp;  nx=1681; ny=2881

            ! ANTARCTICA DOMAINS =======================
            ! EPSG-3031 polar stereographic (lambda=0, phi=-71), centered on
            ! the projection origin (no explicit corner).

            case("ANT-80KM")
                center=.TRUE.; dx=80.0_dp; dy=80.0_dp; nx=79;   ny=74
                lambda=0.0_dp; phi=-71.0_dp; alpha=19.0_dp
            case("ANT-40KM")
                center=.TRUE.; dx=40.0_dp; dy=40.0_dp; nx=157;  ny=147
                lambda=0.0_dp; phi=-71.0_dp; alpha=19.0_dp
            case("ANT-20KM")
                center=.TRUE.; dx=20.0_dp; dy=20.0_dp; nx=313;  ny=293
                lambda=0.0_dp; phi=-71.0_dp; alpha=19.0_dp
            case("ANT-10KM")
                center=.TRUE.; dx=10.0_dp; dy=10.0_dp; nx=625;  ny=585
                lambda=0.0_dp; phi=-71.0_dp; alpha=19.0_dp
            case("ANT-5KM")
                center=.TRUE.; dx=5.0_dp;  dy=5.0_dp;  nx=1249; ny=1169
                lambda=0.0_dp; phi=-71.0_dp; alpha=19.0_dp
            case("ANT-1KM")
                center=.TRUE.; dx=1.0_dp;  dy=1.0_dp;  nx=6241; ny=5841
                lambda=0.0_dp; phi=-71.0_dp; alpha=19.0_dp

            case DEFAULT
                write(*,*) "yelmo_init_grid_fromname:: error: grid name not recognized: "//trim(grid_name)
                stop

        end select

        ! Convert catalog values from km to meters (yelmo grid geometry is in m)
        x0 = x0*1.d3
        y0 = y0*1.d3
        dx = dx*1.d3
        dy = dy*1.d3

        if (center) then
            call grid_init(grd,name=trim(grid_name),mtype=trim(mtype),units="meters", &
                           planet="WGS84",lon180=.TRUE.,dx=dx,nx=nx,dy=dy,ny=ny, &
                           lambda=lambda,phi=phi,alpha=alpha)
        else
            call grid_init(grd,name=trim(grid_name),mtype=trim(mtype),units="meters", &
                           planet="WGS84",lon180=.TRUE.,x0=x0,dx=dx,nx=nx,y0=y0,dy=dy,ny=ny, &
                           lambda=lambda,phi=phi,alpha=alpha)
        end if

        return

    end subroutine yelmo_init_grid_fromname

    subroutine yelmo_init_grid_fromfile(grd,filename,grid_name)
        ! Build a yelmo domain grid from a netcdf grid file. coords grid_init
        ! detects the netcdf file by extension, reads the axes (xc/yc) and the
        ! CF crs projection metadata, and computes lon/lat and cell areas.
        ! The grid geometry is returned in meters.

        implicit none

        type(grid_class), intent(INOUT) :: grd
        character(len=*), intent(IN)    :: filename
        character(len=*), intent(IN)    :: grid_name

        call grid_init(grd,filename)

        ! Use the yelmo grid_name as the grid identifier
        grd%name = trim(grid_name)

        return

    end subroutine yelmo_init_grid_fromfile

    subroutine yelmo_init_grid_fromopt(grd,grid_name,units,x0,dx,nx,y0,dy,ny,lon,lat,area)
        ! Build a cartesian grid from explicit options (used by the idealized
        ! test/benchmark drivers). Inputs are in `units` (km or m); the grid is
        ! built in meters. Optional lon/lat/area override the coords values.

        implicit none

        type(grid_class), intent(INOUT) :: grd
        character(len=*), intent(IN)    :: grid_name
        character(len=*), intent(IN)    :: units

        real(wp), intent(IN), optional :: x0
        real(wp), intent(IN) :: dx
        integer,  intent(IN) :: nx
        real(wp), intent(IN), optional :: y0
        real(wp), intent(IN) :: dy
        integer,  intent(IN) :: ny
        real(wp), intent(IN), optional :: lon(:,:)
        real(wp), intent(IN), optional :: lat(:,:)
        real(wp), intent(IN), optional :: area(:,:)

        ! Local variables
        real(dp) :: conv, x0d, y0d, dxd, dyd

        ! Determine conversion from input units to meters
        select case(trim(units))
            case("kilometers","km")
                conv = 1.d3
            case("meters","m")
                conv = 1.d0
            case DEFAULT
                write(*,*) "yelmo_init_grid_fromopt:: Error: units of input grid parameters &
                            &must be one of: 'kilometers', 'km', 'meters', 'm'."
                write(*,*) "units: ", trim(units)
                stop
        end select

        dxd = dx*conv
        dyd = dy*conv

        if (present(x0) .and. present(y0)) then
            x0d = x0*conv
            y0d = y0*conv
            call grid_init(grd,name=trim(grid_name),mtype="cartesian",units="meters", &
                           x0=x0d,dx=dxd,nx=nx,y0=y0d,dy=dyd,ny=ny)
        else
            call grid_init(grd,name=trim(grid_name),mtype="cartesian",units="meters", &
                           dx=dxd,nx=nx,dy=dyd,ny=ny)
        end if

        ! Optional overrides (values assumed to be in yelmo internal units: [m], [m^2])
        if (present(lon))  grd%lon  = real(lon,dp)
        if (present(lat))  grd%lat  = real(lat,dp)
        if (present(area)) grd%area = real(area,dp)

        return

    end subroutine yelmo_init_grid_fromopt

    subroutine yelmo_init_grid_fromaxes(grd,grid_name,xc,yc,lon,lat,area)
        ! Build a cartesian grid from explicit axis vectors (used by the C API
        ! coupling interface). Axes are assumed to be in meters. Optional
        ! lon/lat/area override the coords values.

        implicit none

        type(grid_class), intent(INOUT) :: grd
        character(len=*), intent(IN)    :: grid_name
        real(wp),         intent(IN)    :: xc(:)
        real(wp),         intent(IN)    :: yc(:)
        real(wp),         intent(IN), optional :: lon(:,:)
        real(wp),         intent(IN), optional :: lat(:,:)
        real(wp),         intent(IN), optional :: area(:,:)

        call grid_init(grd,name=trim(grid_name),mtype="cartesian",units="meters", &
                       x=real(xc,dp),y=real(yc,dp))

        ! Optional overrides (values assumed to be in yelmo internal units: [m], [m^2])
        if (present(lon))  grd%lon  = real(lon,dp)
        if (present(lat))  grd%lat  = real(lat,dp)
        if (present(area)) grd%area = real(area,dp)

        return

    end subroutine yelmo_init_grid_fromaxes

    subroutine yelmo_grid_write(grid, fnm, domain, grid_name, create, irange, jrange)
        ! Write grid info to netcdf file respecting coordinate conventions
        ! for easier interpretation/plotting etc.
        !
        ! Note: this is a yelmo-specific writer (rather than coords grid_write)
        ! because it (a) tags the file with domain/grid_name global attributes
        ! used by the restart machinery, (b) writes km/km^2 for plotting, and
        ! (c) supports regional subsetting via irange/jrange.

        implicit none
        type(grid_class), intent(IN) :: grid
        character(len=*),  intent(IN) :: fnm
        character(len=*),  intent(IN) :: domain
        character(len=*),  intent(IN) :: grid_name
        logical,           intent(IN) :: create
        integer,           intent(IN), optional :: irange(2)
        integer,           intent(IN), optional :: jrange(2)

        ! Local variables
        character(len=16) :: xnm
        character(len=16) :: ynm
        character(len=16) :: grid_mapping_name
        integer :: i1, i2, j1, j2

        xnm = "xc"
        ynm = "yc"
        grid_mapping_name = "crs"

        call get_region_indices(i1,i2,j1,j2,grid%G%nx,grid%G%ny,irange,jrange)

        ! Create the netcdf file if desired
        if (create) then

            ! Create the empty netcdf file
            call nc_create(fnm)

            ! Write some general attributes that can be useful (e.g., for interpolation etc)
            call nc_write_attr(fnm, "domain",    trim(domain))
            call nc_write_attr(fnm, "grid_name", trim(grid_name))


            ! Add grid axis variables to netcdf file
            call nc_write_dim(fnm,xnm,x=grid%G%x(i1:i2)*1d-3,units="km")
            call nc_write_dim(fnm,ynm,x=grid%G%y(j1:j2)*1d-3,units="km")

            if (grid%cs%is_projection) then
                call nc_write_attr(fnm,xnm,"standard_name","projection_x_coordinate")
                call nc_write_attr(fnm,ynm,"standard_name","projection_y_coordinate")
            end if

        end if

        ! Add projection information if needed
        if (grid%cs%is_projection) then
            call nc_write_map(fnm, grid%cs%mtype, grid%cs%proj%lambda, phi=grid%cs%proj%phi, &
                alpha=grid%cs%proj%alpha, x_e=grid%cs%proj%x_e, y_n=grid%cs%proj%y_n, &
                is_sphere=grid%cs%planet%is_sphere, semi_major_axis=grid%cs%planet%a, &
                inverse_flattening=grid%cs%planet%inverse_flattening)
        end if

        call nc_write(fnm, "x2D", grid%x(i1:i2, j1:j2)*1d-3, dim1=xnm, dim2=ynm, &
            units="km", grid_mapping=grid_mapping_name)
        call nc_write(fnm, "y2D", grid%y(i1:i2, j1:j2)*1d-3, dim1=xnm, dim2=ynm, &
            units="km", grid_mapping=grid_mapping_name)

        if (grid%cs%is_projection) then
            call nc_write(fnm, "lon2D", grid%lon(i1:i2, j1:j2), dim1=xnm, dim2=ynm, &
                grid_mapping=grid_mapping_name)
            call nc_write_attr(fnm, "lon2D", "units", "degrees_east")
            call nc_write(fnm, "lat2D", grid%lat(i1:i2, j1:j2), dim1=xnm,dim2=ynm, &
                grid_mapping=grid_mapping_name)
            call nc_write_attr(fnm,"lat2D","units","degrees_north")
        end if

        call nc_write(fnm,"area",  grid%area(i1:i2, j1:j2)*1d-6, dim1=xnm,dim2=ynm, &
            grid_mapping=grid_mapping_name,units="km^2")
        if (grid%cs%is_projection) call nc_write_attr(fnm,"area","coordinates","lat2D lon2D")

        return

    end subroutine yelmo_grid_write

end module yelmo_grid
