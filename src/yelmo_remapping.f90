module yelmo_remapping
    ! Map 2D and 3D fields from another source (another model, another
    ! resolution, data) onto the Yelmo grid.
    !
    !   yelmo_remap       in memory: horizontal remapping with a coords map
    !                     (2D), and vertical interpolation onto a Yelmo axis
    !                     followed by horizontal remapping (3D)
    !   yelmo_load_map    build the map from a file's xc/yc axes onto the
    !                     Yelmo grid (model projection)
    !   yelmo_read_remap  read a 2D or 3D field from a NetCDF file and map it
    !                     onto the Yelmo grid (remapped only if the file's
    !                     grid differs)
    !
    ! Vertical coordinates are normalized heights zeta (0 at the base, 1 at
    ! the top of the layer), as Yelmo's zeta_aa, zeta_ac and the bedrock axes.

    use yelmo_defs, only : sp, dp, wp, io_unit_err
    use ncio
    use coords,  only : grid_class, grid_init
    use mapping, only : map_class, map_read, map_init, map_field

    implicit none

    interface yelmo_remap
        module procedure yelmo_remap_2D
        module procedure yelmo_remap_3D
    end interface

    interface yelmo_read_remap
        module procedure yelmo_read_remap_2D
        module procedure yelmo_read_remap_3D
    end interface

    private
    public :: yelmo_remap
    public :: yelmo_load_map
    public :: yelmo_read_remap

contains

    subroutine yelmo_remap_2D(var,var_in,mp,name)
        ! Map the 2D field var_in from the source grid onto var on the
        ! target grid of the map mp.

        implicit none

        real(wp),          intent(INOUT) :: var(:,:)        ! Target grid
        real(wp),          intent(IN)    :: var_in(:,:)     ! Source grid
        type(map_class),   intent(IN)    :: mp
        character(len=*),  intent(IN), optional :: name     ! Field name (messages only)

        character(len=56) :: vnm

        vnm = "var"
        if (present(name)) vnm = trim(name)

        call map_field(mp,vnm,var_in,var)

        return

    end subroutine yelmo_remap_2D

    subroutine yelmo_remap_3D(var,zeta,var_in,zeta_in,mp,name)
        ! Map the 3D field var_in on levels zeta_in onto var on levels zeta:
        ! linear interpolation in each column (constant beyond the end
        ! levels), then, if mp is given, horizontal remapping of each level.
        ! Without mp, var_in must already be on the target horizontal grid.

        implicit none

        real(wp),          intent(INOUT) :: var(:,:,:)      ! Target grid, levels zeta
        real(wp),          intent(IN)    :: zeta(:)         ! Target levels (ascending)
        real(wp),          intent(IN)    :: var_in(:,:,:)   ! Source grid, levels zeta_in
        real(wp),          intent(IN)    :: zeta_in(:)      ! Source levels (ascending)
        type(map_class),   intent(IN), optional :: mp
        character(len=*),  intent(IN), optional :: name     ! Field name (messages only)

        ! Local variables
        integer :: i, j, k, nx_in, ny_in, nz
        real(wp), allocatable :: tmp(:,:,:)
        character(len=56) :: vnm

        vnm = "var"
        if (present(name)) vnm = trim(name)

        nx_in = size(var_in,1)
        ny_in = size(var_in,2)
        nz    = size(zeta)

        if (size(var,3) .ne. nz .or. size(var_in,3) .ne. size(zeta_in)) then
            write(io_unit_err,*) "yelmo_remap_3D:: Error: the vertical sizes of the fields &
                                 &do not match their axes: ", trim(vnm)
            write(io_unit_err,*) "size(var,3), size(zeta)       = ", size(var,3), nz
            write(io_unit_err,*) "size(var_in,3), size(zeta_in) = ", size(var_in,3), size(zeta_in)
            error stop 1
        end if

        if (any(zeta_in(2:) .le. zeta_in(:size(zeta_in)-1))) then
            write(io_unit_err,*) "yelmo_remap_3D:: Error: zeta_in must be strictly ascending: ", trim(vnm)
            error stop 1
        end if

        ! Vertical interpolation on the source grid
        allocate(tmp(nx_in,ny_in,nz))
        do j = 1, ny_in
        do i = 1, nx_in
            tmp(i,j,:) = interp_column(zeta_in,var_in(i,j,:),zeta)
        end do
        end do

        ! Horizontal remapping of each level
        if (present(mp)) then
            do k = 1, nz
                call map_field(mp,vnm,tmp(:,:,k),var(:,:,k))
            end do
        else
            if (size(var,1) .ne. nx_in .or. size(var,2) .ne. ny_in) then
                write(io_unit_err,*) "yelmo_remap_3D:: Error: without a map, var_in must be &
                                     &on the target horizontal grid: ", trim(vnm)
                error stop 1
            end if
            var = tmp
        end if

        return

    end subroutine yelmo_remap_3D

    function interp_column(z_in,v_in,z) result(v)
        ! Linear interpolation of v_in(z_in) onto z, constant beyond the end
        ! levels. z_in must be strictly ascending.

        implicit none

        real(wp), intent(IN) :: z_in(:)
        real(wp), intent(IN) :: v_in(:)
        real(wp), intent(IN) :: z(:)
        real(wp) :: v(size(z))

        integer  :: k, m, n
        real(wp) :: f

        n = size(z_in)
        m = 1
        do k = 1, size(z)
            if (z(k) .le. z_in(1)) then
                v(k) = v_in(1)
            else if (z(k) .ge. z_in(n)) then
                v(k) = v_in(n)
            else
                do while (z_in(m+1) .lt. z(k))
                    m = m+1
                end do
                f    = (z(k)-z_in(m)) / (z_in(m+1)-z_in(m))
                v(k) = (1.0_wp-f)*v_in(m) + f*v_in(m+1)
            end if
        end do

        return

    end function interp_column

    subroutine yelmo_load_map(mp,grd,filename,src_grid_name,method,gen)
        ! Build the map from the grid of a NetCDF file onto the Yelmo grid grd.
        !   gen = "coords" (default): generate the weights in-package. The
        !         source grid is built from the file's xc/yc axes ([km] or
        !         [m]) and the projection of the Yelmo grid, so source and
        !         target must share the model projection.
        !   gen = "cdo": load a pre-generated cdo SCRIP map from maps/
        !         (src_grid_name => grd%name).
        ! method: mapping method of map_init (default "con", conservative).

        implicit none

        type(map_class),   intent(OUT) :: mp
        type(grid_class),  intent(IN)  :: grd               ! Target (Yelmo) grid
        character(len=*),  intent(IN)  :: filename          ! File holding the source axes
        character(len=*),  intent(IN)  :: src_grid_name     ! Source grid name
        character(len=*),  intent(IN), optional :: method
        character(len=*),  intent(IN), optional :: gen

        ! Local variables
        type(grid_class)      :: grid_src
        real(wp), allocatable :: xc_src(:), yc_src(:)
        character(len=56)     :: mtd, gn

        mtd = "con"
        if (present(method)) mtd = trim(method)
        gn  = "coords"
        if (present(gen)) gn = trim(gen)

        select case(trim(gn))

            case("cdo")

                call map_read(mp,src_grid_name,grd%name,"maps",mtd)
                write(*,*) "yelmo_load_map:: loaded "//trim(mtd)//" SCRIP map (cdo): " &
                            //trim(src_grid_name)//" => "//trim(grd%name)

            case("coords")

                call read_axes_m(filename,xc_src,yc_src)

                ! Source grid inherits the model grid's projection; only the axes differ
                call grid_init(grid_src,grd,name=src_grid_name,x=real(xc_src,dp),y=real(yc_src,dp))

                call map_init(mp,grid_src,grd,method=mtd,gen="coords")
                write(*,*) "yelmo_load_map:: generated "//trim(mtd)//" coords map: " &
                            //trim(src_grid_name)//" => "//trim(grd%name)

            case default

                write(io_unit_err,*) "yelmo_load_map:: Error: unknown gen '"//trim(gn)// &
                                     "'. Expected 'cdo' or 'coords'."
                error stop 1

        end select

        return

    end subroutine yelmo_load_map

    subroutine yelmo_read_remap_2D(var,grd,filename,varname,method)
        ! Read the 2D field varname from filename and map it onto var on the
        ! Yelmo grid grd. If the field has a time dimension, the last record
        ! is read.

        implicit none

        real(wp),          intent(INOUT) :: var(:,:)
        type(grid_class),  intent(IN)    :: grd
        character(len=*),  intent(IN)    :: filename
        character(len=*),  intent(IN)    :: varname
        character(len=*),  intent(IN), optional :: method   ! Mapping method (default "con")

        ! Local variables
        type(map_class)       :: mp
        real(wp), allocatable :: var_in(:,:)
        integer :: nx_in, ny_in

        nx_in = nc_size(filename,"xc")
        ny_in = nc_size(filename,"yc")
        allocate(var_in(nx_in,ny_in))
        call read_field(filename,varname,var_in)

        if (same_grid(filename,grd)) then
            var = var_in
        else
            call yelmo_load_map(mp,grd,filename,trim(grd%name)//"-src",method=method)
            call yelmo_remap(var,var_in,mp,name=varname)
        end if

        return

    end subroutine yelmo_read_remap_2D

    subroutine yelmo_read_remap_3D(var,grd,filename,varname,zeta,zeta_name,method)
        ! Read the 3D field varname on the levels zeta_name from filename and
        ! map it onto var on the Yelmo grid grd and the levels zeta. If the
        ! field has a time dimension, the last record is read.

        implicit none

        real(wp),          intent(INOUT) :: var(:,:,:)
        type(grid_class),  intent(IN)    :: grd
        character(len=*),  intent(IN)    :: filename
        character(len=*),  intent(IN)    :: varname
        real(wp),          intent(IN)    :: zeta(:)         ! Target levels (ascending)
        character(len=*),  intent(IN)    :: zeta_name       ! Name of the source levels in the file
        character(len=*),  intent(IN), optional :: method   ! Mapping method (default "con")

        ! Local variables
        type(map_class)       :: mp
        real(wp), allocatable :: var_in(:,:,:), zeta_in(:)
        integer :: nx_in, ny_in, nz_in

        nx_in = nc_size(filename,"xc")
        ny_in = nc_size(filename,"yc")
        nz_in = nc_size(filename,zeta_name)
        allocate(var_in(nx_in,ny_in,nz_in),zeta_in(nz_in))
        call nc_read(filename,zeta_name,zeta_in)
        call read_field(filename,varname,var_in)

        if (same_grid(filename,grd)) then
            call yelmo_remap(var,zeta,var_in,zeta_in,name=varname)
        else
            call yelmo_load_map(mp,grd,filename,trim(grd%name)//"-src",method=method)
            call yelmo_remap(var,zeta,var_in,zeta_in,mp,name=varname)
        end if

        return

    end subroutine yelmo_read_remap_3D

    subroutine read_field(filename,varname,var)
        ! Read a 2D or 3D field; with one extra (time) dimension, read the
        ! last record.

        implicit none

        character(len=*), intent(IN)  :: filename
        character(len=*), intent(IN)  :: varname
        real(wp),         intent(OUT) :: var(..)

        ! Local variables
        character(len=32), allocatable :: names(:)
        integer, allocatable :: dims(:)
        integer :: nd

        call nc_dims(filename,varname,names,dims)
        nd = size(dims)

        select rank(var)
            rank(2)
                if (nd .eq. 2) then
                    call nc_read(filename,varname,var)
                else if (nd .eq. 3) then
                    call nc_read(filename,varname,var,start=[1,1,dims(3)],count=[dims(1),dims(2),1])
                else
                    call error_rank(varname,nd,2)
                end if
            rank(3)
                if (nd .eq. 3) then
                    call nc_read(filename,varname,var)
                else if (nd .eq. 4) then
                    call nc_read(filename,varname,var,start=[1,1,1,dims(4)],count=[dims(1),dims(2),dims(3),1])
                else
                    call error_rank(varname,nd,3)
                end if
            rank default
                call error_rank(varname,nd,0)
        end select

        return

    end subroutine read_field

    subroutine error_rank(varname,nd,nd_expected)

        implicit none

        character(len=*), intent(IN) :: varname
        integer,          intent(IN) :: nd, nd_expected

        write(io_unit_err,*) "yelmo_read_remap:: Error: ", trim(varname), " has ", nd, &
                             " dimensions; expected ", nd_expected, " (plus an optional time dimension)."
        error stop 1

    end subroutine error_rank

    subroutine read_axes_m(filename,xc,yc)
        ! Read the xc/yc axes of a file in [m] (converted from [km] if needed).

        implicit none

        character(len=*),      intent(IN)  :: filename
        real(wp), allocatable, intent(OUT) :: xc(:), yc(:)

        character(len=56) :: units

        allocate(xc(nc_size(filename,"xc")),yc(nc_size(filename,"yc")))
        call nc_read(filename,"xc",xc)
        call nc_read(filename,"yc",yc)

        units = "m"
        if (nc_exists_attr(filename,"xc","units")) call nc_read_attr(filename,"xc","units",units)
        if (trim(units) .eq. "kilometers" .or. trim(units) .eq. "km") then
            xc = xc*1e3_wp
            yc = yc*1e3_wp
        end if

        return

    end subroutine read_axes_m

    logical function same_grid(filename,grd)
        ! True if the file's xc/yc axes equal the axes of grd (to 1 m).

        implicit none

        character(len=*), intent(IN) :: filename
        type(grid_class), intent(IN) :: grd

        real(wp), allocatable :: xc(:), yc(:)

        call read_axes_m(filename,xc,yc)

        same_grid = size(xc) .eq. grd%G%nx .and. size(yc) .eq. grd%G%ny
        if (same_grid) same_grid = maxval(abs(real(xc,dp)-grd%G%x)) .lt. 1.0_dp .and. &
                                   maxval(abs(real(yc,dp)-grd%G%y)) .lt. 1.0_dp

        return

    end function same_grid

end module yelmo_remapping
