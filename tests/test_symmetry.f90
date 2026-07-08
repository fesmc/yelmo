program test_symmetry
    ! Offline symmetry-regression driver.
    !
    ! Usage:  test_symmetry.x <file.nc> [varname]
    !
    ! Reads the LAST time slice of a 2D field from a NetCDF file and reports its
    ! reflection-symmetry metrics. If no variable name is given, tries "lithk"
    ! (CalvingMIP convention) then "H_ice" (yelmo convention).
    !
    ! Exits with status 1 on symmetry failure (for use as a CI gate).

    use ncio
    use yelmo_defs,     only : wp
    use yelmo_symmetry

    implicit none

    character(len=1024) :: filename
    character(len=256)  :: varname
    integer             :: narg
    integer, allocatable :: dims(:)
    real(wp), allocatable :: H(:,:)
    integer :: nx, ny, nt
    type(sym_metrics_class) :: met
    logical  :: all_pass
    real(wp), parameter :: tol = 1.0e-3_wp

    narg = command_argument_count()
    if (narg < 1) then
        write(*,*) "Usage: test_symmetry.x <file.nc> [varname]"
        stop 1
    end if

    call get_command_argument(1, filename)

    if (narg >= 2) then
        call get_command_argument(2, varname)
    else
        ! Default variable: try lithk (CalvingMIP), else H_ice (yelmo)
        if (nc_exists_var(trim(filename),"lithk")) then
            varname = "lithk"
        else
            varname = "H_ice"
        end if
    end if

    if (.not. nc_exists_var(trim(filename), trim(varname))) then
        write(*,*) "test_symmetry:: Error: variable not found in file: "//trim(varname)
        write(*,*) "  file: "//trim(filename)
        stop 1
    end if

    ! Determine variable dimensions (assumed order: x, y, [time])
    call nc_dims(trim(filename), trim(varname), dims=dims)
    nx = dims(1)
    ny = dims(2)
    if (size(dims) >= 3) then
        nt = dims(size(dims))
    else
        nt = 1
    end if

    allocate(H(nx,ny))

    ! Read the last time slice (or the whole field if no time dimension)
    if (size(dims) >= 3) then
        call nc_read(trim(filename), trim(varname), H, start=[1,1,nt], count=[nx,ny,1])
    else
        call nc_read(trim(filename), trim(varname), H)
    end if

    call calc_symmetry_metrics(H, tol, met, all_pass)
    call report_symmetry("offline: "//trim(varname)//" (last slice)", met, tol, all_pass)

    if (.not. all_pass) then
        write(*,*) "test_symmetry:: symmetry regression check FAILED."
        stop 1
    end if

end program test_symmetry
