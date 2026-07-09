program test_symmetry
    ! Offline symmetry-regression driver.
    !
    ! Usage:  test_symmetry.x <file.nc> [varname] [axis]
    !
    ! Reads the LAST time slice of a 2D field from a NetCDF file and reports its
    ! reflection-symmetry metrics. If no variable name is given, tries "lithk"
    ! (CalvingMIP convention) then "H_ice" (yelmo convention).
    !
    ! The optional axis argument restricts the pass/fail gate to one reflection:
    !   all (default) - require all computed reflections symmetric (fully symmetric case)
    !   lr            - left<->right only
    !   tb            - top<->bottom only  (e.g. N-S symmetric channel/trough setups)
    !   rot           - 180-degree rotation only
    ! All metrics are still reported; only the gate is narrowed.
    !
    ! Exits with status 1 on symmetry failure (for use as a CI gate).

    use ncio
    use yelmo_defs,     only : wp
    use yelmo_symmetry

    implicit none

    character(len=1024) :: filename
    character(len=256)  :: varname
    character(len=16)   :: axis
    integer             :: narg
    integer, allocatable :: dims(:)
    real(wp), allocatable :: H(:,:)
    integer :: nx, ny, nt
    type(sym_metrics_class) :: met
    logical  :: all_pass
    real(wp), parameter :: tol = 1.0e-3_wp

    narg = command_argument_count()
    if (narg < 1) then
        write(*,*) "Usage: test_symmetry.x <file.nc> [varname] [axis=all|lr|tb|rot]"
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

    if (narg >= 3) then
        call get_command_argument(3, axis)
    else
        axis = "all"
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

    ! Narrow the pass/fail gate to the requested reflection axis (metrics for all
    ! reflections are still reported by report_symmetry below).
    select case (trim(axis))
        case ("all")
            ! keep all_pass from calc_symmetry_metrics (all computed reflections)
        case ("lr")
            all_pass = (met%linf_lr  < tol)
        case ("tb")
            all_pass = (met%linf_tb  < tol)
        case ("rot")
            all_pass = (met%linf_rot < tol)
        case default
            write(*,*) "test_symmetry:: Error: unknown axis '"//trim(axis)//"' (use all|lr|tb|rot)"
            stop 1
    end select

    call report_symmetry("offline: "//trim(varname)//" ["//trim(axis)//"] (last slice)", met, tol, all_pass)

    if (.not. all_pass) then
        write(*,*) "test_symmetry:: symmetry regression check FAILED."
        stop 1
    end if

end program test_symmetry
