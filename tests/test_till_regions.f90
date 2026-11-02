program test_till_regions
    ! Regions or basins with their own till cf_ref (the field dyn%now%cf_ref;
    ! ydyn_till_init, set_cf_ref_ids, calc_cb_ref).
    !
    ! 1. ydyn_till_init sets cf_ref = ytill.cf_ref everywhere. set_cf_ref_ids:
    !    no codes leave cf_ref unchanged; a code selects the cells of exactly
    !    that id (not its subregions); other cells keep their value; where
    !    codes repeat, later entries win; a code not in the field changes nothing.
    ! 2. calc_cb_ref with a regional cf_ref field equals, in each region,
    !    calc_cb_ref with a constant field of the value of that region (cf_min,
    !    elevation and sediment scaling apply unchanged on top), for every
    !    scale_zb and scale_sed.
    ! 3. calc_ydyn_cb_ref (ytill parameters, bnd, cf_ref field) equals
    !    calc_cb_ref with that field.

    use yelmo_defs,     only : wp, ydyn_class, ybound_class
    use basal_dragging, only : set_cf_ref_ids, calc_cb_ref
    use yelmo_dynamics, only : calc_ydyn_cb_ref, ydyn_till_init

    implicit none

    integer,  parameter :: nx = 24, ny = 20, n_sd = 5
    real(wp), parameter :: cf_def = 0.1_wp, cf_min = 0.001_wp
    real(wp), parameter :: f_sed = 0.5_wp, sed_min = 5.0_wp, sed_max = 15.0_wp
    real(wp), parameter :: z0 = -100.0_wp, z1 = 100.0_wp

    ! Region codes of the test domain (deepest level): none, Northern
    ! Hemisphere only, North America, its subregions 1.1.1 and 1.1.2,
    ! Greenland, Iceland, Arctic Ocean
    integer,  parameter :: codes(8) = [0, 1, 101, 10101, 10102, 103, 105, 151]

    real(wp) :: regions(nx,ny), z_bed(nx,ny), z_bed_sd(nx,ny), z_sl(nx,ny), H_sed(nx,ny)
    real(wp) :: cf_ref(nx,ny), cf_exp(nx,ny), cb_a(nx,ny), cb_b(nx,ny), cb_c(nx,ny), r(nx,ny)
    integer  :: i, j, k, izb, ised, n_fail
    type(ydyn_class)   :: dyn
    type(ybound_class) :: bnd
    integer, allocatable :: seed(:)

    n_fail = 0

    ! Fixed seed, so that the random fields are the same in every run
    call random_seed(size=k)
    allocate(seed(k))
    seed = 20261031
    call random_seed(put=seed)

    call random_number(r)
    do j = 1, ny
    do i = 1, nx
        regions(i,j) = real(codes(1+mod(int(8.0_wp*r(i,j)),8)),wp)
    end do
    end do
    call random_number(z_bed);    z_bed    = 600.0_wp*(z_bed-0.6_wp)
    call random_number(z_bed_sd); z_bed_sd = 80.0_wp*z_bed_sd
    call random_number(H_sed);    H_sed    = 20.0_wp*H_sed
    z_sl = 0.0_wp

    ! === 1. ydyn_till_init, set_cf_ref_ids ==================================

    allocate(dyn%now%cf_ref(nx,ny))
    dyn%now%cf_ref = -1.0_wp
    dyn%par%till_cf_ref = cf_def
    call ydyn_till_init(dyn,bnd)
    call check("ydyn_till_init: cf_ref everywhere", all(dyn%now%cf_ref .eq. cf_def))

    cf_ref = dyn%now%cf_ref
    call set_cf_ref_ids(cf_ref,regions,[integer ::],[real(wp) ::])
    call check("no codes: cf_ref unchanged", all(cf_ref .eq. cf_def))

    ! Greenland: its cells only
    call set_cf_ref_ids(cf_ref,regions,[103],[0.05_wp])
    call check("Greenland (103): its cells only", &
               all(merge(0.05_wp,cf_def,nint(regions) .eq. 103) .eq. cf_ref))

    ! North America and its subregion 1.1.1: exact ids, 10102 is not selected
    ! by 101; the earlier value of Greenland is kept
    call set_cf_ref_ids(cf_ref,regions,[101, 10101],[0.2_wp, 0.3_wp])
    call check("North America (101) and 10101: exact ids, other cells kept", &
               all(pack(cf_ref,nint(regions) .eq. 101)   .eq. 0.2_wp) .and. &
               all(pack(cf_ref,nint(regions) .eq. 10101) .eq. 0.3_wp) .and. &
               all(pack(cf_ref,nint(regions) .eq. 10102) .eq. cf_def) .and. &
               all(pack(cf_ref,nint(regions) .eq. 103)   .eq. 0.05_wp) .and. &
               all(pack(cf_ref,nint(regions) .eq. 0)     .eq. cf_def) .and. &
               all(pack(cf_ref,nint(regions) .eq. 1)     .eq. cf_def) .and. &
               all(pack(cf_ref,nint(regions) .eq. 105)   .eq. cf_def) .and. &
               all(pack(cf_ref,nint(regions) .eq. 151)   .eq. cf_def))

    ! A repeated code: the later entry wins
    call set_cf_ref_ids(cf_ref,regions,[105, 105],[0.3_wp, 0.4_wp])
    call check("repeated code: later entry wins", all(pack(cf_ref,nint(regions) .eq. 105) .eq. 0.4_wp))

    ! A code not in the field: nothing changes
    cf_exp = cf_ref
    call set_cf_ref_ids(cf_ref,regions,[999],[0.7_wp])
    call check("code not in the field: cf_ref unchanged", all(cf_ref .eq. cf_exp))

    ! === 2. calc_cb_ref with a regional cf_ref field ========================

    cf_ref = cf_def
    call set_cf_ref_ids(cf_ref,regions,[103, 101],[0.05_wp, 0.2_wp])

    do izb = 0, 2
    do ised = 0, 3

        call calc_cb_ref(cb_b,z_bed,z_bed_sd,z_sl,H_sed,f_sed,sed_min,sed_max, &
                         cf_ref,cf_min,z0,z1,n_sd,izb,ised)

        ! Constant fields of the value of each region
        call calc_cb_ref(cb_a,z_bed,z_bed_sd,z_sl,H_sed,f_sed,sed_min,sed_max, &
                         spread(spread(cf_def,1,nx),2,ny),cf_min,z0,z1,n_sd,izb,ised)
        call calc_cb_ref(cb_c,z_bed,z_bed_sd,z_sl,H_sed,f_sed,sed_min,sed_max, &
                         spread(spread(0.05_wp,1,nx),2,ny),cf_min,z0,z1,n_sd,izb,ised)
        where (cf_ref .eq. 0.05_wp) cb_a = cb_c
        call calc_cb_ref(cb_c,z_bed,z_bed_sd,z_sl,H_sed,f_sed,sed_min,sed_max, &
                         spread(spread(0.2_wp,1,nx),2,ny),cf_min,z0,z1,n_sd,izb,ised)
        where (cf_ref .eq. 0.2_wp) cb_a = cb_c

        call check("regional field = constant field per region, "//label(izb,ised), all(cb_a .eq. cb_b))

    end do
    end do

    ! === 3. calc_ydyn_cb_ref =================================================

    dyn%par%till_f_sed     = f_sed
    dyn%par%till_sed_min   = sed_min
    dyn%par%till_sed_max   = sed_max
    dyn%par%till_z0        = z0
    dyn%par%till_z1        = z1
    dyn%par%till_n_sd      = n_sd
    dyn%par%till_cf_min    = cf_min
    dyn%par%till_cf_ref    = cf_def
    dyn%par%till_scale_zb  = 2
    dyn%par%till_scale_sed = 2
    bnd%regions  = regions
    bnd%z_bed    = z_bed
    bnd%z_bed_sd = z_bed_sd
    bnd%z_sl     = z_sl
    bnd%H_sed    = H_sed

    call calc_ydyn_cb_ref(cb_a,dyn%par,bnd,cf_ref)
    call calc_cb_ref(cb_b,z_bed,z_bed_sd,z_sl,H_sed,f_sed,sed_min,sed_max, &
                     cf_ref,cf_min,z0,z1,n_sd,2,2)
    call check("calc_ydyn_cb_ref = calc_cb_ref with the cf_ref field", all(cb_a .eq. cb_b))

    if (n_fail .eq. 0) then
        write(*,*) "test_till_regions: PASSED"
    else
        write(*,*) "test_till_regions: FAILED (", n_fail, " checks)"
        error stop 1
    end if

contains

    function label(izb,ised) result(str)
        integer, intent(IN) :: izb, ised
        character(len=24)   :: str
        write(str,'(a,i0,a,i0)') "scale_zb ", izb, ", scale_sed ", ised
    end function label

    subroutine check(lbl,ok)
        character(len=*), intent(IN) :: lbl
        logical,          intent(IN) :: ok
        if (ok) then
            write(*,'(a,a)') "  ok    ", lbl
        else
            write(*,'(a,a)') "  FAIL  ", lbl
            n_fail = n_fail + 1
        end if
    end subroutine check

end program test_till_regions
