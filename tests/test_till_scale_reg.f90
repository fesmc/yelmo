program test_till_scale_reg
    ! Elevation of the till scaling (ytill.scale_reg):
    !
    ! 1. calc_basin_median on a small field with exact medians: odd and even
    !    numbers of cells, a cell of a basin outside basin_mask (not counted,
    !    gets the median), cells without basin and a basin without cells in
    !    basin_mask (keep their own value); and on random fields against a
    !    reference (insertion sort).
    ! 2. scale_reg = 0, no regions: calc_ydyn_cb_ref is bit-identical to
    !    calc_cb_ref on the current z_bed (n_sd samples over z_bed_sd) with
    !    cf_ref everywhere, for scale_zb = 1, 2, scale_sed = 0..3, n_sd = 1, 10.
    ! 3. scale_reg = 1: ydyn_till_init sets z_bed_reg (nothing with scale_reg = 0);
    !    cb_ref is calc_cb_ref on z_bed_reg with n_sd = 1, does not depend on
    !    z_bed, z_bed_sd or n_sd, is one value per basin and lower in lower basins.
    ! 4. scale_reg = 1 with a region of its own cf_ref (set_cf_ref_ids):
    !    the linear law on the basin median runs from cf_min to the cell's
    !    cf_ref, so one basin has two values where the region cuts it.

    use yelmo_defs,      only : wp, ydyn_class, ybound_class
    use basal_dragging,  only : calc_cb_ref, set_cf_ref_ids, calc_basin_median
    use yelmo_dynamics,  only : calc_ydyn_cb_ref, ydyn_till_init

    implicit none

    integer, parameter :: nx = 40, ny = 30

    type(ydyn_class)   :: dyn
    type(ybound_class) :: bnd
    real(wp) :: z(6,5), bas(6,5), msk(6,5), zm(6,5), zexp(6,5)
    real(wp), allocatable :: zr(:,:), zr_ref(:,:), cb(:,:), cb_ref(:,:), cb2(:,:), r(:,:), cfr(:,:), lam(:,:)
    integer, parameter :: ids6(6) = [1,2,3,5,7,9]
    real(wp) :: zb_b(9), cb_b(9)
    integer  :: i, j, k, b, b2, n_fail, i_zb, i_sed, i_sd, n_sd
    integer, allocatable :: seed(:)
    logical  :: ok

    n_fail = 0

    call random_seed(size=k)
    allocate(seed(k))
    seed = 20261101
    call random_seed(put=seed)

    ! === 1. calc_basin_median ===============================================

    do k = 1, 5
        z(:,k) = 1000.0_wp + 10.0_wp*[1,2,3,4,5,6] + real(k,wp)
    end do
    bas = 0.0_wp
    msk = 0.0_wp

    ! Basin 1: 5 cells in basin_mask (median 50), one outside (9999, not counted)
    bas(1:3,1) = 1.0_wp;  msk(1:3,1) = 1.0_wp;  z(1:3,1) = [-500.0_wp, 100.0_wp, -200.0_wp]
    bas(1:2,2) = 1.0_wp;  msk(1:2,2) = 1.0_wp;  z(1:2,2) = [300.0_wp, 50.0_wp]
    bas(3,2)   = 1.0_wp;                         z(3,2)   = 9999.0_wp
    ! Basin 2: 4 cells (median (10+20)/2 = 15)
    bas(4:5,1) = 2.0_wp;  msk(4:5,1) = 1.0_wp;  z(4:5,1) = [10.0_wp, -30.0_wp]
    bas(4:5,2) = 2.0_wp;  msk(4:5,2) = 1.0_wp;  z(4:5,2) = [70.0_wp, 20.0_wp]
    ! Basin 4 (no basin 3): no cells in basin_mask
    bas(5:6,5) = 4.0_wp

    zexp = z
    where (bas .eq. 1.0_wp) zexp = 50.0_wp
    where (bas .eq. 2.0_wp) zexp = 15.0_wp

    call calc_basin_median(zm,z,bas,msk)
    call check("basin median, small field (odd, even, outside mask, no basin)", all(zm .eq. zexp))

    ! Random fields: 6 basins (ids 1, 2, 3, 5, 7, 9), random basin_mask
    allocate(zr(nx,ny), zr_ref(nx,ny), r(nx,ny))
    allocate(bnd%z_bed(nx,ny), bnd%z_bed_sd(nx,ny), bnd%z_sl(nx,ny), bnd%H_sed(nx,ny))
    allocate(bnd%z_bed_ref(nx,ny), bnd%basins(nx,ny), bnd%basin_mask(nx,ny), bnd%regions(nx,ny))
    bnd%regions = 0.0_wp

    call random_number(r)
    bnd%z_bed_ref = 3000.0_wp*(r - 0.5_wp)
    call random_number(r)
    do j = 1, ny
    do i = 1, nx
        bnd%basins(i,j) = real(ids6(1 + int(6.0_wp*r(i,j))),wp)
    end do
    end do
    call random_number(r)
    where (r .lt. 0.05_wp) bnd%basins = 0.0_wp
    call random_number(r)
    bnd%basin_mask = merge(1.0_wp, 0.0_wp, r .lt. 0.8_wp .and. bnd%basins .gt. 0.0_wp)

    call calc_basin_median(zr,bnd%z_bed_ref,bnd%basins,bnd%basin_mask)
    call median_reference(zr_ref,bnd%z_bed_ref,bnd%basins,bnd%basin_mask)
    call check("basin median, random fields against the reference", all(zr .eq. zr_ref))

    ! === 2. scale_reg = 0: bit-identical to the old call =====================

    call random_number(r)
    bnd%z_bed    = 3000.0_wp*(r - 0.5_wp)
    call random_number(r)
    bnd%z_bed_sd = 300.0_wp*r
    bnd%z_sl     = 0.0_wp
    call random_number(r)
    bnd%H_sed    = 20.0_wp*r

    dyn%par%till_f_sed   = 0.3_wp
    dyn%par%till_sed_min = 5.0_wp
    dyn%par%till_sed_max = 15.0_wp
    dyn%par%till_z0      = -400.0_wp
    dyn%par%till_z1      = 200.0_wp
    dyn%par%till_cf_min  = 0.05_wp
    dyn%par%till_cf_ref  = 0.5_wp

    allocate(cb(nx,ny), cb_ref(nx,ny), cb2(nx,ny), cfr(nx,ny), lam(nx,ny))
    cfr = dyn%par%till_cf_ref
    allocate(dyn%now%cf_ref(nx,ny))
    dyn%now%cf_ref = cfr
    allocate(dyn%now%z_bed_reg(nx,ny))
    dyn%now%z_bed_reg = -1.0e6_wp       ! not used with scale_reg = 0

    ok = .TRUE.
    do i_zb = 1, 2
    do i_sed = 0, 3
    do i_sd = 1, 2
        n_sd = merge(1, 10, i_sd .eq. 1)
        dyn%par%till_scale_zb  = i_zb
        dyn%par%till_scale_sed = i_sed
        dyn%par%till_n_sd      = n_sd
        dyn%par%till_scale_reg = 0

        call calc_cb_ref(cb_ref,bnd%z_bed,bnd%z_bed_sd,bnd%z_sl,bnd%H_sed,dyn%par%till_f_sed, &
                dyn%par%till_sed_min,dyn%par%till_sed_max,cfr,dyn%par%till_cf_min, &
                dyn%par%till_z0,dyn%par%till_z1,dyn%par%till_n_sd,dyn%par%till_scale_zb,dyn%par%till_scale_sed)
        call calc_ydyn_cb_ref(cb,dyn%par,bnd,dyn%now%cf_ref,dyn%now%z_bed_reg)

        ok = ok .and. all(cb .eq. cb_ref)
    end do
    end do
    end do
    call check("scale_reg = 0, no regions: bit-identical to calc_cb_ref on z_bed (scale_zb 1,2; scale_sed 0-3; n_sd 1,10)", ok)

    ! === 3. scale_reg = 1 ====================================================

    dyn%par%till_scale_reg = 0
    dyn%now%z_bed_reg = -1.0e6_wp
    call ydyn_till_init(dyn,bnd)
    call check("ydyn_till_init, scale_reg = 0: z_bed_reg not set", all(dyn%now%z_bed_reg .eq. -1.0e6_wp))

    dyn%par%till_scale_reg = 1
    call ydyn_till_init(dyn,bnd)
    call check("ydyn_till_init, scale_reg = 1: z_bed_reg = basin median of z_bed_ref", all(dyn%now%z_bed_reg .eq. zr))

    ok = .TRUE.
    do i_zb = 1, 2
    do i_sed = 0, 3
        dyn%par%till_scale_zb  = i_zb
        dyn%par%till_scale_sed = i_sed

        call calc_cb_ref(cb_ref,dyn%now%z_bed_reg,bnd%z_bed_sd,bnd%z_sl,bnd%H_sed,dyn%par%till_f_sed, &
                dyn%par%till_sed_min,dyn%par%till_sed_max,cfr,dyn%par%till_cf_min, &
                dyn%par%till_z0,dyn%par%till_z1,1,dyn%par%till_scale_zb,dyn%par%till_scale_sed)

        dyn%par%till_n_sd = 10
        call calc_ydyn_cb_ref(cb,dyn%par,bnd,dyn%now%cf_ref,dyn%now%z_bed_reg)
        ok = ok .and. all(cb .eq. cb_ref)

        ! A depressed bed, other z_bed_sd, n_sd = 1: the same
        bnd%z_bed    = bnd%z_bed - 800.0_wp
        bnd%z_bed_sd = 2.0_wp*bnd%z_bed_sd
        dyn%par%till_n_sd = 1
        call calc_ydyn_cb_ref(cb2,dyn%par,bnd,dyn%now%cf_ref,dyn%now%z_bed_reg)
        ok = ok .and. all(cb2 .eq. cb_ref)
        bnd%z_bed    = bnd%z_bed + 800.0_wp
        bnd%z_bed_sd = 0.5_wp*bnd%z_bed_sd
    end do
    end do
    call check("scale_reg = 1: calc_cb_ref on z_bed_reg, n_sd = 1; independent of z_bed, z_bed_sd, n_sd", ok)

    ! One value per basin (no sediment scaling), lower in lower basins
    dyn%par%till_scale_zb  = 1
    dyn%par%till_scale_sed = 0
    dyn%par%till_n_sd      = 10
    dyn%par%till_z0        = -1500.0_wp      ! keep the basin medians within the linear ramp
    dyn%par%till_z1        =  1500.0_wp
    call calc_ydyn_cb_ref(cb,dyn%par,bnd,dyn%now%cf_ref,dyn%now%z_bed_reg)

    ok = .TRUE.
    do k = 1, 6
        b = ids6(k)
        ok = ok .and. maxval(cb,mask=nint(bnd%basins) .eq. b) .eq. minval(cb,mask=nint(bnd%basins) .eq. b)
        zb_b(k) = maxval(dyn%now%z_bed_reg,mask=nint(bnd%basins) .eq. b)
        cb_b(k) = maxval(cb,mask=nint(bnd%basins) .eq. b)
    end do
    call check("scale_reg = 1: one cb_tgt per basin", ok)

    ok = .TRUE.
    do k = 1, 6
    do b2 = 1, 6
        if (zb_b(k) .lt. zb_b(b2)) ok = ok .and. cb_b(k) .lt. cb_b(b2)
    end do
    end do
    call check("scale_reg = 1: lower basin median, lower cb_tgt", ok)

    ! === 4. scale_reg = 1 with a region of its own cf_ref ====================

    call random_number(r)
    bnd%regions = merge(20101.0_wp, 20102.0_wp, r .lt. 0.5_wp)     ! two regions across the basins
    dyn%now%cf_ref = dyn%par%till_cf_ref
    call set_cf_ref_ids(dyn%now%cf_ref,bnd%regions,[20101],[0.2_wp])
    dyn%par%till_scale_zb      = 1
    dyn%par%till_scale_sed     = 0
    call calc_ydyn_cb_ref(cb,dyn%par,bnd,dyn%now%cf_ref,dyn%now%z_bed_reg)

    cfr = dyn%now%cf_ref
    call calc_cb_ref(cb_ref,dyn%now%z_bed_reg,bnd%z_bed_sd,bnd%z_sl,bnd%H_sed,dyn%par%till_f_sed, &
            dyn%par%till_sed_min,dyn%par%till_sed_max,cfr,dyn%par%till_cf_min, &
            dyn%par%till_z0,dyn%par%till_z1,1,dyn%par%till_scale_zb,dyn%par%till_scale_sed)
    call check("scale_reg = 1, regional cf_ref: calc_cb_ref on z_bed_reg with the regional cf_ref field", all(cb .eq. cb_ref))

    ! The law, explicitly: max(cf_min, cf_ref(cell) * lambda(basin median))
    lam = min(1.0_wp, max(0.0_wp, (dyn%now%z_bed_reg - dyn%par%till_z0)/(dyn%par%till_z1 - dyn%par%till_z0)))
    cb2 = max(dyn%par%till_cf_min, merge(0.2_wp, dyn%par%till_cf_ref, bnd%regions .eq. 20101.0_wp)*lam)
    call check("scale_reg = 1, regional cf_ref: cb_ref = max(cf_min, cf_ref(cell)*lambda(basin median))", &
               all(abs(cb - cb2) .le. 1.0e-6_wp))

    ok = .TRUE.
    do k = 1, 6
        b = ids6(k)
        ok = ok .and. count(nint(bnd%basins) .eq. b .and. bnd%regions .eq. 20101.0_wp) .gt. 0
        ok = ok .and. maxval(cb,mask=nint(bnd%basins) .eq. b .and. bnd%regions .eq. 20101.0_wp) .eq. &
                      minval(cb,mask=nint(bnd%basins) .eq. b .and. bnd%regions .eq. 20101.0_wp)
        ok = ok .and. maxval(cb,mask=nint(bnd%basins) .eq. b .and. bnd%regions .eq. 20102.0_wp) .eq. &
                      minval(cb,mask=nint(bnd%basins) .eq. b .and. bnd%regions .eq. 20102.0_wp)
    end do
    call check("scale_reg = 1, regional cf_ref: one value per basin and region", ok)

    if (n_fail .eq. 0) then
        write(*,*) "test_till_scale_reg: PASSED"
    else
        write(*,*) "test_till_scale_reg: FAILED (", n_fail, " checks)"
        error stop 1
    end if

contains

    subroutine median_reference(zm,z,basins,basin_mask)
        ! Basin medians by insertion sort, as documented for calc_basin_median
        real(wp), intent(OUT) :: zm(:,:)
        real(wp), intent(IN)  :: z(:,:), basins(:,:), basin_mask(:,:)
        real(wp), allocatable :: v(:)
        real(wp) :: t
        integer  :: b, n, i, j

        zm = z
        do b = 1, nint(maxval(basins))
            v = pack(z, nint(basins) .eq. b .and. basin_mask .eq. 1.0_wp)
            n = size(v)
            if (n .eq. 0) cycle
            do i = 2, n
                t = v(i)
                j = i - 1
                do while (j .ge. 1)
                    if (v(j) .le. t) exit
                    v(j+1) = v(j)
                    j = j - 1
                end do
                v(j+1) = t
            end do
            if (mod(n,2) .eq. 1) then
                where (nint(basins) .eq. b) zm = v((n+1)/2)
            else
                where (nint(basins) .eq. b) zm = 0.5_wp*(v(n/2) + v(n/2+1))
            end if
        end do
    end subroutine median_reference

    subroutine check(label,pass)
        character(len=*), intent(IN) :: label
        logical,          intent(IN) :: pass
        if (pass) then
            write(*,'(a,a)') "  ok    ", label
        else
            write(*,'(a,a)') "  FAIL  ", label
            n_fail = n_fail + 1
        end if
    end subroutine check

end program test_till_scale_reg
