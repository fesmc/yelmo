# ----------------------------------------------------------------------
# Margin radii along the axes and diagonals, and the comparison of a run
# with its rotated counterpart (ISLAND4 vs ISLAND4-R, ISLAND4-L vs ISLAND4-L-R).
#
# The ISLAND4 grids have nx = ny even, with the origin on the corner of the
# four central cells. Along the axes, the rays run through the two central
# rows (columns) of cell centres at ±dx/2 from the axis; along the diagonals,
# cell centres lie on the diagonal.
# ----------------------------------------------------------------------

export margin_radii, read_margins, rotation_compare

"Number of contiguous cells with H > 0 from the start of `idx` outward."
function ice_run(H, idx)
    n = 0
    for (i, j) in idx
        H[i, j] > 0 || break
        n += 1
    end
    return n
end

"""
    margin_radii(H, dx) -> NamedTuple

Margin radius [same units as dx] along the axes and the diagonals of an
ice-thickness field on a square grid with the origin on a cell corner. The
radius of a ray is the outer edge of its last contiguous ice-covered cell:
m dx along the axes (8 rays: 4 directions × 2 central rows/columns) and
√2 m dx along the diagonals (4 rays). Returns the mean and the spread
(max − min) over the rays: (axis, axis_spread, diag, diag_spread).
"""
function margin_radii(H::AbstractMatrix, dx::Real)
    nx, ny = size(H)
    (nx == ny && iseven(nx)) || error("margin_radii: needs nx = ny even (origin on a cell corner).")
    c = nx ÷ 2
    out  = c+1:nx
    back = c:-1:1
    axis = Float64[]
    for j0 in (c, c+1)
        push!(axis, ice_run(H, ((i, j0) for i in out)), ice_run(H, ((i, j0) for i in back)))
        push!(axis, ice_run(H, ((j0, j) for j in out)), ice_run(H, ((j0, j) for j in back)))
    end
    diag = Float64[ice_run(H, ((c+m, c+m) for m in 1:c)),   ice_run(H, ((c+1-m, c+m) for m in 1:c)),
                   ice_run(H, ((c+m, c+1-m) for m in 1:c)), ice_run(H, ((c+1-m, c+1-m) for m in 1:c))]
    axis .*= dx
    diag .*= sqrt(2) * dx
    return (axis = sum(axis)/length(axis), axis_spread = maximum(axis) - minimum(axis),
            diag = sum(diag)/length(diag), diag_spread = maximum(diag) - minimum(diag))
end

"Output times, margin radii (yelmo.nc), V and A (yelmo_ts.nc) at the 2D output times of a run."
function read_margins(run::AbstractString)
    t, rad = NCDataset(joinpath(run, "yelmo.nc")) do ds
        xc = Float64.(ds["xc"][:])
        Float64.(ds["time"][:]), [margin_radii(read2d(ds, "H_ice"; k), xc[2] - xc[1]) for k in 1:length(ds["time"])]
    end
    V, A = NCDataset(joinpath(run, "yelmo_ts.nc")) do ds
        ts  = Float64.(ds["time"][:])
        idx = [findfirst(≈(ti), ts) for ti in t]
        any(isnothing, idx) && error("read_margins: 2D output times missing in yelmo_ts.nc of $run.")
        Float64.(ds["V_ice"][idx]), Float64.(ds["A_ice"][idx])
    end
    return (time = t, radii = rad, V = V, A = A)
end

"""
    rotation_compare(run, run_rot) -> NamedTuple

Compare a run (troughs along the diagonals, rot = 0) with its rotated
counterpart (troughs along the axes, rot = 45) at the 2D output times. Axes and
diagonals are exchanged: the trough radius is the diagonal radius of `run` and
the axis radius of `run_rot`, and vice versa for the ridges. V and A [1e6 km³,
1e6 km²] come from yelmo_ts.nc, the radii [km] from H_ice in yelmo.nc.
Returns vectors over time: time, V, V_rot, A, A_rot, r_trough, r_trough_rot,
r_ridge, r_ridge_rot.
"""
function rotation_compare(run::AbstractString, run_rot::AbstractString)
    a = read_margins(run)
    b = read_margins(run_rot)
    a.time == b.time || error("rotation_compare: the runs have different output times.")
    return (time = a.time, V = a.V, V_rot = b.V, A = a.A, A_rot = b.A,
            r_trough = [r.diag for r in a.radii], r_trough_rot = [r.axis for r in b.radii],
            r_ridge  = [r.axis for r in a.radii], r_ridge_rot  = [r.diag for r in b.radii])
end
