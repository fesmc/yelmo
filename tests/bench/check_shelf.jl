# check_shelf.jl
#
# Compare the velocity of an A2 run (radial floating shelf) with the analytical
# solution. The benchmark parameters are read from the fixture attributes.
#
# Usage (from the yelmo root):
#   julia --project=tests/bench tests/bench/check_shelf.jl <fixture.nc> <run-dir or yelmo.nc>

using YelmoBench, NCDatasets

fixture, run = ARGS[1], ARGS[2]
file = isdir(run) ? joinpath(run, "yelmo.nc") : run

b = NCDataset(fixture) do ds
    a = ds.attrib
    xc = ds["xc"][:]
    ShelfRadialBenchmark(; dx_km = a["dx_km"], H = a["H"], R_s = a["R_s"], L = (xc[end] + a["dx_km"]/2) * 1e3,
                           A = a["A"], n = a["n"], beta_reg = a["beta_reg"])
end
println("== ", file, "  (strain rate ", round(shelf_strain_rate(b.H; A = b.A, n = b.n); sigdigits = 4), " 1/yr)")
velocity_errors(b, file)
symmetry_report(file)
