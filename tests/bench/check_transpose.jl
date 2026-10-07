# check_transpose.jl
#
# Compare a run with a run on the transposed grid (x ↔ y), e.g. the B2
# flowline periodic in y with the transposed strip periodic in x
# (make_fixture.jl flowline transpose=true, yelmo.experiment = "periodic-x").
# Prints the maximum relative difference over all output times for the
# standard fields, which should be at round-off.
#
# Usage (from the yelmo root):
#   julia --project=tests/bench tests/bench/check_transpose.jl <run-dir or yelmo.nc> <transposed run-dir or yelmo.nc>

using YelmoBench

length(ARGS) == 2 || error("usage: check_transpose.jl <run> <transposed run>")
files = [isdir(a) ? joinpath(a, "yelmo.nc") : a for a in ARGS]
println("== ", files[1], "  vs transposed  ", files[2])
transpose_report(files...)
