# check_symmetry.jl
#
# Print the maximum D4 symmetry error of the standard fields of Yelmo output files.
#
# Usage (from the yelmo root):
#   julia --project=tests/bench tests/bench/check_symmetry.jl <run-dir or yelmo.nc> [...]

using YelmoBench

for a in ARGS
    file = isdir(a) ? joinpath(a, "yelmo.nc") : a
    println("== ", file)
    symmetry_report(file)
end
