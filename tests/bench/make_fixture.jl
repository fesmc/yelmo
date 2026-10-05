# make_fixture.jl
#
# Write a benchmark fixture for tests/yelmo_bench.f90.
#
# Usage (from the yelmo root):
#   julia --project=tests/bench tests/bench/make_fixture.jl <out.nc> <benchmark> [key=value ...]
# examples:
#   julia --project=tests/bench tests/bench/make_fixture.jl output/bench/island4-16km.nc island4 dx_km=16
#   julia --project=tests/bench tests/bench/make_fixture.jl output/bench/island4-od.nc island4 B_od=700 exp=smb
#   julia --project=tests/bench tests/bench/make_fixture.jl input/bench/shelf-r-10km.nc shelf-r dx_km=10
#
# Keys are the keyword arguments of the benchmark constructor (exp selects the
# experiment, island4 only). Numbers are parsed as Float64, true/false as Bool, the rest as Symbol.

using YelmoBench

function parse_value(s)
    s in ("true", "false") && return s == "true"
    v = tryparse(Float64, s)
    return v === nothing ? Symbol(s) : v
end

function main(args)
    length(args) >= 2 || error("usage: make_fixture.jl <out.nc> <benchmark> [key=value ...]")
    out, name = args[1], lowercase(args[2])
    kw = Dict{Symbol,Any}()
    for a in args[3:end]
        k, v = split(a, "="; limit = 2)
        kw[Symbol(k)] = parse_value(v)
    end
    exp = pop!(kw, :exp, :ctrl)

    b = if name == "island4"
        Island4Benchmark(exp; kw...)
    elseif name == "shelf-r"
        ShelfRadialBenchmark(; kw...)
    else
        error("make_fixture.jl: unknown benchmark $name. Available: island4, shelf-r.")
    end

    write_fixture!(b, out)
    println("Wrote ", out, " (", typeof(b), ", exp = ", exp, ")")
end

main(ARGS)
