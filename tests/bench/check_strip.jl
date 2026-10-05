# check_strip.jl
#
# Compare the ice temperature of an A4 run (thermodynamics strip) with the
# analytic column solution of IceColumnSolutions.jl. The benchmark parameters
# are read from the fixture attributes. For each row (one parameter set) the
# script prints the Péclet number and the maximum error over the compared
# columns of its middle row at each output time, and the spread across these
# columns (which should be at round-off, since the solution is uniform in x).
#
# Usage (from the yelmo root):
#   julia --project=tests/bench tests/bench/check_strip.jl <fixture.nc> <run-dir or yelmo3D.nc> [n_modes]

using YelmoBench, NCDatasets, Printf

fixture, run = ARGS[1], ARGS[2]
n_modes = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10
file = isdir(run) ? joinpath(run, "yelmo3D.nc") : run

b = NCDataset(fixture) do ds
    a = ds.attrib
    StripThermoBenchmark(Symbol(a["exp"]); dx_km = a["dx_km"], nx = Int(a["nx"]),
                         rows = collect(zip(Float64.(a["rows_H"]), Float64.(a["rows_smb"]))),
                         T_srf = a["T_srf"], Q_geo = a["Q_geo"], dT = a["dT"], fsmb = a["fsmb"],
                         z_bed = a["z_bed"])
end

ds    = NCDataset(file)
t     = Float64.(ds["time"][:])
zeta  = Float64.(ds["zeta"][:])
T_ice = Float64.(ds["T_ice"][:, :, :, :])
close(ds)
nx, ny = size(T_ice, 1), size(T_ice, 2)
ic = 3:nx-2                                   # compared columns (see strip.jl)
# Times printed: all, or a selection when the run has many outputs (the
# maximum error covers all times)
tsel = length(t) <= 8 ? collect(eachindex(t)) :
       unique([argmin(abs.(t .- x)) for x in (0.0, 1e3, 5e3, 1e4, 2e4, 5e4) if x <= t[end]])

println("== ", file, "  (exp = ", b.exp, ", nz = ", length(zeta), ", n_modes = ", n_modes, ")")
@printf("%4s %7s %6s %6s  %s\n", "row", "H [m]", "SMB", "Pe", join([@sprintf("%9s", "t=$(Int(round(t[n])))") for n in tsel], " "))
emax = 0.0; spread = 0.0
for k in eachindex(b.rows)
    j = strip_check_row(b, k)
    H, smb = b.rows[k]
    p  = strip_params(b, k)
    Ta = strip_solution(b, k, zeta, t; n_modes)
    err = [maximum(abs.(T_ice[i, j, :, n] .- Ta[:, n])) for i in ic, n in eachindex(t)]
    e   = vec(maximum(err; dims = 1))
    global emax   = max(emax, maximum(e))
    global spread = max(spread, maximum(maximum(T_ice[ic, j, :, :]; dims = 1) .- minimum(T_ice[ic, j, :, :]; dims = 1)))
    @printf("%4d %7.0f %6.2f %6.2f  %s\n", k, H, p.w0, p.Pe, join([@sprintf("%9.3g", e[n]) for n in tsel], " "))
end
@printf("max |T - T_analytic| = %.3g K,  max spread across columns = %.3g K\n", emax, spread)
