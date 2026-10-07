# check_flowline.jl
#
# Compare the grounding line of a B2 run (marine flowline) with the steady
# grounding line of the boundary-layer theory (Schoof, 2007). The benchmark
# parameters are read from the fixture attributes. The grounding line is the
# zero of the thickness above flotation, H_grnd = H_ice − (ρ_w/ρ_i)(z_sl − z_bed),
# interpolated linearly between the last grounded and the first floating cell
# centre (as f_grnd_acx in Yelmo, ytopo.gl_sep = 1), on both sides of the
# divide. The script prints x_g over time, the final x_g against the
# reference, dx_g/dt over the last output interval and the asymmetry between
# the two sides and across the rows.
#
# Usage (from the yelmo root):
#   julia --project=tests/bench tests/bench/check_flowline.jl <fixture.nc> <run-dir or yelmo.nc>

using YelmoBench, NCDatasets, Printf

fixture, run = ARGS[1], ARGS[2]
file = isdir(run) ? joinpath(run, "yelmo.nc") : run

b = NCDataset(fixture) do ds
    a = ds.attrib
    FlowlineBenchmark(Symbol(a["exp"]); A = a["A"], dx_km = a["dx_km"], L = a["L"], x_cf = a["x_cf"],
                      ny = a["ny"], transpose = a["transpose"] == 1, C = a["C"], m = a["m"], n = a["n"],
                      smb = a["smb"], rho_ice = a["rho_ice"], rho_sw = a["rho_sw"], g = a["g"],
                      sec_year = a["sec_year"])
end
p    = YelmoBench.flowline_yelmo_params(b)
xg_r = flowline_xg(b)

# The run must use the rate factor and friction coefficient of the fixture
if isdir(run)
    for (key, val) in (("rf_const", p.rf_const), ("cf_ref", p.cf_ref))
        for f in filter(endswith(".nml"), readdir(run; join = true))
            mt = match(Regex("^\\s*" * key * "\\s*=\\s*([-+0-9.eEdD]+)", "m"), read(f, String))
            mt === nothing && continue
            v = parse(Float64, replace(mt[1], r"[dD]" => "e"))
            isapprox(v, val; rtol = 1e-6) ||
                @warn "$(basename(f)): $key = $v, but the fixture needs $key = $val"
        end
    end
end

ds = NCDataset(file)
t     = Float64.(ds["time"][:])
H     = Float64.(coalesce.(ds["H_ice"][:, :, :], 0.0))
z_bed = Float64.(coalesce.(ds["z_bed"][:, :, :], 0.0))
z_sl  = Float64.(coalesce.(ds["z_sl"][:, :, :], 0.0))
close(ds)
if b.transpose
    H, z_bed, z_sl = (permutedims(v, (2, 1, 3)) for v in (H, z_bed, z_sl))
end
x  = b.xc
dx = b.dx_km * 1e3
H_grnd = H .- (b.rho_sw / b.rho_ice) .* (z_sl .- z_bed)
i0 = findfirst(==(0.0), x)

"Grounding line [m] along row j at time index k on the side s = +1 (x > 0) or −1 (x < 0)."
function grounding_line(Hg, k, j, s)
    i = i0
    while 1 <= i + s <= length(x) && Hg[i+s, j, k] > 0
        i += s
    end
    Hg[i, j, k] > 0 || return NaN
    1 <= i + s <= length(x) || return NaN
    f = Hg[i, j, k] / (Hg[i, j, k] - Hg[i+s, j, k])
    return abs(x[i]) + f * dx
end

ny  = size(H, 2)
xgp = [grounding_line(H_grnd, k, j, +1) for k in eachindex(t), j in 1:ny]
xgm = [grounding_line(H_grnd, k, j, -1) for k in eachindex(t), j in 1:ny]
xg  = 0.5 .* (xgp .+ xgm)
xgm_t = vec(sum(xg; dims = 2)) ./ ny
asym  = maximum(abs.(xgp .- xgm))
rows  = maximum(maximum(xg; dims = 2) .- minimum(xg; dims = 2))

println("== ", file)
@printf("A = %.4g Pa^-3 s^-1 (rf_const = %.6g Pa^-3 a^-1), dx = %g km, x_cf = %g km\n",
        b.A, p.rf_const, b.dx_km, b.x_cf / 1e3)
@printf("%10s %10s %12s\n", "t [yr]", "x_g [km]", "dx_g/dt [m/a]")
for k in eachindex(t)
    dxdt = k == 1 ? NaN : (xgm_t[k] - xgm_t[k-1]) / (t[k] - t[k-1])
    @printf("%10.0f %10.2f %12.3g\n", t[k], xgm_t[k] / 1e3, dxdt)
end
dxdt_end = length(t) > 1 ? (xgm_t[end] - xgm_t[end-1]) / (t[end] - t[end-1]) : NaN
@printf("x_g (model, t = %.0f yr) = %.2f km, x_g (Schoof 2007) = %.2f km, difference = %.2f km (%.2f %%)\n",
        t[end], xgm_t[end] / 1e3, xg_r / 1e3, (xgm_t[end] - xg_r) / 1e3, 100 * (xgm_t[end] - xg_r) / xg_r)
@printf("dx_g/dt over the last output interval = %.3g m/a\n", dxdt_end)
@printf("max |x_g(+x) − x_g(−x)| = %.3g m, max spread across rows = %.3g m\n", asym, rows)
