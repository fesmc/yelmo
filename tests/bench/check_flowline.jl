# check_flowline.jl
#
# Compare the grounding line of a B2 run (marine flowline) with the steady
# grounding line of the boundary-layer theory (Schoof, 2007). The benchmark
# parameters (domain, A, start) are read from the fixture attributes. The
# grounding line is the zero of the thickness above flotation, interpolated
# linearly between the last grounded and the first floating cell centre (as
# f_grnd_acx in Yelmo, ytopo.gl_sep = 1), see `flowline_grounding_line`. The
# script prints x_g over time, the final x_g against the reference, dx_g/dt over
# the last output interval and over the last 2 kyr, the model flux across the
# grounding line against the boundary-layer flux q_g(h_g), and the spread across
# the rows (and, for the symmetric domain, between the two sides of the divide).
#
# Usage (from the yelmo root):
#   julia --project=tests/bench tests/bench/check_flowline.jl <fixture.nc> <run-dir or yelmo.nc>

using YelmoBench, Printf

fixture, run = ARGS[1], ARGS[2]

b    = flowline_from_fixture(fixture)
p    = YelmoBench.flowline_yelmo_params(b)
xg_r = flowline_xg(b)

# The run must use the rate factor, friction coefficient and boundary treatment of the fixture
if isdir(run)
    for f in filter(endswith(".nml"), readdir(run; join = true))
        txt = read(f, String)
        for (key, val) in (("rf_const", p.rf_const), ("cf_ref", p.cf_ref))
            mt = match(Regex("^\\s*" * key * "\\s*=\\s*([-+0-9.eEdD]+)", "m"), txt)
            mt === nothing && continue
            v = parse(Float64, replace(mt[1], r"[dD]" => "e"))
            isapprox(v, val; rtol = 1e-6) ||
                @warn "$(basename(f)): $key = $v, but the fixture needs $key = $val"
        end
        mt = match(r"^\s*experiment\s*=\s*\"([^\"]*)\""m, txt)
        mt === nothing || mt[1] == flowline_experiment(b) ||
            @warn "$(basename(f)): experiment = $(mt[1]), but the fixture needs experiment = $(flowline_experiment(b))"
    end
end

gl = flowline_grounding_line(b, run)
t, xg = gl.time, gl.xg
a  = b.smb                                       # [m/a]
sy = b.sec_year

println("== ", run)
@printf("A = %.4g Pa^-3 s^-1 (rf_const = %.6g Pa^-3 a^-1), dx = %g km, domain = %s, x_cf = %g km\n",
        b.A, p.rf_const, b.dx_km, b.domain, b.x_cf / 1e3)
@printf("start: x_g = %.2f km (reference %+.1f km)\n", (xg_r + b.dxg_start) / 1e3, b.dxg_start / 1e3)
@printf("%10s %10s %12s %14s %14s\n", "t [yr]", "x_g [km]", "dx_g/dt [m/a]", "q model [m2/a]", "q_g(h_g) [m2/a]")
for k in eachindex(t)
    dxdt = k == 1 ? NaN : (xg[k] - xg[k-1]) / (t[k] - t[k-1])
    @printf("%10.0f %10.2f %12.3g %14.4g %14.4g\n", t[k], xg[k] / 1e3, dxdt, gl.qg[k], flowline_qg(b, gl.hg[k]) * sy)
end
dxdt_end = length(t) > 1 ? (xg[end] - xg[end-1]) / (t[end] - t[end-1]) : NaN
k2 = findfirst(>=(t[end] - 2000.0), t)
dxdt_2k = t[end] > t[k2] ? (xg[end] - xg[k2]) / (t[end] - t[k2]) : NaN
@printf("x_g (model, t = %.0f yr) = %.2f km, x_g (Schoof 2007) = %.2f km, difference = %.2f km (%.2f %%)\n",
        t[end], xg[end] / 1e3, xg_r / 1e3, (xg[end] - xg_r) / 1e3, 100 * (xg[end] - xg_r) / xg_r)
@printf("dx_g/dt: last output interval = %.3g m/a, last %.0f yr = %.3g m/a\n", dxdt_end, t[end] - t[k2], dxdt_2k)
@printf("flux at x_g: model = %.4g, balance a x_g = %.4g, boundary layer q_g(h_g) = %.4g m2/a (h_g = %.1f m)\n",
        gl.qg[end], a * xg[end], flowline_qg(b, gl.hg[end]) * sy, gl.hg[end])
@printf("max spread across rows = %.3g m", gl.rows)
b.domain == :symmetric ? @printf(", max |x_g(+x) − x_g(−x)| = %.3g m\n", gl.asym) : println()
