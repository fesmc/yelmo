# check_island4l.jl
#
# B1 ISLAND4-L checks (docs/dev/benchmark-protocol): D4 symmetry error over
# time (H_ice, velocity, enthalpy), mass budget residual per output interval,
# and the comparison with the rotated geometry ISLAND4-L-R.
#
# Usage (from the yelmo root):
#   julia --project=tests/bench tests/bench/check_island4l.jl <run> [<run>...]
#   julia --project=tests/bench tests/bench/check_island4l.jl <run-L> --rot <run-L-R>
# Each run directory holds yelmo.nc, yelmo_ts.nc and optionally yelmo3D.nc.
# With --rot, the second run (troughs along the axes) is also compared with the first.

using YelmoBench, Printf

fmt(x) = @sprintf("%.1e", x)

"Value of `v` at the time closest to `t`."
at(time, v, t) = isempty(v) ? NaN : v[argmin(abs.(time .- t))]

function report_run(run)
    s = symmetry_series(run)
    b = mass_budget(run)
    println("== ", run)
    println("  D4 symmetry error  ", rpad("t [yr]", 8), rpad("H_ice", 10), rpad("vel", 10), "enth")
    for t in (100.0, 500.0, 1000.0, 2000.0, 5000.0)
        t <= s.time[end] || continue
        println("                     ", rpad(Int(t), 8), rpad(fmt(at(s.time, s.H_ice, t)), 10),
                rpad(fmt(at(s.time, s.vel, t)), 10), fmt(at(s.time, s.enth, t)))
    end
    k = argmax(abs.(b.r_M))
    @printf("  mass budget        max |r_M| = %.1e at t = %.0f yr (interval %.0f yr), mean |r_M| = %.1e\n",
            abs(b.r_M[k]), b.time[k], b.time[2] - b.time[1], sum(abs, b.r_M[2:end]) / (length(b.r_M) - 1))
    c = budget_closure(run)
    kc = argmax(abs.(c.r_C))
    @printf("  closed budget      max |r_C| = %.1e at t = %.0f yr, max |r_resid| (mb_resid + mb_relax) = %.1e  (%d intervals ending at 2D times)\n",
            abs(c.r_C[kc]), c.time[kc], maximum(abs, c.r_resid), length(c.time))
    @printf("  V: %.4f → %.4f 1e6 km³\n", b.V[1] * 1e-15, b.V[end] * 1e-15)
end

function report_rotation(run, run_rot)
    c = rotation_compare(run, run_rot)
    println("== rotation: ", run, " vs ", run_rot, " (relative differences; radii in km, rot − ref)")
    println("  ", rpad("t [yr]", 8), rpad("dV/V", 10), rpad("dA/A", 10), rpad("r_trough", 18), "r_ridge")
    for t in (0.0, 1000.0, 2000.0, 5000.0)
        t <= c.time[end] || continue
        k = argmin(abs.(c.time .- t))
        @printf("  %-8d%-10s%-10s%4.0f/%4.0f (%+4.0f)  %4.0f/%4.0f (%+4.0f)\n", Int(t),
                fmt((c.V_rot[k] - c.V[k]) / c.V[k]), fmt((c.A_rot[k] - c.A[k]) / c.A[k]),
                c.r_trough[k], c.r_trough_rot[k], c.r_trough_rot[k] - c.r_trough[k],
                c.r_ridge[k], c.r_ridge_rot[k], c.r_ridge_rot[k] - c.r_ridge[k])
    end
end

function main(args)
    i = findfirst(==("--rot"), args)
    runs = filter(!=("--rot"), args)
    isempty(runs) && error("usage: check_island4l.jl <run> [<run>...] | <run-L> --rot <run-L-R>")
    foreach(report_run, runs)
    if i !== nothing
        (i == 2 && length(args) == 3) || error("check_island4l.jl: --rot takes exactly one run before and one after it.")
        report_rotation(args[1], args[3])
    end
end

main(ARGS)
