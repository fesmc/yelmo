# Centreline time series (yc=0) of velocity and ice thickness at x = X0 km (default 300) for the
# clip and drag ts_* runs (0-6 kyr, 2D output every 10 yr; see README.md).
# Usage: ./vl_ts_fetch.sh; julia --project=.. vl_centreline.jl [x_km]
using NCDatasets, CairoMakie, Printf
const D = joinpath(@__DIR__, "data")
const RUNS = [("clip, u_max 5000", "ts_clip", :black, 2.0), ("drag, u_max 5000", "ts_drag5k", :red, 1.2),
              ("drag, u_max 6000", "ts_drag6k", :blue, 1.2), ("drag, u_max 8000", "ts_drag8k", :green, 1.2),
              ("drag, u_max 10000", "ts_drag10k", :purple, 1.2), ("drag, u_max 15000", "ts_drag15k", :orange, 1.2), ("drag, u_max 20000", "ts_drag20k", :cyan3, 1.2), ("drag, u_max 30000", "ts_drag30k", :brown, 1.2), ("drag, u_max 50000", "ts_drag50k", :magenta, 1.2)]
const X0 = parse(Float64, get(ARGS, 1, "300"))
function load(k)
    ds = NCDataset(joinpath(D, k, "centreline.nc"))
    t = Float64.(ds["time"][:]); x = Float64.(ds["xc"][:]); i = argmin(abs.(x .- X0))
    sq(v) = Float64.(dropdims(ds[v][i:i, :, :]; dims = (1, 2)))   # (xc, yc=1, time)
    o = (; t, x = x[i], u = sq("uxy_bar"), ub = sq("uxy_b"), H = sq("H_ice"), fg = sq("f_grnd"))
    close(ds); o
end
f = Figure(size = (1200, 800))
axu = Axis(f[1, 1], ylabel = "uxy_bar [m/yr]", title = @sprintf("centreline, x = %.0f km", X0))
axH = Axis(f[2, 1], ylabel = "H_ice [m]", xlabel = "time [yr]")
azu = Axis(f[1, 2], title = "activation 1 (zoom)"); azH = Axis(f[2, 2], xlabel = "time [yr]")
for (l, k, c, lw) in RUNS
    isfile(joinpath(D, k, "centreline.nc")) || continue
    o = load(k)
    for (a, v) in ((axu, o.u), (azu, o.u), (axH, o.H), (azH, o.H)); lines!(a, o.t, v; color = c, linewidth = lw, label = l); end
    # peaks per activation
    pk = [(w, (j = findall(t -> a <= t < b, o.t); isempty(j) ? (NaN, NaN, NaN) : (o.t[j[argmax(o.u[j])]], maximum(o.u[j]), minimum(o.H[j]))))
          for (w, a, b) in (("act1", 3200, 4200), ("act2", 5200, 6000))]
    println(@sprintf("%-18s x=%.0f ", l, o.x), join([@sprintf("%s: peak u %5.0f at %5.0f, min H %4.0f", w, p[2], p[1], p[3]) for (w, p) in pk], " | "))
end
for a in (azu, azH); xlims!(a, 3250, 3750); end
for a in (axu, axH); xlims!(a, 0, 6000); end
axislegend(axu; position = :lt, framevisible = false)
save(joinpath(@__DIR__, "plots", @sprintf("vl_centreline_x%.0f.png", X0)), f)
