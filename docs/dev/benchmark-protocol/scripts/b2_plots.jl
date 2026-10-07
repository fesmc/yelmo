# B2 marine flowline: grounding line over time, final grounding line vs resolution,
# final profiles and grounding-line flux, for the runs <dir>/os-<dx>km-a<k>-<start>
# (dx = 8, 4, 2 km; k = 1, 4, 7 for A = 4.6416e-24, -25, -26; start = ref, adv, ret)
# with the fixtures <fixdir>/flowline-<dx>km-a<k>-<start>.nc.
# Usage: julia --project=docs/dev/benchmark-protocol/scripts b2_plots.jl <dir> <fixdir> [outdir]
# Prints the results table (markdown) and writes b2_xg_time.png, b2_xg_resolution.png,
# b2_profiles.png and b2_flux.png to outdir (default: current directory).
using YelmoBench, NCDatasets, CairoMakie, Printf

dir, fixdir = ARGS[1], ARGS[2]
outdir = length(ARGS) > 2 ? ARGS[3] : "."

dxs    = (8, 4, 2)
akeys  = (1, 4, 7)
starts = ("ref", "adv", "ret")
colors = Dict(8 => :darkorange, 4 => :seagreen, 2 => :navy)
styles = Dict("ref" => :solid, "adv" => :dash, "ret" => :dot)
sname  = Dict("ref" => "reference", "adv" => "advanced", "ret" => "retreated")
rundir(dx, k, st)  = joinpath(dir, "os-$(dx)km-a$k-$st")
fixture(dx, k, st) = joinpath(fixdir, "flowline-$(dx)km-a$k-$st.nc")
done(dx, k, st)    = isfile(joinpath(rundir(dx, k, st), "yelmo.nc")) && isfile(fixture(dx, k, st))

"Wall time [min] of a run, from the line `Time = … min` the driver prints at the end of out.out (NaN if absent)."
function walltime(r)
    p = joinpath(r, "out.out")
    isfile(p) || return NaN
    mt = match(r"^Time\s+=\s*([0-9.]+)\s*min"m, read(p, String))
    return mt === nothing ? NaN : parse(Float64, mt[1])
end

# Grounding line of all runs
gl = Dict{Tuple{Int,Int,String},Any}()
bm = Dict{Tuple{Int,Int,String},FlowlineBenchmark}()
for dx in dxs, k in akeys, st in starts
    done(dx, k, st) || continue
    b = flowline_from_fixture(fixture(dx, k, st))
    bm[(dx, k, st)] = b
    gl[(dx, k, st)] = flowline_grounding_line(b, rundir(dx, k, st))
end
isempty(gl) && error("b2_plots.jl: no runs found in $dir")
bref(k) = first(b for ((_, kk, _), b) in bm if kk == k)

# Results table
println("| A [Pa⁻³ s⁻¹] | dx [km] | Start | x_g final [km] | Difference [km] | Difference [%] | dx_g/dt, last 2 kyr [m a⁻¹] | t_end [kyr] | Wall time [min] |")
println("|---|---|---|---|---|---|---|---|---|")
for k in akeys, dx in dxs, st in starts
    haskey(gl, (dx, k, st)) || continue
    g, b = gl[(dx, k, st)], bm[(dx, k, st)]
    xr = flowline_xg(b)
    k2 = findfirst(>=(g.time[end] - 2000.0), g.time)
    v  = (g.xg[end] - g.xg[k2]) / (g.time[end] - g.time[k2])
    @printf("| %.4g | %d | %s | %.1f | %+.1f | %+.2f | %.3g | %.0f | %.0f |\n", b.A, dx, st, g.xg[end] / 1e3,
            (g.xg[end] - xr) / 1e3, 100 * (g.xg[end] - xr) / xr, v, g.time[end] / 1e3, walltime(rundir(dx, k, st)))
end

Atitle(k) = @sprintf("A = %.4g Pa⁻³ s⁻¹", bref(k).A)

# 1. x_g(t)
fig = Figure(size = (1300, 420))
for (c, k) in enumerate(akeys)
    ax = Axis(fig[1, c]; title = Atitle(k), xlabel = "time [kyr]", ylabel = c == 1 ? "x_g [km]" : "")
    hlines!(ax, [flowline_xg(bref(k)) / 1e3]; color = :black, linewidth = 1, label = "Schoof (2007)")
    for dx in dxs, st in starts
        g = get(gl, (dx, k, st), nothing)
        g === nothing && continue
        lines!(ax, g.time ./ 1e3, g.xg ./ 1e3; color = colors[dx], linestyle = styles[st],
               label = "$dx km, $(sname[st])")
    end
    c == 3 && Legend(fig[1, 4], ax; unique = true, framevisible = false)
end
save(joinpath(outdir, "b2_xg_time.png"), fig)

# 2. Final x_g vs resolution
fig = Figure(size = (1200, 380))
markers = Dict("ref" => :circle, "adv" => :utriangle, "ret" => :dtriangle)
for (c, k) in enumerate(akeys)
    ax = Axis(fig[1, c]; title = Atitle(k), xlabel = "dx [km]", ylabel = c == 1 ? "x_g final − x_g Schoof [km]" : "",
              xscale = log2, xticks = [2, 4, 8])
    hlines!(ax, [0.0]; color = :black, linewidth = 1, label = "Schoof (2007)")
    for st in starts
        pts = [(dx, gl[(dx, k, st)].xg[end]) for dx in dxs if haskey(gl, (dx, k, st))]
        isempty(pts) && continue
        xr = flowline_xg(bref(k))
        scatterlines!(ax, first.(pts), (last.(pts) .- xr) ./ 1e3; marker = markers[st], markersize = 12,
                      color = :black, linestyle = styles[st], label = sname[st])
    end
    c == 3 && Legend(fig[1, 4], ax; unique = true, framevisible = false)
end
save(joinpath(outdir, "b2_xg_resolution.png"), fig)

# 3. Final profiles of the reference starts against the semi-analytic profile
fig = Figure(size = (1300, 700))
for (c, k) in enumerate(akeys)
    b0 = bref(k)
    xr = flowline_xg(b0)
    xs = collect(0.0:500.0:b0.x_cf - 1.0)
    Hs = flowline_thickness(b0, xs)
    zb = flowline_bed.(xs)
    zs = [x <= xr ? zb[i] + Hs[i] : (1 - b0.rho_ice / b0.rho_sw) * Hs[i] for (i, x) in enumerate(xs)]
    ax1 = Axis(fig[1, c]; title = Atitle(k), xlabel = "x [km]", ylabel = c == 1 ? "elevation [m]" : "")
    ax2 = Axis(fig[2, c]; title = "near the grounding line", xlabel = "x [km]", ylabel = c == 1 ? "elevation [m]" : "")
    for ax in (ax1, ax2)
        lines!(ax, xs ./ 1e3, zb; color = :saddlebrown, linewidth = 1)
        lines!(ax, xs ./ 1e3, zs; color = :black, linestyle = :dash, label = "semi-analytic")
        lines!(ax, xs ./ 1e3, zs .- Hs; color = :black, linestyle = :dash)
        vlines!(ax, [xr / 1e3]; color = :black, linewidth = 0.8)
        for dx in dxs
            haskey(gl, (dx, k, "ref")) || continue
            b = bm[(dx, k, "ref")]
            H, zsrf = NCDataset(joinpath(rundir(dx, k, "ref"), "yelmo.nc")) do ds
                Float64.(coalesce.(ds["H_ice"][:, 2, end], 0.0)), Float64.(coalesce.(ds["z_srf"][:, 2, end], 0.0))
            end
            i = (b.xc .>= 0) .& (H .> 0)
            lines!(ax, b.xc[i] ./ 1e3, zsrf[i]; color = colors[dx], label = "$dx km")
            lines!(ax, b.xc[i] ./ 1e3, zsrf[i] .- H[i]; color = colors[dx])
            vlines!(ax, [gl[(dx, k, "ref")].xg[end] / 1e3]; color = colors[dx], linewidth = 0.8)
        end
    end
    xlims!(ax2, xr / 1e3 - 100, xr / 1e3 + 100)
    ylims!(ax2, -900, 300)
    c == 3 && Legend(fig[1:2, 4], ax1; unique = true, framevisible = false)
end
save(joinpath(outdir, "b2_profiles.png"), fig)

# 4. Flux across the grounding line vs the grounding-line thickness
fig = Figure(size = (1300, 420))
for (c, k) in enumerate(akeys)
    b0 = bref(k)
    sy = b0.sec_year
    xr = flowline_xg(b0)
    hr = -(b0.rho_sw / b0.rho_ice) * flowline_bed(xr)
    hh = collect(range(0.5hr, 1.6hr; length = 200))
    xh = 750e3 .* (720.0 .+ (b0.rho_ice / b0.rho_sw) .* hh) ./ 778.5     # x where h_f(x) = h
    ax = Axis(fig[1, c]; title = Atitle(k), xlabel = "h_g [m]", ylabel = c == 1 ? "flux [m² a⁻¹]" : "",
              yscale = log10)
    lines!(ax, hh, flowline_qg.(Ref(b0), hh) .* sy; color = :black, label = "q_g(h_g), Schoof (2007)")
    lines!(ax, hh, b0.smb .* xh; color = :grey, linestyle = :dash, label = "balance a x_g")
    for dx in dxs, st in starts
        g = get(gl, (dx, k, st), nothing)
        g === nothing && continue
        ok = .!isnan.(g.qg) .& (g.qg .> 0)
        lines!(ax, g.hg[ok], g.qg[ok]; color = colors[dx], linestyle = styles[st], linewidth = 1,
               label = "$dx km, $(sname[st])")
        scatter!(ax, [g.hg[end]], [g.qg[end]]; color = colors[dx], markersize = 8)
    end
    c == 3 && Legend(fig[1, 4], ax; unique = true, framevisible = false)
end
save(joinpath(outdir, "b2_flux.png"), fig)
println("wrote b2_xg_time.png, b2_xg_resolution.png, b2_profiles.png, b2_flux.png to ", outdir)
