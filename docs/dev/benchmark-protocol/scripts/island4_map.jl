# island4_map.jl
#
# Map and profiles of the ISLAND4 benchmark geometry and its radial forcing
# (docs/dev/benchmark-protocol/index.md, Sect. "ISLAND4 domain").
#
# Usage (from the yelmo root):
#   julia --project=docs/dev/benchmark-protocol/scripts docs/dev/benchmark-protocol/scripts/island4_map.jl

using CairoMakie

const FIGDIR = normpath(joinpath(@__DIR__, "..", "figures"))
mkpath(FIGDIR)

# Geometry and forcing from YelmoBench (tests/bench/src/island4.jl)
using YelmoBench: island4_bed, island4_base, island4_smb

const B_od  = 700.0   # [m]  overdeepening used in panels (b), (d), (e)
const r_ela = 800e3   # [m]  radius of the equilibrium line (SMB = 0), ISLAND4
const r_lim = 750e3   # [m]  no ice for r >= r_lim

z_bed(x, y; B_od = B_od, rot = 0.0) = island4_bed(x, y; B_od, rot)
smb(r) = island4_smb(r; r_ela)
bed_base(r) = island4_base(r)

function main()
    xs = range(-800e3, 800e3, length = 401)
    km = xs ./ 1e3

    circle(r) = (r / 1e3 .* cos.(range(0, 2π, length = 361)),
                 r / 1e3 .* sin.(range(0, 2π, length = 361)))

    cases = [("(a) ISLAND4",                 (x, y) -> z_bed(x, y; B_od = 0.0)),
             ("(b) ISLAND4, overdeepened",   (x, y) -> z_bed(x, y)),
             ("(c) ISLAND4-R (rotated 45°)", (x, y) -> z_bed(x, y; B_od = 0.0, rot = 45.0))]

    fig = Figure(size = (1500, 1050), fontsize = 15)
    cmap = :oleron
    crange = (-2500, 2500)
    hm = nothing
    for (k, (title, f)) in enumerate(cases)
        z = [f(x, y) for x in xs, y in xs]
        ax = Axis(fig[1, k], aspect = 1, title = title, xlabel = "x (km)",
                  ylabel = k == 1 ? "y (km)" : "")
        hm = heatmap!(ax, km, km, z, colormap = cmap, colorrange = crange)
        contour!(ax, km, km, z, levels = -2000:500:-500, color = (:white, 0.5), linewidth = 0.7)
        contour!(ax, km, km, z, levels = [0.0], color = :black, linewidth = 1.5)
        lines!(ax, circle(r_ela)..., color = :firebrick, linestyle = :dash, linewidth = 1.5)
        lines!(ax, circle(r_lim)..., color = :black, linestyle = :dot, linewidth = 1.5)
        if k < 3
            lines!(ax, [0, 800], [0, 0], color = :gold, linewidth = 2.5)
            lines!(ax, [0, 800 / √2], [0, 800 / √2], color = :darkorange3, linewidth = 2.5)
        end
    end
    Colorbar(fig[1, 4], hm, label = "Bed elevation (m)")

    # (d) profiles along the axis and the diagonal
    rs = range(0, 800e3, length = 401)
    ax = Axis(fig[2, 1:2], title = "(d) Profiles", xlabel = "r (km)",
              ylabel = "Bed elevation (m)")
    hlines!(ax, [0.0], color = :steelblue, linewidth = 1)
    lines!(ax, rs ./ 1e3, [z_bed(r, 0.0; B_od = 0.0) for r in rs], color = :gold,
           linewidth = 2.5, label = "axis (θ = 0°)")
    lines!(ax, rs ./ 1e3, [z_bed(r / √2, r / √2; B_od = 0.0) for r in rs],
           color = :darkorange3, linewidth = 2.5, label = "diagonal (θ = 45°)")
    lines!(ax, rs ./ 1e3, [z_bed(r / √2, r / √2) for r in rs], color = :darkorange3,
           linewidth = 2.5, linestyle = :dash, label = "diagonal, overdeepened")
    vlines!(ax, [r_ela / 1e3], color = :firebrick, linestyle = :dash)
    vlines!(ax, [r_lim / 1e3], color = :black, linestyle = :dot)
    axislegend(ax, position = :lb, framevisible = false)
    axb = Axis(fig[2, 1:2], yaxisposition = :right, ylabel = "SMB (m/a)",
               ylabelcolor = :firebrick, yticklabelcolor = :firebrick)
    hidespines!(axb); hidexdecorations!(axb); hideydecorations!(axb, label = false,
               ticklabels = false, ticks = false)
    lines!(axb, rs ./ 1e3, smb.(rs), color = :firebrick, linewidth = 1.5)
    ylims!(axb, -0.5, 1.0)

    # (e) cross-section through the trough
    ηs = range(-300e3, 300e3, length = 301)
    ax2 = Axis(fig[2, 3:4], title = "(e) Across the trough at r = 450 km",
               xlabel = "Distance from trough axis (km)", ylabel = "Bed elevation (m)")
    hlines!(ax2, [0.0], color = :steelblue, linewidth = 1)
    xc, yc = 450e3 / √2, 450e3 / √2
    lines!(ax2, ηs ./ 1e3, [z_bed(xc - η / √2, yc + η / √2; B_od = 0.0) for η in ηs],
           color = :darkorange3, linewidth = 2.5)
    lines!(ax2, ηs ./ 1e3, [z_bed(xc - η / √2, yc + η / √2) for η in ηs],
           color = :darkorange3, linewidth = 2.5, linestyle = :dash)

    out = joinpath(FIGDIR, "island4.png")
    save(out, fig, px_per_unit = 2)
    println("Saved ", out)

    # Coastline radii and trough-bed slopes for the protocol text
    coast(f) = rs[findfirst(r -> f(r) < 0, rs)] / 1e3
    println("Coast along axis:     ", coast(r -> z_bed(r, 0.0; B_od = 0.0)), " km")
    println("Coast along diagonal: ", coast(r -> z_bed(r / √2, r / √2; B_od = 0.0)), " km")
    zd = [z_bed(r / √2, r / √2) for r in rs]
    retro = rs[2:end][diff(zd) .> 0] ./ 1e3
    isempty(retro) || println("Retrograde (overdeepened): ", first(retro), "–", last(retro), " km")
end

main()
