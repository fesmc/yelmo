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

# Island base profile
const Bc = 900.0      # [m]  bed elevation at the centre
const Bl = -2000.0    # [m]  bed elevation at r = R0
const R0 = 1000e3     # [m]  radius where the base profile reaches Bl

# Troughs (one per diagonal)
const r_h   = 250e3   # [m]  radius of the trough head (apex of the V)
const alpha = 120.0   # [deg] opening angle of the trough
const ell   = 50e3    # [m]  width of the trough walls and head
const D0    = 1500.0  # [m]  trough depth below the base profile
const B_od  = 700.0   # [m]  depth of the overdeepening (0: off)
const r_od  = 350e3   # [m]  radius of the deepest point of the overdeepening
const w_od  = 60e3    # [m]  half-width of the overdeepening

# Forcing
const smb0  = 0.5     # [m/a] SMB at the centre
const r_ela = 450e3   # [m]   radius of the equilibrium line (SMB = 0)
const r_lim = 750e3   # [m]   calving mask: no ice for r >= r_lim

"Radial base profile of the island."
bed_base(r) = Bc - (Bc - Bl) * r^2 / R0^2

"""
Depth of one V-shaped trough in local coordinates (ξ along, η across the trough
axis). The walls are straight lines from the apex at ξ = r_h, smoothed over ell.
"""
function trough(ξ, η; alpha = alpha, B_od = B_od)
    φ = deg2rad(alpha / 2)
    T = 0.5 * (1 + tanh(((ξ - r_h) * tan(φ) - sqrt(η^2 + ell^2)) / ell))
    D = D0 + B_od * exp(-((ξ - r_od) / w_od)^2)
    return D * T
end

"""
ISLAND4 bed elevation [m]. Troughs lie along the diagonals (rot = 0) or along
the axes (rot = 45, ISLAND4-R). Invariant under the D4 group of the grid.
"""
function z_bed(x, y; alpha = alpha, B_od = B_od, rot = 0.0)
    z = bed_base(hypot(x, y))
    for ψ in deg2rad.((45.0, 135.0, 225.0, 315.0) .+ rot)
        ξ =  x * cos(ψ) + y * sin(ψ)
        η = -x * sin(ψ) + y * cos(ψ)
        z -= trough(ξ, η; alpha, B_od)
    end
    return z
end

"Radial surface mass balance [m/a]."
smb(r) = smb0 * (1 - r / r_ela)

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
