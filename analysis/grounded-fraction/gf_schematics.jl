# Schematics for docs/physics/grounded-fraction.md (gl_sep = 3, calc_f_grnd_subgrid_area).
# Usage (from the yelmo root): julia --project=analysis analysis/grounded-fraction/gf_schematics.jl
# Writes docs/img/grounded-fraction-{grid,quarter,island}.png
using CairoMakie

const OUT = joinpath(@__DIR__, "..", "..", "docs", "img")

# Bilinear function on the unit square with corner values h00, h10, h01, h11
bilin(h00, h10, h01, h11, x, y) = h00*(1-x)*(1-y) + h10*x*(1-y) + h01*(1-x)*y + h11*x*y

# Piecewise-bilinear interpolant between cell centres, cell centres at integer
# coordinates (1..n), evaluated at (x, y) inside the block (centres 1..n)
function interp_centres(H, x, y)
    i = clamp(floor(Int, x), 1, size(H,1)-1); j = clamp(floor(Int, y), 1, size(H,2)-1)
    bilin(H[i,j], H[i+1,j], H[i,j+1], H[i+1,j+1], x-i, y-j)
end

# Interpolant of the removed gl_sep = 2: per cell, bilinear between the four corner means
function interp_corner_means(H, x, y)
    i = clamp(round(Int, x), 2, size(H,1)-1); j = clamp(round(Int, y), 2, size(H,2)-1)
    c(di, dj) = 0.25*(H[i,j] + H[i+di,j] + H[i,j+dj] + H[i+di,j+dj])
    bilin(c(-1,-1), c(1,-1), c(-1,1), c(1,1), x-(i-0.5), y-(j-0.5))
end

# ---------------------------------------------------------------------------
# 1. Grid: cell, quarters, interpolation nodes and the control volumes
# ---------------------------------------------------------------------------
function fig_grid()
    fig = Figure(size = (1000, 480), fontsize = 15)
    qcol = Dict(:SW => (:steelblue, 0.30), :SE => (:darkorange, 0.30), :NW => (:seagreen, 0.30), :NE => (:purple, 0.30))

    # (a) one cell, its quarters and the nine nodes of the interpolant
    ax = Axis(fig[1,1], aspect = DataAspect(), title = "(a) Quarters of cell (i,j)")
    hidedecorations!(ax); hidespines!(ax)
    for (q, (x0, y0)) in ((:SW, (-0.5, -0.5)), (:SE, (0.0, -0.5)), (:NW, (-0.5, 0.0)), (:NE, (0.0, 0.0)))
        poly!(ax, Rect(x0, y0, 0.5, 0.5), color = qcol[q], strokewidth = 0)
        text!(ax, x0+0.25, y0+0.25, text = string(q), align = (:center, :center), fontsize = 15)
    end
    # neighbouring cells and centres
    for x in -1.5:1:1.5; lines!(ax, [x, x], [-1.6, 1.6], color = :gray60, linewidth = 1); end
    for y in -1.5:1:1.5; lines!(ax, [-1.6, 1.6], [y, y], color = :gray60, linewidth = 1); end
    lines!(ax, [-0.5, 0.5, 0.5, -0.5, -0.5], [-0.5, -0.5, 0.5, 0.5, -0.5], color = :black, linewidth = 2)
    lines!(ax, [0, 0], [-0.5, 0.5], color = :black, linestyle = :dash, linewidth = 1)
    lines!(ax, [-0.5, 0.5], [0, 0], color = :black, linestyle = :dash, linewidth = 1)
    cx = [x for x in -1:1, y in -1:1][:]; cy = [y for x in -1:1, y in -1:1][:]
    scatter!(ax, cx, cy, color = :black, markersize = 11)
    scatter!(ax, [0.5, -0.5, 0, 0], [0, 0, 0.5, -0.5], color = :white, strokecolor = :black, strokewidth = 1.5,
             marker = :rect, markersize = 11)
    scatter!(ax, [0.5, -0.5, 0.5, -0.5], [0.5, 0.5, -0.5, -0.5], color = :white, strokecolor = :black,
             strokewidth = 1.5, marker = :diamond, markersize = 13)
    # dual cell (bilinear between the four centres) of the NE quarter
    lines!(ax, [0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color = :firebrick, linewidth = 2)
    text!(ax, 1.02, 1.02, text = "bilinear between\nfour centres", color = :firebrick, fontsize = 13,
          align = (:left, :bottom))
    text!(ax, 0.0, -0.62, text = "cell (i,j)", fontsize = 13, align = (:center, :top))
    limits!(ax, -1.6, 2.1, -1.6, 1.9)

    # legend for the node types
    Legend(fig[2,1], [MarkerElement(marker = :circle, color = :black, markersize = 11),
                      MarkerElement(marker = :rect, color = :white, strokecolor = :black, strokewidth = 1.5, markersize = 11),
                      MarkerElement(marker = :diamond, color = :white, strokecolor = :black, strokewidth = 1.5, markersize = 13)],
           ["cell centre: H_grnd", "face midpoint: mean of 2", "cell corner: mean of 4"],
           orientation = :horizontal, framevisible = false, labelsize = 13)

    # (b) control volumes as unions of quarters
    ax2 = Axis(fig[1,2], aspect = DataAspect(), title = "(b) Cell, face and corner: four quarters each")
    hidedecorations!(ax2); hidespines!(ax2)
    for x in -1.5:1:2.5; lines!(ax2, [x, x], [-1.6, 1.6], color = :gray60, linewidth = 1); end
    for y in -1.5:1:1.5; lines!(ax2, [-1.6, 2.6], [y, y], color = :gray60, linewidth = 1); end
    for x in -1:1:2, y in -1:1; scatter!(ax2, [x], [y], color = :black, markersize = 8); end
    # aa: cell (i,j) = its four quarters
    poly!(ax2, Rect(-1.5, -1.5, 1.0, 1.0), color = (:steelblue, 0.35))
    lines!(ax2, [-1, -1], [-1.5, -0.5], color = :black, linestyle = :dash, linewidth = 1)
    lines!(ax2, [-1.5, -0.5], [-1, -1], color = :black, linestyle = :dash, linewidth = 1)
    text!(ax2, -1.0, -0.38, text = "f_grnd (aa)", align = (:center, :center), fontsize = 13)
    # acx: right face of cell = E quarters of cell (i) + W quarters of cell (i+1)
    poly!(ax2, Rect(0.0, -0.5, 1.0, 1.0), color = (:darkorange, 0.35))
    lines!(ax2, [0.5, 0.5], [-0.5, 0.5], color = :darkorange, linewidth = 3)
    lines!(ax2, [0.0, 1.0], [0.0, 0.0], color = :black, linestyle = :dash, linewidth = 1)
    text!(ax2, 0.5, 0.62, text = "f_grnd_acx", align = (:center, :center), fontsize = 13)
    # ab: corner = one quarter of each of the four cells around it
    poly!(ax2, Rect(1.0, 0.0, 1.0, 1.0), color = (:seagreen, 0.35))
    scatter!(ax2, [1.5], [0.5], color = :seagreen, marker = :diamond, markersize = 14)
    text!(ax2, 1.5, 1.12, text = "f_grnd_ab", align = (:center, :center), fontsize = 13)
    limits!(ax2, -1.6, 2.6, -1.6, 1.6)

    rowsize!(fig.layout, 2, Auto(0.08))
    save(joinpath(OUT, "grounded-fraction-grid.png"), fig, px_per_unit = 2)
end

# ---------------------------------------------------------------------------
# 2. One quarter: grounded length on vertical lines, pieces between the roots
# ---------------------------------------------------------------------------
function fig_quarter()
    h00, h10, h01, h11 = -40.0, 20.0, -10.0, 80.0          # B: -40 -> 20, T: -10 -> 80
    rB = -h00/(h10-h00); rT = -h01/(h11-h01)
    fig = Figure(size = (900, 430), fontsize = 15)

    ax = Axis(fig[1,1], aspect = DataAspect(), xlabel = "x", ylabel = "y",
              title = "(a) Quarter: H_grnd ≥ 0 (grey) and the pieces")
    xs = range(0, 1, length = 401)
    z = [bilin(h00, h10, h01, h11, x, y) for x in xs, y in xs]
    contourf!(ax, xs, xs, z, levels = [-1e9, 0, 1e9], colormap = [:white, (:gray70)])
    contour!(ax, xs, xs, z, levels = [0], color = :black, linewidth = 2)
    for r in (rT, rB); lines!(ax, [r, r], [0, 1], color = :firebrick, linestyle = :dash); end
    text!(ax, rT, 1.02, text = "root of T", color = :firebrick, align = (:center, :bottom), fontsize = 12)
    text!(ax, rB, 1.02, text = "root of B", color = :firebrick, align = (:center, :bottom), fontsize = 12)
    for (xm, lbl) in ((rT/2, "0"), ((rT+rB)/2, "P/D"), ((rB+1)/2, "1"))
        text!(ax, xm, -0.09, text = lbl, align = (:center, :center), fontsize = 14, color = :firebrick)
    end
    # one vertical line in the mixed piece
    xl = 0.4; yl = -bilin(h00, h10, 0.0, 0.0, xl, 0.0) / (bilin(0.0, 0.0, h01, h11, xl, 1.0) - bilin(h00, h10, 0.0, 0.0, xl, 0.0))
    lines!(ax, [xl, xl], [0, 1], color = :steelblue, linewidth = 1.5)
    lines!(ax, [xl, xl], [yl, 1], color = :steelblue, linewidth = 5)
    text!(ax, xl+0.02, 0.0, text = "B(x)", color = :steelblue, align = (:left, :bottom), fontsize = 13)
    text!(ax, xl+0.02, 1.0, text = "T(x)", color = :steelblue, align = (:left, :top), fontsize = 13)
    for (x, y, v) in ((0, 0, h00), (1, 0, h10), (0, 1, h01), (1, 1, h11))
        text!(ax, x, y, text = string(Int(v)), align = (x == 0 ? :right : :left, y == 0 ? :top : :bottom),
              offset = (x == 0 ? -4 : 4, y == 0 ? -4 : 4), fontsize = 12)
    end
    limits!(ax, -0.15, 1.15, -0.17, 1.12)

    # (b) the grounded length along x
    ax2 = Axis(fig[1,2], xlabel = "x", ylabel = "grounded length of the line x",
               title = "(b) Integrand: 0, P/D or 1")
    len(x) = (B = h00 + (h10-h00)*x; T = h01 + (h11-h01)*x;
              B >= 0 && T >= 0 ? 1.0 : (B < 0 && T < 0 ? 0.0 : max(B, T)/abs(B-T)))
    band!(ax2, xs, zeros(length(xs)), len.(xs), color = (:gray70, 0.8))
    lines!(ax2, xs, len.(xs), color = :black, linewidth = 2)
    for r in (rT, rB); vlines!(ax2, r, color = :firebrick, linestyle = :dash); end
    scatter!(ax2, [xl], [len(xl)], color = :steelblue, markersize = 10)
    ylims!(ax2, 0, 1.08); xlims!(ax2, 0, 1)
    save(joinpath(OUT, "grounded-fraction-quarter.png"), fig, px_per_unit = 2)
end

# ---------------------------------------------------------------------------
# 3. Grounded centre next to deep ocean: corner means (removed gl_sep = 2) vs centres
# ---------------------------------------------------------------------------
function fig_island()
    H = fill(-150.0, 5, 5); H[3,3] = 40.0                  # 5x5 cells, centres at 1..5
    fig = Figure(size = (900, 420), fontsize = 15)
    xs = range(1.5, 4.5, length = 301)
    for (k, (f, ttl)) in enumerate(((interp_corner_means, "(a) Between corner means (gl_sep = 2, removed)"),
                                    (interp_centres, "(b) Between cell centres (gl_sep = 3)")))
        ax = Axis(fig[1,k], aspect = DataAspect(), title = ttl, titlesize = 14)
        hidedecorations!(ax)
        z = [f(H, x, y) for x in xs, y in xs]
        hm = heatmap!(ax, xs, xs, z, colormap = :RdBu, colorrange = (-150, 150))
        contour!(ax, xs, xs, z, levels = [0], color = :black, linewidth = 2)
        for x in 1.5:1:4.5; lines!(ax, [x, x], [1.5, 4.5], color = :gray40, linewidth = 0.8); end
        for y in 1.5:1:4.5; lines!(ax, [1.5, 4.5], [y, y], color = :gray40, linewidth = 0.8); end
        scatter!(ax, [3.0], [3.0], color = :black, markersize = 8)
        text!(ax, 3.0, 3.42, text = "H_grnd = 40 m", align = (:center, :center), fontsize = 12)
        text!(ax, 2.0, 1.9, text = "-150 m", align = (:center, :center), fontsize = 12, color = :white)
        k == 2 && Colorbar(fig[1,3], hm, label = "H_grnd [m]")
    end
    save(joinpath(OUT, "grounded-fraction-island.png"), fig, px_per_unit = 2)
end

fig_grid()
fig_quarter()
fig_island()
