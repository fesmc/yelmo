# B1 ISLAND4-L: D4 symmetry error, mass budget residual and ISLAND4-L vs ISLAND4-L-R
# for the runs <dir>/{L,L-R}-{32,16}km-{iso,coup}.
# Usage: julia --project=docs/dev/benchmark-protocol/scripts b1_plots.jl <dir> [out.png]
# Prints the maximum symmetry error per kyr and writes the figure (default b1.png).
using YelmoBench, CairoMakie, Printf

dir = ARGS[1]
out = length(ARGS) > 1 ? ARGS[2] : "b1.png"

cases  = [(th, dx) for th in ("iso", "coup") for dx in (32, 16)]
colors = Dict(("iso", 32) => :steelblue, ("iso", 16) => :navy, ("coup", 32) => :orange, ("coup", 16) => :firebrick)
run(v, th, dx) = joinpath(dir, "$v-$(dx)km-$th")

fig = Figure(size = (1300, 750))
logax = ((1, 1), (1, 2), (1, 3), (2, 1))
ax = [Axis(fig[i, j]; xlabel = "time [kyr]", yscale = (i, j) in logax ? log10 : identity) for i in 1:2, j in 1:3]
ax[1, 1].title = "D4 error H_ice"; ax[1, 2].title = "D4 error velocity"; ax[1, 3].title = "D4 error enthalpy (coupled)"
ax[2, 1].title = "|r_M| per 10 yr (dots: |r_C|, closed)"; ax[2, 2].title = "V [1e6 km³] (dashed: L-R)"
ax[2, 3].title = "margin radius [km]: trough (thick), ridge (thin); dashed: L-R"
floor(v) = max.(v, 1e-9)

println(rpad("run", 18), "max D4 error per kyr (H / vel / enth): 0-1, 1-2, 2-3, 3-4, 4-5 kyr")
for (th, dx) in cases, v in ("L", "L-R")
    r = run(v, th, dx)
    isdir(r) || continue
    c  = colors[(th, dx)]
    ls = v == "L" ? :solid : :dash
    s = symmetry_series(r)
    t = s.time ./ 1e3
    lines!(ax[1, 1], t, floor(s.H_ice); color = c, linestyle = ls, label = "$th $(dx) km")
    lines!(ax[1, 2], t, floor(s.vel);   color = c, linestyle = ls)
    th == "coup" && lines!(ax[1, 3], t, floor(s.enth); color = c, linestyle = ls)
    win(e) = join([@sprintf("%.0e", maximum(e[(s.time .> 1e3(k-1)) .& (s.time .<= 1e3k)])) for k in 1:5], " ")
    println(rpad(basename(r), 18), win(s.H_ice), " / ", win(s.vel), th == "coup" ? " / " * win(s.enth) : "")

    b = mass_budget(r)
    lines!(ax[2, 1], b.time[2:end] ./ 1e3, floor(abs.(b.r_M[2:end])); color = c, linestyle = ls, linewidth = 0.8)
    cl = budget_closure(r)
    scatter!(ax[2, 1], cl.time ./ 1e3, floor(abs.(cl.r_C)); color = c, markersize = 4)
    lines!(ax[2, 2], b.time ./ 1e3, b.V .* 1e-15; color = c, linestyle = ls)

    m = read_margins(r)
    trough = v == "L" ? [x.diag for x in m.radii] : [x.axis for x in m.radii]
    ridge  = v == "L" ? [x.axis for x in m.radii] : [x.diag for x in m.radii]
    lines!(ax[2, 3], m.time ./ 1e3, trough; color = c, linestyle = ls, linewidth = 2.5)
    lines!(ax[2, 3], m.time ./ 1e3, ridge;  color = c, linestyle = ls, linewidth = 1)
end
hlines!(ax[2, 3], [750.0]; color = :grey, linestyle = :dot)
axislegend(ax[1, 1]; position = :rb, unique = true, framevisible = false)
save(out, fig)
println("wrote ", out)
