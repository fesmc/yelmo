# Kleiner (2015) Exp B: CTS of the old (CTS diffusivity override) and new (split
# sensible/latent face flux) enthalpy solvers. Reads <dir>/{old,new}/kb_nz<nz>.nc,
# test_enthalpy kleiner-b output (cr = 1e-4) with every step written for the last
# 1 kyr (see ../index.md), and the analytic profile in tests/data/Kleiner2015.
# Usage: julia --project=docs/dev/cts-subgrid/scripts cts_plots.jl <dir> [out.png]
using NCDatasets, CairoMakie, Statistics, Printf

dir  = ARGS[1]
out  = length(ARGS) > 1 ? ARGS[2] : "cts_kleiner_b.png"
nzs  = [51, 101, 201, 401, 801]
cp, L, H = 2009.0, 3.35e5, 200.0

# Analytic profile: columns z, enth, T' [C], omega [%]
an  = filter(l -> !isempty(strip(l)), readlines(joinpath(@__DIR__, "../../../../tests/data/Kleiner2015/Kleiner2015_EXPB_analytic_nz401_z.dat"))[2:end])
A   = reduce(hcat, [parse.(Float64, split(l)) for l in an])'
z_an, ep_an = A[:, 1], cp .* A[:, 3] .+ L .* A[:, 4] ./ 100
k = findlast(A[:, 4] .> 0)
cts_an = z_an[k] + A[k, 4] / (A[k, 4] - A[k+1, 4]) * (z_an[k+1] - z_an[k])   # omega -> 0

function load(v, nz)
    NCDataset(joinpath(dir, v, "kb_nz$nz.nc")) do ds
        t = Float64.(ds["time"][:]); i = findall(t .> 49000)
        (t = t[i], H_cts = Float64.(ds["H_cts"][i]), z = H .* Float64.(ds["zeta"][:]),
         ep = Float64.(ds["enth"][:, i[end]]) .- cp .* Float64.(ds["T_pmp"][:, i[end]]),
         om0 = Float64(ds["omega"][1, i[end]]))
    end
end
# CTS by linear interpolation of E - E_pmp (the old diagnostic)
function cts_lin(z, ep)
    k = findfirst(ep .< 0) - 1
    z[k] + ep[k] / (ep[k] - ep[k+1]) * (z[k+1] - z[k])
end

R = Dict((v, nz) => load(v, nz) for v in ("old", "new"), nz in nzs)
col = Dict("old" => :firebrick, "new" => :steelblue)

fig = Figure(size = (1200, 800))
ax1 = Axis(fig[1, 1]; xlabel = "time [a]", ylabel = "H_cts [m]", title = "(a) CTS height, every step (nz = 201)")
ax2 = Axis(fig[1, 2]; xlabel = "E - E_pmp [J/kg]", ylabel = "z [m]", title = "(b) profile at 50 ka (nz = 201)")
ax3 = Axis(fig[2, 1]; xlabel = "nz", ylabel = "H_cts [m]", xscale = log2, xticks = nzs, title = "(c) CTS height vs resolution")
ax4 = Axis(fig[2, 2]; xlabel = "nz", ylabel = "base omega [%]", xscale = log2, xticks = nzs, title = "(d) basal water content")

for v in ("old", "new")
    r = R[(v, 201)]; i = findall(r.t .> r.t[end] - 60)
    scatterlines!(ax1, r.t[i], r.H_cts[i]; color = col[v], markersize = 5, label = v)
    scatterlines!(ax2, r.ep, r.z; color = col[v], markersize = 5, label = v)
end
lines!(ax2, ep_an, z_an; color = :black, linestyle = :dash, label = "analytic")
xlims!(ax2, -60, 400); ylims!(ax2, 15, 26)
hlines!(ax1, [cts_an]; color = :black, linestyle = :dash)

o = [R[("old", nz)] for nz in nzs]; n = [R[("new", nz)] for nz in nzs]
rangebars!(ax3, nzs, minimum.(getfield.(o, :H_cts)), maximum.(getfield.(o, :H_cts)); color = col["old"], whiskerwidth = 8)
scatter!(ax3, nzs, mean.(getfield.(o, :H_cts)); color = col["old"], label = "old (mean, range)")
scatterlines!(ax3, nzs, [cts_lin(r.z, r.ep) for r in n]; color = col["new"], linestyle = :dot, label = "new, linear interpolation")
scatterlines!(ax3, nzs, [r.H_cts[end] for r in n]; color = col["new"], label = "new, H_cts")
hlines!(ax3, [cts_an]; color = :black, linestyle = :dash, label = "analytic")
scatterlines!(ax4, nzs, 100 .* getfield.(o, :om0); color = col["old"], label = "old")
scatterlines!(ax4, nzs, 100 .* getfield.(n, :om0); color = col["new"], label = "new")
hlines!(ax4, [A[1, 4]]; color = :black, linestyle = :dash, label = "analytic")
axislegend(ax1; position = :rt); axislegend(ax2; position = :rt); axislegend(ax3; position = :rt); axislegend(ax4; position = :rb)
save(out, fig)

@printf("analytic CTS %.2f m, base omega %.3f %%\n", cts_an, A[1, 4])
println(rpad("nz", 6), "old H_cts mean [min, max]      new H_cts  new lin   old om0  new om0")
for (nz, a, b) in zip(nzs, o, n)
    @printf("%-6d %6.2f [%6.2f, %6.2f]   %8.2f  %8.2f  %7.4f  %7.4f\n", nz, mean(a.H_cts), minimum(a.H_cts),
            maximum(a.H_cts), b.H_cts[end], cts_lin(b.z, b.ep), a.om0, b.om0)
end
