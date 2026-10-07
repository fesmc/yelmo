# C0 final states (H with grounding line and coastline, basal homologous temperature) from yelmo_restart.nc.
# Usage: julia --project=docs/dev/benchmark-protocol/scripts c0_maps.jl <run dirs...>  (writes c0_maps.png)
using NCDatasets, CairoMakie
runs = ARGS
fig = Figure(size=(350*length(runs), 700))
for (k,run) in enumerate(runs)
    d = NCDataset(joinpath(run,"yelmo_restart.nc"))
    x = Float64.(d["xc"][:]); y = Float64.(d["yc"][:])
    H = Float64.(d["H_ice"][:,:,end]); fg = Float64.(d["f_grnd"][:,:,end]); zb = Float64.(d["z_bed"][:,:,end]); Tb = Float64.(d["T_prime_b"][:,:,end])
    ax = Axis(fig[1,k], title="$(run): H", aspect=DataAspect()); hidedecorations!(ax)
    hm = heatmap!(ax, x, y, ifelse.(H .> 0, H, NaN), colorrange=(0,3500), nan_color=:grey90)
    contour!(ax, x, y, fg, levels=[0.5], color=:red, linewidth=1); contour!(ax, x, y, zb, levels=[0.0], color=:black, linewidth=0.6)
    k == length(runs) && Colorbar(fig[1,k+1], hm)
    ax2 = Axis(fig[2,k], title="T'_b", aspect=DataAspect()); hidedecorations!(ax2)
    hm2 = heatmap!(ax2, x, y, ifelse.((H .> 0) .& (fg .> 0), Tb, NaN), colormap=:thermal, colorrange=(-20,0), nan_color=:grey90)
    contour!(ax2, x, y, fg, levels=[0.5], color=:red, linewidth=1)
    k == length(runs) && Colorbar(fig[2,k+1], hm2)
    close(d)
end
save("c0_maps.png", fig)
