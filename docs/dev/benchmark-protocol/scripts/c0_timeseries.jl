# C0 time series (V, A_g, A_f, drift, grounded f_pmp, speed) from yelmo_ts.nc of each run directory.
# Usage: julia --project=docs/dev/benchmark-protocol/scripts c0_timeseries.jl <run dirs...>  (writes c0_timeseries.png)
using NCDatasets, Printf, Statistics, CairoMakie
fig = Figure(size=(1300,700)); axs = [Axis(fig[i,j]) for i in 1:2, j in 1:3]
titles = ["V_ice [1e6 km3]","A_ice_g [1e6 km2]","A_ice_f [1e6 km2]","dV/dt / V [% per kyr]","grounded f_pmp","uxy_s_g [m/yr]"]
for (k,(i,j)) in enumerate(((1,1),(1,2),(1,3),(2,1),(2,2),(2,3))); axs[i,j].title = titles[k]; axs[i,j].xlabel = "time [kyr]"; end
for run in ARGS
    d = NCDataset(joinpath(run,"yelmo_ts.nc"))
    t = Float64.(d["time"][:]); V = Float64.(d["V_ice"][:]); Ag = Float64.(d["A_ice_g"][:]); Af = Float64.(d["A_ice_f"][:])
    fp = Float64.(d["f_pmp"][:]); ug = Float64.(d["uxy_s_g"][:])
    # drift: linear fit of V over last 5 kyr
    m = t .>= t[end]-5e3
    tt = t[m] .- mean(t[m]); slope = sum(tt .* (V[m] .- mean(V[m]))) / sum(tt.^2)   # per yr
    drift = 100*slope*1e3/mean(V[m])
    # running drift over 1-kyr windows
    tk = Float64[]; dk = Float64[]
    for t0 in 1e3:1e3:t[end]
        w = (t .>= t0-1e3) .& (t .<= t0); length(findall(w)) < 3 && continue
        x = t[w] .- mean(t[w]); s = sum(x .* (V[w] .- mean(V[w])))/sum(x.^2)
        push!(tk, t0/1e3); push!(dk, 100*s*1e3/mean(V[w]))
    end
    @printf("%-18s t=%5.0f kyr | V=%.4f (units as file) Ag=%.4f Af=%.4f | drift(last 5 kyr) = %+.3f %%/kyr | f_pmp_g=%.3f | uxy_s_g=%.1f\n",
        run, t[end]/1e3, V[end], Ag[end], Af[end], drift, fp[end], ug[end])
    lbl = basename(run)
    lines!(axs[1,1], t./1e3, V, label=lbl); lines!(axs[1,2], t./1e3, Ag); lines!(axs[1,3], t./1e3, Af)
    lines!(axs[2,1], tk, dk); lines!(axs[2,2], t./1e3, fp); lines!(axs[2,3], t./1e3, ug)
    close(d)
end
hlines!(axs[2,1], [-0.1, 0.1], color=:grey, linestyle=:dash); ylims!(axs[2,1], -2, 2)
axislegend(axs[1,1], position=:rt)
save("c0_timeseries.png", fig)
