# Mean time step and model speed vs grid resolution for 1-kyr present-day
# initmip runs (GRL and ANT), in the style of Robinson et al. (2022), Fig. 3a,b.
#
# Data: timing_t1k4.csv, from the Levante runs ~/models/yelmo/output/t1k4/<dom><dx>
# (dev 6d0b30d6, par/yelmo_initmip.nml defaults with ctrl.time_end=1000
# ctrl.time_equil=0; 16 OpenMP threads, shared queue). loop_min is the main-loop
# wall time printed at the end of out.out ("Time = ... min"); mean dt = 1000 yr / steps.
#
# Reference: DIVA (black dots) digitised from Robinson et al. (2022), Fig. 3
# (robinson2022_fig3.png): circle centres and axis ticks were located in the
# image pixels, and the values follow from a log-linear fit of the tick rows.
# The 2022 test was Greenland, so the reference is shown in the Greenland figure only.
# Their speed is per hour on one processor; ours is wall time with 16 OpenMP threads.
#
# Usage: julia plot_timing_resolution.jl   (writes timing_resolution_grl.png, timing_resolution_ant.png)

using CairoMakie, DelimitedFiles, Printf

const DIR = @__DIR__

# --- Robinson et al. (2022), Fig. 3, DIVA ---------------------------------
# Pixel rows of the y-axis ticks (left frame) and their values
const TICK_A = ([84.5, 130.5, 166.0, 201.0, 247.5, 282.5, 317.5, 364.0, 399.0, 434.0],
                [5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005])
const TICK_B = ([95.0, 117.0, 138.5, 167.5, 189.5, 211.5, 240.5, 262.5, 284.0, 313.0, 334.5, 356.5, 385.5, 407.5, 429.5],
                [200, 100, 50, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005])
# Pixel rows of the DIVA circle centres at dx = 4, 8, 16, 32 km. Panel (b) at
# 32 km is partly covered by other symbols: row from its widest visible extent.
const ROW_A = [270.7, 202.4, 136.9, 85.5]
const ROW_B = [336.0, 235.4, 145.8, 68.5]
const DX_2022 = [4.0, 8.0, 16.0, 32.0]

function row2val(rows, tick)
    # least-squares fit log10(value) = c0 + c1*row
    r, v = tick
    A = hcat(ones(length(r)), r)
    c = A \ log10.(Float64.(v))
    return 10 .^ (c[1] .+ c[2] .* rows)
end

dt_2022    = row2val(ROW_A, TICK_A)     # [yr]
speed_2022 = row2val(ROW_B, TICK_B)     # [kyr/hr], one processor

# --- Current runs -----------------------------------------------------------
d, hdr = readdlm(joinpath(DIR, "timing_t1k4.csv"), ',', header=true)
col(name) = d[:, findfirst(==(name), vec(hdr))]
dom   = String.(col("domain"))
dx    = Float64.(col("dx_km"))
steps = Float64.(col("steps"))
tend  = Float64.(col("time_end_yr"))
loop  = Float64.(col("loop_min"))
nthr  = Float64.(col("omp_threads"))
hash  = String(col("git_hash")[1])

dt_mean = tend ./ steps                      # [yr]
speed   = (tend ./ 1e3) ./ (loop ./ 60)      # [kyr/hr], wall time

# Slope p of dt ~ dx^p (least squares in log-log)
slope(x, y) = (hcat(ones(length(x)), log10.(x)) \ log10.(y))[2]

xt  = ([4, 8, 16, 32], ["4", "8", "16", "32"])
ytd = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]
yts = 10.0 .^ (-2:3)
fitline(x, y) = (c = hcat(ones(length(x)), log10.(x)) \ log10.(y); xx = [3.5, 36.0]; (xx, 10 .^ (c[1] .+ c[2] .* log10.(xx))))

for (name, tag) in (("Greenland", "grl"), ("Antarctica", "ant"))
    m = dom .== name
    fig = Figure(size = (900, 470), fontsize = 14)
    axa = Axis(fig[1, 1], xscale = log2, yscale = log10, xticks = xt, yticks = (ytd, string.(ytd)),
               xlabel = "Grid resolution (km)", ylabel = "Mean time step (yr)", title = "(a) $name")
    axb = Axis(fig[1, 2], xscale = log2, yscale = log10, xticks = xt,
               yticks = (yts, [@sprintf("%g", v) for v in yts]),
               xlabel = "Grid resolution (km)", ylabel = "Model speed (kyr/hr)", title = "(b) $name")

    # 2022 reference (DIVA, Greenland)
    ref = name == "Greenland"
    if ref
        p22 = slope(DX_2022, dt_2022)
        lines!(axa, fitline(DX_2022, dt_2022)..., color = :grey55, linewidth = 1)
        s22 = scatter!(axa, DX_2022, dt_2022, color = :grey55, markersize = 13)
        scatter!(axb, DX_2022, speed_2022, color = :grey55, markersize = 13)
        text!(axa, 4.2, 3.0, text = @sprintf("p = %.2f", p22), color = :grey45, fontsize = 13)
    end

    # Current dev
    x, y, s = dx[m], dt_mean[m], speed[m]
    p = slope(x, y)
    lines!(axa, fitline(x, y)..., color = :black, linewidth = 1)
    snow = scatter!(axa, x, y, color = :black, markersize = 13)
    scatter!(axb, x, s, color = :black, markersize = 13)
    text!(axa, 4.2, 1.6, text = @sprintf("p = %.2f", p), color = :black, fontsize = 13)

    for ax in (axa, axb); xlims!(ax, 3.4, 37); end
    ylims!(axa, 0.004, 8); ylims!(axb, 3e-3, 1.2e3)
    nt = Int(nthr[m][1])
    els  = ref ? [s22, snow] : [snow]
    labs = ["DIVA, dev $hash: wall time, $nt OpenMP threads"]
    ref && pushfirst!(labs, "DIVA, Robinson et al. (2022): 1 processor")
    Legend(fig[2, 1:2], els, labs, orientation = :horizontal, framevisible = false, labelsize = 12, nbanks = 1)
    save(joinpath(DIR, "timing_resolution_$(tag).png"), fig, px_per_unit = 2)
end

@printf("2022 DIVA: dt = %s yr, speed = %s kyr/hr\n", round.(dt_2022, sigdigits = 3), round.(speed_2022, sigdigits = 3))
for i in eachindex(dx)
    @printf("%-10s %4.0f km  dt %.2f yr  speed %.2f kyr/hr (wall, %d thr)\n", dom[i], dx[i], dt_mean[i], speed[i], nthr[i])
end
