# vel-lim: compare clip (dev, branch) and drag runs of TROUGH-F17, 0-8 kyr (see README.md).
# Solver statistics per activation window, clip dev-vs-branch bit-identity, max speeds, surge timing/volume.
# Usage: ./vl_fetch.sh; julia --project=.. vl_analyse.jl
using NCDatasets, Printf, Statistics, CairoMakie
const D = joinpath(@__DIR__, "data")
const RUNS = [("dev clip", "dev_clip", :black), ("branch clip", "br_clip", :gray60), ("branch drag", "br_drag", :red), ("branch drag 6000", "br_drag6k", :blue), ("branch drag 8000", "br_drag8k", :green), ("branch drag 10000", "br_drag10k", :purple)]
const WIN = [("pre 0-3300", 0, 3300), ("act1 3350-3700", 3350, 3700), ("3700-5450", 3700, 5450),
             ("act2 5450-5900", 5450, 5900), ("5900-7650", 5900, 7650), ("act3 7650-8000", 7650, 8000)]
function steps(k)
    ds = NCDataset(joinpath(D, k, "timesteps.nc"))
    o = (; t = Float64.(ds["time"][:]), dt = Float64.(ds["dt_now"][:]), eta = Float64.(ds["pc_eta"][:]),
          ssa = Int.(ds["ssa_iter"][:]), lfail = Int.(ds["ssa_lin_fail"][:]), redo = Int.(ds["iter_redo"][:]))
    close(ds)
    # files copied from running jobs: truncate at the first non-increasing time after the initial record
    i = findfirst(k -> o.t[k] <= o.t[k-1], 3:length(o.t)); n = i === nothing ? length(o.t) : i + 1
    map(v -> v[1:n], o)
end
function speeds(k)
    ds = NCDataset(joinpath(D, k, "yelmo.nc"))
    t = Float64.(ds["time"][:]); u = ds["uxy_bar"][:, :, :]; fg = ds["f_grnd"][:, :, :]; H = ds["H_ice"][:, :, :]
    close(ds)
    g = [(@views m = (fg[:, :, n] .== 1) .& (H[:, :, n] .> 0); u[:, :, n][m]) for n in eachindex(t)]
    f = [(@views m = (fg[:, :, n] .< 1) .& (H[:, :, n] .> 0); u[:, :, n][m]) for n in eachindex(t)]
    (; t, gmax = maximum.(g; init = 0.0), fmax = maximum.(f; init = 0.0), n4 = count.(>(4000), g), n5 = count.(>(4999), g), H, u)
end
have(k) = isfile(joinpath(D, k, "timesteps.nc"))
for (l, k, _) in RUNS
    have(k) || continue
    T = steps(k)
    @printf("\n%-12s end t=%7.1f steps %6d  dt_min steps %5d  redos %d\n", l, T.t[end], length(T.t), count(T.dt .< 0.0101), sum(T.redo))
    @printf("  %-15s %7s %7s %9s %9s %9s %9s %9s\n", "window", "steps", "dtmin", "dt med", "picard", "at 20", "linfail", "eta med")
    for (w, a, b) in WIN
        i = findall(x -> a <= x < b, T.t); isempty(i) && continue
        @printf("  %-15s %7d %7d %9.3g %9.2f %8.0f%% %9.2f %9.3g\n", w, length(i), count(T.dt[i] .< 0.0101), median(T.dt[i]),
                mean(T.ssa[i]), 100 * count(T.ssa[i] .>= 20) / length(i), mean(T.lfail[i] ./ max.(T.ssa[i], 1)), median(T.eta[i]))
    end
end
# bit-identity of branch clip vs dev clip
if have("dev_clip") && have("br_clip")
    A = steps("dev_clip"); B = steps("br_clip"); n = min(length(A.t), length(B.t))
    i = findfirst(j -> A.t[j] != B.t[j] || A.dt[j] != B.dt[j] || A.eta[j] != B.eta[j], 1:n)
    println("\nclip dev vs branch, timesteps: ", i === nothing ? "identical over $n steps (t <= $(A.t[n]))" : "first difference at step $i, t=$(A.t[i])")
    if isfile(joinpath(D, "dev_clip", "yelmo.nc")) && isfile(joinpath(D, "br_clip", "yelmo.nc"))
        a = speeds("dev_clip"); b = speeds("br_clip"); m = min(length(a.t), length(b.t))
        println("clip dev vs branch, yelmo.nc H_ice/uxy_bar identical over $m records: ", a.H[:, :, 1:m] == b.H[:, :, 1:m] && a.u[:, :, 1:m] == b.u[:, :, 1:m])
    end
end
println("\nmax grounded speed / cells >4000 / cells >=5000 (2D output):")
for (l, k, _) in RUNS
    isfile(joinpath(D, k, "yelmo.nc")) || continue
    S = speeds(k); j = findall(t -> any(a <= t < b for (_, a, b) in WIN[[2, 4, 6]]), S.t)
    @printf("  %-12s max grounded %7.0f  max floating %7.0f  | activation records: ", l, maximum(S.gmax), maximum(S.fmax))
    println(join([@sprintf("%g:%.0f/%d/%d", S.t[n], S.gmax[n], S.n4[n], S.n5[n]) for n in j], "  "))
end
f = Figure(size = (1100, 800))
ax1 = Axis(f[1, 1], ylabel = "dt [yr]", yscale = log10)
ax2 = Axis(f[2, 1], ylabel = "Picard iterations (20-step mean)")
ax3 = Axis(f[3, 1], ylabel = "pc_eta", yscale = log10, xlabel = "time [yr]")
for (l, k, c) in RUNS
    have(k) || continue
    T = steps(k); s = [mean(T.ssa[max(1, i - 19):i]) for i in eachindex(T.ssa)]
    lines!(ax1, T.t, max.(T.dt, 1e-3); color = c, linewidth = 0.6, label = l)
    lines!(ax2, T.t, s; color = c, linewidth = 0.6)
    lines!(ax3, T.t, clamp.(T.eta, 1e-5, 10); color = c, linewidth = 0.4)
end
for ax in (ax1, ax2, ax3); xlims!(ax, 0, 8000); end
axislegend(ax1; position = :lb, framevisible = false)
save(joinpath(@__DIR__, "plots", "vl_runs.png"), f)
# surge characterisation per activation: onset/end = first/last step with dt below the CFL plateau,
# ice volume (yelmo_ts, 100-yr records) before and after
println("\nsurges: onset - end (duration) [yr] | V_ice before -> min after [1e6 km3] | V_ice_g")
for (l, k, _) in RUNS
    have(k) || continue
    T = steps(k); ds = NCDataset(joinpath(D, k, "yelmo_ts.nc")); tt = Float64.(ds["time"][:]); V = Float64.(ds["V_ice"][:]); Vg = Float64.(ds["V_ice_g"][:]); close(ds)
    s = String[]
    for (w, a, b) in WIN[[2, 4, 6]]
        js = findall(j -> a - 100 <= T.t[j] < b && T.dt[j] < 2.4, eachindex(T.t)); isempty(js) && continue
        t0, t1 = T.t[js[1]], T.t[js[end]]
        jb = findall(x -> t0 - 300 <= x <= t0, tt); ja = findall(x -> t0 <= x <= t1 + 200, tt)
        (isempty(jb) || isempty(ja)) && continue
        push!(s, @sprintf("%5.0f-%5.0f (%3.0f) V %.4f->%.4f Vg %.4f->%.4f", t0, t1, t1 - t0, maximum(V[jb]), minimum(V[ja]), maximum(Vg[jb]), minimum(Vg[ja])))
    end
    println(@sprintf("  %-18s ", l), join(s, " | "))
end
