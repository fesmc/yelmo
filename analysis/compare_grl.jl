# compare_grl.jl
#
# Compare the enthalpy (method="enth") and temperature (method="temp")
# thermodynamics solvers on the realistic initmip-Greenland domain.
#
# Reads the yelmo1D.nc (time series) and yelmo2D.nc (maps) of an enth and a
# temp run and writes comparison figures (Makie/CairoMakie, styled after
# Robinson et al., 2020) into analysis/figures/. This is the 2D real-domain
# counterpart to compare_enth_temp.jl (EISMINT); it is the validation that
# the enth solver reproduces temp on a real ice sheet after the cp_ref (A1)
# enthalpy fix and the dzsdt->uz fix.
#
# Usage (from the worktree root):
#   julia analysis/compare_grl.jl
# Run directories are configured in EXPERIMENTS below.
#
# analysis/figures/ is intentionally untracked (see analysis/.gitignore).

using CairoMakie
using NCDatasets
using Statistics
using Printf

const HERE   = @__DIR__
const ROOT   = normpath(joinpath(HERE, ".."))
const FIGDIR = joinpath(HERE, "figures")
mkpath(FIGDIR)

# Each experiment: enth/temp run dirs. Point these at the runs you want to
# compare (default: the initmip-grl-16km 1 ka runs used for the A1/dzsdt fix).
const EXPERIMENTS = [
    (name = "initmip-grl-16km", key = "grl16",
     enth = joinpath(ROOT, "output/therm-dev/grl16-enth-final"),
     temp = joinpath(ROOT, "output/therm-dev/grl16-temp-final")),
]

CairoMakie.activate!(type = "png")

# --- Robinson et al. (2020) style theme ------------------------------------
const COL_ENTH = RGBf(0.75, 0.22, 0.17)   # deep red
const COL_TEMP = RGBf(0.17, 0.44, 0.71)   # steel blue

set_theme!(Theme(
    fontsize = 13,
    figure_padding = 14,
    Axis = (
        xgridvisible = false, ygridvisible = false,
        xtickalign = 1, ytickalign = 1, xticksize = 4, yticksize = 4,
        spinewidth = 1.0, titlesize = 14, titlegap = 5,
        titlefont = :bold, xlabelpadding = 3, ylabelpadding = 3,
    ),
    Colorbar = (spinewidth = 1.0, tickalign = 1, ticksize = 3, width = 12),
    Lines = (linewidth = 2.2,),
))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

tof64(a) = map(x -> ismissing(x) ? NaN : Float64(x), a)

function ncget(ds, name)
    haskey(ds, name) || error("variable '$name' not found")
    return tof64(ds[name][:])
end

"Load 1D time-series diagnostics."
function load_1d(dir)
    ds = NCDataset(joinpath(dir, "yelmo1D.nc"))
    out = Dict{String,Any}()
    out["time"] = ncget(ds, "time") ./ 1e3        # kyr
    for v in ("V_ice", "A_ice", "V_sle", "uxy_s", "uxy_bar", "dzsdt")
        haskey(ds, v) && (out[v] = ncget(ds, v))
    end
    close(ds)
    return out
end

"Load final-time 2D maps. Basal fields are taken from the base level (k=1) of
the 3D T_prime; f_pmp/uxy_s are 2D."
function load_2d(dir)
    ds = NCDataset(joinpath(dir, "yelmo2D.nc"))
    out = Dict{String,Any}()
    out["xc"] = ncget(ds, "xc")
    out["yc"] = ncget(ds, "yc")
    nt = ds.dim["time"]
    out["time_end"] = ncget(ds, "time")[nt] / 1e3   # kyr
    for v in ("H_ice", "f_pmp", "uxy_s", "uxy_b", "z_srf")
        haskey(ds, v) && (out[v] = tof64(ds[v][:, :, nt]))
    end
    # basal homologous temperature = T_prime at k=1 (base)
    haskey(ds, "T_prime") && (out["T_prime_b"] = tof64(ds["T_prime"][:, :, 1, nt]))
    close(ds)
    return out
end

maskice(field, H; thr = 10.0) = [H[i] > thr ? field[i] : NaN for i in eachindex(field)]
rs(x, y, z) = reshape(x, length(y), length(z))
noaxis!(a) = hidedecorations!(a, label = false)

# ===========================================================================
# Figures
# ===========================================================================

"Maps: temp | enth | (enth − temp) for the key thermo/geometry fields."
function fig_maps(ex, e2, t2)
    # (key, label, colormap, log10?)
    specs = [("H_ice",     "Ice thickness [m]",           :viridis,  false),
             ("T_prime_b", "Basal homol. T′_b [°C]",       :inferno,  false),
             ("f_pmp",     "Basal temperate frac [-]",     :dense,    false),
             ("uxy_s",     "Surface speed [m/yr]",         :speed,    true)]
    xc, yc = e2["xc"] ./ 1e3, e2["yc"] ./ 1e3   # km
    He, Ht = e2["H_ice"], t2["H_ice"]
    fig = Figure(size = (1150, 1500))
    for (r, (k, lab, cmap, logscale)) in enumerate(specs)
        (haskey(e2, k) && haskey(t2, k)) || continue
        fe = maskice(e2[k], He); ft = maskice(t2[k], Ht)
        if logscale
            fe = [x > 0 ? log10(x) : NaN for x in fe]
            ft = [x > 0 ? log10(x) : NaN for x in ft]
            lab = "log₁₀ " * lab
        end
        dd = fe .- ft
        vals = filter(isfinite, vcat(fe, ft))
        lo, hi = isempty(vals) ? (0.0, 1.0) : (minimum(vals), maximum(vals))
        dvals = filter(isfinite, dd)
        dm = isempty(dvals) ? 1e-12 : max(maximum(abs, dvals), 1e-12)
        ax1 = Axis(fig[r, 1]; aspect = DataAspect(), title = r == 1 ? "temp" : "", ylabel = lab)
        ax2 = Axis(fig[r, 2]; aspect = DataAspect(), title = r == 1 ? "enth" : "")
        ax3 = Axis(fig[r, 3]; aspect = DataAspect(), title = r == 1 ? "enth − temp" : "")
        for a in (ax1, ax2, ax3); noaxis!(a); end
        h1 = heatmap!(ax1, xc, yc, rs(ft, xc, yc); colormap = cmap, colorrange = (lo, hi))
        heatmap!(ax2, xc, yc, rs(fe, xc, yc); colormap = cmap, colorrange = (lo, hi))
        h3 = heatmap!(ax3, xc, yc, rs(dd, xc, yc); colormap = :balance, colorrange = (-dm, dm))
        for (a, H) in ((ax1, Ht), (ax2, He), (ax3, He))
            contour!(a, xc, yc, rs(H, xc, yc); levels = [10.0], color = :black, linewidth = 0.5)
        end
        Colorbar(fig[r, 4], h1); Colorbar(fig[r, 5], h3)
        # annotate the RMS difference over ice
        rms = isempty(dvals) ? 0.0 : sqrt(mean(dvals .^ 2))
        text!(ax3, 0.03, 0.02; text = @sprintf("RMS = %.3g", rms),
              space = :relative, align = (:left, :bottom), fontsize = 11)
    end
    Label(fig[0, :], @sprintf("%s at %.1f kyr — temperature | enthalpy | difference",
          ex.name, e2["time_end"]), fontsize = 17, font = :bold)
    p = joinpath(FIGDIR, "$(ex.key)_maps.png"); save(p, fig); println("wrote ", p)
end

"Time series: integrated volume/area/velocity for temp vs enth."
function fig_timeseries(ex, e1, t1)
    specs = [("V_ice", "Ice volume [10⁶ km³]"),
             ("A_ice", "Ice area [10⁶ km²]"),
             ("uxy_s", "Mean surface speed [m/yr]"),
             ("V_sle", "Sea-level equiv. [m]")]
    fig = Figure(size = (1000, 720))
    for (n, (k, ylab)) in enumerate(specs)
        ax = Axis(fig[(n-1)÷2+1, (n-1)%2+1]; title = ylab, xlabel = "time [kyr]")
        haskey(t1, k) && lines!(ax, t1["time"], t1[k]; color = COL_TEMP, label = "temp")
        haskey(e1, k) && lines!(ax, e1["time"], e1[k]; color = COL_ENTH, label = "enth")
        n == 1 && axislegend(ax; position = :rb, framevisible = false, rowgap = 0)
    end
    Label(fig[0, :], "$(ex.name): enthalpy vs temperature solver — time series",
          fontsize = 17, font = :bold)
    p = joinpath(FIGDIR, "$(ex.key)_timeseries.png"); save(p, fig); println("wrote ", p)
end

# ===========================================================================
# Main
# ===========================================================================

for ex in EXPERIMENTS
    isdir(ex.enth) && isdir(ex.temp) || (@warn "missing run dirs for $(ex.name), skipping"; continue)
    println("\n=== $(ex.name) ===")
    e2, t2 = load_2d(ex.enth), load_2d(ex.temp)
    fig_maps(ex, e2, t2)
    e1, t1 = load_1d(ex.enth), load_1d(ex.temp)
    fig_timeseries(ex, e1, t1)

    # quick numeric summary to stdout (maskice already NaNs non-ice cells)
    He, Ht = e2["H_ice"], t2["H_ice"]
    dTb = filter(isfinite, maskice(e2["T_prime_b"], He) .- maskice(t2["T_prime_b"], Ht))
    @printf("basal T′_b (enth−temp): mean %.2f °C, RMS %.2f °C\n",
            mean(dTb), sqrt(mean(dTb .^ 2)))
    @printf("mean thickness: temp %.1f m, enth %.1f m\n",
            mean(Ht[Ht .> 10]), mean(He[He .> 10]))
end
