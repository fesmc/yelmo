# compare_enth_temp.jl
#
# Compare the enthalpy (method="enth") and temperature (method="temp")
# thermodynamics solvers on the 2D EISMINT-A benchmark.
#
# Reads the yelmo1D.nc (time series) and yelmo2D.nc (maps) output from two runs
# and writes comparison figures (Makie/CairoMakie) into analysis/figures/.
#
# Usage (from the worktree root):
#   julia analysis/compare_enth_temp.jl [ENTH_RUNDIR] [TEMP_RUNDIR]
# Defaults: ENTH_RUNDIR=RUNDIR, TEMP_RUNDIR=RUNDIR_temp
#
# The figures/ folder is intentionally untracked (see analysis/.gitignore).

using CairoMakie
using NCDatasets
using Statistics
using Printf

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

const HERE     = @__DIR__
const ROOT     = normpath(joinpath(HERE, ".."))
enth_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(ROOT, "RUNDIR")
temp_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(ROOT, "RUNDIR_temp")
const FIGDIR   = joinpath(HERE, "figures")
mkpath(FIGDIR)

CairoMakie.activate!(type = "png")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"Convert an array (possibly holding missing) to Float64 with missing -> NaN."
tof64(a) = map(x -> ismissing(x) ? NaN : Float64(x), a)

"Read a variable from an NCDataset as a plain Float64 array, missing -> NaN."
function ncget(ds, name)
    haskey(ds, name) || error("variable '$name' not found")
    return tof64(ds[name][:])
end

"Load the 1D time-series diagnostics we compare."
function load_1d(dir)
    ds = NCDataset(joinpath(dir, "yelmo1D.nc"))
    out = Dict{String,Any}()
    out["time"] = ncget(ds, "time") ./ 1e3          # kyr
    for v in ("V_ice", "A_ice", "f_pmp", "bmb", "H_ice", "W_til")
        haskey(ds, v) && (out[v] = ncget(ds, v))
    end
    close(ds)
    return out
end

"Load the final-time 2D maps and axes we compare, plus a mid-row cross-section."
function load_2d(dir)
    ds = NCDataset(joinpath(dir, "yelmo2D.nc"))
    out = Dict{String,Any}()
    out["xc"] = ncget(ds, "xc")
    out["yc"] = ncget(ds, "yc")
    out["zeta"] = ncget(ds, "zeta")
    nt = ds.dim["time"]
    out["time_end"] = ncget(ds, "time")[nt] / 1e3    # kyr
    # 2D surface/basal fields: NCDatasets gives (x, y, time)
    for v in ("H_ice", "T_prime_b", "f_pmp", "bmb_grnd", "H_cts")
        haskey(ds, v) && (out[v] = tof64(ds[v][:, :, nt]))
    end
    # 3D fields: (x, y, zeta, time)
    jc = (length(out["yc"]) + 1) ÷ 2
    for v in ("omega", "T_prime")
        if haskey(ds, v)
            out[v*"_base"] = tof64(ds[v][:, :, 1, nt])
            out[v*"_xsec"] = tof64(ds[v][:, jc, :, nt])
        end
    end
    out["jc"] = (length(out["yc"]) + 1) ÷ 2
    close(ds)
    return out
end

"Mask a field to ice-covered cells (H_ice > thr), else NaN."
maskice(field, H; thr = 1.0) = [H[i] > thr ? field[i] : NaN for i in eachindex(field)]

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

println("Reading enth from: ", enth_dir)
println("Reading temp from: ", temp_dir)

e1 = load_1d(enth_dir); t1 = load_1d(temp_dir)
e2 = load_2d(enth_dir); t2 = load_2d(temp_dir)

const COL_ENTH = :firebrick
const COL_TEMP = :steelblue

# ===========================================================================
# Figure 1 — time series
# ===========================================================================

series_specs = [
    ("V_ice", "Ice volume", "V_ice"),
    ("A_ice", "Ice area",   "A_ice"),
    ("f_pmp", "Basal temperate fraction f_pmp", "f_pmp [-]"),
    ("W_til", "Mean basal water thickness W_til", "W_til [m]"),
]

fig1 = Figure(size = (1000, 720))
for (n, (key, title, ylab)) in enumerate(series_specs)
    row = (n - 1) ÷ 2 + 1
    col = (n - 1) % 2 + 1
    ax = Axis(fig1[row, col]; title = title, xlabel = "time [kyr]", ylabel = ylab)
    if haskey(t1, key)
        lines!(ax, t1["time"], t1[key]; color = COL_TEMP, label = "temp", linewidth = 2)
    end
    if haskey(e1, key)
        lines!(ax, e1["time"], e1[key]; color = COL_ENTH, label = "enth", linewidth = 2)
    end
    n == 1 && axislegend(ax; position = :rb)
end
Label(fig1[0, :], "EISMINT-A: enthalpy vs temperature solver — time series",
      fontsize = 18, font = :bold)
save(joinpath(FIGDIR, "timeseries_enth_vs_temp.png"), fig1)
println("wrote ", joinpath(FIGDIR, "timeseries_enth_vs_temp.png"))

# ===========================================================================
# Figure 2 — final-time maps: temp | enth | (enth - temp)
# ===========================================================================

map_specs = [
    # key,        label,                         cmap,       symmetric-diff?
    ("H_ice",     "Ice thickness [m]",           :viridis,   false),
    ("T_prime_b", "Basal homol. temp T'_b [K]",  :thermal,   false),
    ("f_pmp",     "Basal temperate frac f_pmp",  :viridis,   false),
    ("bmb_grnd",  "Grounded bmb [m/yr]",         :balance,   true),
]

xc = e2["xc"]; yc = e2["yc"]
Hmask_e = e2["H_ice"]; Hmask_t = t2["H_ice"]

fig2 = Figure(size = (1180, 1180))
for (r, (key, label, cmap, symdiff)) in enumerate(map_specs)
    haskey(e2, key) && haskey(t2, key) || continue
    fe = maskice(e2[key], Hmask_e)
    ft = maskice(t2[key], Hmask_t)
    dd = fe .- ft
    # shared color range for temp/enth
    vals = filter(isfinite, vcat(fe, ft))
    lo, hi = isempty(vals) ? (0.0, 1.0) : (minimum(vals), maximum(vals))
    if symdiff
        m = maximum(abs, filter(isfinite, vcat(fe, ft)); init = 1.0)
        lo, hi = -m, m
    end
    dvals = filter(isfinite, dd)
    dm = isempty(dvals) ? 1.0 : maximum(abs, dvals; init = 1e-12)
    dm = dm == 0 ? 1e-12 : dm

    ax1 = Axis(fig2[r, 1]; aspect = 1, title = r == 1 ? "temp" : "", ylabel = label)
    ax2 = Axis(fig2[r, 2]; aspect = 1, title = r == 1 ? "enth" : "")
    ax3 = Axis(fig2[r, 3]; aspect = 1, title = r == 1 ? "enth − temp" : "")
    for a in (ax1, ax2, ax3); hidedecorations!(a, label = false); end

    hm1 = heatmap!(ax1, xc, yc, reshape(ft, length(xc), length(yc));
                   colormap = cmap, colorrange = (lo, hi))
    heatmap!(ax2, xc, yc, reshape(fe, length(xc), length(yc));
             colormap = cmap, colorrange = (lo, hi))
    hm3 = heatmap!(ax3, xc, yc, reshape(dd, length(xc), length(yc));
                   colormap = :balance, colorrange = (-dm, dm))
    Colorbar(fig2[r, 4], hm1)
    Colorbar(fig2[r, 5], hm3)
end
Label(fig2[0, :], @sprintf("EISMINT-A maps at %.0f kyr — temp | enth | difference", e2["time_end"]),
      fontsize = 18, font = :bold)
save(joinpath(FIGDIR, "maps_enth_vs_temp_final.png"), fig2)
println("wrote ", joinpath(FIGDIR, "maps_enth_vs_temp_final.png"))

# ===========================================================================
# Figure 3 — enthalpy-only polythermal structure
# ===========================================================================

fig3 = Figure(size = (1150, 460))

# basal water content map
if haskey(e2, "omega_base")
    ax = Axis(fig3[1, 1]; aspect = 1, title = "Basal water content ω (enth)",
              xlabel = "x [km]", ylabel = "y [km]")
    hm = heatmap!(ax, xc, yc, reshape(maskice(e2["omega_base"], Hmask_e), length(xc), length(yc));
                  colormap = :dense)
    Colorbar(fig3[1, 2], hm)
end

# CTS height map
if haskey(e2, "H_cts")
    ax = Axis(fig3[1, 3]; aspect = 1, title = "CTS height H_cts [m] (enth)",
              xlabel = "x [km]", ylabel = "y [km]")
    hcts = maskice(e2["H_cts"], Hmask_e)
    hm = heatmap!(ax, xc, yc, reshape(hcts, length(xc), length(yc)); colormap = :haline)
    Colorbar(fig3[1, 4], hm)
end

# vertical cross-section of omega through the dome centre
if haskey(e2, "omega_xsec")
    ax = Axis(fig3[1, 5]; title = "ω cross-section (y = 0)",
              xlabel = "x [km]", ylabel = "ζ (0=base, 1=surf)")
    om = e2["omega_xsec"]                  # (x, zeta)
    # mask columns with no ice
    jc = e2["jc"]
    Hrow = e2["H_ice"][:, jc]
    for i in axes(om, 1), k in axes(om, 2)
        Hrow[i] > 1.0 || (om[i, k] = NaN)
    end
    hm = heatmap!(ax, xc, e2["zeta"], om; colormap = :dense)
    Colorbar(fig3[1, 6], hm; label = "ω [-]")
end
Label(fig3[0, :], @sprintf("EISMINT-A polythermal structure (enth) at %.0f kyr", e2["time_end"]),
      fontsize = 18, font = :bold)
save(joinpath(FIGDIR, "polythermal_structure_enth.png"), fig3)
println("wrote ", joinpath(FIGDIR, "polythermal_structure_enth.png"))

# ===========================================================================
# Quantitative summary to stdout
# ===========================================================================

println("\n================ quantitative summary (final time) ================")
finalval(d, k) = haskey(d, k) ? d[k][end] : NaN
@printf("V_ice   : temp %.4g   enth %.4g   Δ %.3g%%\n",
        finalval(t1,"V_ice"), finalval(e1,"V_ice"),
        100*(finalval(e1,"V_ice")-finalval(t1,"V_ice"))/finalval(t1,"V_ice"))
@printf("A_ice   : temp %.4g   enth %.4g   Δ %.3g%%\n",
        finalval(t1,"A_ice"), finalval(e1,"A_ice"),
        100*(finalval(e1,"A_ice")-finalval(t1,"A_ice"))/finalval(t1,"A_ice"))
@printf("f_pmp   : temp %.4g   enth %.4g\n", finalval(t1,"f_pmp"), finalval(e1,"f_pmp"))
if haskey(e2,"T_prime_b") && haskey(t2,"T_prime_b")
    de = maskice(e2["T_prime_b"], Hmask_e); dt = maskice(t2["T_prime_b"], Hmask_t)
    diff = filter(isfinite, de .- dt)
    @printf("T'_b    : mean|Δ| %.4g K   max|Δ| %.4g K   (ice cells)\n",
            mean(abs, diff), maximum(abs, diff))
end
if haskey(e2,"omega_base")
    ob = filter(isfinite, maskice(e2["omega_base"], Hmask_e))
    @printf("ω_base  : max %.4g   mean(>0) %.4g   (enth)\n",
            maximum(ob; init=0.0), (any(ob .> 0) ? mean(ob[ob .> 0]) : 0.0))
end
println("===================================================================")
