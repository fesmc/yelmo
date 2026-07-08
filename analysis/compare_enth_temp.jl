# compare_enth_temp.jl
#
# Compare the enthalpy (method="enth") and temperature (method="temp")
# thermodynamics solvers on the 2D EISMINT-2 benchmark experiments A and F.
#
# For each experiment it reads the yelmo1D.nc (time series) and yelmo2D.nc
# (maps) output of an enth and a temp run and writes publication-style figures
# (Makie/CairoMakie, styled after Robinson et al., 2020) into analysis/figures/,
# plus a Robinson-2020-style summary table (analysis/results_eismint.{md,csv}).
#
# References:
#   Robinson et al. (2020), GMD 13:2805-2823   - Yelmo description/validation
#   Greve & Blatter (2016), Polar Sci. 10:11-23 - enthalpy schemes on EISMINT-2
#   Payne et al. (2000), J. Glaciol. 46:227-238 - EISMINT-2 intercomparison
#
# Usage (from the worktree root):
#   julia analysis/compare_enth_temp.jl
# Run directories are configured in EXPERIMENTS below.
#
# analysis/figures/ is intentionally untracked (see analysis/.gitignore).

using CairoMakie
using NCDatasets
using Statistics
using Printf

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

const HERE   = @__DIR__
const ROOT   = normpath(joinpath(HERE, ".."))
const FIGDIR = joinpath(HERE, "figures")
mkpath(FIGDIR)

# Each experiment: enth/temp run dirs and (optional) EISMINT-2 reference values.
# Reference = Payne et al. (2000) EISMINT-2 intercomparison means. Exp A is
# well established; Exp F numeric references are not reproduced here (nothing).
# NOTE: the Yelmo EISMINT par uses hybrid dynamics *with sliding* whereas the
# Payne reference is pure-SIA no-slip, so absolute values are not expected to
# match exactly - the temp<->enth difference is the informative quantity.
const EXPERIMENTS = [
    (name = "EISMINT-A", key = "A",
     enth = joinpath(ROOT, "RUNDIR"),
     temp = joinpath(ROOT, "RUNDIR_temp"),
     ref  = (H0 = 3990.0, Tpb = -13.34, V = 2.128, A = 1.034, fmelt = 0.718)),
    (name = "EISMINT-F", key = "F",
     enth = joinpath(ROOT, "RUNDIR_expf_enth"),
     temp = joinpath(ROOT, "RUNDIR_expf_temp"),
     ref  = nothing),
]

CairoMakie.activate!(type = "png")

# --- Robinson et al. (2020) style theme -----------------------------------
const COL_ENTH = RGBf(0.75, 0.22, 0.17)   # deep red
const COL_TEMP = RGBf(0.17, 0.44, 0.71)   # steel blue
const COL_REF  = RGBf(0.35, 0.35, 0.35)   # grey (EISMINT reference)

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
    for v in ("V_ice", "A_ice", "f_pmp", "H_ice", "W_til")
        haskey(ds, v) && (out[v] = ncget(ds, v))
    end
    close(ds)
    return out
end

"Load final-time 2D maps, a mid-row cross-section, and scalar summit/integrated diagnostics."
function load_2d(dir)
    ds = NCDataset(joinpath(dir, "yelmo2D.nc"))
    out = Dict{String,Any}()
    out["xc"] = ncget(ds, "xc")
    out["yc"] = ncget(ds, "yc")
    out["zeta"] = ncget(ds, "zeta")
    nt = ds.dim["time"]
    out["time_end"] = ncget(ds, "time")[nt] / 1e3   # kyr
    nx, ny = length(out["xc"]), length(out["yc"])
    for v in ("H_ice", "T_prime_b", "f_pmp", "bmb_grnd", "H_cts")
        haskey(ds, v) && (out[v] = tof64(ds[v][:, :, nt]))
    end
    jc = (ny + 1) ÷ 2
    for v in ("omega", "T_prime")
        if haskey(ds, v)
            out[v*"_base"] = tof64(ds[v][:, :, 1, nt])
            out[v*"_xsec"] = tof64(ds[v][:, jc, :, nt])
        end
    end
    out["jc"] = jc
    # --- scalar diagnostics ---
    ic = (nx + 1) ÷ 2
    H = out["H_ice"]
    ice = H .> 1.0
    out["H0"]  = H[ic, jc]                             # summit thickness [m]
    out["Tpb"] = get(out, "T_prime_b", fill(NaN,nx,ny))[ic, jc]  # summit homol. temp
    if haskey(out, "H_cts")
        hc = out["H_cts"]
        out["Vtemp_frac"] = sum(hc[ice]) / sum(H[ice]) * 100      # temperate vol %
        out["Hcts_max"]   = maximum(hc[ice]; init = 0.0)          # max CTS [m]
    end
    close(ds)
    return out
end

maskice(field, H; thr = 1.0) = [H[i] > thr ? field[i] : NaN for i in eachindex(field)]

"Deviation of a 2D field from its azimuthal (radial-ring) mean over ice cells;
reveals grid-aligned spokes. Returns (asym_map, rms_asym)."
function azimuthal_asym(field, xc, yc, H; thr = 1.0)
    nx, ny = length(xc), length(yc)
    dr = abs(xc[2] - xc[1])
    R  = [hypot(xc[i], yc[j]) for i in 1:nx, j in 1:ny]
    nb = ceil(Int, maximum(R) / dr) + 1
    binof(r) = clamp(floor(Int, r / dr) + 1, 1, nb)
    rsum = zeros(nb); rcnt = zeros(Int, nb)
    for i in 1:nx, j in 1:ny
        (H[i, j] > thr && isfinite(field[i, j])) || continue
        b = binof(R[i, j]); rsum[b] += field[i, j]; rcnt[b] += 1
    end
    rmean = [rcnt[b] > 0 ? rsum[b] / rcnt[b] : NaN for b in 1:nb]
    asym = fill(NaN, nx, ny)
    for i in 1:nx, j in 1:ny
        (H[i, j] > thr && isfinite(field[i, j])) || continue
        asym[i, j] = field[i, j] - rmean[binof(R[i, j])]
    end
    vals = filter(isfinite, asym)
    return asym, (isempty(vals) ? 0.0 : sqrt(mean(vals .^ 2)))
end

rs(x, y, z) = reshape(x, length(y), length(z))
noaxis!(a) = hidedecorations!(a, label = false)

# ===========================================================================
# Figures
# ===========================================================================

function fig_timeseries(ex, e1, t1)
    specs = [("V_ice", "Ice volume [10⁶ km³]", ex.ref === nothing ? nothing : ex.ref.V),
             ("A_ice", "Ice area [10⁶ km²]",   ex.ref === nothing ? nothing : ex.ref.A),
             ("f_pmp", "Basal temperate fraction [-]", ex.ref === nothing ? nothing : ex.ref.fmelt),
             ("W_til", "Mean basal water W_til [m]", nothing)]
    fig = Figure(size = (1000, 720))
    for (n, (k, ylab, refv)) in enumerate(specs)
        ax = Axis(fig[(n-1)÷2+1, (n-1)%2+1]; title = ylab, xlabel = "time [kyr]")
        if refv !== nothing
            hlines!(ax, [refv]; color = COL_REF, linestyle = :dash, linewidth = 1.5,
                    label = "EISMINT ref")
        end
        haskey(t1, k) && lines!(ax, t1["time"], t1[k]; color = COL_TEMP, label = "temp")
        haskey(e1, k) && lines!(ax, e1["time"], e1[k]; color = COL_ENTH, label = "enth")
        n == 1 && axislegend(ax; position = :rb, framevisible = false, rowgap = 0)
    end
    Label(fig[0, :], "$(ex.name): enthalpy vs temperature solver — time series",
          fontsize = 17, font = :bold)
    p = joinpath(FIGDIR, "$(ex.key)_timeseries.png"); save(p, fig); println("wrote ", p)
end

function fig_maps(ex, e2, t2)
    specs = [("H_ice", "Ice thickness [m]", :ice, false),
             ("T_prime_b", "Basal homol. T′_b [°C]", :thermal, false),
             ("f_pmp", "Basal temperate frac [-]", :dense, false),
             ("bmb_grnd", "Grounded bmb [m/yr]", :curl, true)]
    xc, yc = e2["xc"], e2["yc"]
    He, Ht = e2["H_ice"], t2["H_ice"]
    fig = Figure(size = (1180, 1200))
    for (r, (k, lab, cmap, symrange)) in enumerate(specs)
        (haskey(e2, k) && haskey(t2, k)) || continue
        fe = maskice(e2[k], He); ft = maskice(t2[k], Ht); dd = fe .- ft
        vals = filter(isfinite, vcat(fe, ft))
        lo, hi = isempty(vals) ? (0.0, 1.0) : (minimum(vals), maximum(vals))
        symrange && (m = maximum(abs, vals; init = 1.0); (lo, hi) = (-m, m))
        dvals = filter(isfinite, dd)
        dm = isempty(dvals) ? 1e-12 : max(maximum(abs, dvals), 1e-12)
        ax1 = Axis(fig[r, 1]; aspect = 1, title = r == 1 ? "temp" : "", ylabel = lab)
        ax2 = Axis(fig[r, 2]; aspect = 1, title = r == 1 ? "enth" : "")
        ax3 = Axis(fig[r, 3]; aspect = 1, title = r == 1 ? "enth − temp" : "")
        for a in (ax1, ax2, ax3); noaxis!(a); end
        h1 = heatmap!(ax1, xc, yc, rs(ft, xc, yc); colormap = cmap, colorrange = (lo, hi))
        heatmap!(ax2, xc, yc, rs(fe, xc, yc); colormap = cmap, colorrange = (lo, hi))
        h3 = heatmap!(ax3, xc, yc, rs(dd, xc, yc); colormap = :balance, colorrange = (-dm, dm))
        for (a, H) in ((ax1, Ht), (ax2, He), (ax3, He))
            contour!(a, xc, yc, rs(H, xc, yc); levels = [1.0], color = :black, linewidth = 0.6)
        end
        Colorbar(fig[r, 4], h1); Colorbar(fig[r, 5], h3)
    end
    Label(fig[0, :], @sprintf("%s maps at %.0f kyr — temp | enth | difference", ex.name, e2["time_end"]),
          fontsize = 17, font = :bold)
    p = joinpath(FIGDIR, "$(ex.key)_maps.png"); save(p, fig); println("wrote ", p)
end

function fig_polythermal(ex, e2)
    xc, yc = e2["xc"], e2["yc"]; He = e2["H_ice"]
    fig = Figure(size = (1150, 470))
    if haskey(e2, "omega_base")
        ax = Axis(fig[1, 1]; aspect = 1, title = "Basal water ω (enth)", xlabel = "x [km]", ylabel = "y [km]")
        h = heatmap!(ax, xc, yc, rs(maskice(e2["omega_base"], He), xc, yc); colormap = :dense)
        contour!(ax, xc, yc, rs(He, xc, yc); levels=[1.0], color=:black, linewidth=0.6)
        Colorbar(fig[1, 2], h)
    end
    if haskey(e2, "H_cts")
        ax = Axis(fig[1, 3]; aspect = 1, title = "CTS height H_cts [m] (enth)", xlabel = "x [km]")
        h = heatmap!(ax, xc, yc, rs(maskice(e2["H_cts"], He), xc, yc); colormap = :haline)
        contour!(ax, xc, yc, rs(He, xc, yc); levels=[1.0], color=:black, linewidth=0.6)
        Colorbar(fig[1, 4], h)
    end
    if haskey(e2, "omega_xsec")
        ax = Axis(fig[1, 5]; title = "ω cross-section (y=0)", xlabel = "x [km]", ylabel = "ζ (0=base,1=surf)")
        om = copy(e2["omega_xsec"]); Hrow = He[:, e2["jc"]]
        for i in axes(om, 1), k in axes(om, 2); Hrow[i] > 1.0 || (om[i, k] = NaN); end
        h = heatmap!(ax, xc, e2["zeta"], om; colormap = :dense)
        Colorbar(fig[1, 6], h; label = "ω [-]")
    end
    Label(fig[0, :], @sprintf("%s polythermal structure (enth) at %.0f kyr", ex.name, e2["time_end"]),
          fontsize = 17, font = :bold)
    p = joinpath(FIGDIR, "$(ex.key)_polythermal_enth.png"); save(p, fig); println("wrote ", p)
end

"Symmetry / spoke diagnostic: basal homologous temperature and its deviation
from the azimuthal (radial) mean, for temp and enth."
function fig_symmetry(ex, e2, t2)
    xc, yc = e2["xc"], e2["yc"]; He, Ht = e2["H_ice"], t2["H_ice"]
    fe, rmse = azimuthal_asym(e2["T_prime_b"], xc, yc, He)
    ft, rmst = azimuthal_asym(t2["T_prime_b"], xc, yc, Ht)
    am = max(maximum(abs, filter(isfinite, vcat(vec(fe), vec(ft)))), 1e-6)
    tb_e = maskice(e2["T_prime_b"], He); tb_t = maskice(t2["T_prime_b"], Ht)
    vals = filter(isfinite, vcat(tb_e, tb_t)); lo, hi = minimum(vals), maximum(vals)
    fig = Figure(size = (900, 900))
    rows = (("temp", tb_t, ft, rmst, Ht), ("enth", tb_e, fe, rmse, He))
    for (r, (nm, tb, asym, rms, H)) in enumerate(rows)
        ax1 = Axis(fig[r, 1]; aspect = 1, ylabel = nm, title = r == 1 ? "Basal homol. T′_b [°C]" : "")
        ax2 = Axis(fig[r, 2]; aspect = 1, title = r == 1 ? "T′_b − azimuthal mean [°C]" : "")
        for a in (ax1, ax2); noaxis!(a); end
        h1 = heatmap!(ax1, xc, yc, rs(tb, xc, yc); colormap = :thermal, colorrange = (lo, hi))
        h2 = heatmap!(ax2, xc, yc, rs(asym, xc, yc); colormap = :balance, colorrange = (-am, am))
        contour!(ax1, xc, yc, rs(H, xc, yc); levels=[1.0], color=:black, linewidth=0.6)
        contour!(ax2, xc, yc, rs(H, xc, yc); levels=[1.0], color=:black, linewidth=0.6)
        text!(ax2, 0.02, 0.02; text = @sprintf("RMS asym = %.3f °C", rms),
              space = :relative, align = (:left, :bottom), fontsize = 11)
        r == 1 && Colorbar(fig[1, 3], h1)
        r == 2 && Colorbar(fig[2, 3], h2)
    end
    Label(fig[0, :], @sprintf("%s basal-temperature symmetry at %.0f kyr", ex.name, e2["time_end"]),
          fontsize = 17, font = :bold)
    p = joinpath(FIGDIR, "$(ex.key)_symmetry.png"); save(p, fig); println("wrote ", p)
    return (temp = rmst, enth = rmse)
end

# ===========================================================================
# Main: per-experiment figures + collected diagnostics
# ===========================================================================

rows = Vector{NamedTuple}()   # table rows

for ex in EXPERIMENTS
    isdir(ex.enth) && isdir(ex.temp) || (@warn "missing run dirs for $(ex.name), skipping"; continue)
    println("\n=== $(ex.name) ===")
    e1, t1 = load_1d(ex.enth), load_1d(ex.temp)
    e2, t2 = load_2d(ex.enth), load_2d(ex.temp)
    fig_timeseries(ex, e1, t1)
    fig_maps(ex, e2, t2)
    fig_polythermal(ex, e2)
    sym = fig_symmetry(ex, e2, t2)
    for (nm, d1, d2, r) in (("temp", t1, t2, 0.0), ("enth", e1, e2, 0.0))
        push!(rows, (exp = ex.name, solver = nm,
                     H0 = d2["H0"], Tpb = d2["Tpb"],
                     V = d1["V_ice"][end], A = d1["A_ice"][end],
                     fmelt = d1["f_pmp"][end],
                     Vtemp = get(d2, "Vtemp_frac", 0.0), Hcts = get(d2, "Hcts_max", 0.0),
                     symrms = nm == "temp" ? sym.temp : sym.enth))
    end
    if ex.ref !== nothing
        push!(rows, (exp = ex.name, solver = "EISMINT ref", H0 = ex.ref.H0, Tpb = ex.ref.Tpb,
                     V = ex.ref.V, A = ex.ref.A, fmelt = ex.ref.fmelt,
                     Vtemp = NaN, Hcts = NaN, symrms = NaN))
    end
end

# ===========================================================================
# Table (Robinson et al. 2020 style) -> markdown + csv + stdout
# ===========================================================================

fmt(x; d = 2) = (x isa Number && !isnan(x)) ? string(round(x; digits = d)) : "—"

hdr = ["Experiment", "Solver", "Summit H₀ [m]", "Summit T′_b [°C]",
       "Volume [10⁶ km³]", "Area [10⁶ km²]", "Melt frac [-]",
       "Temp. vol [%]", "Max CTS [m]", "T′_b asym RMS [°C]"]

rowstrs = [[r.exp, r.solver, fmt(r.H0; d=0), fmt(r.Tpb; d=2), fmt(r.V; d=3),
            fmt(r.A; d=3), fmt(r.fmelt; d=3), fmt(r.Vtemp; d=2), fmt(r.Hcts; d=1),
            fmt(r.symrms; d=3)] for r in rows]

open(joinpath(HERE, "results_eismint.md"), "w") do io
    println(io, "# EISMINT-2 enthalpy vs temperature — summary\n")
    println(io, "Yelmo `therm-dev`, 100 ka, dx=25 km, hybrid dynamics (with sliding).")
    println(io, "Reference: Payne et al. (2000) EISMINT-2 means (pure-SIA no-slip; Exp A only).")
    println(io, "Style after Robinson et al. (2020); enthalpy context: Greve & Blatter (2016).\n")
    println(io, "| " * join(hdr, " | ") * " |")
    println(io, "|" * repeat(" --- |", length(hdr)))
    for rr in rowstrs; println(io, "| " * join(rr, " | ") * " |"); end
end
open(joinpath(HERE, "results_eismint.csv"), "w") do io
    println(io, join(hdr, ","))
    for rr in rowstrs; println(io, join(rr, ",")); end
end

println("\n================ EISMINT-2 summary table ================")
println(join(hdr, " | "))
for rr in rowstrs; println(join(rr, " | ")); end
println("wrote ", joinpath(HERE, "results_eismint.md"), " (+ .csv)")
