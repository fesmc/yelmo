# compare_column_analytic.jl
#
# Referee the Yelmo temperature ("temp") and enthalpy ("enth") column solvers
# against the exact 1D advection-diffusion solutions of IceColumnSolutions.jl
# (Moreno-Parada et al., 2024) in the constant-property limit.
#
# Companion to docs/physics/enth-temp-discretization.md. In the constant-cp/kt
# regime both Yelmo solvers *should* reproduce the analytic column exactly; this
# script confirms that, measures each scheme's order of accuracy on an nz
# refinement sweep, and repeats across a (Pe, gamma) parameter grid. (The
# variable-cp(T) divergence discussed in the note lies outside the analytic's
# linear/constant-coefficient validity domain.)
#
# The driver orchestrates the Fortran harness directly: for each configuration
# it runs the "ics" experiment of tests/test_icetemp.x (constant properties,
# linear vertical velocity, cold column, geothermal flux prescribed at the base,
# uniform zeta grid) and reads back the steady profile.
#
# Usage (from the worktree root, after `make icetemp`):
#   julia --project=analysis analysis/compare_column_analytic.jl

using NCDatasets
using IceColumnSolutions
using Printf
using Statistics
using CairoMakie

const HERE   = @__DIR__
const ROOT   = normpath(joinpath(HERE, ".."))
const OUTDIR = joinpath(ROOT, "output")
const FIGDIR = joinpath(HERE, "figures")
const BIN    = joinpath(ROOT, "libyelmo", "bin", "test_icetemp.x")
mkpath(FIGDIR)

# EISMINT reference values (input/yelmo_phys_const.nml, &EISMINT). The ics run
# uses constant cp/kt; rho_ice is not written to output, so pin it here.
const RHO_ICE  = 910.0        # [kg m-3]
const SEC_YEAR = 31556926.0   # [s a-1]

const NZSWEEP = [26, 51, 101, 201, 401]
const SOLVERS = ["temp", "enth"]

const COL = Dict("temp" => RGBf(0.17, 0.44, 0.71),   # steel blue
                 "enth" => RGBf(0.75, 0.22, 0.17))    # deep red

# --- run the ics harness and read back the steady profile + parameters --------
function run_ics(solver, nz; smb, qgeo)
    isfile(BIN) || error("binary not found: $BIN (run `make icetemp` first)")
    run(pipeline(`$BIN $nz 1e-3 $solver ics $smb $qgeo`; stdout = devnull, stderr = devnull))
    fn = joinpath(OUTDIR, "test_ics_$(solver)_nz$(nz).nc")
    ds = NCDataset(fn)
    zeta = Float64.(ds["zeta"][:])
    T    = Float64.(ds["T_ice"][:, end])     # [K] final (steady) time
    Tpmp = Float64.(ds["T_pmp"][:, end])     # [K]
    H    = Float64(ds["H_ice"][end])
    Tsrf = Float64(ds["T_srf"][end])
    smb_ = Float64(ds["smb"][end])
    Qgeo = Float64(ds["Q_geo"][end])
    cp   = Float64(ds["cp"][1, end])
    kt   = Float64(ds["kt"][1, end])
    close(ds)
    return (; zeta, T, Tpmp, H, Tsrf, smb = smb_, Qgeo, cp, kt)
end

# --- map Yelmo parameters -> IceColumnPar and evaluate the exact steady T(zeta) ---
# Constant properties, Dirichlet surface (beta=0). Accumulation is Pe > 0 in the
# IceColumnSolutions convention (real-erf branch), so w0 = +smb. The analytic
# grid range(0,1,nz) coincides node-for-node with Yelmo's linear ics zeta grid.
function analytic_solution(d)
    k     = d.kt / SEC_YEAR              # [W m-1 K-1]
    kappa = k / (RHO_ICE * d.cp)         # [m2 s-1]
    w0    = d.smb                        # [m a-1]
    G     = d.Qgeo * 1e-3               # [mW m-2] -> [W m-2]
    par   = IceColumnPar(d.H, d.Tsrf, kappa, k, 0.0, G; w0 = w0)
    return solve_stationary(par; nz = length(d.zeta)).T_eq, par
end

# base at least `tol` K below pressure melting => analytic (no pmp cap) applies
is_cold(d; tol = 0.5) = (d.Tpmp[1] - d.T[1]) > tol

# error ~ C*dz^p, dz ∝ 1/(nz-1); slope of log(err) vs log(dz)
function fit_order(nz, err)
    logdz  = log.(1.0 ./ (nz .- 1))
    logerr = log.(err)
    n = length(nz)
    (n * sum(logdz .* logerr) - sum(logdz) * sum(logerr)) /
        (n * sum(logdz .^ 2) - sum(logdz)^2)
end

rms(e) = sqrt(mean(e .^ 2))

# ===========================================================================
# Part A: convergence study at the default configuration
# ===========================================================================
println("\n########## Part A: convergence at default config (smb=0.1, Qgeo=42) ##########")

convA = Dict{String,Any}()
local parA
for solver in SOLVERS
    rmss = Float64[]; maxs = Float64[]
    for nz in NZSWEEP
        d = run_ics(solver, nz; smb = 0.1, qgeo = 42.0)
        Ta, par = analytic_solution(d)
        global parA = par
        e = d.T .- Ta
        push!(rmss, rms(e)); push!(maxs, maximum(abs, e))
    end
    convA[solver] = (rms = rmss, maxerr = maxs)
end

@printf("Pe = %.3f   gamma = %.3f\n", parA.Pe, parA.gamma)
@printf("%-6s | %4s | %12s | %12s\n", "solver", "nz", "RMS err [K]", "max err [K]")
println("-"^46)
for solver in SOLVERS, i in eachindex(NZSWEEP)
    @printf("%-6s | %4d | %12.3e | %12.3e\n",
            solver, NZSWEEP[i], convA[solver].rms[i], convA[solver].maxerr[i])
end
println("-"^46)
for solver in SOLVERS
    @printf("%-6s : convergence order (RMS) = %.2f\n", solver, fit_order(NZSWEEP, convA[solver].rms))
end

# ===========================================================================
# Part B: (Pe, gamma) parameter sweep  [smb sets Pe, Qgeo sets gamma]
# ===========================================================================
println("\n########## Part B: (Pe, gamma) parameter sweep ##########")

const SMB_SWEEP  = [0.0, 0.05, 0.1, 0.2, 0.4]
const QGEO_SWEEP = [15.0, 30.0, 45.0]

sweeprows = NamedTuple[]
skipped   = String[]
for smb in SMB_SWEEP, qgeo in QGEO_SWEEP
    res = Dict{String,Any}()
    Pe = NaN; gam = NaN; ok = true
    for solver in SOLVERS
        rmss = Float64[]
        for nz in NZSWEEP
            d = run_ics(solver, nz; smb = smb, qgeo = qgeo)
            if !is_cold(d)
                ok = false; break
            end
            Ta, par = analytic_solution(d)
            Pe = par.Pe; gam = par.gamma
            push!(rmss, rms(d.T .- Ta))
        end
        ok || break
        res[solver] = rmss
    end
    if ok
        push!(sweeprows, (smb = smb, qgeo = qgeo, Pe = Pe, gamma = gam,
                          rms_temp = res["temp"][end], rms_enth = res["enth"][end],
                          ord_temp = fit_order(NZSWEEP, res["temp"]),
                          ord_enth = fit_order(NZSWEEP, res["enth"])))
    else
        push!(skipped, @sprintf("smb=%.2f Qgeo=%.0f (base reached pmp)", smb, qgeo))
    end
end

@printf("%6s | %6s | %6s | %11s | %11s | %8s | %8s\n",
        "smb", "Qgeo", "Pe", "RMS temp", "RMS enth", "ord temp", "ord enth")
println("-"^72)
for r in sweeprows
    @printf("%6.2f | %6.0f | %6.2f | %11.3e | %11.3e | %8.2f | %8.2f\n",
            r.smb, r.qgeo, r.Pe, r.rms_temp, r.rms_enth, r.ord_temp, r.ord_enth)
end
println("-"^72)
if !isempty(skipped)
    println("skipped (outside analytic validity — base at pressure melting):")
    for s in skipped; println("  - ", s); end
end

# ===========================================================================
# Figures
# ===========================================================================
fig = Figure(size = (1050, 430))

ax1 = Axis(fig[1, 1]; title = "Convergence at default config",
           xlabel = "nz", ylabel = "RMS |T − analytic| [K]",
           xscale = log10, yscale = log10)
for solver in SOLVERS
    scatterlines!(ax1, Float64.(NZSWEEP), convA[solver].rms; color = COL[solver],
                  linewidth = 2, label = @sprintf("%s (order %.2f)", solver,
                  fit_order(NZSWEEP, convA[solver].rms)))
end
# 2nd-order reference slope
let n = Float64.(NZSWEEP)
    ref = convA["temp"].rms[1] .* (n[1] ./ n) .^ 2
    lines!(ax1, n, ref; color = :gray, linestyle = :dash, label = "order 2 ref")
end
axislegend(ax1; position = :lb, framevisible = false)

ax2 = Axis(fig[1, 2]; title = "Accuracy across the sweep",
           xlabel = "Péclet number Pe", ylabel = "RMS |T − analytic| [K] (nz=$(NZSWEEP[end]))",
           yscale = log10)
scatter!(ax2, [r.Pe for r in sweeprows], [r.rms_temp for r in sweeprows];
         color = COL["temp"], label = "temp", markersize = 10)
scatter!(ax2, [r.Pe for r in sweeprows], [r.rms_enth for r in sweeprows];
         color = COL["enth"], label = "enth", markersize = 10, marker = :diamond)
axislegend(ax2; position = :lt, framevisible = false)

Label(fig[0, :], "Yelmo temp/enth column solvers vs analytic (constant properties)",
      fontsize = 16, font = :bold)
p = joinpath(FIGDIR, "column_analytic.png")
save(p, fig)
println("\nwrote ", p)
