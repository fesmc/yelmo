#!/usr/bin/env julia
# Parameterize and visualize the ytill cb_ref elevation-scaling methods.
#
# Mirrors the online calculation in src/physics/basal_dragging.f90
# (subroutine calc_cb_ref + calc_lambda_bed_lin / calc_lambda_bed_exp).
# With z_bed_sd ignored (single sample, n_sd=1) the per-cell value is:
#
#     scale_zb=1 (:lin)  lambda = clamp((z_bed - z0)/(z1 - z0), 0, 1)
#     scale_zb=2 (:exp)  lambda = min(exp((z_bed - z1)/(z1 - z0)), 1)
#     cb_ref            = max(cf_ref * lambda, cf_min)            (both)
#
# Note z_rel = z_bed in the Fortran (relative to present-day sea level),
# so this is purely a function of bedrock elevation.
#
# Main use case: pick exp (scale_zb=2) parameters, then find the lin
# (scale_zb=1) z0/z1 that best reproduce the same cb_ref(z_bed) curve.
#
# Usage:
#     julia --project=@. scripts/plot_till_scaling.jl [<outfile.png>]

using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Core scaling functions (faithful to basal_dragging.f90)
# ----------------------------------------------------------------------

lambda_lin(z_bed, z0, z1) = clamp((z_bed - z0) / (z1 - z0), 0.0, 1.0)

lambda_exp(z_bed, z0, z1) = min(exp((z_bed - z1) / (z1 - z0)), 1.0)

"cb_ref as a function of bedrock elevation for a given scaling method."
function cb_ref(z_bed; method::Symbol, z0, z1, cf_ref, cf_min)
    lambda = method === :lin ? lambda_lin(z_bed, z0, z1) :
             method === :exp ? lambda_exp(z_bed, z0, z1) :
             error("unknown method $method (use :lin or :exp)")
    return max(cf_ref * lambda, cf_min)
end

# ----------------------------------------------------------------------
# Fit lin z0/z1 to an exp reference curve (brute-force grid search)
# ----------------------------------------------------------------------

"""
    fit_lin_to_exp(; z0_exp, z1_exp, cf_ref, cf_min, zrange)

Find the (z0, z1) for the :lin method that best matches the :exp curve
defined by (z0_exp, z1_exp) over the bedrock-elevation grid `zrange`,
minimizing the sum of squared differences in cb_ref. Returns
(z0_best, z1_best, rmse).
"""
function fit_lin_to_exp(; z0_exp, z1_exp, cf_ref, cf_min, zrange,
                          z0_grid = range(-1000.0, 0.0; length = 201),
                          z1_grid = range(0.0, 1000.0; length = 201))
    target = [cb_ref(z; method = :exp, z0 = z0_exp, z1 = z1_exp,
                      cf_ref = cf_ref, cf_min = cf_min) for z in zrange]

    best = (z0 = NaN, z1 = NaN, sse = Inf)
    for z0 in z0_grid, z1 in z1_grid
        z1 <= z0 && continue
        sse = 0.0
        @inbounds for (k, z) in enumerate(zrange)
            d = cb_ref(z; method = :lin, z0 = z0, z1 = z1,
                       cf_ref = cf_ref, cf_min = cf_min) - target[k]
            sse += d * d
        end
        if sse < best.sse
            best = (z0 = z0, z1 = z1, sse = sse)
        end
    end
    rmse = sqrt(best.sse / length(zrange))
    return (best.z0, best.z1, rmse)
end

# ----------------------------------------------------------------------
# Configuration: edit these to explore parameter sets
# ----------------------------------------------------------------------

# Bedrock-elevation axis [m rel. to present-day sea level]
const ZRANGE = range(-1000.0, 1000.0; length = 801)

# Reference friction limits (shared across cases below)
const CF_REF = 0.1
const CF_MIN = 0.001

# Exp case to emulate (your current scale_zb=2 setup)
const EXP_Z0 = -100.0
const EXP_Z1 =  200.0

# Hand-picked lin cases to overlay for comparison
# (label, z0, z1)
const LIN_CASES = [
    ("lin z0=-100 z1=200", -100.0, 200.0),
    ("lin z0=0 z1=500",       0.0, 500.0),
]

# ----------------------------------------------------------------------
# Build plot
# ----------------------------------------------------------------------

function main(outfile)
    z = collect(ZRANGE)

    # Best-fit lin to the exp curve
    z0_fit, z1_fit, rmse = fit_lin_to_exp(; z0_exp = EXP_Z0, z1_exp = EXP_Z1,
                                            cf_ref = CF_REF, cf_min = CF_MIN,
                                            zrange = ZRANGE)
    @printf("Best-fit :lin for :exp(z0=%.1f, z1=%.1f): z0=%.1f, z1=%.1f (RMSE=%.4g)\n",
            EXP_Z0, EXP_Z1, z0_fit, z1_fit, rmse)

    fig = Figure(size = (900, 600))
    ax = Axis(fig[1, 1];
              xlabel = "bedrock elevation z_bed [m rel. PD sea level]",
              ylabel = "cb_ref [--]",
              title  = @sprintf("till cb_ref scaling  (cf_ref=%.3g, cf_min=%.3g)",
                                CF_REF, CF_MIN))

    # exp reference (thick black)
    cb_exp = [cb_ref(zb; method = :exp, z0 = EXP_Z0, z1 = EXP_Z1,
                     cf_ref = CF_REF, cf_min = CF_MIN) for zb in z]
    lines!(ax, z, cb_exp;
           color = :black, linewidth = 3,
           label = @sprintf(":exp z0=%.0f z1=%.0f (reference)", EXP_Z0, EXP_Z1))

    # best-fit lin (dashed red)
    cb_fit = [cb_ref(zb; method = :lin, z0 = z0_fit, z1 = z1_fit,
                     cf_ref = CF_REF, cf_min = CF_MIN) for zb in z]
    lines!(ax, z, cb_fit;
           color = :red, linewidth = 2.5, linestyle = :dash,
           label = @sprintf(":lin best-fit z0=%.0f z1=%.0f", z0_fit, z1_fit))

    # hand-picked lin overlays
    for (label, z0, z1) in LIN_CASES
        cb = [cb_ref(zb; method = :lin, z0 = z0, z1 = z1,
                     cf_ref = CF_REF, cf_min = CF_MIN) for zb in z]
        lines!(ax, z, cb; linewidth = 1.5, label = label)
    end

    axislegend(ax; position = :rb, framevisible = true)
    mkpath(dirname(abspath(outfile)))
    save(outfile, fig)
    println("Wrote $outfile")
    return nothing
end

outfile = length(ARGS) >= 1 ? ARGS[1] : "till_scaling.png"
main(outfile)
