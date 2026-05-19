#!/usr/bin/env julia
# Plot CalvingMIP front diagnostics (position, speed, thickness) along
# selected profiles for exp2 and exp4, using the standard CalvingMIP
# output file written by yelmo_calving.f90:
#
#     <expdir>/CalvingMIP_EXP{N}_YELMO_AWI.nc
#
# Per-profile variables read: Time1, xcf<P>, ycf<P>, lithkcf<P>,
# xvelmeancf<P>, yvelmeancf<P>.
#
# Usage:
#     julia --project=@. scripts/plot_calvingmip.jl [<output-folder>]
#
# Default <output-folder> is tmp/yelmo-calvingmip-2026-05-19.

using CairoMakie
using NCDatasets

# Profile origins per CalvingMIP wiki spec. Front "position" is the
# Euclidean distance from this origin to (xcf, ycf).
const PROFILE_ORIGIN = Dict(
    "A" => (0.0, 0.0), "B" => (0.0, 0.0), "C" => (0.0, 0.0), "D" => (0.0, 0.0),
    "E" => (0.0, 0.0), "F" => (0.0, 0.0), "G" => (0.0, 0.0), "H" => (0.0, 0.0),
    "CapA" => (-390e3, 0.0), "CapB" => ( 390e3, 0.0),
    "CapC" => (-390e3, 0.0), "CapD" => ( 390e3, 0.0),
    "HalA" => (-150e3, 0.0), "HalB" => ( 150e3, 0.0),
    "HalC" => (-150e3, 0.0), "HalD" => ( 150e3, 0.0),
)

function load_front_series(ds, prof)
    t     = Array(ds["Time1"])
    xcf   = Array(ds["xcf"*prof])
    ycf   = Array(ds["ycf"*prof])
    lithk = Array(ds["lithkcf"*prof])
    ux    = Array(ds["xvelmeancf"*prof])
    uy    = Array(ds["yvelmeancf"*prof])

    x0, y0 = PROFILE_ORIGIN[prof]
    dist   = @. sqrt((xcf - x0)^2 + (ycf - y0)^2) / 1e3   # km
    speed  = @. sqrt(ux^2 + uy^2)                          # m/yr

    return (; t, dist, speed, thk = lithk)
end

# CalvingMIP exp2/exp4 forcing imposes a front-normal velocity anomaly
#     wv(t) = -300 sin(2π t / 1000)   [m yr⁻¹]
# on top of the steady-state front (cf. calvmip_exp2 in calving_ac.f90).
# Time-integrating wv from t=0 gives a closed-form front displacement,
# returned here in km. r0 anchors the curve at the simulated t=0 front.
function exp2_front_analytical(t, r0_km)
    return r0_km .+ (150.0 / π) .* (cos.(2π .* t ./ 1000.0) .- 1.0)
end

function plot_experiment(ncpath, profiles, outpath;
                        suptitle = "", with_analytical = false)
    isfile(ncpath) || error("CalvingMIP file not found: $ncpath")

    NCDataset(ncpath, "r") do ds
        series = [load_front_series(ds, p) for p in profiles]

        fig = Figure(size = (900, 800))
        if !isempty(suptitle)
            Label(fig[0, :], suptitle; fontsize = 18, font = :bold)
        end

        rowlabels = ("Front position [km]",
                     "Front speed [m yr⁻¹]",
                     "Front thickness [m]")
        fields    = (:dist, :speed, :thk)

        for (col, (prof, s)) in enumerate(zip(profiles, series))
            for (row, (ylabel, field)) in enumerate(zip(rowlabels, fields))
                ax = Axis(fig[row, col];
                    xlabel = row == length(fields) ? "Time [yr]" : "",
                    ylabel = col == 1 ? ylabel : "",
                    title  = row == 1 ? "Profile $prof" : "")

                if row == 1 && with_analytical
                    r_an = exp2_front_analytical(s.t, first(s.dist))
                    band!(ax, s.t, r_an .- 5.0, r_an .+ 5.0;
                          color = (:tomato, 0.25), label = "Analytical ±5 km")
                    lines!(ax, s.t, r_an;
                           color = :tomato, linestyle = :dash,
                           label = "Analytical")
                    lines!(ax, s.t, s.dist;
                           color = :black, label = "Yelmo")
                    axislegend(ax; position = :rb, framevisible = false,
                               labelsize = 10, padding = (4, 4, 4, 4))
                else
                    lines!(ax, s.t, getfield(s, field); color = :black)
                end
            end
        end

        # Force equal column widths regardless of axis-decoration size.
        for col in 1:length(profiles)
            colsize!(fig.layout, col, Relative(1.0 / length(profiles)))
        end

        save(outpath, fig)
        @info "Wrote $outpath"
    end
    return outpath
end

function main()
    fldr = length(ARGS) >= 1 ? ARGS[1] : "tmp/yelmo-calvingmip-2026-05-19"

    plot_experiment(
        joinpath(fldr, "exp2", "CalvingMIP_EXP2_YELMO_AWI.nc"),
        ["A", "B"],
        joinpath(fldr, "exp2", "calvingmip_front_exp2.pdf");
        suptitle = "CalvingMIP EXP2 — front diagnostics",
        with_analytical = true,
    )

    plot_experiment(
        joinpath(fldr, "exp4", "CalvingMIP_EXP4_YELMO_AWI.nc"),
        ["CapA", "HalA"],
        joinpath(fldr, "exp4", "calvingmip_front_exp4.pdf");
        suptitle = "CalvingMIP EXP4 — front diagnostics",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
