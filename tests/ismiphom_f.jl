# ISMIP-HOM Experiment F: Yelmo steady state vs the models of Pattyn et al. (2008)
# along the central flowline (y = 0).
#
# usage: julia tests/ismiphom_f.jl <out.png> <exp>:<label>=<run> [<exp>:<label>=<run> ...]
#   <exp>   F1 or F2
#   <run>   a Yelmo run directory (containing yelmo.nc) or a NetCDF file with
#           xc, yc, z_srf and ux_s, uy_s; the last time record is used
# example:
#   julia tests/ismiphom_f.jl ismiphom-f.png F1:DIVA=output/ismiphom-f1 F2:DIVA=output/ismiphom-f2
#
# Prints the hump/trough of z_s and the extremes of |u_s| with their positions,
# and the rms difference to the ensemble mean and to full Stokes (cma1).
# Reference data: tests/data/ISMIPHOM-F (see its README.md).

using NCDatasets, CairoMakie, DelimitedFiles, Statistics, Printf

const REF    = joinpath(@__DIR__, "data", "ISMIPHOM-F")
const MODELS = ["cma1", "cma2", "fpa1", "fsa1", "mbr1", "mtk1", "oga1"]
const FS     = ["cma1", "oga1"]     # full-Stokes models (Pattyn et al., 2008, Table 3)
const L      = 100.0                # [km] domain size (period)

"Reference profile along y = 0: x [km], z_s [m], |v_s| [m/a] (horizontal)."
function ref_profile(model, exp)
    d = readdlm(joinpath(REF, "$(model)f00$(exp - 1).txt"))
    x, y = d[:, 1], d[:, 2]
    if model in ("cma1", "cma2")    # normalised [0,1] coordinates, bump at 0.5
        x = (x .- 0.5) .* L
        y = (y .- 0.5) .* L
    end
    s  = sqrt.(d[:, 4] .^ 2 .+ d[:, 5] .^ 2)
    yr = round.(y; digits=2)        # mtk1 coordinates carry ~1e-3 km jitter
    yu = unique(yr)
    y0 = minimum(abs.(yu))
    # Rows nearest to y = 0 (fpa1 has y = +-1.28 km only): average them
    rows = [findall(yr .== v) for v in yu[abs.(abs.(yu) .- y0) .< 1e-6]]
    ord  = [r[sortperm(x[r])] for r in rows]
    xp = mean(hcat([x[o] for o in ord]...), dims=2)[:, 1]
    zp = mean(hcat([d[o, 3] for o in ord]...), dims=2)[:, 1]
    sp = mean(hcat([s[o] for o in ord]...), dims=2)[:, 1]
    return (x=xp, zs=zp, us=sp)
end

"Yelmo profile along y = 0 from the last time record."
function yelmo_profile(run)
    f  = isdir(run) ? joinpath(run, "yelmo.nc") : run
    ds = NCDataset(f)
    x  = Float64.(ds["xc"][:])
    y  = Float64.(ds["yc"][:])
    zs = Float64.(ds["z_srf"][:, :, end])
    ux = Float64.(ds["ux_s"][:, :, end])
    uy = Float64.(ds["uy_s"][:, :, end])
    close(ds)
    # ux_s on acx-nodes (i+1/2), uy_s on acy-nodes (j+1/2): to aa-nodes (periodic)
    uxa = 0.5 .* (ux .+ circshift(ux, (1, 0)))
    uya = 0.5 .* (uy .+ circshift(uy, (0, 1)))
    j0  = argmin(abs.(y))
    return (x=x, zs=zs[:, j0], us=sqrt.(uxa[:, j0] .^ 2 .+ uya[:, j0] .^ 2))
end

"Periodic linear interpolation of (x,f) to xi."
function pinterp(x, f, xi)
    xm = mod.(x .+ L / 2, L) .- L / 2
    p  = sortperm(xm)
    xm = xm[p]; fm = f[p]
    keep = [true; diff(xm) .> 1e-6]      # x = -50 and +50 are the same point
    xm = xm[keep]; fm = fm[keep]
    xe = [xm[end] - L; xm; xm[1] + L]
    fe = [fm[end]; fm; fm[1]]
    map(xi) do xv
        xv = mod(xv + L / 2, L) - L / 2
        k  = clamp(searchsortedlast(xe, xv), 1, length(xe) - 1)
        w  = (xv - xe[k]) / (xe[k+1] - xe[k])
        (1 - w) * fe[k] + w * fe[k+1]
    end
end

"Extremum and its position, refined by a parabola through the neighbours."
function extremum(x, f, fn)
    k = fn == :max ? argmax(f) : argmin(f)
    if 1 < k < length(f)
        f0, f1, f2 = f[k-1], f[k], f[k+1]
        den = f0 - 2f1 + f2
        if den != 0
            d = 0.5 * (f0 - f2) / den
            return (f1 - 0.25 * (f0 - f2) * d, x[k] + d * (x[k+1] - x[k]))
        end
    end
    return (f[k], x[k])
end

function main(args)
    length(args) >= 2 || error("usage: julia tests/ismiphom_f.jl <out.png> <exp>:<label>=<run> ...")
    fout = args[1]
    runs = map(args[2:end]) do a
        m = match(r"^F([12]):([^=]+)=(.+)$", a)
        m === nothing && error("argument not of the form F1:label=run: $a")
        (exp=parse(Int, m[1]), label=m[2], run=m[3])
    end

    xg   = collect(-50.0:0.5:49.5)
    rms(a) = sqrt(mean(a .^ 2))
    cols = Dict(zip(MODELS, Makie.wong_colors()[1:7]))
    styles = [:solid, :dash, :dot, :dashdot]

    fig = Figure(size=(1000, 760), fontsize=13)
    for exp in (1, 2)
        ax1 = Axis(fig[exp, 1], xlabel="x (km)", ylabel="z_s perturbation (m)",
                   title="F$exp ($(exp == 1 ? "no slip" : "slip ratio 1")): surface elevation")
        ax2 = Axis(fig[exp, 2], xlabel="x (km)", ylabel="|u_s| (m/a)",
                   title="F$exp: surface speed")

        refs = Dict(m => ref_profile(m, exp) for m in MODELS)
        for m in MODELS
            lw = m in FS ? 2.0 : 1.2
            lines!(ax1, refs[m].x, refs[m].zs, color=cols[m], linewidth=lw, label=m * (m in FS ? " (FS)" : ""))
            lines!(ax2, refs[m].x, refs[m].us, color=cols[m], linewidth=lw)
        end

        ens_zs = mean(hcat([pinterp(refs[m].x, refs[m].zs, xg) for m in MODELS]...), dims=2)[:, 1]
        ens_us = mean(hcat([pinterp(refs[m].x, refs[m].us, xg) for m in MODELS]...), dims=2)[:, 1]
        fs_zs  = pinterp(refs["cma1"].x, refs["cma1"].zs, xg)
        fs_us  = pinterp(refs["cma1"].x, refs["cma1"].us, xg)

        println("\nF$exp, central flowline y = 0")
        @printf("%-16s %16s %16s %16s %16s %9s %9s %9s %9s\n", "", "z_s max @ x", "z_s min @ x",
                "|u_s| max @ x", "|u_s| min @ x", "dzs ens", "dus ens", "dzs cma1", "dus cma1")
        function row(name, x, zs, us)
            zi = pinterp(x, zs, xg); ui = pinterp(x, us, xg)
            zx, xzx = extremum(x, zs, :max); zn, xzn = extremum(x, zs, :min)
            ux, xux = extremum(x, us, :max); un, xun = extremum(x, us, :min)
            @printf("%-16s %8.2f @ %5.1f %8.2f @ %5.1f %8.2f @ %5.1f %8.2f @ %5.1f %9.2f %9.2f %9.2f %9.2f\n",
                    name, zx, xzx, zn, xzn, ux, xux, un, xun,
                    rms(zi .- ens_zs), rms(ui .- ens_us), rms(zi .- fs_zs), rms(ui .- fs_us))
        end
        for m in MODELS
            row(m, refs[m].x, refs[m].zs, refs[m].us)
        end

        k = 0
        for r in runs
            r.exp == exp || continue
            k += 1
            P = yelmo_profile(r.run)
            row("Yelmo " * r.label, P.x, P.zs, P.us)
            ls = styles[mod1(k, length(styles))]
            lines!(ax1, P.x, P.zs, color=:black, linewidth=2.5, linestyle=ls, label="Yelmo " * r.label)
            scatter!(ax1, P.x, P.zs, color=:black, markersize=4)
            lines!(ax2, P.x, P.us, color=:black, linewidth=2.5, linestyle=ls)
            scatter!(ax2, P.x, P.us, color=:black, markersize=4)
        end
        xlims!(ax1, -50, 50); xlims!(ax2, -50, 50)
        Legend(fig[exp, 3], ax1, framevisible=false)
    end

    save(fout, fig, px_per_unit=2)
    println("\nsaved ", fout)
end

main(ARGS)
