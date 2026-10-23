# ----------------------------------------------------------------------
# Comparison of Yelmo output with analytical velocity solutions.
#
# Yelmo writes ux_bar(i, j) on the right face of cell i and uy_bar(i, j) on
# the top face of cell j, both with shape (Nx, Ny). analytical_velocity
# returns all faces, (Nx+1, Ny) and (Nx, Ny+1), so Yelmo face i corresponds
# to analytical face i+1.
# ----------------------------------------------------------------------

export velocity_errors

"""
    velocity_errors(b, file; k = nothing, io = stdout) -> NamedTuple

Error of the depth-averaged velocity in a Yelmo output file against
`analytical_velocity(b, 0)`, at time index `k` (default: last). Faces between
two ice-covered cells are interior faces; faces between an ice-covered and an
ice-free cell are front faces. Errors are relative to the maximum analytical
speed on the interior faces.
"""
function velocity_errors(b::AbstractBenchmark, file::AbstractString; k = nothing, io = stdout)
    uxa, uya = analytical_velocity(b, 0.0)
    ux, uy, H = NCDataset(file) do ds
        read2d(ds, "ux_bar"; k), read2d(ds, "uy_bar"; k), read2d(ds, "H_ice"; k)
    end
    Nx, Ny = size(H)
    ice = H .> 0

    errs = Dict(:interior => Float64[], :front => Float64[])
    ref  = Float64[]
    for j in 1:Ny, i in 1:Nx-1
        n_ice = ice[i, j] + ice[i+1, j]
        n_ice == 0 && continue
        push!(errs[n_ice == 2 ? :interior : :front], ux[i, j] - uxa[i+1, j])
        n_ice == 2 && push!(ref, abs(uxa[i+1, j]))
    end
    for j in 1:Ny-1, i in 1:Nx
        n_ice = ice[i, j] + ice[i, j+1]
        n_ice == 0 && continue
        push!(errs[n_ice == 2 ? :interior : :front], uy[i, j] - uya[i, j+1])
        n_ice == 2 && push!(ref, abs(uya[i, j+1]))
    end

    umax = maximum(ref)
    rms(v) = sqrt(sum(abs2, v) / length(v))
    out = (umax = umax,
           interior_max = maximum(abs, errs[:interior]) / umax,
           interior_rms = rms(errs[:interior]) / umax,
           front_max    = maximum(abs, errs[:front]) / umax,
           front_rms    = rms(errs[:front]) / umax)
    println(io, "max analytical speed (interior faces) = ", round(umax; sigdigits = 4), " m/yr")
    println(io, "interior faces: max rel. error = ", round(out.interior_max; sigdigits = 3),
                ", rms = ", round(out.interior_rms; sigdigits = 3))
    println(io, "front faces:    max rel. error = ", round(out.front_max; sigdigits = 3),
                ", rms = ", round(out.front_rms; sigdigits = 3))
    return out
end
