# ----------------------------------------------------------------------
# Mass budget residual from the Yelmo time series (yelmo_ts.nc).
#
# The driver writes yelmo_ts.nc after every outer step dtt. The budget terms
# (smb_tot, bmb_tot, fmb_tot, dmb, cmb) are the rates applied by Yelmo,
# averaged over that step (calc_ytopo_rates), so their integral over the
# output interval (t[n-1], t[n]] is (t[n] - t[n-1]) times the value at t[n].
# ----------------------------------------------------------------------

export mass_budget, budget_closure

"""
    mass_budget(run) -> NamedTuple

Mass budget residual per output interval of yelmo_ts.nc,

    r_M[n] = [ΔV − Δt (SMB + BMB + FMB + DMB + CMB)] / V[n],

with all terms integrated over the region of the time series (the calving
mask of the benchmark). CMB is negative for calving. Returns (time, r_M, V,
dV, flux), with dV and flux in m³ per interval; r_M[1] (initial state) is 0.
"""
function mass_budget(run::AbstractString)
    file = isdir(run) ? joinpath(run, "yelmo_ts.nc") : run
    NCDataset(file) do ds
        t = Float64.(ds["time"][:])
        V = Float64.(ds["V_ice"][:]) .* 1e15          # [1e6 km³] → [m³]
        rate = zeros(length(t))
        for name in ("smb_tot", "bmb_tot", "fmb_tot", "dmb", "cmb")   # [m³/yr]
            rate .+= Float64.(coalesce.(ds[name][:], 0.0))
        end
        dt   = [0.0; diff(t)]
        dV   = [0.0; diff(V)]
        flux = dt .* rate
        flux[1] = 0.0
        r_M  = [V[n] > 0 ? (dV[n] - flux[n]) / V[n] : 0.0 for n in eachindex(t)]
        return (time = t, r_M = r_M, V = V, dV = dV, flux = flux)
    end
end

"""
    budget_closure(run) -> NamedTuple

Closed mass budget for the output intervals of yelmo_ts.nc that end at a 2D
output time (yelmo.nc):

    r_C = [ΔV − Δt Σ (mb_net + cmb) dx²] / V,

where mb_net (2D, averaged over the interval like the time-series terms) also
holds the margin terms that the time series does not resolve (mb_resid: thin
and isolated margin ice, the mask_ice constraint; mb_relax). A nonzero r_C is
mass that leaves without any budget term, e.g. transport into cells with
mask_ice = none, which the advection solver holds at zero thickness. Returns
(time, r_C, r_resid), with r_resid = Δt Σ (mb_net − smb − bmb − fmb) dx² / V
the part of r_M explained by mb_resid and mb_relax (fmb when written; the 2D
output has no dmb).
"""
function budget_closure(run::AbstractString)
    b = mass_budget(run)
    NCDataset(joinpath(run, "yelmo.nc")) do ds
        xc = Float64.(ds["xc"][:])
        a  = ((xc[2] - xc[1]) * 1e3)^2                 # [km] → cell area [m²]
        t2 = Float64.(ds["time"][:])
        time, r_C, r_resid = Float64[], Float64[], Float64[]
        for k in eachindex(t2)
            n = findfirst(≈(t2[k]), b.time)
            (n === nothing || n == 1) && continue
            dt  = b.time[n] - b.time[n-1]
            sum2(name) = sum(Float64.(coalesce.(ds[name][:, :, k], 0.0)))
            mb  = sum2("mb_net") * a * dt
            cmb = sum2("cmb") * a * dt
            ext = haskey(ds, "fmb") ? sum2("fmb") : 0.0
            push!(time, t2[k])
            push!(r_C, (b.dV[n] - mb - cmb) / b.V[n])
            push!(r_resid, (mb - (sum2("smb") + sum2("bmb") + ext) * a * dt) / b.V[n])
        end
        return (time = time, r_C = r_C, r_resid = r_resid)
    end
end
