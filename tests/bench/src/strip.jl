# ----------------------------------------------------------------------
# StripThermoBenchmark (A4) — 3D thermodynamics on a strip with prescribed
# plug flow, where every column satisfies the assumptions of the analytic
# column solution of IceColumnSolutions.jl (Moreno-Parada et al., 2024).
#
# Reference: docs/dev/benchmark-protocol/index.md, Sect. "A4 3D
# thermodynamics strip". Fortran driver: yelmo/tests/yelmo_bench.f90 with
# ydyn.solver = "fixed" and ytopo.topo_fixed = True.
#
# Domain: x ∈ [−nx dx/2, nx dx/2] (cell-centred, x = 0 on a cell centre),
# three identical rows per parameter set along y, plus one border row at each
# end that repeats its neighbour. Only the middle row of each set and the
# columns two or more cells from the x-borders are compared: the default
# vertical velocity (ydyn.uz_method = 3) averages the divergence over the
# neighbouring cells (weights 1/4, 1/2, 1/4), which mixes parameter sets in
# adjacent rows and reaches one cell in from the borders.
#
# Velocity (plug flow, uniform with depth):
#   u = c(y) x,  v = 0,  c = SMB/H.
# The horizontal divergence is SMB/H, so continuity gives the linear vertical
# velocity uz = −SMB ζ, and the thickness is steady. The temperature is
# uniform in x and v = 0, so horizontal advection vanishes, also in the
# discrete equations. With solver = "fixed" the basal stress is zero, so
# there is no frictional heating.
#
# Experiments:
#   :stationary  T_ice and forcing from parameter set P1 (drift test, A4a)
#   :tsrf        T_ice from P1, forcing P2 = P1 with T_srf + dT   (A4b)
#   :smb         T_ice from P1, forcing P2 = P1 with SMB × fsmb  (A4b)
# ----------------------------------------------------------------------

export StripThermoBenchmark, strip_params, strip_column, strip_solution, strip_check_row

# Default rows (H [m], SMB [m/yr]): Péclet numbers ~0.5–33, all with a basal
# temperature at least ~5 K below the pressure-melting point for
# T_srf = −40 °C and Q_geo = 50 mW m⁻².
const STRIP_ROWS = [(1000.0, 0.02), (1000.0, 0.05), (1000.0, 0.1), (1000.0, 0.2), (1000.0, 0.4),
                    (2000.0, 0.05), (2000.0, 0.1),  (2000.0, 0.2), (2000.0, 0.4),
                    (3000.0, 0.1),  (3000.0, 0.2),  (3000.0, 0.4)]

"""
    StripThermoBenchmark(exp::Symbol = :stationary; dx_km = 10.0, nx = 7,
                         rows = STRIP_ROWS, T_srf = 233.15, Q_geo = 50.0,
                         dT = -10.0, fsmb = 2.0, z_bed = 500.0)

Thermodynamics strip (A4). `rows` holds one (H, SMB) pair per row. `T_srf` [K]
and `Q_geo` [mW m⁻²] are uniform. `dT` [K] and `fsmb` define the step change
of the transient experiments `:tsrf` and `:smb`.
"""
struct StripThermoBenchmark <: AbstractBenchmark
    exp    ::Symbol
    xc     ::Vector{Float64}
    yc     ::Vector{Float64}
    dx_km  ::Float64
    rows   ::Vector{Tuple{Float64,Float64}}
    T_srf  ::Float64
    Q_geo  ::Float64
    dT     ::Float64
    fsmb   ::Float64
    z_bed  ::Float64
end

function StripThermoBenchmark(exp::Symbol = :stationary;
                              dx_km::Real = 10.0, nx::Integer = 7,
                              rows = STRIP_ROWS, T_srf::Real = 233.15, Q_geo::Real = 50.0,
                              dT::Real = -10.0, fsmb::Real = 2.0, z_bed::Real = 500.0)
    exp in (:stationary, :tsrf, :smb) ||
        error("StripThermoBenchmark: unsupported exp = $exp. Supported: :stationary, :tsrf, :smb.")
    isodd(nx) || error("StripThermoBenchmark: nx must be odd, so that x = 0 is a cell centre (got $nx).")
    dx = Float64(dx_km) * 1e3
    xc = collect((-(nx ÷ 2):(nx ÷ 2)) .* dx)
    yc = collect(((0:STRIP_NREP*length(rows)+1) .+ 0.5) .* dx)
    b = StripThermoBenchmark(exp, xc, yc, Float64(dx_km), [Tuple(Float64.(r)) for r in rows],
                             Float64(T_srf), Float64(Q_geo), Float64(dT), Float64(fsmb), Float64(z_bed))
    for (k, (H, smb)) in enumerate(b.rows)
        Tb   = strip_column(b, k, [0.0])[1]
        Tpmp = 273.15 - 9.8e-8 * 910.0 * 9.81 * H
        Tb < Tpmp - 1.0 ||
            error("StripThermoBenchmark: row $k (H = $H, SMB = $smb) has a temperate base " *
                  "(T_b − T_pmp = $(round(Tb - Tpmp; digits = 2)) K); the analytic solution needs a cold base.")
    end
    return b
end

"Rows per parameter set (the middle one is compared)."
const STRIP_NREP = 3

"Parameter set of grid row j (border rows repeat their neighbour)."
strip_row(b::StripThermoBenchmark, j) = clamp(cld(j - 1, STRIP_NREP), 1, length(b.rows))

"Grid row compared for parameter set k (the middle row of its block)."
strip_check_row(b::StripThermoBenchmark, k) = STRIP_NREP * k

"""
    strip_params(b, k; forcing = :P2) -> IceColumnPar

Column parameters of row `k`, for the initial parameter set (`:P1`) or the
forcing of the experiment (`:P2`, equal to P1 for `:stationary`).
"""
function strip_params(b::StripThermoBenchmark, k; forcing::Symbol = :P2)
    H, smb = b.rows[k]
    T_srf  = b.T_srf
    if forcing == :P2
        b.exp == :tsrf && (T_srf += b.dT)
        b.exp == :smb  && (smb *= b.fsmb)
    end
    # Downward surface velocity w0 = −SMB (w0 < 0 is downward in IceColumnSolutions.jl)
    return IceColumnPar(H, T_srf, COLUMN_KAPPA, COLUMN_K, 0.0, b.Q_geo * 1e-3; w0 = -smb)
end

"""
    strip_column(b, k, zeta) -> Vector

Initial (stationary P1) temperature [K] of row `k` on the levels `zeta`.
"""
function strip_column(b::StripThermoBenchmark, k, zeta)
    return solve_stationary(strip_params(b, k; forcing = :P1); zeta = collect(Float64, zeta)).T_eq
end

"""
    strip_solution(b, k, zeta, times; n_modes = 10) -> Matrix

Analytic temperature [K] of row `k` on the levels `zeta` (rows) at the times
`times` [yr] (columns): the transient solution under the forcing P2 that starts
from the stationary P1 profile (Moreno-Parada et al., 2024, Appendix A). For
`:stationary` the P1 profile at every time. Ten modes reproduce 30 to ~1e-3 K
for t ≥ 1 kyr (each further mode adds tens of seconds per column).
"""
function strip_solution(b::StripThermoBenchmark, k, zeta, times; n_modes::Int = 10)
    z  = collect(Float64, zeta)
    b.exp == :stationary && return repeat(strip_column(b, k, z), 1, length(times))
    # Initial condition: the P1 profile, interpolated from a fine grid
    zf = collect(range(0.0, 1.0; length = 2001))
    Tf = strip_column(b, k, zf)
    θ0(ξ, p) = _interp_linear(zf, Tf, ξ) / p.T_air
    sol = solve(strip_params(b, k; forcing = :P2), collect(Float64, times); init = θ0, n_modes, zeta = z)
    return sol.T
end

function _interp_linear(x, y, xi)
    i = clamp(searchsortedlast(x, xi), 1, length(x) - 1)
    w = (xi - x[i]) / (x[i+1] - x[i])
    return (1 - w) * y[i] + w * y[i+1]
end

"""
    state(b::StripThermoBenchmark, t) -> NamedTuple

Geometry, forcing (P2), prescribed velocity and initial temperature (P1) at
`t = 0`. `ux_bar` and `uy_bar` are on the Yelmo C-grid faces: ux_bar[i, j] at
x = xc[i] + dx/2, uy_bar[i, j] at y = yc[j] + dx/2.
"""
function state(b::StripThermoBenchmark, t::Real)
    Float64(t) == 0.0 || error("StripThermoBenchmark.state: only t = 0 is supported (got t = $t).")
    Nx, Ny = length(b.xc), length(b.yc)
    dx   = b.dx_km * 1e3
    zeta = collect(range(0.0, 1.0; length = 201))

    H_ice  = zeros(Nx, Ny); smb = zeros(Nx, Ny); T_srf = zeros(Nx, Ny)
    ux_bar = zeros(Nx, Ny); T_ice = zeros(Nx, Ny, length(zeta))
    for j in 1:Ny
        k  = strip_row(b, j)
        p2 = strip_params(b, k; forcing = :P2)
        T1 = strip_column(b, k, zeta)
        for i in 1:Nx
            H_ice[i, j]  = p2.L
            smb[i, j]    = -p2.w0
            T_srf[i, j]  = p2.T_air
            ux_bar[i, j] = -p2.w0 / p2.L * (b.xc[i] + dx / 2)
            T_ice[i, j, :] = T1
        end
    end

    return (xc = b.xc, yc = b.yc, zeta = zeta,
            H_ice = H_ice, z_bed = fill(b.z_bed, Nx, Ny), z_sl = zeros(Nx, Ny),
            smb_ref = smb, T_srf = T_srf, Q_geo = fill(b.Q_geo, Nx, Ny),
            bmb_shlf = zeros(Nx, Ny), T_shlf = fill(271.15, Nx, Ny), H_sed = zeros(Nx, Ny),
            mask_ice = fill(MASK_ICE_DYNAMIC, Nx, Ny),
            ux_bar = ux_bar, uy_bar = zeros(Nx, Ny), T_ice = T_ice)
end

function write_fixture!(b::StripThermoBenchmark, path::AbstractString;
                        times::AbstractVector{<:Real} = [0.0])
    (length(times) == 1 && first(times) == 0.0) ||
        error("write_fixture!(StripThermoBenchmark, …): only times = [0.0] is supported.")
    fields = (FIELDS_2D...,
              ("ux_bar", "m/yr", "Prescribed velocity, x-component (x-faces, plug flow)"),
              ("uy_bar", "m/yr", "Prescribed velocity, y-component (y-faces, plug flow)"))
    attrs = Dict("benchmark" => "STRIP-A4", "exp" => String(b.exp), "dx_km" => b.dx_km,
                 "nx" => length(b.xc), "rows_H" => first.(b.rows), "rows_smb" => last.(b.rows),
                 "T_srf" => b.T_srf, "Q_geo" => b.Q_geo, "dT" => b.dT, "fsmb" => b.fsmb,
                 "z_bed" => b.z_bed)
    write_fixture_nc(path, state(b, 0.0), fields; fields3D = FIELDS_3D, attrs)
    return [path]
end
