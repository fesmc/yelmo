# ----------------------------------------------------------------------
# Island4Benchmark — closed island with land-terminating, marine-terminating
# and floating margins, invariant under the D4 group of the grid.
#
# Reference: docs/dev/benchmark-protocol/index.md, Sect. "ISLAND4 domain".
# Fortran driver: yelmo/tests/yelmo_bench.f90 (online forcing in
# yelmo/tests/bench_forcing.f90, mirrored here by island4_tsrf and
# island4_bmb).
#
# Domain:
#   x, y ∈ [-800, 800] km, cell-centred axes (origin on a cell corner).
#   Bed:  radial base profile minus four V-shaped troughs along the
#         diagonals (rot = 0) or the axes (rot = 45, ISLAND4-R).
#
# Forcing:
#   smb_ref   radial, SMB0 (1 − r/r_ela)                      (fixed)
#   Q_geo     uniform                                          (fixed)
#   mask_ice  no ice for r ≥ r_lim                             (fixed)
#   T_srf     lapse rate on the evolving surface               (online)
#   bmb_shlf  MISMIP+ Ice1 (Asay-Davis et al., 2016)           (online)
# ----------------------------------------------------------------------

export Island4Benchmark
export island4_bed, island4_smb, island4_tsrf, island4_bmb, island4_vialov

# Values for mask_ice (as in yelmo/src/yelmo_defs.f90)
const MASK_ICE_NONE    = 0
const MASK_ICE_DYNAMIC = 2

# -----------------------------------------------------------------------
# Geometry
# -----------------------------------------------------------------------

"Radial base profile of the island [m]."
island4_base(r; Bc = 900.0, Bl = -2000.0, R0 = 1000e3) = Bc - (Bc - Bl) * r^2 / R0^2

"""
    island4_trough(ξ, η; r_h, alpha, ell, D0, B_od, r_od, w_od) -> Float64

Depth [m] of one V-shaped trough in local coordinates (ξ along, η across the
trough axis). The walls are straight lines from the apex at ξ = r_h with
opening angle `alpha` [deg], smoothed over `ell`. `B_od` > 0 adds an
overdeepening centred at ξ = r_od.
"""
function island4_trough(ξ, η; r_h = 250e3, alpha = 120.0, ell = 50e3, D0 = 1500.0,
                        B_od = 0.0, r_od = 350e3, w_od = 60e3)
    φ = deg2rad(alpha / 2)
    T = 0.5 * (1 + tanh(((ξ - r_h) * tan(φ) - sqrt(η^2 + ell^2)) / ell))
    D = D0 + B_od * exp(-((ξ - r_od) / w_od)^2)
    return D * T
end

"""
    island4_bed(x, y; B_od = 0.0, rot = 0.0, dz = 0.0, kw...) -> Float64

ISLAND4 bed elevation [m] at (x, y) [m]. Troughs lie along the diagonals for
rot = 0 and along the axes for rot = 45 (ISLAND4-R). `dz` raises the whole bed
(ISLAND4-L uses dz = 2500 m). Further keywords go to `island4_trough`.
"""
function island4_bed(x, y; B_od = 0.0, rot = 0.0, dz = 0.0, kw...)
    z = island4_base(hypot(x, y)) + dz
    for ψ in deg2rad.((45.0, 135.0, 225.0, 315.0) .+ rot)
        ξ =  x * cos(ψ) + y * sin(ψ)
        η = -x * sin(ψ) + y * cos(ψ)
        z -= island4_trough(ξ, η; B_od, kw...)
    end
    return z
end

"Radial surface mass balance [m/a ice eq.]."
island4_smb(r; smb0 = 0.5, r_ela = 450e3) = smb0 * (1 - r / r_ela)

"""
    island4_vialov(r; H0 = 3500.0, R_i = 650e3, n = 3) -> Float64

Vialov profile [m], used as initial thickness (A3, B1, C0).
"""
function island4_vialov(r; H0 = 3500.0, R_i = 650e3, n = 3)
    r >= R_i && return 0.0
    return H0 * (1 - (r / R_i)^((n + 1) / n))^(n / (2n + 2))
end

# -----------------------------------------------------------------------
# Online forcing (Julia counterparts of tests/bench_forcing.f90)
# -----------------------------------------------------------------------

"""
    island4_tsrf(z_srf; T_sl = -10.0, lapse = 8e-3, T0 = 273.15) -> Float64

Surface temperature [K] from a lapse rate on the surface elevation [m].
Mirrors `bench_tsrf_lapse` in tests/bench_forcing.f90.
"""
island4_tsrf(z_srf; T_sl = -10.0, lapse = 8e-3, T0 = 273.15) = T0 + T_sl - lapse * z_srf

"""
    island4_bmb(H_ice, z_bed, z_sl; rho_ice, rho_sw, Omega = 0.2, Hc0 = 75.0, z0 = -100.0)

Basal mass balance below floating ice [m/a ice eq., negative for melt], from the
MISMIP+ Ice1 parameterization m = Ω tanh(H_c/H_c0) max(z0 − z_d, 0), where z_d
is the depth of the ice base and H_c the water-column thickness. Zero for
grounded ice and ice-free points. Mirrors `bench_bmb_mismipplus` in
tests/bench_forcing.f90.
"""
function island4_bmb(H_ice, z_bed, z_sl; rho_ice = 910.0, rho_sw = 1028.0,
                     Omega = 0.2, Hc0 = 75.0, z0 = -100.0)
    z_d = max(z_bed, z_sl - H_ice * rho_ice / rho_sw)
    H_c = z_d - z_bed
    return -Omega * tanh(H_c / Hc0) * max(z0 - z_d, 0.0)
end

# -----------------------------------------------------------------------
# Island4Benchmark struct
# -----------------------------------------------------------------------

"""
    Island4Benchmark(exp::Symbol = :ctrl; dx_km = 16.0, B_od = 0.0, rot = 0.0,
                     land = false, init = :vialov, ...)

ISLAND4 benchmark (docs/dev/benchmark-protocol). `exp` selects the fixed
forcing: `:ctrl` (control) or `:smb` (SMB − 0.1 m/a, C2c). The ocean
experiments (C2a, C2b) change the online forcing in the parameter file and
use the `:ctrl` fixture.

Keywords:
  - `dx_km`  grid resolution [km]; the axes are cell-centred on [-800, 800] km.
  - `B_od`   depth of the trough overdeepening [m] (0: off).
  - `rot`    rotation of the troughs [deg]: 0 (diagonals) or 45 (axes, ISLAND4-R).
  - `land`   true: ISLAND4-L, bed raised by 2500 m (B1).
  - `init`   initial thickness: `:vialov` or `:zero`.
"""
struct Island4Benchmark <: AbstractBenchmark
    exp      ::Symbol
    xc       ::Vector{Float64}
    yc       ::Vector{Float64}
    dx_km    ::Float64
    B_od     ::Float64
    rot      ::Float64
    land     ::Bool
    init     ::Symbol
    smb0     ::Float64
    r_ela    ::Float64
    dsmb     ::Float64
    r_lim    ::Float64
    Q_geo    ::Float64
    T_shlf   ::Float64
end

function Island4Benchmark(exp::Symbol = :ctrl;
                          dx_km::Real  = 16.0,
                          B_od::Real   = 0.0,
                          rot::Real    = 0.0,
                          land::Bool   = false,
                          init::Symbol = :vialov,
                          smb0::Real   = 0.5,
                          r_ela::Real  = 450e3,
                          r_lim::Real  = 750e3,
                          Q_geo::Real  = 50.0,
                          T_shlf::Real = 271.15)
    exp in (:ctrl, :smb) ||
        error("Island4Benchmark: unsupported exp = $exp. Supported: :ctrl, :smb.")
    init in (:vialov, :zero) ||
        error("Island4Benchmark: unsupported init = $init. Supported: :vialov, :zero.")
    rot in (0.0, 45.0) ||
        error("Island4Benchmark: rot must be 0 or 45 to keep the D4 symmetry (got $rot).")

    dx_m     = Float64(dx_km) * 1e3
    extent_m = 1600e3
    N  = Int(round(extent_m / dx_m))
    N * dx_m ≈ extent_m || error("Island4Benchmark: dx_km must divide 1600 km (got $dx_km).")
    xc = collect(range(-extent_m/2 + dx_m/2, extent_m/2 - dx_m/2; length = N))

    dsmb = exp == :smb ? -0.1 : 0.0

    return Island4Benchmark(exp, xc, copy(xc), Float64(dx_km), Float64(B_od), Float64(rot),
                            land, init, Float64(smb0), Float64(r_ela), dsmb,
                            Float64(r_lim), Float64(Q_geo), Float64(T_shlf))
end

"Short name of the variant, e.g. \"ISLAND4-L-R-OD\"."
function island4_name(b::Island4Benchmark)
    name = "ISLAND4"
    b.land      && (name *= "-L")
    b.rot != 0  && (name *= "-R")
    b.B_od != 0 && (name *= "-OD")
    return name
end

# -----------------------------------------------------------------------
# Analytical t = 0 state
# -----------------------------------------------------------------------

"""
    state(b::Island4Benchmark, t) -> NamedTuple

Initial state and fixed forcing at `t = 0`. T_srf and bmb_shlf are the online
forcing evaluated on the initial geometry (z_sl = 0); the driver recomputes them
during the run.
"""
function state(b::Island4Benchmark, t::Real)
    Float64(t) == 0.0 || error(
        "Island4Benchmark.state: only t = 0 is supported (got t = $t).")

    Nx, Ny = length(b.xc), length(b.yc)
    # Evaluate every field at the image of (x, y) in the octant 0 ≤ y ≤ x, so
    # that the discrete fields are exactly invariant under D4.
    oct(i, j) = (max(abs(b.xc[i]), abs(b.yc[j])), min(abs(b.xc[i]), abs(b.yc[j])))
    r(i, j)   = hypot(oct(i, j)...)
    dz = b.land ? 2500.0 : 0.0

    z_bed = [island4_bed(oct(i, j)...; B_od = b.B_od, rot = b.rot, dz) for i in 1:Nx, j in 1:Ny]
    smb   = [island4_smb(r(i, j); smb0 = b.smb0, r_ela = b.r_ela) + b.dsmb for i in 1:Nx, j in 1:Ny]
    mask_ice = [r(i, j) < b.r_lim ? MASK_ICE_DYNAMIC : MASK_ICE_NONE for i in 1:Nx, j in 1:Ny]

    H_ice = b.init == :vialov ? [island4_vialov(r(i, j)) for i in 1:Nx, j in 1:Ny] : zeros(Nx, Ny)
    H_ice[mask_ice .== MASK_ICE_NONE] .= 0.0

    z_sl  = zeros(Nx, Ny)
    z_srf = max.(z_bed .+ H_ice, z_sl .+ H_ice .* (1 - 910.0 / 1028.0))
    T_srf    = island4_tsrf.(z_srf)
    bmb_shlf = island4_bmb.(H_ice, z_bed, z_sl)

    # Initial ice temperature (stationary column solution)
    zeta  = collect(range(0.0, 1.0; length = 51))
    flt   = H_ice .* 910.0 ./ 1028.0 .< z_sl .- z_bed
    T_ice = zeros(Nx, Ny, length(zeta))
    for j in 1:Ny, i in 1:Nx
        T_ice[i, j, :] = column_temperature(zeta, H_ice[i, j], T_srf[i, j], smb[i, j], b.Q_geo;
                                            floating = flt[i, j], T_shlf = b.T_shlf)
    end

    return (xc = b.xc, yc = b.yc, zeta = zeta,
            H_ice = H_ice, z_bed = z_bed, z_sl = z_sl,
            smb_ref = smb, T_srf = T_srf, Q_geo = fill(b.Q_geo, Nx, Ny),
            bmb_shlf = bmb_shlf, T_shlf = fill(b.T_shlf, Nx, Ny),
            H_sed = zeros(Nx, Ny), mask_ice = mask_ice, T_ice = T_ice)
end

"""
    write_fixture!(b::Island4Benchmark, path; times = [0.0]) -> Vector{String}

Write `state(b, 0)` to the NetCDF fixture read by tests/yelmo_bench.f90.
"""
function write_fixture!(b::Island4Benchmark, path::AbstractString;
                        times::AbstractVector{<:Real} = [0.0])
    (length(times) == 1 && first(times) == 0.0) ||
        error("write_fixture!(Island4Benchmark, …): only times = [0.0] is supported.")

    s = state(b, 0.0)
    attrs = Dict("benchmark" => island4_name(b), "exp" => String(b.exp),
                 "dx_km" => b.dx_km, "B_od" => b.B_od, "rot" => b.rot,
                 "init" => String(b.init))
    write_fixture_nc(path, s, FIELDS_2D; fields3D = FIELDS_3D, attrs)
    return [path]
end
