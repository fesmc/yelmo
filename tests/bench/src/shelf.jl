# ----------------------------------------------------------------------
# ShelfRadialBenchmark (A2) — circular floating shelf of constant thickness.
#
# Reference: docs/dev/benchmark-protocol/index.md, Sect. "A2 Radial floating
# shelf". Fortran driver: yelmo/tests/yelmo_bench.f90 (diagnostic run).
#
# For constant thickness H and rate factor A, the exact solution of the SSA
# with a calving-front stress condition is isotropic spreading,
#   u = ε̇ x,  v = ε̇ y,  ε̇ = A S^n 3^(−(n+1)/2),  S = ½ ρ_i g H (1 − ρ_i/ρ_w).
# The membrane stresses are uniform, so the interior momentum balance holds
# trivially and the front condition sets ε̇. A small uniform friction β_reg
# removes the rigid-body nullspace of a shelf that is floating everywhere; it
# changes the solution by a relative amount of order β_reg R_s² / (η H).
# ----------------------------------------------------------------------

export ShelfRadialBenchmark, shelf_strain_rate

"""
    shelf_strain_rate(H; A, n = 3, rho_ice = 910.0, rho_sw = 1028.0, g = 9.81)

Spreading rate ε̇ [1/yr] of a free-floating shelf of constant thickness H [m]
spreading isotropically in two dimensions, for the rate factor A [Pa^-n yr^-1].
"""
function shelf_strain_rate(H; A, n = 3, rho_ice = 910.0, rho_sw = 1028.0, g = 9.81)
    S = 0.5 * rho_ice * g * H * (1 - rho_ice / rho_sw)
    return A * S^n * 3.0^(-(n + 1) / 2)
end

"""
    ShelfRadialBenchmark(; dx_km = 10.0, H = 400.0, R_s = 300e3, L = 400e3,
                           A = 1e-18, n = 3.0, beta_reg = 1e-3)

Circular floating shelf (A2) of thickness `H` [m] and radius `R_s` [m] on a
square domain [−L, L]², with uniform friction `beta_reg` [Pa yr m^-1] and
rate factor `A` [Pa^-n yr^-1] (Glen exponent `n`).
"""
struct ShelfRadialBenchmark <: AbstractBenchmark
    xc       ::Vector{Float64}
    yc       ::Vector{Float64}
    dx_km    ::Float64
    H        ::Float64
    R_s      ::Float64
    A        ::Float64
    n        ::Float64
    beta_reg ::Float64
    z_bed    ::Float64
end

function ShelfRadialBenchmark(; dx_km::Real = 10.0, H::Real = 400.0, R_s::Real = 300e3,
                                L::Real = 400e3, A::Real = 1e-18, n::Real = 3.0,
                                beta_reg::Real = 1e-3, z_bed::Real = -2000.0)
    dx_m = Float64(dx_km) * 1e3
    N = Int(round(2L / dx_m))
    N * dx_m ≈ 2L || error("ShelfRadialBenchmark: dx_km must divide 2L (got dx_km = $dx_km, L = $L).")
    xc = collect(range(-L + dx_m/2, L - dx_m/2; length = N))
    return ShelfRadialBenchmark(xc, copy(xc), Float64(dx_km), Float64(H), Float64(R_s),
                                Float64(A), Float64(n), Float64(beta_reg), Float64(z_bed))
end

"""
    state(b::ShelfRadialBenchmark, t) -> NamedTuple

Shelf geometry and boundary fields at `t = 0`, including the friction field
`beta`. The shelf is the set of cells with centre at r < R_s.
"""
function state(b::ShelfRadialBenchmark, t::Real)
    Float64(t) == 0.0 || error("ShelfRadialBenchmark.state: only t = 0 is supported (got t = $t).")
    Nx, Ny = length(b.xc), length(b.yc)
    H_ice = [hypot(b.xc[i], b.yc[j]) < b.R_s ? b.H : 0.0 for i in 1:Nx, j in 1:Ny]
    return (xc = b.xc, yc = b.yc,
            H_ice = H_ice, z_bed = fill(b.z_bed, Nx, Ny), z_sl = zeros(Nx, Ny),
            smb_ref = zeros(Nx, Ny), T_srf = fill(263.15, Nx, Ny), Q_geo = fill(50.0, Nx, Ny),
            bmb_shlf = zeros(Nx, Ny), T_shlf = fill(271.15, Nx, Ny), H_sed = zeros(Nx, Ny),
            mask_ice = fill(MASK_ICE_DYNAMIC, Nx, Ny), beta = fill(b.beta_reg, Nx, Ny))
end

"""
    analytical_velocity(b::ShelfRadialBenchmark, t) -> (ux_bar, uy_bar)

Exact depth-averaged velocity on the cell faces: `ux_bar` of shape (Nx+1, Ny)
at x-faces, `uy_bar` of shape (Nx, Ny+1) at y-faces (IceSheetBenchmarks
convention). Defined everywhere; the comparison selects the ice-covered faces.
"""
function analytical_velocity(b::ShelfRadialBenchmark, t::Real)
    ε = shelf_strain_rate(b.H; A = b.A, n = b.n)
    dx = b.dx_km * 1e3
    xf = vcat(b.xc .- dx/2, b.xc[end] + dx/2)
    yf = vcat(b.yc .- dx/2, b.yc[end] + dx/2)
    ux = [ε * xf[i] for i in eachindex(xf), j in eachindex(b.yc)]
    uy = [ε * yf[j] for i in eachindex(b.xc), j in eachindex(yf)]
    return ux, uy
end

function write_fixture!(b::ShelfRadialBenchmark, path::AbstractString;
                        times::AbstractVector{<:Real} = [0.0])
    (length(times) == 1 && first(times) == 0.0) ||
        error("write_fixture!(ShelfRadialBenchmark, …): only times = [0.0] is supported.")
    fields = (FIELDS_2D..., ("beta", "Pa yr m^-1", "Basal friction coefficient (imposed)"))
    attrs = Dict("benchmark" => "SHELF-R", "dx_km" => b.dx_km, "H" => b.H, "R_s" => b.R_s,
                 "A" => b.A, "n" => b.n, "beta_reg" => b.beta_reg,
                 "strain_rate" => shelf_strain_rate(b.H; A = b.A, n = b.n))
    write_fixture_nc(path, state(b, 0.0), fields; attrs)
    return [path]
end
