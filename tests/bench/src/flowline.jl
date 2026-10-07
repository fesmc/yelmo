# ----------------------------------------------------------------------
# FlowlineBenchmark (B2) — marine ice sheet on a linear prograde bed
# (MISMIP experiment 1, Pattyn et al., 2012) in a strip that is periodic
# across the flow, with the steady grounding line of the boundary-layer
# theory of Schoof (2007) as reference.
#
# Reference: docs/dev/benchmark-protocol/index.md, Sect. "B2 Marine
# flowline". Fortran driver: yelmo/tests/yelmo_bench.f90 with
# yelmo.experiment = "periodic-y" (or "periodic-x" for the transposed
# strip), constant rate factor and the power-law friction
# tau_b = C |u|^(m-1) u (ydyn.beta_method = 2, beta_q = m, beta_u0 = 1 m/a,
# c_bed = cf_ref = C sec_year^(-m), N_eff = 1 Pa).
#
# Domain:
#   x ∈ [−L, L] with the ice divide at x = 0 on a cell centre (nx odd), so
#   the strip has no divide or lateral boundary condition; ny identical rows
#   across the flow (periodic). Ice is allowed for |x| < x_cf (fixed
#   calving front at a cell edge, mask_ice), the borders are ice-free.
#   Bed:  b(x) = 720 − 778.5 |x| / 750 km   [m]
#   SMB:  a = 0.3 m/a everywhere, no basal melt.
#
# Initial state (semi-analytic steady profile):
#   grounding line x_g from q_g(x_g) = a x_g with the boundary-layer flux
#   q_g(h) = [A (ρ_i g)^(n+1) (1 − ρ_i/ρ_w)^n / (4^n C)]^(1/(m+1)) h^((m+n+3)/(m+1))
#   (Schoof, 2007, Eq. 29), h_g = −(ρ_w/ρ_i) b(x_g);
#   grounded ice: the outer (sliding) solution of Schoof (2007),
#   C (a x/h)^m = −ρ_i g h d(h+b)/dx, integrated upstream from h(x_g) = h_g;
#   floating ice: the unconfined shelf, d(uh)/dx = a and
#   du/dx = A [ρ_i g (1 − ρ_i/ρ_w) h/4]^n, integrated downstream from x_g.
# ----------------------------------------------------------------------

export FlowlineBenchmark, flowline_bed, flowline_qg, flowline_xg, flowline_thickness
export MISMIP_A, MISMIP_SEC_YEAR

"Seconds per year of MISMIP (Pattyn et al., 2012)."
const MISMIP_SEC_YEAR = 31556926.0

"Rate factors A [Pa⁻³ s⁻¹] of MISMIP experiment 1 (steps 1–9)."
const MISMIP_A = [4.6416e-24, 2.1544e-24, 1.0e-24, 4.6416e-25, 2.1544e-25, 1.0e-25,
                  4.6416e-26, 2.1544e-26, 1.0e-26]

"""
    FlowlineBenchmark(; A = 4.6416e-24, dx_km = 8.0, L = 1800e3, x_cf = 1700e3, ny = 3,
                      transpose = false, C = 7.624e6, m = 1/3, n = 3.0, smb = 0.3,
                      rho_ice = 900.0, rho_sw = 1000.0, g = 9.8, sec_year = MISMIP_SEC_YEAR)

Marine flowline (B2). `A` [Pa⁻³ s⁻¹] and `C` [Pa m^(−m) s^m] are in SI units as
in MISMIP; the parameter file needs ymat.rf_const = A sec_year [Pa⁻³ a⁻¹] and
ytill.cf_ref = C sec_year^(−m) [Pa (m/a)^(−m)] (see `flowline_yelmo_params`).
`L` is the half-width of the domain and `x_cf` the position of the fixed
calving front [m], rounded to the nearest cell edge. `ny` is the number of
identical rows across the flow. `transpose = true` swaps x and y (flow along
y, periodic in x), for the test of periodic-y against periodic-x. The
densities and g are those of the MISMIP3D group of yelmo_phys_const.nml
(set yelmo.phys_const = "MISMIP3D").
"""
struct FlowlineBenchmark <: AbstractBenchmark
    exp       ::Symbol
    xc        ::Vector{Float64}
    yc        ::Vector{Float64}
    dx_km     ::Float64
    A         ::Float64
    C         ::Float64
    m         ::Float64
    n         ::Float64
    smb       ::Float64
    rho_ice   ::Float64
    rho_sw    ::Float64
    g         ::Float64
    sec_year  ::Float64
    L         ::Float64
    x_cf      ::Float64
    ny        ::Int
    transpose ::Bool
end

function FlowlineBenchmark(exp::Symbol = :ctrl; A::Real = 4.6416e-24, dx_km::Real = 8.0,
                           L::Real = 1800e3, x_cf::Real = 1700e3, ny::Real = 3,
                           transpose::Bool = false, C::Real = 7.624e6, m::Real = 1/3,
                           n::Real = 3.0, smb::Real = 0.3, rho_ice::Real = 900.0,
                           rho_sw::Real = 1000.0, g::Real = 9.8, sec_year::Real = MISMIP_SEC_YEAR)
    exp == :ctrl || error("FlowlineBenchmark: unsupported exp = $exp. Supported: :ctrl.")
    (isinteger(ny) && ny >= 1) || error("FlowlineBenchmark: ny must be a positive integer (got $ny).")
    dx  = Float64(dx_km) * 1e3
    nh  = Int(round(L / dx))
    xc  = collect((-nh:nh) .* dx)
    ncf = Int(round(x_cf / dx - 0.5))
    x_cf_edge = (ncf + 0.5) * dx
    x_cf_edge < nh * dx || error("FlowlineBenchmark: x_cf must lie inside the domain (x_cf = $x_cf, L = $L).")
    yc  = collect(((1:Int(ny)) .- (Int(ny) + 1) / 2) .* dx)
    b = FlowlineBenchmark(exp, xc, yc, Float64(dx_km), Float64(A), Float64(C), Float64(m), Float64(n),
                          Float64(smb), Float64(rho_ice), Float64(rho_sw), Float64(g), Float64(sec_year),
                          nh * dx, x_cf_edge, Int(ny), transpose)
    xg = flowline_xg(b)
    xg < b.x_cf - 2dx ||
        error("FlowlineBenchmark: the steady grounding line (x_g = $(round(xg/1e3; digits = 1)) km) " *
              "lies beyond the calving front (x_cf = $(b.x_cf/1e3) km).")
    return b
end

"Parameters of the Yelmo parameter file that depend on the benchmark (Yelmo units)."
flowline_yelmo_params(b::FlowlineBenchmark) =
    (rf_const = b.A * b.sec_year, cf_ref = b.C * b.sec_year^(-b.m), beta_q = b.m, beta_u0 = 1.0)

# -----------------------------------------------------------------------
# Geometry and reference
# -----------------------------------------------------------------------

"Bed elevation [m] of MISMIP experiment 1 at distance x [m] from the divide."
flowline_bed(x) = 720.0 - 778.5 * abs(x) / 750e3

"Flotation thickness [m] at x [m]."
flowline_hf(b::FlowlineBenchmark, x) = max(-(b.rho_sw / b.rho_ice) * flowline_bed(x), 0.0)

"""
    flowline_qg(b, h) -> Float64

Ice flux [m² s⁻¹] across a grounding line with thickness h [m], from the
boundary-layer theory of Schoof (2007, Eq. 29) for the power-law friction
τ_b = C |u|^(m−1) u.
"""
function flowline_qg(b::FlowlineBenchmark, h)
    (; A, C, m, n, rho_ice, rho_sw, g) = b
    return (A * (rho_ice * g)^(n + 1) * (1 - rho_ice / rho_sw)^n / (4^n * C))^(1 / (m + 1)) *
           h^((m + n + 3) / (m + 1))
end

"""
    flowline_xg(b) -> Float64

Steady grounding-line position [m] of the boundary-layer theory: the root of
q_g(h_f(x)) = a x on the marine part of the bed (bisection).
"""
function flowline_xg(b::FlowlineBenchmark)
    a  = b.smb / b.sec_year
    f(x) = flowline_qg(b, flowline_hf(b, x)) - a * x
    lo = 750e3 * 720.0 / 778.5 + 1.0       # Bed at sea level
    hi = 10 * b.L
    f(lo) < 0 && f(hi) > 0 || error("flowline_xg: no grounding line on the marine bed.")
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        f(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return 0.5 * (lo + hi)
end

"RK4 integration of dh/dx = F(x, h) from (x0, h0) to x1 in steps of about ds [m]; returns (x, h)."
function _rk4(F, x0, h0, x1; ds = 100.0)
    N  = max(1, ceil(Int, abs(x1 - x0) / ds))
    dx = (x1 - x0) / N
    xs = collect(range(x0, x1; length = N + 1))
    hs = zeros(N + 1); hs[1] = h0
    for k in 1:N
        x, h = xs[k], hs[k]
        k1 = F(x, h)
        k2 = F(x + dx / 2, h + dx / 2 * k1)
        k3 = F(x + dx / 2, h + dx / 2 * k2)
        k4 = F(x + dx, h + dx * k3)
        hs[k+1] = h + dx / 6 * (k1 + 2k2 + 2k3 + k4)
    end
    return xs, hs
end

"""
    flowline_thickness(b, x; ds = 100.0) -> Vector

Semi-analytic steady thickness [m] at the distances `x` [m] from the divide
(|x| is used): the outer sliding solution of Schoof (2007) upstream of the
boundary-layer grounding line, and the unconfined shelf downstream of it, up
to the calving front (zero beyond).
"""
function flowline_thickness(b::FlowlineBenchmark, x::AbstractVector; ds = 100.0)
    (; A, C, m, n, rho_ice, rho_sw, g) = b
    a    = b.smb / b.sec_year
    xg   = flowline_xg(b)
    hg   = flowline_hf(b, xg)
    dbdx = -778.5 / 750e3
    # Grounded: ρ_i g h (dh/dx + db/dx) = −C (a x/h)^m
    Fg(x, h) = -dbdx - C * (a * x / h)^m / (rho_ice * g * h)
    xgr, hgr = _rk4(Fg, xg, hg, 0.0; ds)
    reverse!(xgr); reverse!(hgr)
    # Floating: q = a x, u = q/h, du/dx = A (Cs h)^n
    Cs = rho_ice * g * (1 - rho_ice / rho_sw) / 4
    Ff(x, h) = h * (a - A * (Cs * h)^n * h) / (a * x)
    xfl, hfl = _rk4(Ff, xg, hg, b.x_cf; ds)
    H = zeros(length(x))
    for (k, xk) in enumerate(abs.(x))
        if xk <= xg
            H[k] = _interp_linear(xgr, hgr, xk)
        elseif xk < b.x_cf
            H[k] = _interp_linear(xfl, hfl, xk)
        end
    end
    return H
end

# -----------------------------------------------------------------------
# State and fixture
# -----------------------------------------------------------------------

"""
    state(b::FlowlineBenchmark, t) -> NamedTuple

Semi-analytic steady state and fixed forcing at `t = 0`. With `transpose`,
x and y are swapped (flow along y).
"""
function state(b::FlowlineBenchmark, t::Real)
    Float64(t) == 0.0 || error("FlowlineBenchmark.state: only t = 0 is supported (got t = $t).")
    Nx, Ny = length(b.xc), length(b.yc)
    H1 = flowline_thickness(b, b.xc)
    z1 = flowline_bed.(b.xc)
    m1 = [abs(x) < b.x_cf ? MASK_ICE_DYNAMIC : MASK_ICE_NONE for x in b.xc]
    rows(v) = repeat(v, 1, Ny)
    s = (H_ice = rows(H1), z_bed = rows(z1), z_sl = zeros(Nx, Ny),
         smb_ref = fill(b.smb, Nx, Ny), T_srf = fill(253.15, Nx, Ny), Q_geo = fill(50.0, Nx, Ny),
         bmb_shlf = zeros(Nx, Ny), T_shlf = fill(271.15, Nx, Ny), H_sed = zeros(Nx, Ny),
         mask_ice = rows(m1))
    if b.transpose
        return (; xc = b.yc, yc = b.xc, map(permutedims, s)...)
    else
        return (; xc = b.xc, yc = b.yc, s...)
    end
end

function write_fixture!(b::FlowlineBenchmark, path::AbstractString;
                        times::AbstractVector{<:Real} = [0.0])
    (length(times) == 1 && first(times) == 0.0) ||
        error("write_fixture!(FlowlineBenchmark, …): only times = [0.0] is supported.")
    p = flowline_yelmo_params(b)
    attrs = Dict("benchmark" => b.transpose ? "FLOWLINE-T" : "FLOWLINE", "exp" => String(b.exp),
                 "dx_km" => b.dx_km, "A" => b.A, "C" => b.C, "m" => b.m, "n" => b.n, "smb" => b.smb,
                 "rho_ice" => b.rho_ice, "rho_sw" => b.rho_sw, "g" => b.g, "sec_year" => b.sec_year,
                 "L" => b.L, "x_cf" => b.x_cf, "ny" => b.ny, "transpose" => Int(b.transpose),
                 "x_g" => flowline_xg(b),
                 "rf_const" => p.rf_const, "cf_ref" => p.cf_ref)
    write_fixture_nc(path, state(b, 0.0), FIELDS_2D; attrs)
    return [path]
end
