# ----------------------------------------------------------------------
# FlowlineBenchmark (B2) — marine ice sheet on a linear prograde bed
# (MISMIP experiment 1, Pattyn et al., 2012) in a strip that is periodic
# across the flow, with the steady grounding line of the boundary-layer
# theory of Schoof (2007) as reference.
#
# Reference: docs/dev/benchmark-protocol/index.md, Sect. "B2 Marine
# flowline". Fortran driver: yelmo/tests/yelmo_bench.f90 with constant rate
# factor and the power-law friction tau_b = C |u|^(m-1) u
# (ydyn.beta_method = 2, beta_q = m, beta_u0 = 1 m/a,
# c_bed = cf_ref = C sec_year^(-m), N_eff = 1 Pa).
#
# Domain (keyword `domain`):
#   :onesided (default) x ∈ [0, L] with the ice divide at x = 0 on the face
#     between the first two cells (cell centres −dx/2, dx/2, 3dx/2, …), with
#     yelmo.experiment = "MISMIP3D": the first cell mirrors the second (H, β),
#     ux = 0 on the divide face (no-slip low-x border of the SSA), the strip
#     is periodic in y.
#   :symmetric x ∈ [−L, L] with the divide at x = 0 on the centre cell (nx
#     odd), so the strip has no divide boundary condition, with
#     yelmo.experiment = "periodic-y" ("periodic-x" with transpose = true).
#   ny identical rows across the flow (periodic). Ice is allowed for
#   |x| < x_cf (fixed calving front at a cell edge, mask_ice), the far
#   x-borders are ice-free.
#   Bed:  b(x) = 720 − 778.5 |x| / 750 km   [m]
#   SMB:  a = 0.3 m/a everywhere, no basal melt.
#
# Initial state (semi-analytic profile):
#   grounding line x_g from q_g(x_g) = a x_g with the boundary-layer flux
#   q_g(h) = [A (ρ_i g)^(n+1) (1 − ρ_i/ρ_w)^n / (4^n C)]^(1/(m+1)) h^((m+n+3)/(m+1))
#   (Schoof, 2007, Eq. 29), h_g = −(ρ_w/ρ_i) b(x_g);
#   grounded ice: the outer (sliding) solution of Schoof (2007),
#   C (a x/h)^m = −ρ_i g h d(h+b)/dx, integrated upstream from h(x_g) = h_g;
#   floating ice: the unconfined shelf, d(uh)/dx = a and
#   du/dx = A [ρ_i g (1 − ρ_i/ρ_w) h/4]^n, integrated downstream from x_g.
#   Perturbed starts (keyword dxg_start): the same construction with the
#   grounding line at x_g + dxg_start (advanced for dxg_start > 0,
#   retreated for dxg_start < 0).
# ----------------------------------------------------------------------

export FlowlineBenchmark, flowline_from_fixture, flowline_experiment, flowline_bed, flowline_qg,
       flowline_xg, flowline_thickness, flowline_grounding_line
export MISMIP_A, MISMIP_SEC_YEAR

"Seconds per year of MISMIP (Pattyn et al., 2012)."
const MISMIP_SEC_YEAR = 31556926.0

"Rate factors A [Pa⁻³ s⁻¹] of MISMIP experiment 1 (steps 1–9)."
const MISMIP_A = [4.6416e-24, 2.1544e-24, 1.0e-24, 4.6416e-25, 2.1544e-25, 1.0e-25,
                  4.6416e-26, 2.1544e-26, 1.0e-26]

"Distance [m] from the divide where the bed of MISMIP experiment 1 crosses sea level."
const FLOWLINE_X_SL = 750e3 * 720.0 / 778.5

"""
    FlowlineBenchmark(; A = 4.6416e-24, dx_km = 8.0, domain = :onesided, L = 1800e3,
                      x_cf = 1700e3, ny = 3, dxg_start = 0.0, transpose = false,
                      C = 7.624e6, m = 1/3, n = 3.0, smb = 0.3, rho_ice = 900.0,
                      rho_sw = 1000.0, g = 9.8, sec_year = MISMIP_SEC_YEAR)

Marine flowline (B2). `A` [Pa⁻³ s⁻¹] and `C` [Pa m^(−m) s^m] are in SI units as
in MISMIP; the parameter file needs ymat.rf_const = A sec_year [Pa⁻³ a⁻¹] and
ytill.cf_ref = C sec_year^(−m) [Pa (m/a)^(−m)] (see `flowline_yelmo_params`).
`domain` is `:onesided` (x ∈ [0, L], divide at the low-x border,
yelmo.experiment = "MISMIP3D") or `:symmetric` (x ∈ [−L, L], divide on the
centre cell, yelmo.experiment = "periodic-y"); see `flowline_experiment`. `L`
is the extent of the domain from the divide and `x_cf` the position of the
fixed calving front [m], rounded to the nearest cell edge. `ny` is the number
of identical rows across the flow. `dxg_start` [m] shifts the grounding line of
the initial profile from the reference (perturbed starts). `transpose = true`
swaps x and y (flow along y, periodic in x; symmetric domain only), for the
test of periodic-y against periodic-x. The densities and g are those of the
MISMIP3D group of yelmo_phys_const.nml (set yelmo.phys_const = "MISMIP3D").
"""
struct FlowlineBenchmark <: AbstractBenchmark
    exp       ::Symbol
    domain    ::Symbol
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
    dxg_start ::Float64
    transpose ::Bool
end

function FlowlineBenchmark(exp::Symbol = :ctrl; A::Real = 4.6416e-24, dx_km::Real = 8.0,
                           domain::Symbol = :onesided, L::Real = 1800e3, x_cf::Real = 1700e3,
                           ny::Real = 3, dxg_start::Real = 0.0, transpose::Bool = false,
                           C::Real = 7.624e6, m::Real = 1/3, n::Real = 3.0, smb::Real = 0.3,
                           rho_ice::Real = 900.0, rho_sw::Real = 1000.0, g::Real = 9.8,
                           sec_year::Real = MISMIP_SEC_YEAR)
    exp == :ctrl || error("FlowlineBenchmark: unsupported exp = $exp. Supported: :ctrl.")
    (isinteger(ny) && ny >= 1) || error("FlowlineBenchmark: ny must be a positive integer (got $ny).")
    dx = Float64(dx_km) * 1e3
    nh = Int(round(L / dx))
    if domain == :symmetric
        # Divide on the centre cell, cell edges at (k + 1/2) dx
        xc = collect((-nh:nh) .* dx)
        x_cf_edge = (round(x_cf / dx - 0.5) + 0.5) * dx
    elseif domain == :onesided
        # Divide on the face between the first two cells, cell edges at k dx.
        # The first cell (x = −dx/2) mirrors the second (yelmo.experiment = "MISMIP3D").
        transpose && error("FlowlineBenchmark: transpose = true needs domain = :symmetric " *
                           "(Yelmo has no transposed MISMIP3D boundary treatment).")
        xc = collect(((0:nh) .- 0.5) .* dx)
        x_cf_edge = round(x_cf / dx) * dx
    else
        error("FlowlineBenchmark: unsupported domain = $domain. Supported: :onesided, :symmetric.")
    end
    x_cf_edge < maximum(xc) || error("FlowlineBenchmark: x_cf must lie inside the domain (x_cf = $x_cf, L = $L).")
    yc = collect(((1:Int(ny)) .- (Int(ny) + 1) / 2) .* dx)
    b = FlowlineBenchmark(exp, domain, xc, yc, Float64(dx_km), Float64(A), Float64(C), Float64(m),
                          Float64(n), Float64(smb), Float64(rho_ice), Float64(rho_sw), Float64(g),
                          Float64(sec_year), nh * dx, x_cf_edge, Int(ny), Float64(dxg_start), transpose)
    xg = flowline_xg(b) + b.dxg_start
    xg < b.x_cf - 2dx ||
        error("FlowlineBenchmark: the initial grounding line (x_g = $(round(xg/1e3; digits = 1)) km) " *
              "lies beyond the calving front (x_cf = $(b.x_cf/1e3) km).")
    xg > FLOWLINE_X_SL + 2dx ||
        error("FlowlineBenchmark: the initial grounding line (x_g = $(round(xg/1e3; digits = 1)) km) " *
              "does not lie on the marine part of the bed.")
    return b
end

"""
    flowline_from_fixture(path) -> FlowlineBenchmark

The benchmark of a B2 fixture, from its attributes. Fixtures written before the
`domain` and `dxg_start` attributes existed are symmetric reference starts.
"""
function flowline_from_fixture(path::AbstractString)
    NCDataset(path) do ds
        a = ds.attrib
        domain    = haskey(a, "domain") ? Symbol(a["domain"]) : :symmetric
        dxg_start = haskey(a, "dxg_start") ? a["dxg_start"] : 0.0
        FlowlineBenchmark(Symbol(a["exp"]); A = a["A"], dx_km = a["dx_km"], domain, L = a["L"],
                          x_cf = a["x_cf"], ny = a["ny"], dxg_start, transpose = a["transpose"] == 1,
                          C = a["C"], m = a["m"], n = a["n"], smb = a["smb"], rho_ice = a["rho_ice"],
                          rho_sw = a["rho_sw"], g = a["g"], sec_year = a["sec_year"])
    end
end

"Value of yelmo.experiment for the domain of the benchmark."
flowline_experiment(b::FlowlineBenchmark) =
    b.domain == :onesided ? "MISMIP3D" : (b.transpose ? "periodic-x" : "periodic-y")

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
    lo = FLOWLINE_X_SL + 1.0               # Bed at sea level
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
    flowline_thickness(b, x; xg = flowline_xg(b), ds = 100.0) -> Vector

Semi-analytic thickness [m] at the distances `x` [m] from the divide (|x| is
used): the outer sliding solution of Schoof (2007) upstream of the grounding
line `xg`, and the unconfined shelf downstream of it, up to the calving front
(zero beyond). With the default `xg` (the boundary-layer grounding line) this
is the steady profile; another `xg` gives the perturbed starts.
"""
function flowline_thickness(b::FlowlineBenchmark, x::AbstractVector; xg::Real = flowline_xg(b), ds = 100.0)
    (; A, C, m, n, rho_ice, rho_sw, g) = b
    a    = b.smb / b.sec_year
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

Initial state (semi-analytic profile with the grounding line at
x_g + dxg_start) and fixed forcing at `t = 0`. With `transpose`, x and y are
swapped (flow along y).
"""
function state(b::FlowlineBenchmark, t::Real)
    Float64(t) == 0.0 || error("FlowlineBenchmark.state: only t = 0 is supported (got t = $t).")
    Nx, Ny = length(b.xc), length(b.yc)
    H1 = flowline_thickness(b, b.xc; xg = flowline_xg(b) + b.dxg_start)
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
                 "domain" => String(b.domain), "experiment" => flowline_experiment(b),
                 "dx_km" => b.dx_km, "A" => b.A, "C" => b.C, "m" => b.m, "n" => b.n, "smb" => b.smb,
                 "rho_ice" => b.rho_ice, "rho_sw" => b.rho_sw, "g" => b.g, "sec_year" => b.sec_year,
                 "L" => b.L, "x_cf" => b.x_cf, "ny" => b.ny, "transpose" => Int(b.transpose),
                 "x_g" => flowline_xg(b), "dxg_start" => b.dxg_start,
                 "x_g_start" => flowline_xg(b) + b.dxg_start,
                 "rf_const" => p.rf_const, "cf_ref" => p.cf_ref)
    write_fixture_nc(path, state(b, 0.0), FIELDS_2D; attrs)
    return [path]
end

# -----------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------

"""
    flowline_grounding_line(b, run) -> NamedTuple

Grounding line of a B2 run (`run` is the run directory or its yelmo.nc) at
each output time. The grounding line is the zero of the thickness above
flotation, H_grnd = H_ice − (ρ_w/ρ_i)(z_sl − z_bed), interpolated linearly
between the last grounded and the first floating cell centre of each row (as
f_grnd_acx in Yelmo, ytopo.gl_sep = 1), searching outward from the divide; for
the symmetric domain on both sides. Returns

- `time` [yr]
- `xg` [m]: mean over the rows (and sides)
- `xg_side` [m]: per time, row and side (side 2 is NaN for the one-sided domain)
- `hg` [m]: thickness at x_g, interpolated like H_grnd (mean over rows and sides)
- `qg` [m² a⁻¹]: model flux |ux_bar| H across the face between the last grounded
  and the first floating cell, with the upwind (grounded) thickness as in the
  thickness advection (mean over rows and sides)
- `asym` [m]: max |x_g(+x) − x_g(−x)| (symmetric domain, else NaN)
- `rows` [m]: max spread of x_g across the rows
"""
function flowline_grounding_line(b::FlowlineBenchmark, run::AbstractString)
    file = isdir(run) ? joinpath(run, "yelmo.nc") : run
    t, H, z_bed, z_sl, u = NCDataset(file) do ds
        rd(v) = Float64.(coalesce.(ds[v][:, :, :], 0.0))
        (Float64.(ds["time"][:]), rd("H_ice"), rd("z_bed"), rd("z_sl"),
         b.transpose ? rd("uy_bar") : rd("ux_bar"))
    end
    if b.transpose
        H, z_bed, z_sl, u = (permutedims(v, (2, 1, 3)) for v in (H, z_bed, z_sl, u))
    end
    x  = b.xc
    dx = b.dx_km * 1e3
    Hg = H .- (b.rho_sw / b.rho_ice) .* (z_sl .- z_bed)
    i0, sides = b.domain == :onesided ? (2, (1,)) : (findfirst(==(0.0), x), (1, -1))
    nt, ny = length(t), size(H, 2)
    xg = fill(NaN, nt, ny, 2)
    hg = fill(NaN, nt, ny, 2)
    qg = fill(NaN, nt, ny, 2)
    for k in 1:nt, j in 1:ny, (is, s) in enumerate(sides)
        i = i0
        while 1 <= i + s <= length(x) && Hg[i+s, j, k] > 0
            i += s
        end
        (Hg[i, j, k] > 0 && 1 <= i + s <= length(x)) || continue
        f = Hg[i, j, k] / (Hg[i, j, k] - Hg[i+s, j, k])
        xg[k, j, is] = abs(x[i]) + f * dx
        hg[k, j, is] = (1 - f) * H[i, j, k] + f * H[i+s, j, k]
        # ux_bar(i) lies on the face between cells i and i+1; upwind (grounded-cell)
        # thickness, as in the thickness advection, so that q = a x in steady state
        iface = s > 0 ? i : i - 1
        qg[k, j, is] = abs(u[iface, j, k]) * H[i, j, k]
    end
    ns = length(sides)
    mean_rs(v) = vec(sum(v[:, :, 1:ns]; dims = (2, 3))) ./ (ny * ns)
    asym = ns == 2 ? maximum(abs.(xg[:, :, 1] .- xg[:, :, 2])) : NaN
    rows = maximum(maximum(xg[:, :, 1:ns]; dims = 2) .- minimum(xg[:, :, 1:ns]; dims = 2))
    return (; time = t, xg = mean_rs(xg), xg_side = xg, hg = mean_rs(hg), qg = mean_rs(qg), asym, rows)
end
