# ----------------------------------------------------------------------
# Initial ice temperature from the stationary column solution of
# IceColumnSolutions.jl (Moreno-Parada et al., 2024): constant properties,
# linear vertical velocity, prescribed surface temperature and geothermal
# heat flux at the base. Cold-base solution only; the model caps the
# temperature at the pressure melting point (thrm_method = "prescribed").
# ----------------------------------------------------------------------

export column_temperature

const COLUMN_K     = 2.1                    # [W m-1 K-1] Heat conductivity
const COLUMN_KAPPA = 2.1 / (910.0 * 2009.0) # [m2 s-1] Heat diffusivity (rho_ice = 910, cp = 2009)

"""
    column_temperature(zeta, H, T_srf, smb, Q_geo; floating = false, T_shlf = 271.15) -> Vector

Temperature [K] on the levels `zeta` of a column of thickness H [m], surface
temperature T_srf [K], surface mass balance smb [m/yr] and geothermal heat flux
Q_geo [mW m-2]. Grounded columns use the stationary solution with the vertical
velocity w0 = max(smb, 0) (no upward flow in ablation zones). Floating columns
use a linear profile from T_shlf at the base to T_srf. Columns thinner than
1 m are at T_srf.
"""
function column_temperature(zeta, H, T_srf, smb, Q_geo; floating = false, T_shlf = 271.15)
    H < 1.0 && return fill(T_srf, length(zeta))
    floating && return T_shlf .+ (T_srf - T_shlf) .* zeta
    par = IceColumnPar(H, T_srf, COLUMN_KAPPA, COLUMN_K, 0.0, Q_geo * 1e-3; w0 = max(smb, 0.0))
    sol = solve_stationary(par; nz = length(zeta))
    sol.zeta ≈ zeta || error("column_temperature: zeta must be uniform on [0, 1] (got $(zeta)).")
    return sol.T_eq
end
