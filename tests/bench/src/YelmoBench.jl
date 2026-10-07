module YelmoBench

# Fixtures, online-forcing mirrors and diagnostics for the Yelmo benchmark
# protocol (docs/dev/benchmark-protocol). The benchmark types follow the
# conventions of IceSheetBenchmarks.jl (subtype AbstractBenchmark, implement
# `state` and `write_fixture!`), so that they can move there unchanged.

using NCDatasets
using IceColumnSolutions: IceColumnPar, solve_stationary, solve
using IceSheetBenchmarks: AbstractBenchmark
import IceSheetBenchmarks: state, write_fixture!, analytical_velocity

export state, write_fixture!, analytical_velocity

# Fixture fields: (name, units, long name). xc and yc are written in km.
const FIELDS_2D = (
    ("H_ice",    "m",       "Ice thickness"),
    ("z_bed",    "m",       "Bedrock elevation"),
    ("z_sl",     "m",       "Sea level"),
    ("smb_ref",  "m/yr",    "Surface mass balance"),
    ("T_srf",    "K",       "Surface temperature"),
    ("Q_geo",    "mW m^-2", "Geothermal heat flux"),
    ("bmb_shlf", "m/yr",    "Basal mass balance below floating ice"),
    ("T_shlf",   "K",       "Temperature at the base of floating ice"),
    ("H_sed",    "m",       "Sediment thickness"),
    ("mask_ice", "1",       "Ice mask (0: no ice, 1: fixed, 2: dynamic)"),
)

# Optional 3D fixture fields on the levels zeta (0 at the base, 1 at the surface)
const FIELDS_3D = (
    ("T_ice", "K", "Ice temperature"),
)

"""
    write_fixture_nc(path, s, fields; fields3D = (), attrs = Dict())

Write the 2D fields `fields` of the state NamedTuple `s` to a NetCDF fixture
with axes xc, yc in km. Integer fields (mask_ice) are written as Int32. The 3D
fields `fields3D` are written on the levels `s.zeta`.
"""
function write_fixture_nc(path::AbstractString, s, fields; fields3D = (), attrs = Dict())
    mkpath(dirname(abspath(path)))
    isfile(path) && rm(path)
    NCDataset(path, "c") do ds
        defDim(ds, "xc", length(s.xc))
        defDim(ds, "yc", length(s.yc))
        xv = defVar(ds, "xc", Float64, ("xc",)); xv[:] = s.xc ./ 1e3; xv.attrib["units"] = "km"
        yv = defVar(ds, "yc", Float64, ("yc",)); yv[:] = s.yc ./ 1e3; yv.attrib["units"] = "km"
        for (name, units, longname) in fields
            data = getproperty(s, Symbol(name))
            T = eltype(data) <: Integer ? Int32 : Float64
            v = defVar(ds, name, T, ("xc", "yc"))
            v[:, :] = T.(data)
            v.attrib["units"]     = units
            v.attrib["long_name"] = longname
        end
        if !isempty(fields3D)
            defDim(ds, "zeta", length(s.zeta))
            zv = defVar(ds, "zeta", Float64, ("zeta",)); zv[:] = s.zeta; zv.attrib["units"] = "1"
            zv.attrib["long_name"] = "Normalized height (0: base, 1: surface)"
            for (name, units, longname) in fields3D
                v = defVar(ds, name, Float64, ("xc", "yc", "zeta"))
                v[:, :, :] = getproperty(s, Symbol(name))
                v.attrib["units"]     = units
                v.attrib["long_name"] = longname
            end
        end
        for (k, v) in attrs
            ds.attrib[k] = v
        end
    end
    return path
end

include("column.jl")
include("island4.jl")
include("shelf.jl")
include("strip.jl")
include("flowline.jl")
include("symmetry.jl")
include("compare.jl")
include("budget.jl")
include("margins.jl")

end # module
