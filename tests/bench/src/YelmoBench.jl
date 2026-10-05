module YelmoBench

# Fixtures, online-forcing mirrors and diagnostics for the Yelmo benchmark
# protocol (docs/dev/benchmark-protocol). The benchmark types follow the
# conventions of IceSheetBenchmarks.jl (subtype AbstractBenchmark, implement
# `state` and `write_fixture!`), so that they can move there unchanged.

using NCDatasets
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

"""
    write_fixture_nc(path, s, fields; attrs = Dict())

Write the 2D fields `fields` of the state NamedTuple `s` to a NetCDF fixture
with axes xc, yc in km. Integer fields (mask_ice) are written as Int32.
"""
function write_fixture_nc(path::AbstractString, s, fields; attrs = Dict())
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
        for (k, v) in attrs
            ds.attrib[k] = v
        end
    end
    return path
end

include("island4.jl")
include("shelf.jl")
include("symmetry.jl")
include("compare.jl")

end # module
