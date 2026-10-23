# Yelmo variable tables

Here are tables containing all variables available within Yelmo with dimensions and units.
These tables are used directly in the code for output writing routines, to select which variables to write.
The pages are copied from `input/yelmo-variables-*.md` when the site is rendered (`_scripts/sync-variable-tables.sh`).

The `prec` column sets the precision of the values written to the output files (ncio bit rounding,
which makes the compressed files much smaller):

- `-`: lossless
- `1e-3`: absolute precision, in the units of the field
- `1e-4r`: relative precision
- `1e-3,1e-4r`: both; the larger of the two errors is allowed

The values written carry the attributes `precision_abs`/`precision_rel`. Integer fields are never
rounded. Restart files (`yelmo_restart_write`) are always written lossless.

- [Yelmo topography](yelmo-variables-ytopo.md)
- [Yelmo dynamics](yelmo-variables-ydyn.md)
- [Yelmo material](yelmo-variables-ymat.md)
- [Yelmo tracers](yelmo-variables-ytrc.md)
- [Yelmo thermodynamics](yelmo-variables-ytherm.md)
- [Yelmo hydrology](yelmo-variables-yhyd.md)
- [Yelmo boundaries](yelmo-variables-ybound.md)
- [Yelmo data](yelmo-variables-ydata.md)
