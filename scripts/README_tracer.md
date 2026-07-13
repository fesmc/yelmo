# Running the %trc passive-tracer backend

The `%trc` subsystem runs three interchangeable englacial tracing backends —
`euler` (in-tree), `trc` (Lagrangian particles, `tracer` library) and `elsa`
(layer stack, `elsa` library). This note covers activating them from the command
line without bespoke namelist files.

## Defaults files (the schema)

`tracer` and `elsa` each own a defaults file listing every parameter of their
namelist group, mirroring yelmo's `input/yelmo_defaults.nml`:

- `input/tracer_defaults.nml` — group `&trc` (synced copy of `tracer/input/tracer_defaults.nml`)
- `input/elsa_defaults.nml`   — group `&elsa` (synced copy of `elsa/input/elsa_defaults.nml`)

A run's own namelist need only list the entries it overrides; anything omitted
falls back to these defaults, and `nml_validate` rejects unknown keys (typos).

> These copies are kept in sync by hand for now. If you edit a default in the
> `tracer`/`elsa` repo, copy it back into `yelmo/input/`.

## Activating the tracer with `runme -p`

`scripts/run_trc_eismint.sh` runs EISMINT-moving with the particle tracer and its
gridded stats, set up purely with `runme -p`:

```bash
scripts/run_trc_eismint.sh [OUTDIR]      # default OUTDIR: output/trc-eismint
```

Under the hood:

```bash
runme -r -e benchmarks -o OUTDIR -n par/yelmo_EISMINT_moving.nml \
    -p ctrl.time_end=2000 ctrl.dt2D_out=1000 \
       ytrc.use_euler=True ytrc.use_tracer=True ytrc.calc_age=True \
       ytrc.tracer_nml=yelmo_EISMINT_moving.nml \
       trc.stats=True
```

- `ytrc.tracer_nml=yelmo_EISMINT_moving.nml` self-references the staged par
  (runme copies it into OUTDIR under the same basename and runs there), so the
  `&trc` group and the `&ytrc` switches live in one file.
- All other `&trc` parameters (`n`, `dt_dep`, `n_depth`, `time_iso`, …) come from
  `input/tracer_defaults.nml`.

### Two `runme -p` constraints

1. **It overrides existing keys, it does not create them.** The base par carries
   a dormant `&trc` group (just `stats = False`) so `-p trc.stats=True` has a key
   to hit. To override another `&trc` entry via `-p`, add it to that group first.
2. **Commas mean an ensemble.** Array parameters like `time_iso` cannot be set
   with `-p` (`a=1,2,3` defines an ensemble dimension). Set arrays in the par or
   in the defaults file instead — hence keeping `time_iso` as a default.

## Output

The tracer's gridded stats are written to `OUTDIR/yelmo_restart.nc`:

- `trc_count(time, depth_norm, yc, xc)` — particles per normalized-depth band
- `trc_depth_iso(time, time_iso, yc, xc)` — isochrone depth from the particle cloud

The benchmark driver's 2D writer emits a fixed variable list, so these do not
appear in `yelmo2D.nc` unless a driver requests them by name. The `t=0` restart
is the intended product for comparison with data.

Set `trc.stats=False` (the default) to disable stats entirely for production.
