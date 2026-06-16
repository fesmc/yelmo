# yelmo-config

A command-line tool to **discover, manage, compare and validate** Yelmo
ice-sheet model parameter files.

Yelmo now reads a canonical default parameter set from
[`input/yelmo_defaults.nml`](../../input/yelmo_defaults.nml); a user par file
only needs to list the parameters it wants to override. `yelmo-config` makes
that workflow easy: it knows every parameter, its default, its documentation and
its constraints, and can build, diff and check par files accordingly.

## Install

System-wide (from a Yelmo checkout):

```bash
pip install -e tools/yelmo-config      # editable: always tracks the repo's defaults
# or
pip install tools/yelmo-config
```

This installs the `yelmo-config` console command. Python ≥ 3.11, no third-party
dependencies.

### How it finds the model files

The tool reads `input/yelmo_defaults.nml` (parameters, defaults, docs) and
`src/*.f90` (enum constraints). It locates them, in order, via:

1. `--defaults PATH` / `--src PATH`
2. `$YELMO_DEFAULTS`, or `$YELMO_ROOT` (expects `$YELMO_ROOT/input/yelmo_defaults.nml`)
3. walking up from the current directory
4. walking up from the installed package (works for `-e` installs)

So run it from inside a Yelmo checkout, or set `YELMO_ROOT`, or pass `--defaults`.

## Commands

| Command | Purpose |
|---|---|
| `groups` | List the canonical library groups and the `&yelmo` component→group mapping. |
| `list [group]` | List parameters with defaults and docs, grouped by section. |
| `show group.name` | Full detail for one parameter: default, type, doc, allowed values, constraints. |
| `search PATTERN` | Regex search over parameter names and documentation. |
| `write [user.nml] [-o out]` | Emit a **complete**, documented par file from the defaults, with any overrides from `user.nml` applied. |
| `diff A [B] [--raw]` | Compare two par files, or one against the defaults. Supports per-group selectors. |
| `check FILE [--files]` | Validate a user par file against the schema and all constraints. |
| `completion {bash,zsh}` | Print a shell completion script. |

### Output formats

`list`, `show`, `search`, `write`, `diff` and `check` accept `--format`
(after the subcommand):

- `text` (default) — human-readable, coloured.
- `compact` — one `group.par=value` per line; greppable and diff-friendly.
- `json` — structured output for scripting (`jq`, etc.).

```bash
yelmo-config list ydyn --format compact         # ydyn.solver="diva" ...
yelmo-config write par/run.nml --format compact # full effective config, flat
yelmo-config check par/run.nml --format json | jq '.issues'
yelmo-config diff a.nml b.nml --format json
```

### Examples

```bash
yelmo-config list ydyn                       # browse the dynamics parameters
yelmo-config show ytill.z0                    # one parameter in detail
yelmo-config search 'calv'                    # find calving-related parameters

yelmo-config write -o par/complete.nml        # full default par file
yelmo-config write par/yelmo_initmip.nml      # defaults + initmip overrides -> stdout

yelmo-config check par/yelmo_initmip.nml      # validate a par file

# diff the effective config of a run against the defaults:
yelmo-config diff par/yelmo_initmip.nml

# diff two runs (effective config = defaults + each file's overrides):
yelmo-config diff runA.nml runB.nml
yelmo-config diff --raw runA.nml runB.nml     # only params literally written in each file
```

### Per-group diff and renamed groups

A par file may rename a component group via `&yelmo`, e.g.
`nml_ydyn = "ydyn_north"` means the `&ydyn_north` group is read in place of
`&ydyn`. `yelmo-config` resolves this everywhere (validation, diff, merge).

You can diff a single group, including across differently-named groups, with the
`file:group` selector:

```bash
# compare the dynamics block of two runs even if one renamed it:
yelmo-config diff runA.nml:ydyn_north runB.nml:ydyn

# compare one file's group against the defaults for that group:
yelmo-config diff runA.nml:ydyn_north
```

### Shell completion

Completion covers subcommands, options, `--format` values, and — when run inside
a checkout — group names (`list`) and `group.name` parameters (`show`).

```bash
# bash
yelmo-config completion bash | sudo tee /etc/bash_completion.d/yelmo-config
#   or, per-shell:  echo 'source <(yelmo-config completion bash)' >> ~/.bashrc

# zsh  (into a directory on your $fpath)
yelmo-config completion zsh > "${fpath[1]}/_yelmo-config"
#   or, per-shell:  echo 'source <(yelmo-config completion zsh)' >> ~/.zshrc
```

## What `check` validates

Only the canonical library groups are validated; driver groups such as `&ctrl`
or `&set_*` are reported once as informational pass-throughs.

- **unknown parameters** — names not declared in the defaults for that group
- **enum violations** — values outside the `yelmo_check_enum` allowed set
- **range violations** — scalars outside their numeric range
- **ordering violations** — inter-parameter ordering (e.g. `sd_min < sd_max`)
- **type mismatches** — e.g. a string where the default is numeric (warning)
- **missing `domain` / `grid_name`** — still set to the `"None"` dummy
- **missing file inputs** (`--files`) — a `*_path` whose `*_load` is true but
  whose file is absent (after `{domain}` / `{grid_name}` substitution)

## Where constraints come from

- **Enum** and **file** constraints are extracted automatically from the
  `yelmo_check_enum` / `yelmo_check_file` calls in `src/*.f90` — the model itself
  is the single source of truth.
- **Range** and **ordering** constraints are curated in
  [`yelmo_config/data/constraints.toml`](yelmo_config/data/constraints.toml),
  because the Fortran refers to these parameters by struct-member name (which
  differs from the namelist name) inside free-form `if` blocks. Keep that file
  in sync if those checks change in `src/`.
