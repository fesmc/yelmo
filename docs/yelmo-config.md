# yelmo-config

`yelmo-config` is a command-line tool to **discover, manage, compare and
validate** Yelmo parameter files. It lives in the Yelmo repository under
[`tools/yelmo-config/`](https://github.com/fesmc/yelmo/tree/main/tools/yelmo-config).

## The defaults + overrides model

Yelmo reads a canonical default parameter set from `input/yelmo_defaults.nml`,
which declares **every** parameter the library uses together with its default
value and an inline comment documenting its meaning. A user parameter file only
needs to list the parameters it wants to **override** — anything missing is taken
from the defaults. Unknown parameter names (not declared in the defaults) are
reported as errors at load time.

Each component group can also be renamed via the `&yelmo` block — for example
`nml_ydyn = "ydyn_north"` tells Yelmo to read the dynamics parameters from a
`&ydyn_north` group instead of `&ydyn`. `yelmo-config` understands this mapping
everywhere it inspects a file.

In addition to the defaults, the model enforces consistency checks at load time:

- **enum** checks — a parameter whose value must be one of a fixed set
  (e.g. `ydyn.solver`), including checks that depend on another parameter
  (e.g. `ycalv.calv_flt_method` allows a different set depending on `use_lsf`);
- **range** checks — e.g. `yelmo.cfl_max` must be in `(0, 1]`;
- **ordering** checks — e.g. `ycalv.sd_min < sd_max`.

`yelmo-config` mirrors all of these so you can catch problems before launching a
run.

## Installation

Install system-wide straight from the repository:

```bash
pip install -U "git+https://github.com/fesmc/yelmo#subdirectory=tools/yelmo-config"
```

Or from a local checkout:

```bash
pip install -e tools/yelmo-config     # editable: tracks the checkout's files
```

This installs the `yelmo-config` command. Python ≥ 3.11 is required; there are no
third-party dependencies. To update later:

```bash
yelmo-config update            # latest from the default branch
yelmo-config update dev        # a specific branch / tag / commit SHA
```

### Works with or without a checkout

The tool needs the defaults (`input/yelmo_defaults.nml`) and the enum
constraints (extracted from `src/*.f90`). A snapshot of both is **bundled into
the package**, so `yelmo-config` works even with no copy of Yelmo present. When
run inside a checkout it uses that copy's live files instead; if the local files
differ from the bundled snapshot, it prints a `warning:` so you know the
installed tool is out of date relative to your checkout (or vice-versa).

Resolution order: `--defaults` / `--src` → `$YELMO_DEFAULTS` / `$YELMO_ROOT` →
the current directory (walking up) → the installed package → the bundled
snapshot.

## Commands

| Command | Purpose |
|---|---|
| `groups` | List the canonical groups and the `&yelmo` component→group mapping. |
| `list [group]` | List parameters with defaults and docs, grouped by section. |
| `show group.name` | Full detail for one parameter: default, type, doc, allowed values, constraints. |
| `search PATTERN` | Regex search over parameter names and documentation. |
| `write [user.nml]` | Emit a complete, documented parameter file from the defaults, applying any overrides. |
| `diff A [B]` | Compare two parameter files, or one against the defaults. |
| `check FILE` | Validate a parameter file against the schema and all constraints. |
| `update [ref]` | Self-update the tool from its git repo. |
| `snapshot` | Refresh the bundled defaults/enums snapshot from a checkout (maintainers). |
| `completion {bash,zsh}` | Print a shell completion script. |

### Discovering parameters

```bash
yelmo-config groups                  # overview of all groups
yelmo-config list ydyn               # the dynamics parameters, with docs
yelmo-config show ytill.z0           # one parameter in detail
yelmo-config search 'calv'           # find calving-related parameters
```

### Building a parameter file

```bash
yelmo-config write -o par/complete.nml        # full default parameter file
yelmo-config write par/yelmo_initmip.nml      # defaults + initmip overrides
```

### Comparing parameter sets

`diff` compares the *effective* configuration (defaults + each file's overrides),
so it shows only meaningful differences:

```bash
yelmo-config diff par/yelmo_initmip.nml       # a file vs the defaults
yelmo-config diff runA.nml runB.nml           # two runs
yelmo-config diff --raw runA.nml runB.nml     # only params literally written
```

You can also compare a single group, even across differently-named groups:

```bash
yelmo-config diff runA.nml:ydyn_north runB.nml:ydyn
```

### Validating a parameter file

```bash
yelmo-config check par/yelmo_initmip.nml
```

`check` reports unknown parameters, enum / range / ordering violations, type
mismatches, an unset `domain` / `grid_name`, and (with `--files`) missing input
files. Only the Yelmo library groups are validated; driver groups such as
`&ctrl` are reported once as informational pass-throughs.

### Output formats

`list`, `show`, `search`, `write`, `diff` and `check` accept `--format`
(after the subcommand):

- `text` (default) — human-readable, coloured;
- `compact` — one `group.par=value` per line;
- `json` — structured output for scripting.

```bash
yelmo-config check par/run.nml --format json | jq '.issues'
yelmo-config write par/run.nml --format compact
```

### Shell completion

```bash
# bash
yelmo-config completion bash | sudo tee /etc/bash_completion.d/yelmo-config
# zsh (into a directory on your $fpath)
yelmo-config completion zsh > "${fpath[1]}/_yelmo-config"
```

## Keeping the bundled snapshot current

The bundled defaults/enums snapshot is committed in the repository. Maintainers
should regenerate it from a checkout whenever the defaults or the enum checks in
`src/` change, then commit the result so installs from GitHub ship current data:

```bash
yelmo-config snapshot
```
