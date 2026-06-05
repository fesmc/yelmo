- Production code of Yelmo ice-sheet model. Be very careful about what is modified in here.
- Initialized from tagged version v1.15 from github.com/palma-ice/yelmo repository. New repository here is designed to support development towards new version v2.0.

## New worktree setup

A fresh worktree has no `Makefile`, no `.runme_config`, and no data symlinks. Before building or running anything, set these up from inside the worktree:

```bash
ln -s /Users/alrobi001/models/ice_data ice_data
ln -s /Users/alrobi001/models/fesm-utils fesm-utils
ln -s /Users/alrobi001/models/fasthydrology FastHydrology
runme --config                              # creates .runme_config from .runme/runme_config
python config.py config/macbook_gfortran   # or whichever host config applies
```

After that, `make <target>` and `runme ...` work the same as in the main tree.
`runme` is installed system-wide via pip (`pip install git+https://github.com/fesmc/runme`), not a local script.

The `FastHydrology` symlink (capitalized to match the upstream repo name and
`FASTHYDROROOT` in `config/common.mk`) is expected by every host config in
`config/`. `libfasthydro.a` must be built separately by running
`make fasthydro-static` inside the fasthydrology checkout — yelmo's Makefile
does not recurse into it (same convention as `fesm-utils`).
