# The validation documentation moved

**2026-08-10.** Everything that used to live in `athenapk/docs/validation/` is now at:

```
/beegfs/u/bbg6470/validation/docs/work_packages/
```

This stub exists because ~40 production input decks and submit scripts carry comments
of the form `see docs/validation/WP13b_restart_gpu_amr.md`. Those pointers were left
alone on purpose — rewriting 24 live ensemble decks to fix a comment is churn with a
non-zero chance of touching something that matters. Read them as
`/beegfs/u/bbg6470/validation/docs/work_packages/<file>` instead.

| old path | new path |
|---|---|
| `docs/validation/WP*.md` | `/beegfs/u/bbg6470/validation/docs/work_packages/WP*.md` |
| `docs/validation/B*.md`, `D1_*.md` | `/beegfs/u/bbg6470/validation/docs/work_packages/` |
| `docs/validation/PRODUCTION_SWITCH_2026-08-03.md` | `/beegfs/u/bbg6470/validation/docs/work_packages/` |
| `docs/validation/scripts/*.py` | `/beegfs/u/bbg6470/validation/docs/work_packages/scripts/` |

See `/beegfs/u/bbg6470/validation/README.md` for the whole layout.
