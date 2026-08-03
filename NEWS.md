# NEWS

Notable changes to the **Heat-Stress (WBGT) Forecast Tools**.
Newest entries first. Dates are release/commit dates (EAT).

---

## 2026-08-03 — Weekly post-processing: day+night in one run, redesigned category legend

All changes in `postproc_WBGT_weekly_forecast.py`. No other pipeline stage is
affected, and output file names are unchanged.

### Added

- **Daytime and nighttime are post-processed in a single run.** The `VAR`
  setting is replaced by `VARIABLES = ["WBGT_tmax", "WBGT_tmin"]`; the script
  loops over both and writes all four figures per forecast date. Remove an
  entry to process one variable alone. A variable missing from the input file
  is reported and skipped instead of raising.
- **Single source of truth for the heat-stress bands.** A new `CATEGORIES`
  table holds `(class name, WBGT range, advisory, colour)` per band, with an
  assertion that it stays in step with `CAT_BOUNDS`.

### Changed

- **Category map legend redesigned for briefing use.** The vertical colorbar
  now carries the class names only (*Low → Extreme*), and a horizontal legend
  below the panels pairs each colour with its WBGT range and recommended
  action (`< 25 °C - normal activity` … `≥ 32 °C - avoid outdoor work`),
  arranged on two rows with a bold title.
- **Equal-height category blocks on the colorbar.** The bar is drawn from an
  index-scale `BoundaryNorm` proxy rather than from the contour set, so each
  class gets the same height; previously the wide `0–25` and `32–60` bands
  took ~90% of the bar between them and squashed the operationally important
  classes into a sliver.
- **Colorbar geometry follows the maps.** Its box is read back from the panel
  bounding box after a first draw, so it spans exactly the panel height and
  tracks the right-hand panel's edge. It stays correct when `EXTENT`, the
  figure size, or the number of weeks changes.
- Figure code refactored into `plot_summary()` and `plot_categories()`, using
  explicit `fig.*` calls and closing each figure after saving.

### Fixed

- **Week 1 spanned eight days instead of seven.** `week1_end` used
  `timedelta(days=7)` against an inclusive date slice. Week 1 is now day 0–6
  and Week 2 day 7–13, so both weeks are seven days and no day falls between
  them. *Behaviour change:* Week 1's mean-of-daily-maxima now averages seven
  days, and its weekly maximum no longer includes the day that moved into
  Week 2.
- Removed the superseded legend block on the category figure: a duplicated
  colorbar, a sixth phantom advisory label with no matching colour, and a
  second legend of transparent patches that never rendered.

### Notes

- The category maps classify the **weekly maximum** WBGT — the single hottest
  6-hourly window of the week per grid cell — not the weekly average. This is
  deliberately worst-case for advisory use; the summary figure's *(maximum)*
  column shows the same field on a continuous scale.
- Both variables still share the 23–33 °C colour scale in the summary figure
  so daytime and nighttime panels can be compared directly.
