# Native QC JSON export

CEQC can export the native `+qc` result as a chart-friendly JSON document:

```bash
ceqc +qc +qc_json qc.json obs.26o
ceqc +qc +nav nav.26p +qc_json qc.json obs.26o
ceqc +qc +qc_json qc.json raw.rtcm3
ceqc +qc +qc_json qc.json raw.ubx
```

The JSON export belongs to the native CEQC QC branch.  It is not part of the
`+teqc` text compatibility renderer.  If `+qc_json` is combined with `+teqc`, CEQC
prints a warning and does not write the JSON file.

## Main top-level sections

- `source`: input path, RINEX kind/version, marker, receiver, antenna and NAV input files.
- `time_window`: first/last epoch and estimated observation interval.
- `counts`: epoch, observation, navigation, system and satellite counts.
- `qc`: active options, epoch-SV statistics, observation completeness counters and legacy GPS/GLO counters.
- `snr`: SNR statistics by observation code and a chart-ready array.
- `ionosphere`: ionospheric diagnostic statistics and histograms.
- `multipath`: MP RMS values, samples and chart-ready multipath rows.
- `timeplots`: native compact timeplot strings, observation bin counts and per-satellite rows.
- `skyplot`: per-satellite azimuth/elevation points sampled at maximum elevation in the QC window.
- `rise_set`: per-satellite first/last observation, duration and max elevation.
- `position`: position QC counters, average/header XYZ and optional every-epoch position rows.
- `residuals`: residual QC totals and per-system RMS rows.
- `warnings`: threshold, position, residual and raw translator warnings.
- `charts`: flattened arrays intended for frontend plotting.
