# Anubis QC comparison coverage

CEQC release validation now covers three output families with Anubis 3.11:

1. **RINEX2 legacy GPS/GLO subset** — intended for teqc-era checks.  CEQC writes fixed 80-column RINEX2 OBS/NAV records; the validation script creates CRLF copies for Anubis 3.11 on Linux because this Anubis build aborts on LF-only RINEX2 fixed-column records before QC starts.  This does not alter numeric data.
2. **RINEX3 modern mixed OBS/NAV** — GPS/GLO/BDS/GAL/QZS observations where available.
3. **RINEX4 strict typed NAV/OBS** — RINEX4 typed NAV records are kept for strict tools; RTKLIB/RTKPLOT should still use `+rtklib` compatibility output.

Expected interpretation:

- RINEX2 should be compared only for old GPS/GLONASS L1/L2-compatible metrics.
- RINEX3/4 Anubis modern checks confirm that CEQC OBS for GPS/GLO/BDS are parsed and that epoch/sample counts close.
- Anubis 3.11 may not use all broadcast NAV constellations for elevation/position in the same way CEQC/gfzrnx do.  If an XTR line reports `woElev` for GLO/BDS while OBS counts still close, treat it as a navigation-usage limitation or constellation-specific NAV support gap, not an OBS parser failure.
- Filename/marker mismatch warnings are not RINEX syntax errors; rename the file to match the four-character marker or set `-O.mo` when such warnings need to be eliminated.

Run example:

```bash
tools/anubis_compare.sh \
  --anubis /path/to/anubis-3.11-lin-static-64b \
  --ceqc ./ceqc-linux-amd64 \
  --input sample.rtcm3 \
  --tr rtcm3 \
  --outdir anubis_compare_out
```

For an existing RINEX OBS/NAV pair:

```bash
tools/anubis_compare.sh \
  --anubis /path/to/anubis-3.11-lin-static-64b \
  --ceqc ./ceqc-linux-amd64 \
  --input sample.26o \
  --nav sample.26p \
  --tr rinex \
  --outdir anubis_compare_out
```
