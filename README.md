# CEQC

CEQC is a C++ GNSS preprocessing and quality-checking tool with a teqc-compatible command style.  This build is the full project line, not the earlier UBX-only simplification: it includes RINEX 2/3/4 parsing and synthesis, RTCM3, u-blox UBX RAWX/SFRBX, RINEX editing/splicing/decimation, QC, teqc-style reporting, Linux and Windows builds.

## Project status and compatibility notice

CEQC is an early-stage open-source project.  Please use it with caution and independently validate any generated OBS/NAV/QC output before relying on it in production, survey, operational, or safety-relevant workflows.

Current development is focused on standard RINEX, RTCM3, and u-blox UBX RAWX/SFRBX workflows.  Some older RTK receivers, receiver-private legacy formats, and vendor-specific historical behaviors will not be supported, or may remain only partially implemented.

## Important 0.0.1 fixes

- `+help` / `-help` now prints CEQC's own help text.  It intentionally uses CEQC names and CEQC implementation status instead of embedding the original teqc help text verbatim.
- RINEX target versions are now honored consistently in QC preprocessing and output: `+v2` -> RINEX 2.11, `+v3` -> RINEX 3.05, `+v4` -> RINEX 4.02.
- `rinex:` in QC output reflects the actual parsed/generated RINEX version instead of always showing `3.05`.
- RINEX 4 NAV typed records keep the `> EPH ...` marker, while RINEX 3 NAV omits it and RINEX 2 NAV is restricted to GPS legacy records where possible.

## Build

```bash
cmake -S . -B build -DCEQC_CXX_STANDARD=21 -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Windows cross build:

```bash
cmake -S . -B build-win \
  -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/mingw-w64-x86_64.cmake \
  -DCEQC_CXX_STANDARD=21 \
  -DCEQC_BUILD_TESTS=OFF \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-win -j
```

## Examples

```bash
ceqc +help
ceqc +verify testdata/minimal_v3.23o
ceqc +v2 +qc testdata/minimal_v3.23o
ceqc +v4 +obs out.26o testdata/minimal_v3.23o
ceqc -tr rtcm3 +v3 +obs obs.26o +nav nav.26p 26042313.2026156binRTCM3
ceqc -tr ubx +v2 +G +obs gps.26o raw.ubx
```


## RINEX 2/3/4 subversions

ceqc 0.0.1 preserves parsed RINEX subversions and accepts explicit output selectors such as `+v2.10`, `+v3.03`, `+v4.00` and `+v4.02`.  OBS parsing now handles RINEX 2 mixed-system satellite identifiers and RINEX 3/4 `SYS / # / OBS TYPES` continuation lines.  See `docs/RINEX_SUPPORT.md`.

## 0.0.1 update

- RINEX 4 NAV field mapping is now constellation-specific for GPS/QZSS LNAV/CNAV, GLONASS FDMA, Galileo INAV/FNAV, BeiDou D1/D2 and NavIC-like records.
- RINEX 4 system-data records such as `> STO GPUT ...` are grouped as non-EPH `NavigationRecord` entries rather than being confused with satellite ephemerides.
- RINEX header label recovery handles short or slightly shifted lines where known labels appear before column 61.
- Additional regression tests cover Galileo/BeiDou RINEX 4 NAV field names and non-EPH RINEX 4 system-data records.
- OBS `INTERVAL` is not a fixed constant.  During synthesis/merge/window/decimation CEQC now derives the nominal interval from the actual unique observation epoch sequence, using the dominant adjacent epoch delta in a teqc-like way.
- Stale input `INTERVAL` header lines are not preserved when a new OBS file is synthesized; the output interval is recomputed from the output epochs.

## Questions or Issues?

For any questions or inquiries, please contact me via email at root@cyqsd.cn.