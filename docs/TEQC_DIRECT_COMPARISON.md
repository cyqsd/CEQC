# teqc comparison scope

`teqc 2019Feb25` is useful as a regression reference only for the subset that it can represent: mainly legacy RINEX 2 GPS/GLONASS L1/L2 observation QC and the corresponding old MP12/MP21 style statistics.  It must not be used as a golden reference for modern features that teqc never implemented.

## Comparable with teqc

Use teqc as a reference for these metrics only when the OBS/NAV input pair and time window are the same:

- GPS and GLONASS legacy satellite timelines for the common satellites.
- Observation interval, epoch count, repeated epochs, and duplicate-satellite checks.
- Legacy `Moving average MP12` / `Moving average MP21` for GPS/GLONASS L1/L2.
- Mean `S1` / `S2` for legacy L1/L2 observations.
- Basic GPS unhealthy-SV and GPS/GLONASS missing OBS/NAV lists when the same NAV support files are available.
- Position summary only as a sanity check; differences can appear because teqc ignores modern systems and may reject NAV records whose Toe/Tow difference exceeds its legacy limits.

## Not comparable with teqc

These CEQC outputs are intentionally modern-only and must be evaluated separately:

- BeiDou, Galileo, QZSS, SBAS, and mixed RINEX 3/4 system handling.
- BDS-3 new signal codes and bands such as B1C/B1A, B2a, B2b, B2ab, B3.
- Modern multipath combinations such as BDS `MP26`, `MP27`, `MP62`, `MP67`, GPS/QZSS/Galileo `MP15`, `MP51`, `MP52`, and any combination beyond teqc's L1/L2 MP12/MP21.
- RINEX 3/4 typed NAV records and UBX SFRBX page inventory / guarded ICD decode diagnostics.
- S5/S6/S7/S8 signal-strength statistics and any non-L1/L2 frequency summary.
- RTCM3 MSM native decoding behavior, especially modern cell masks, lock indicators, half-cycle flags, and fine phase-rate fields.

## Practical workflow

1. Generate a legacy RINEX 2 GPS/GLONASS subset for teqc.
2. Run teqc on that legacy OBS/NAV pair.
3. Run CEQC on the same legacy pair or on the original stream.
4. Compare only the common subset listed above.
5. Review CEQC's modern-only section separately; do not treat differences there as teqc regressions.

If the teqc report and the CEQC report show different calendar dates, different OBS/NAV inputs, or different enabled systems, the numeric summary is not a valid direct comparison.  First align the input pair and time window, then compare the common metrics.
