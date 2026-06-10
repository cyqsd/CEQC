# RINEX strict validation notes for ceqc 0.0.1

This release separates standards-oriented output from RTKLIB compatibility output.

## Fixes in 0.0.1

- RINEX 2 OBS no longer copies RINEX 3/4-only `GLONASS SLOT / FRQ #` or `GLONASS COD/PHS/BIS` header records.  Old readers such as teqc treat those lines as malformed in RINEX 2.
- RINEX 4 typed NAV `EPH` records no longer expose CEQC translator source strings such as `RTCM1019`, `RTCM1042`, `RTCM1020`, `UBX-SFRBX`, `RAW`, or `SFRBX-ASSEMBLED` in the typed marker.  They are written with the neutral source `0`, for example:

```text
> EPH C01 D1 0
> EPH G01 LNAV 0
> EPH R05 FDMA 0
```

This avoids strict readers mapping the source token to `XXXX`.

## Validation performed

The sample `26042313.2026156binRTCM3` was translated into RINEX 2.11, 3.05 and 4.02 OBS/NAV.  Internal strict checks verified:

- all header lines are exactly 80 characters with labels in columns 61-80;
- RINEX 2 OBS does not contain RINEX 3/4-only GLONASS header labels;
- RINEX 3/4 epoch satellite counts match actual following records;
- RINEX 4 typed NAV markers do not contain internal diagnostic source tags;
- CEQC `+verify` reports OK for all generated OBS/NAV files.

Anubis 3.11 successfully read RINEX 3.05 and 4.02 OBS/NAV.  The only remaining Anubis warning is the expected filename/site-name mismatch when the test filename does not match the marker name.

## gfzrnx

gfzrnx is not bundled with CEQC and was not available in the build container.  The changes above target the exact strict-format issues that gfzrnx commonly reports after the 0.0.1 changes: RINEX2-only header leakage and unsupported RINEX4 typed NAV source tokens.
