# CEQC 0.0.1 Anubis/RINEX OBS calibration

This release tightens Anubis validation behavior without hard-coding any sample station, satellite, or coordinate.

## RINEX2 transport

Anubis 3.11 can abort on LF-only RINEX2 files in one fixed-column code path even when the records are 80 columns. CEQC now writes RINEX2 output with CRLF line endings. The numeric records and column positions are unchanged. RINEX3/4 continue using LF.

## RINEX4 OBS TYPES

`SYS / # / OBS TYPES` is now derived only from observables that have at least one finite numeric value in the output data. Empty carrier slots used only to preserve LLI/SSI no longer create a declared observation type. This keeps RINEX3 and RINEX4 OBS type declarations aligned for the same source observations.

## GAL/QZS validity

CEQC emits Galileo/QZSS observation systems only when the source contains valid 109x/111x MSM observations with finite values. If no `GALSUM`/`QZSSUM` appears in Anubis for the current RTCM3 sample, it means the sample did not produce valid E/J observations; CEQC does not synthesize them.
