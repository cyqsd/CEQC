# RINEX support notes for ceqc 0.0.1

`ceqc` parses and writes RINEX 2, 3 and 4 families as separate format families rather than treating the version line as a cosmetic header.

## OBS

- RINEX 2.x: legacy epoch records, global `# / TYPES OF OBSERV`, GPS-only numeric PRNs and mixed-system `G/R/E/C/J/I/S` satellite identifiers in the epoch satellite list.
- RINEX 3.x: `>` epoch records, `SYS / # / OBS TYPES`, multi-GNSS satellite identifiers, continuation lines for long observation-type lists.
- RINEX 4.x: OBS files use the RINEX 3-style epoch/observation body with 4.x headers and are parsed through the same multi-GNSS observation model.

## NAV

- RINEX 2.x: GPS legacy NAV by default; GLONASS/Galileo/BeiDou/QZSS/NavIC/SBAS system hints in the version/type line are preserved when present.
- RINEX 3.x: mixed-system NAV records are parsed by satellite prefix (`G`, `R`, `E`, `C`, `J`, `I`, `S`).
- RINEX 4.x: typed records such as `> EPH G01 LNAV ...` are preserved with record type, satellite, message type and subtype.

## Explicit subversions

The command line accepts both major selectors and explicit subversions:

- `+v2`, `+v2.10`, `+v2.11`
- `+v3`, `+v3.00`, `+v3.03`, `+v3.04`, `+v3.05`
- `+v4`, `+v4.00`, `+v4.01`, `+v4.02`

A bare major selector uses ceqc's default subversion for that family; an explicit `.00/.03/.05` selector is kept as requested.


## 0.0.1 additions

- Header label recovery now handles slightly short RINEX header lines that shift known labels left of column 61, while preserving exact fixed-column behavior for normal files.
- RINEX 4 NAV field models are constellation-specific for GPS/QZSS LNAV/CNAV, GLONASS FDMA, Galileo INAV/FNAV, BeiDou D1/D2, and NavIC-like records.
- RINEX 4 system-data records such as `> STO GPUT ...` are grouped as non-EPH NavigationRecord entries instead of being confused with satellite ephemerides.
