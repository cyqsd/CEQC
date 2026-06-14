# RINEX version-profile implementation notes

CEQC treats `+v2.10`, `+v2.11`, `+v3.00` ... `+v4.02` as versioned output
profiles, not as a cosmetic rewrite of the first header line.  The profiles are
implemented from the public IGS/RTCM RINEX specifications and validated with
`gfzrnx` where test data is available.

This document is intentionally conservative: when a sub-version requires a
feature that CEQC cannot yet prove with a real RTCM/UBX/RINEX test sample, the
feature is listed under **not yet fully covered** rather than silently being
claimed as complete.

## Implemented profile gates

| Requested version | CEQC output family | Constellation gate used by CEQC | Main implemented differences |
|---|---|---|---|
| `+v2.10` | RINEX 2 legacy | OBS: GPS/GLONASS; NAV: GPS LNAV only | RINEX2 epoch records, global `# / TYPES OF OBSERV`, RINEX2 NAV continuation indent, CRLF output. |
| `+v2.11` | RINEX 2 legacy | OBS: GPS/GLONASS; NAV: GPS LNAV only | Same conservative standard profile as 2.10 for CEQC output. CEQC does not use non-IGS RINEX2 mixed extensions for BDS/QZSS. |
| `+v3.00` | RINEX 3 early multi-GNSS | GPS/GLONASS/Galileo | `>` epoch records, `SYS / # / OBS TYPES`, RINEX3 mixed NAV body; BDS/QZSS/SBAS/NavIC are filtered from output. |
| `+v3.01` | RINEX 3 early multi-GNSS | GPS/GLONASS/Galileo | Same CEQC profile gate as 3.00. SBAS/GEO-specific 3.01 details are not emitted because CEQC has no validated SBAS test vector yet. |
| `+v3.02` | RINEX 3 current-GNSS profile | GPS/GLONASS/Galileo/BDS/QZSS/SBAS | Enables BDS/QZSS/SBAS output families and the RINEX3 mixed NAV body. NavIC is kept out of 3.02 until a validated 3.02 NavIC sample is added. |
| `+v3.03` | RINEX 3 current-GNSS profile | GPS/GLONASS/Galileo/BDS/QZSS/SBAS/NavIC | Same writer as 3.02 plus NavIC/IRNSS gate when CEQC has records. |
| `+v3.04` | RINEX 3 updated signal profile | GPS/GLONASS/Galileo/BDS/QZSS/SBAS/NavIC | Same structural writer as 3.03; CEQC keeps strict observation-code validation delegated to gfzrnx. |
| `+v3.05` | RINEX 3 final profile | GPS/GLONASS/Galileo/BDS/QZSS/SBAS/NavIC | Same structural writer as 3.04, with BDS-2/BDS-3 and extended GLONASS NAV fields supported by CEQC's decoded field model. |
| `+v4.00` | RINEX 4 typed NAV | GPS/GLONASS/Galileo/BDS/QZSS/SBAS/NavIC | RINEX4 `> EPH` typed NAV markers, source token normalization to `0`, RINEX3-style OBS body. |
| `+v4.01` | RINEX 4 typed NAV | GPS/GLONASS/Galileo/BDS/QZSS/SBAS/NavIC | Same CEQC structural writer as 4.00. Version-specific 4.01 system-data additions are documented as not yet fully covered. |
| `+v4.02` | RINEX 4 typed NAV | GPS/GLONASS/Galileo/BDS/QZSS/SBAS/NavIC | Same CEQC structural writer as 4.01, validated against gfzrnx for CEQC's RTCM3/UBX broadcast-ephemeris records. |

## Important limitations and open items

### RINEX 2.11 multi-GNSS extensions

Some ecosystems used RINEX 2.11 or extended RINEX2 conventions for Galileo/SBAS
and later for other constellations.  CEQC does **not** emit those extensions by
default.  It keeps RINEX2 output to the conservative GPS/GLONASS OBS profile and
GPS LNAV NAV profile because that is what old `teqc`/legacy RINEX2 paths can
handle most reliably.

### RINEX 3.00/3.01 BDS/QZSS exclusion

RINEX 3.00/3.01 output is intentionally restricted to GPS/GLONASS/Galileo.
BDS/QZSS/SBAS output begins at the CEQC 3.02 profile.  This prevents files such
as `+v3.00` from containing modern `Cxx/Jxx` records that validators may accept
syntactically but which do not match the early RINEX3 profile.

### RINEX 3.04/3.05 signal-code granularity

CEQC implements BDS/Galileo/GPS/GLONASS/QZSS observation codes that occur in the
current RTCM3 and UBX test data.  It does not yet contain a hard-coded complete
per-subversion allow-list for every signal code introduced in every RINEX 3.04
or 3.05 appendix table.  Unsupported or unvalidated signal-code combinations are
therefore documented here rather than claimed complete.  The regression gate is
`gfzrnx -meta full:txt` over generated files and per-source QC residual checks.


### Galileo week-number normalization

Galileo broadcast pages carry GST week numbering relative to the GPS week-1024
rollover epoch.  RINEX navigation files require the Galileo week field to be the
continuous GPS-aligned week.  CEQC therefore writes Galileo NAV week as
`broadcast_GST_week + 1024 + rollover_adjustment`, selected near the record
epoch.  This removes `gfzrnx` `wrong sysweek` warnings for UBX SFRBX-derived
Galileo ephemerides.

### RINEX 4.01/4.02 system data records

CEQC can preserve and write non-EPH RINEX4 records such as `STO`, `EOP`, and
`ION` when they are present in parsed input.  Current RTCM3/UBX translators are
focused on broadcast ephemerides and observations, so version-specific system
records that are not present in the source stream are not synthesized.

### RTKLIB/RTKPlot compatibility branch

`+rtklib` / `+rtkplot` is intentionally not a strict RINEX profile.  It mirrors
selected `RTKCONV-EX 2.5.0` formatting quirks only in that explicit branch.  The
mainline `+v2/+v3/+v4` writer keeps gfzrnx/Anubis-oriented strict output.

## Validation policy

For every supported sub-version CEQC tests:

1. The requested header version is written exactly.
2. Unsupported constellation families are filtered according to the profile gate.
3. OBS and NAV are parsed by `gfzrnx -meta full:txt` when records are present.
4. CEQC's own QC/residual path verifies the generated OBS/NAV pair for the most
   complete available profile, currently `+v4.02`.

Passing this matrix means the profile is **accepted and validated for the current
CEQC test data**.  It does not mean every optional record type from every RINEX
appendix has a dedicated source translator and sample yet.
