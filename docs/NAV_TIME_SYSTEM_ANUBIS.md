# NAV time-system correction for Anubis usable GLO/BDS ephemerides

CEQC 0.0.1 fixes RTCM3 -> RINEX NAV time tags so that Anubis counts GLONASS and BeiDou ephemerides as usable rather than only listing them in the navigation inventory.

## Root cause

RTCM observation epochs are commonly represented in a receiver/GPS-aligned timeline, but RINEX NAV ephemeris epochs remain constellation-specific:

- GPS/QZSS/Galileo/BDS ephemeris first-line epoch is the clock epoch (Toc) in that constellation time scale.
- BeiDou RTCM 1042 Toe/Toc are BDT seconds and must not be shifted by +14 s when written as BeiDou RINEX NAV.
- GLONASS RTCM 1020 time tags are GLONASS/UTC(SU)-based; RINEX GLONASS NAV is UTC-based, so CEQC must remove the +3 h GLONASS offset but must not add leap seconds to GPST.

The previous output converted BDS NAV epochs by +14 s and GLONASS NAV epochs by +18 s. Anubis could parse the records, but it did not count them as usable navigation for QC.

## Implemented behavior

- GPS/QZSS/Galileo/BDS NAV line epoch uses decoded Toc where available.
- BeiDou 1042 writes BDT Toc directly, without GPST +14 s.
- GLONASS 1020 writes UTC ephemeris epoch `tb*900 - 3h` and UTC message frame time `tk - 3h`, both normalized to the local day.
- No station, satellite, coordinate, receiver ID, or sample filename is hard-coded.

## Anubis 3.11 regression evidence

Using the same CEQC sample data set, Anubis `navi` changed from inventory-only to usable NAV:

| Scope | GPS have | GLO have | BDS have | Result |
|---|---:|---:|---:|---|
| RINEX2 legacy | 12 | 0 | n/a | expected: RINEX2 GPS NAV only |
| RINEX3 mixed | 12 | 6 | 16 | GLO/BDS usable |
| RINEX4 typed | 12 | 15 | 16 | GLO/BDS usable |

The only remaining Anubis warning in these validation logs is marker name vs. test filename, which is not a NAV parse or usability error.
