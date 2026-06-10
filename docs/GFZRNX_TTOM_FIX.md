# gfzrnx TtoM warning fix in ceqc 0.0.1

gfzrnx can warn on GPS/QZSS LNAV records when the RINEX NAV `TransmissionTime` / `TtoM` field is written as zero:

```text
correcting zero nav record 25/TtoM
correcting zero TtoM nav. record to 0.999900000000e+09
```

In RTCM 1019/1044 streams the RINEX TtoM value is not always available as a complete decoded field. CEQC now treats a decoded zero/non-finite/out-of-range `TransmissionTime` as missing and derives a replacement dynamically from the same navigation record:

1. Prefer decoded `TransmissionTimeRaw` if present and valid.
2. Prefer `Toe` if present and valid.
3. Prefer `Toc` if present and valid.
4. Fall back to the record epoch seconds-of-week in the corresponding system time scale.

No sample-specific values, satellite IDs, file names, receiver IDs, or coordinates are used. The fallback is based only on the current navigation record fields and epoch.

Regression check for the supplied RTCM3 sample:

```text
nav.26p : G/J LNAV records = 24, zero TtoM = 0
nav4.26p: G/J LNAV records = 24, zero TtoM = 0
```
