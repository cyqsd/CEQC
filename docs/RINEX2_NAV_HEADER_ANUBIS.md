# RINEX2 NAV header compatibility for Anubis

CEQC 0.0.1 writes RINEX 2 navigation files with an explicit legacy GPS
navigation type on the first header line:

```text
     2.11           N: GPS NAV DATA                         RINEX VERSION / TYPE
```

CEQC's RINEX2 NAV output contains only GPS LNAV records because RINEX2 mixed
navigation support is not used for modern multi-GNSS data.  Writing the generic
`NAVIGATION DATA` text caused Anubis 3.11 to warn:

```text
warning - RINEXN system not defined, used GPS
```

The fix is based on the actual RINEX2 output scope, not on any sample file name,
station coordinate, receiver model, or satellite list.
