# UBX Anubis/gfzrnx validation notes for ceqc 0.0.1

The 2026-01-24 UBX RAWX/SFRBX sample exposes two different classes of messages:

1. Real CEQC output problems fixed in 0.0.1 validation helpers:
   - The generated OBS now carries a non-zero `APPROX POSITION XYZ` when GPS/QZSS LNAV SFRBX records are available for an SPP estimate.
   - The `tools/anubis_compare.sh` helper derives the Anubis `:gen:sys` list from the actual OBS records, so absent constellations are not requested by the test command.
   - The `tools/gfzrnx_compare.sh` helper skips a NAV file only if it truly contains no ephemeris records; this avoids treating a deliberately empty NAV header as a RINEX syntax failure.

2. Anubis 3.11 behavior that should not be worked around by fabricating data:
   - For RINEX3/4 mixed OBS files, Anubis may warn `GLO SLOT/FREQ not available` and `GLO BIASES not available` even when the OBS file contains no `R` records.
   - CEQC does not write a dummy `GLONASS SLOT / FRQ #` line for no-GLONASS UBX data. A previous attempt to write a zero-count slot line caused Anubis to report `GLONASS SLOT / FRQ zero obs codes`.
   - CEQC only writes GLONASS slot/frequency metadata when real GLONASS ephemerides or observations provide it.

This keeps the main RINEX output truthful and avoids hard-coded or fabricated GLONASS slot/channel values.
