# ceqc compatibility status

Implemented core translators: RINEX, RTCM3, UBX RAWX/SFRBX.

CEQC supports both a teqc-compatible legacy QC subset and modern GNSS QC extensions.  The report now explicitly separates those two scopes so old teqc comparisons do not incorrectly penalize modern BDS-3/Galileo/QZSS/RINEX3/4 features.

Comparable with teqc 2019 only for common legacy RINEX 2 GPS/GLONASS L1/L2 items: epoch/window statistics, duplicate-SV checks, legacy MP12/MP21, S1/S2 summaries, and basic missing OBS/NAV lists.  Modern-only CEQC metrics, including BDS-3 multi-frequency MP combinations, RINEX3/4 typed NAV, UBX SFRBX page inventory, and S5/S6/S7/S8 summaries, must be analyzed independently.

Not yet byte-identical without golden tuning: teqc compact3 exact layout, every historical warning line, full Galileo/BDS SFRBX ICD field decoders, and exact BINEX dialect variations used by specific teqc versions.
