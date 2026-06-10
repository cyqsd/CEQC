# CEQC 0.0.1 RTKLIB_EX / RTKPLOT compatibility notes

The strict CEQC RINEX path is aimed at gfzrnx and Anubis.  In strict RINEX 4
mode it writes typed NAV markers such as `> EPH G01 LNAV 0`.  RTKLIB_EX 2.5.0 /
RTKPLOT should use CEQC's explicit compatibility projection instead:

```bash
ceqc -tr rtcm3 +rtklib +v4.02 +obs sample.obs +nav sample.nav sample.rtcm3
```

or:

```bash
tools/rtklib_plot_export.sh sample.rtcm3 sample rtcm3
```

## Verified behavior

The following points were checked against the current implementation:

- In `+rtklib` / `+rtkplot` mode, an explicit `+v3.05` or `+v4.02` is
  preserved.  If no `+v` is supplied, CEQC defaults the compatibility export to
  3.05.
- The compatibility branch keeps a legacy mixed NAV body and suppresses RINEX 4
  typed `> EPH` markers.
- NAV numeric fields use RTKCONV-style `D` exponents.
- OBS output uses RTKLIB signal aliases where needed, for example GPS
  `C2L/L2L/S2L` -> `C2W/L2W/S2W` and BDS `C7X/L7X/S7X` -> `C7I/L7I/S7I`.
- OBS TYPES are ordered as signal triplets (`C/L/S`) per frequency, matching
  RTKCONV-EX style.
- LLI/SSI columns are blank in the compatibility branch to match RTKCONV-EX /
  RTKPLOT import behavior.
- `+rtklib` / `+rtkplot` writes CRLF for OBS and NAV.  `writeFile()` opens the
  file in binary mode so Windows/MinGW text-mode translation does not turn CRLF
  into `\r\r\n`.
- For RINEX 3/4 compatibility output, CEQC writes an RTKCONV-EX-like
  `PGM / RUN BY / DATE` line and inserts immediate `COMMENT` records:
  `format: <input-format>` and `log: <input-path>`.
- For OBS compatibility output, CEQC normalizes stream-style station metadata to
  RTKCONV-like defaults: `MARKER NAME` = `0000`, blank `MARKER NUMBER`, blank
  `MARKER TYPE`, and blank `OBSERVER / AGENCY`.

## Separation From Strict Output

This mode is intentionally separate from the mainline RINEX writer.  Without
`+rtklib` / `+rtkplot`, CEQC keeps the strict gfzrnx/Anubis-oriented RINEX2/3/4
output unchanged:

- strict RINEX 4 still emits typed NAV markers such as `> EPH G01 LNAV 0`;
- strict RINEX 3/4 still use LF line endings;
- RINEX 2 continues using CRLF for the Anubis compatibility reason documented
  elsewhere.
