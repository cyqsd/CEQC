# ceqc architecture

## Why the layout changed

Earlier C++ iterations used an MVC-style layout (`model`, `view`, `controller`) because that matched the first migration plan.  For a GNSS tool that must support RINEX, RTCM3, UBX, BINEX, QC, editing and synthesis, a domain-oriented layout is easier to maintain than forcing all code into MVC buckets.

## Current module boundaries

- `apps/ceqc`: executable only.  It does not own parsing, QC or format logic.
- `include/ceqc/core`: small shared models such as command options and time helpers.
- `include/ceqc/rinex` and `src/rinex`: RINEX data structures, parser, writer, edits, decimation, windowing and synthesis.
- `include/ceqc/translate` and `src/translate`: native receiver or stream translators.  Implemented: RINEX, RTCM3, UBX.  Registered placeholders exist for remaining teqc receiver families.
- `include/ceqc/qc` and `src/qc`: QC statistics, derived observables, timeplot, slips, residual screening.
- `include/ceqc/cli` and `src/cli`: teqc-style command line and config-file parsing.
- `include/ceqc/app` and `src/app`: high-level workflow that wires translators, RINEX services, QC and output.
- `include/ceqc/io` and `src/io`: help text, issue printing and report rendering.

## Compatibility strategy

1. Preserve teqc-style option spellings in `src/cli/CommandLine.cpp`.
2. Keep options in `ceqc::model::Operation` and `QCOptions`, even when the deep algorithm still needs further work.
3. Add executable behavior for every parsed option wherever possible; unsupported receiver translators fail clearly rather than silently producing false data.
4. Use `+teqc -teqc_golden` for byte-for-byte report comparison against real teqc output.
5. Keep regression tests in `tests/` focused on behavior rather than implementation details.
