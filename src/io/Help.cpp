#include "ceqc/io/Help.hpp"
namespace ceqc::view {
std::string_view ceqcHelpText(){ return R"CEQC_HELP(Usage: ceqc [opts] file1 [file2 [...]]
   or: ceqc [opts] < stdin

CEQC is a clean-room C++ GNSS preprocessing and quality-checking tool with a
teqc-compatible command style.  It supports RINEX 2/3/4, RTCM3, u-blox UBX,
RINEX editing, splicing, decimation, QC reports, and teqc-style golden-diff
workflows.  Some historical teqc receiver-private translators are registered
but still report a clear "not implemented" error until their decoders are added.
Author: cyqsd
Email: root@cyqsd.cn

General options:
	-ver[sion] or +ver[sion]  write CEQC program version and build to stderr
	-id or +id                write CEQC program id to stderr
	-help or +help            output this CEQC on-line help to stderr
	+formats                  list registered input translators and implementation status
	+rtklib or +rtkplot       emit RTKLIB_EX/RTKPLOT-friendly OBS/NAV; teqc-style plus option only; use default mode by simply omitting it
	+err name                 write stderr directly to file 'name'
	++err name                append stderr directly to file 'name'
	+out name                 write stdout directly to file 'name'
	++out name                append stdout directly to file 'name'
	-config name              read file(s) 'name' as configuration file(s)
	+config                   output currently set CEQC parameters as configuration
	++config                  output a superset CEQC configuration template
	+bcf                      output BINEX-style configuration template
	-delim#                   change delimiter to # for separating file names (default = ,)
	+teqc                    render a teqc-style QC report
	-teqc_golden name         compare rendered report byte-for-byte with file 'name'
	-teqc_diff name           write first byte-diff context to file 'name'
	-teqc_eol lf|crlf         force LF or CRLF before golden comparison

Input and output:
	-tr fmt                   force input translator: rinex, rtcm3, ubx, nmea, binex, etc.
	+obs[file(s)] name        output OBS records in RINEX file 'name'
	+nav[file(s)] name        output NAV records in RINEX file 'name'
	+met[file(s)] name        output MET records in RINEX file 'name'
	+binex name               write BINEX 0x00 site metadata file 'name'
	+v2                       output RINEX 2.11 where possible
	+v3                       output RINEX 3.05 where possible
	+v4                       output RINEX 4.02 where possible
	+v2.10/+v2.11             request a versioned RINEX 2 output profile
	+v3.00..+v3.05            request a versioned RINEX 3 output profile
	+v4.00..+v4.02            request a versioned RINEX 4 output profile
	+nav_summary              output a compact NAV/SFRBX summary when available
	+nav_assembled            output only assembled UBX SFRBX navigation records when available

RINEX parser compatibility:
	+extend                   allow extended RINEX version 2
	-extend                   strict requirements on RINEX version 2 (default)
	+relax                    allow relaxed requirements on some RINEX fields
	-relax                    strict requirements on RINEX fields (default)
	+reformat                 allow reading of selected misformatted RINEX data fields
	-reformat                 strict requirements on RINEX data fields (default)
	+sv_dup[licates]          allow duplicate SVs in RINEX output
	-sv_dup[licates]          remove duplicate SVs before RINEX output (default)
	+svo                      order SVs by PRN or slot number
	-svo                      preserve SV detection order (default)

Time windowing and decimation:
	-st[art_window] str       set windowing start time to [[[[[[YY]YY]MM]DD]hh]mm]ss[.sssss]
	-e[nd_window] str         set windowing end time to [[[[[[YY]YY]MM]DD]hh]mm]ss[.sssss]
	-O.dec[imate] #[:#]       decimate OBS epochs by time interval and optional offset
	-N.dec[imate] #[:#]       decimate NAV epochs by time interval and optional offset
	-M.dec[imate] #[:#]       decimate MET epochs by time interval and optional offset

Satellite/system selection:
	+G / -G / +G<list> / -G<list>       include/exclude GPS PRNs
	+R / -R / +R<list> / -R<list>       include/exclude GLONASS slots
	+S / -S / +S<list> / -S<list>       include/exclude SBAS PRNs
	+E / -E / +E<list> / -E<list>       include/exclude Galileo PRNs
	+C / -C / +C<list> / -C<list>       include/exclude BeiDou PRNs
	+J / -J / +J<list> / -J<list>       include/exclude QZSS PRNs
	+I / -I / +I<list> / -I<list>       include/exclude IRNSS/NavIC PRNs
	-max_rx_ch[annels] #       set maximum receiver channels
	-max_rx_SVs #              set maximum SVs allowed per OBS epoch
	-n_GPS #                   set maximum expected GPS PRN
	-n_GLONASS #               set maximum expected GLONASS slot
	-n_SBAS #                  set maximum expected SBAS PRN-119
	-n_Galileo #               set maximum expected Galileo PRN
	-n_Beidou #                set maximum expected BeiDou PRN
	-n_QZSS #                  set maximum expected QZSS PRN-192
	-n_IRNSS #                 set maximum expected IRNSS/NavIC PRN
	-eph <list>                exclude RINEX4 NAV message types, e.g. G:CNAV,E:FNAV
	+eph <list>                include RINEX4 NAV message types, e.g. G:LNAV,E:INAV

Observation header editing (-O.*):
	-O.mo str                  MARKER NAME
	-O.mn str                  MARKER NUMBER
	-O.mt str                  MARKER TYPE
	-O.oi str                  OBSERVER / AGENCY
	-O.rn str                  receiver serial / number in REC # / TYPE / VERS
	-O.rt str                  receiver type in REC # / TYPE / VERS
	-O.rv str                  receiver version in REC # / TYPE / VERS
	-O.an str                  antenna serial / number in ANT # / TYPE
	-O.at str                  antenna type in ANT # / TYPE
	-O.px x y z                APPROX POSITION XYZ in meters
	-O.pg lat lon h            geodetic WGS84 latitude/longitude/height to XYZ
	-O.pe e n u                ANTENNA: DELTA H/E/N
	-O.sl[ant] dh s d          antenna slant height conversion to H/E/N
	-O.start yyyy mm dd hh mm ss   TIME OF FIRST OBS
	-O.obs list                keep OBS types in list
	-O.-obs list               remove OBS types in list
	-O.rename_obs list         rename OBS types in order
	-O.def_wf i j              add default WAVELENGTH FACT L1/2
	-O.mod_wf i j n SV...      add SV-specific WAVELENGTH FACT L1/2
	+O.comment str             append COMMENT
	+O.summary name|e          write observable summary to file or append to OBS output

Navigation header editing (-N.*):
	-N.ionA sys a0 a1 a2 a3    ionospheric alpha/correction terms
	-N.ionB sys b0 b1 b2 b3    ionospheric beta/correction terms
	-N.dUTC sys a0 a1 tot week TIME SYSTEM CORR
	-N.leap #                  LEAP SECONDS
	+N.comment str             append COMMENT to NAV header
	-N.r str                   run-by/date program text helper

Meteorological editing (-M.*):
	-M.obs list                keep MET observation types in list
	-M.-obs list               remove MET observation types in list
	-M.rename_obs list         rename MET observation types in order
	+M.model str               add MET model COMMENT
	+M.position str            add MET position COMMENT
	+M.comment str             append COMMENT to MET header

BINEX 0x00 metadata (-B.*):
	-B.marker str              site marker name
	-B.number str              site marker number
	-B.name str                site name
	-B.receiver str            receiver metadata
	-B.antenna str             antenna metadata
	-B.agency str              agency/operator metadata
	-B.comment str             metadata comment

Quality check modes:
	+qc                        teqc-style CEQC quality-check report
	+qcq                       compact CEQC diagnostic summary
	-qc / -qcq                 disable quality-check output
	-nav file1[,file2...]      use auxiliary RINEX NAV files for QC/residuals
	+nav file1[,file2...]      teqc-compatible QC alias for -nav when +qc/+qcq is used
	-no_orbit <systems>        do not use listed systems for orbit residual QC
	-no_pos <systems>          do not use listed systems for point-position QC
	+ap / -ap                  average position report on/off
	+pos / -pos                point-position report on/off
	+eep / +eepx / +eepg / +eepd  every-epoch position output modes
	+cl / -cl                  clock-slip candidates on/off
	+ion / -ion                ionospheric combination diagnostics on/off
	+iod / -iod                ionospheric delay-rate diagnostics on/off
	+lli / -lli                loss-of-lock indicator diagnostics on/off
	+mp / -mp                  multipath combination diagnostics on/off
	+mp_raw                    report raw multipath series where supported
	+pl / -pl                  pseudorange-minus-phase diagnostics on/off
	+sn/+snr / -sn/-snr        signal-to-noise diagnostics on/off
	+slips name                write slip/event candidates to file name
	++slips name               append slip/event candidates to file name
	-slips                     disable slip/event file output
	+plot / -plot              write compact plot file on/off
	-root str                  root name for plot/slip auxiliary files
	+sym / ++sym               print plot symbol legend / all symbols
	+w # or -w #               report width / compact plot width
	-set_horizon #             set QC horizon angle in degrees
	-set_mask #                set elevation mask angle in degrees
	-set_comp[arison] #        set comparison elevation angle in degrees
	-bins #                    set elevation/time histogram bins
	-sn_bins #                 set SNR bins
	-mp_bins #                 set multipath bins
	-gap_mn #                  report gaps longer than # minutes
	-gap_mx #                  no-NAV gap threshold in minutes
	-min_SVs #                 minimum SVs per epoch for position QC
	-min_L1 # ... -min_L8 #    minimum RINEX S/N code for frequency band
	-mp_sigmas #               multipath outlier sigma threshold
	-mp_win #                  multipath moving-average window
	-mp12_rms #, -mp21_rms #   multipath RMS thresholds in cm
	-ion_jump #                ionospheric jump threshold in cm
	-iod_jump #                IOD jump threshold in cm/min
	-msec_tol #                millisecond clock-slip tolerance
	-code_sigmas #             code sigma threshold
	+rs / -rs                  rise/set report on/off
	+ssv                       show per-SV summaries
	+svpr                      show SV pseudorange summaries
	+Y-code / -Y-code          Y-code report control
	+ceqc_ext / -ceqc_ext      CEQC extended QC details on/off

Examples:
	ceqc +verify file.23o
	ceqc +qc file.23o
	ceqc -tr rtcm3 +v3 +obs out.26o +nav out.26p raw.rtcm3
	ceqc -tr ubx +v2 +G +obs gps.26o raw.ubx
	ceqc +qc +teqc -teqc_golden teqc.qc -teqc_diff diff.txt file.26o

Implementation note:
	CEQC intentionally exposes teqc-style switches, but this help describes CEQC's
	current behavior.  Some receiver-private formats are registered but not yet
	fully decoded; RINEX, RTCM3 and u-blox UBX RAWX/SFRBX are the primary implemented
	paths in this build.
)CEQC_HELP"; }
} // namespace ceqc::view
