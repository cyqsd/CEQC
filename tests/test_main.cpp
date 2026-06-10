#include "ceqc/rinex/RinexService.hpp"
#include "ceqc/qc/QCService.hpp"
#include "ceqc/cli/CommandLine.hpp"
#include <cassert>
#include <iostream>
#include <sstream>
#include <set>
int main(){
  auto rf = ceqc::service::rinex::readFile("testdata/minimal_v3.23o");
  assert(rf.header.kind == ceqc::model::RinexKind::Obs);
  assert(rf.data.observationEpochs.size() == 2);
  assert(rf.data.observationRecords.size() == 3);
  auto q = ceqc::service::qc::analyze(rf);
  assert(q.observationRecords == 3);

  ceqc::model::QCOptions qoff;
  qoff.ion=false; qoff.iod=false; qoff.multipath=false; qoff.snr=false; qoff.lli=false; qoff.pseudorangePhase=false; qoff.clockSlips=false; qoff.width=20;
  auto q2 = ceqc::service::qc::analyze(rf, qoff);
  assert(q2.derived);
  assert(q2.derived->optionsActive.empty());
  assert(q2.derived->snrStats.empty());
  assert(q2.derived->multipathStats.empty());
  assert(q2.derived->ionStats.empty());
  assert(q2.derived->iodStats.empty());
  assert(q2.derived->pseudorangePhase.empty());
  assert(q2.derived->timeplot.size() == 22); // bars plus two delimiters
  ceqc::model::QCOptions qdeep; qdeep.svpr=true; qdeep.yCode=true; qdeep.dataIndicators=true; qdeep.riseSet=true;
  auto qdeepSum = ceqc::service::qc::analyze(rf, qdeep);
  assert(qdeepSum.derived);
  assert(!qdeepSum.derived->riseSetEvents.empty());
  assert(qdeepSum.derived->dataCompleteness.completeRecords + qdeepSum.derived->dataCompleteness.partialRecords == qdeepSum.observationRecords);
  assert(!qdeepSum.derived->svPseudorangeStats.empty());
  assert(qdeepSum.derived->yCodeEnabled);

  auto qop = ceqc::cli::parseArgs({"+qcq", "-ion", "-iod", "-mp", "-sn", "-lli", "-cl", "+w", "20", "-no_orbit", "GPS+BDS", "-no_pos", "GAL+QZS", "testdata/minimal_v3.23o"});
  assert(qop.qc && qop.quietQC);
  assert(!qop.qcOptions.ion && !qop.qcOptions.iod && !qop.qcOptions.multipath && !qop.qcOptions.snr && !qop.qcOptions.lli && !qop.qcOptions.clockSlips);
  assert(qop.qcOptions.width == 20);
  assert(qop.qcOptions.noOrbitSystems["G"] && qop.qcOptions.noOrbitSystems["C"]);
  assert(qop.qcOptions.noPositionSystems["E"] && qop.qcOptions.noPositionSystems["J"]);

  // Structured teqc-style OBS edits should preserve REC/ANT columns and observable filtering.
  ceqc::service::rinex::applyHeaderEdits(rf, {{"-O.rn","RX001"},{"-O.rt","UM960"},{"-O.rv","1.0"},{"-O.an","ANT001"},{"-O.at","ADVNULL"},{"-O.px","1 2 3"}});
  bool recOk=false, antOk=false, posOk=false;
  for(auto& h:rf.header.lines){
    if(h.label=="REC # / TYPE / VERS" && h.raw.find("RX001")!=std::string::npos && h.raw.find("UM960")!=std::string::npos) recOk=true;
    if(h.label=="ANT # / TYPE" && h.raw.find("ANT001")!=std::string::npos && h.raw.find("ADVNULL")!=std::string::npos) antOk=true;
    if(h.label=="APPROX POSITION XYZ" && h.raw.find("1.0000")!=std::string::npos) posOk=true;
  }
  assert(recOk && antOk && posOk);
  ceqc::service::rinex::applyObsTypeFilter(rf, "C1C L1C", false, false);
  for(auto& r:rf.data.observationRecords) assert(r.values.size()==2);
  auto dspec = ceqc::model::DecimationSpec{true, std::chrono::seconds(30), std::chrono::seconds(0), "30"};
  auto decimated = ceqc::service::rinex::decimate(rf, dspec);
  assert(decimated.data.observationEpochs.size()==2);

  auto merged = ceqc::service::rinex::merge({rf}, ceqc::model::RinexKind::Obs, 3.05);
  assert(!merged.body.empty());
  ceqc::model::NavigationRecord nr;
  nr.system="G"; nr.satellite="G01"; nr.recordType="EPH"; nr.messageType="LNAV"; nr.messageSubtype="TEST";
  nr.epoch=ceqc::model::makeUTC(2026,6,5,13,0,0);
  for(int i=0;i<10;++i) nr.values.push_back(0.001*i);
  nr.fields["SqrtA"]={"SqrtA","sqrt(m)",5153.7,0};
  ceqc::model::RinexFile nav; nav.header.kind=ceqc::model::RinexKind::Nav; nav.header.version=4.02; nav.data.navigationRecords.push_back(nr);
  auto cfgop = ceqc::cli::parseArgs({"+formats"});
  assert(cfgop.showFormats);
  auto optop = ceqc::cli::parseArgs({"+qc", "+mask", "-set_comp", "25", "-teqc_eol", "crlf", "testdata/minimal_v3.23o"});
  assert(optop.qc && optop.qcOptions.mask && optop.qcOptions.setComparisonDeg == 25 && optop.teqcEOL == "crlf");


  // Additional teqc-style metadata edits: geodetic conversion, slant eccentricity,
  // wavelength factors, formatted TIME OF FIRST OBS and NAV ion/time-system records.
  auto rf2 = ceqc::service::rinex::readFile("testdata/minimal_v3.23o");
  ceqc::service::rinex::applyHeaderEdits(rf2, {{"-O.pg","0 0 0"},{"-O.sl","5 6 1"},{"-O.start","2026 6 5 13 0 30"},{"-O.def_wf","1 1"},{"-O.mod_wf","1 0 2 G01 G02"}});
  bool pgOk=false, slOk=false, stOk=false, wfOk=false, mwfOk=false;
  for(auto& h: rf2.header.lines){
    if(h.label=="APPROX POSITION XYZ" && h.raw.find("6378137.0000")!=std::string::npos) pgOk=true;
    if(h.label=="ANTENNA: DELTA H/E/N" && h.raw.find("5.0000")!=std::string::npos) slOk=true;
    if(h.label=="TIME OF FIRST OBS" && h.raw.find("2026")!=std::string::npos && h.raw.find("30.0000000")!=std::string::npos) stOk=true;
    if(h.label=="WAVELENGTH FACT L1/2" && h.raw.find("     1     1")!=std::string::npos) wfOk=true;
    if(h.label=="WAVELENGTH FACT L1/2" && h.raw.find("G01")!=std::string::npos && h.raw.find("G02")!=std::string::npos) mwfOk=true;
  }
  assert(pgOk && slOk && stOk && wfOk && mwfOk);

  ceqc::model::RinexFile navEdit; navEdit.header.kind=ceqc::model::RinexKind::Nav; navEdit.header.version=4.02;
  navEdit.header.lines.push_back({"     4.02           NAVIGATION DATA     M                   RINEX VERSION / TYPE","     4.02           NAVIGATION DATA     M","RINEX VERSION / TYPE"});
  navEdit.header.lines.push_back({"                                                            END OF HEADER","","END OF HEADER"});
  navEdit.header.byLabel["RINEX VERSION / TYPE"].push_back(0); navEdit.header.byLabel["END OF HEADER"].push_back(1);
  ceqc::service::rinex::applyHeaderEdits(navEdit, {{"-N.ionA","G 1 2 3 4"},{"-N.ionB","C 5 6 7 8"},{"-N.dUTC","G 1e-9 2e-15 345600 2420"}});
  bool ionOk=false, tsOk=false;
  for(auto& h: navEdit.header.lines){ if(h.label=="IONOSPHERIC CORR" && (h.raw.find("GPSA")!=std::string::npos || h.raw.find("BDSB")!=std::string::npos)) ionOk=true; if(h.label=="TIME SYSTEM CORR" && h.raw.find("GPUT")!=std::string::npos) tsOk=true; }
  assert(ionOk && tsOk);

  auto nmerged=ceqc::service::rinex::merge({nav}, ceqc::model::RinexKind::Nav, 4.02);
  assert(!nmerged.body.empty());
  bool hasTyped=false, hasNumeric=false;
  for(auto& l:nmerged.body){ if(l.find("> EPH G01 LNAV")!=std::string::npos) hasTyped=true; if(l.find("1.000000000000E-03")!=std::string::npos) hasNumeric=true; }
  assert(hasTyped && hasNumeric);

  auto sysop = ceqc::cli::parseArgs({"+G1,2", "-R", "-max_rx_SVs", "2", "+svo", "testdata/minimal_v3.23o"});
  assert(sysop.svSelections["G"].specified && sysop.svSelections["G"].includeMode && sysop.svSelections["G"].prns.size()==2);
  assert(sysop.svSelections["R"].specified && !sysop.svSelections["R"].includeMode && sysop.svSelections["R"].all);
  assert(sysop.maxRxSVs == 2 && sysop.orderSVByPRN);
  auto compactTime = ceqc::cli::parseArgs({"-st", "260605130030", "-e", "20260605130045", "testdata/minimal_v3.23o"});
  assert(compactTime.windowStart && compactTime.windowEnd);
  auto cfgOut = ceqc::cli::parseArgs({"+config"});
  assert(cfgOut.showConfig);
  auto bcfOut = ceqc::cli::parseArgs({"+bcf"});
  assert(bcfOut.showBCF);



  // RINEX 2 mixed-system OBS: satellite symbols in the epoch list must be
  // respected instead of defaulting every numeric/system token to GPS.
  std::string v2mix;
  v2mix += "     2.11           OBSERVATION DATA    M                   RINEX VERSION / TYPE\n";
  v2mix += "CEQC                TEST                                    PGM / RUN BY / DATE\n";
  v2mix += "     3    C1    L1    S1                                    # / TYPES OF OBSERV\n";
  v2mix += "                                                            END OF HEADER\n";
  v2mix += " 26  1 24  4  0  4.0000000  0  4G01R02E11C03\n";
  v2mix += "  22500000.000    118000.000        45.000  \n";
  v2mix += "  22600000.000    119000.000        44.000  \n";
  v2mix += "  22700000.000    120000.000        43.000  \n";
  v2mix += "  22800000.000    121000.000        42.000  \n";
  std::istringstream v2in(v2mix);
  auto v2rf = ceqc::service::rinex::readStream("v2mix.26o", v2in);
  assert(v2rf.header.version == 2.11);
  assert(v2rf.data.observationRecords.size() == 4);
  std::set<std::string> sys2;
  for (auto& r : v2rf.data.observationRecords) sys2.insert(r.system);
  assert(sys2.count("G") && sys2.count("R") && sys2.count("E") && sys2.count("C"));

  // RINEX 3.03/3.04/3.05 continuation lines for SYS / # / OBS TYPES.
  std::string v3cont;
  v3cont += "     3.03           OBSERVATION DATA    G                   RINEX VERSION / TYPE\n";
  v3cont += "G   16 C1C L1C D1C S1C C2W L2W D2W S2W C5X L5X D5X S5X C1W SYS / # / OBS TYPES\n";
  v3cont += "      L1W D1W S1W                                      SYS / # / OBS TYPES\n";
  v3cont += "                                                            END OF HEADER\n";
  v3cont += "> 2026 01 24 04 00  4.0000000  0   1\n";
  v3cont += "G01  1.000  2.000  3.000  4.000  5.000  6.000  7.000  8.000  9.000 10.000 11.000 12.000 13.000\n";
  v3cont += "     14.000 15.000 16.000\n";
  std::istringstream v3in(v3cont);
  auto v3rf = ceqc::service::rinex::readStream("v3cont.26o", v3in);
  assert(v3rf.header.version == 3.03);
  assert(v3rf.header.obsTypes["G"].size() == 16);
  assert(v3rf.data.observationRecords.size() == 1);
  assert(v3rf.data.observationRecords[0].values.size() == 16);


  // OBS INTERVAL must be derived from actual output epochs, not a fixed value
  // and not a stale copied header.  This mirrors teqc's nominal interval logic.
  std::string vint;
  vint += "     3.05           OBSERVATION DATA    G                   RINEX VERSION / TYPE\n";
  vint += "    30.000                                                  INTERVAL\n";
  vint += "G    1 C1C                                                SYS / # / OBS TYPES\n";
  vint += "                                                            END OF HEADER\n";
  vint += "> 2026 01 24 04 00  0.0000000  0   1\n";
  vint += "G01  1.000\n";
  vint += "> 2026 01 24 04 00  5.0000000  0   1\n";
  vint += "G01  2.000\n";
  vint += "> 2026 01 24 04 00 10.0000000  0   1\n";
  vint += "G01  3.000\n";
  std::istringstream vintin(vint);
  auto vintrf = ceqc::service::rinex::readStream("vint.26o", vintin);
  auto vintmerged = ceqc::service::rinex::merge({vintrf}, ceqc::model::RinexKind::Obs, 3.05);
  int intervalLines = 0;
  bool intervalFive = false;
  for (const auto& h : vintmerged.header.lines) {
    if (h.label == "INTERVAL") {
      ++intervalLines;
      intervalFive = intervalFive || h.raw.find("5.000") != std::string::npos;
    }
  }
  assert(intervalLines == 1 && intervalFive);

  // RINEX 4.00 typed NAV records must keep message type/subtype and numeric payload.
  std::string v4nav;
  v4nav += "     4.00           NAVIGATION DATA     M                   RINEX VERSION / TYPE\n";
  v4nav += "                                                            END OF HEADER\n";
  v4nav += "> EPH G01 LNAV TEST\n";
  v4nav += "G01 2026 01 24 04 00 00 1.0D-04 2.0D-12 0.0D+00\n";
  v4nav += "    1.0D+00 2.0D+00 3.0D+00 4.0D+00\n";
  std::istringstream v4in(v4nav);
  auto v4rf = ceqc::service::rinex::readStream("v4nav.26p", v4in);
  assert(v4rf.header.version == 4.00);
  assert(v4rf.data.navigationRecords.size() == 1);
  assert(v4rf.data.navigationRecords[0].messageType == "LNAV");
  assert(v4rf.data.navigationRecords[0].messageSubtype == "TEST");

  auto vop = ceqc::cli::parseArgs({"+v4.00", "+qcq", "testdata/minimal_v3.23o"});
  assert(vop.targetVersion > 4.0 && vop.targetVersion < 4.001);



  // RINEX 4 typed Galileo and BeiDou NAV records should use constellation-specific
  // field names, not the GPS LNAV field order.  This is essential for RINEX 4.01/4.02
  // mixed NAV QC and residual screening.
  std::string v4ec;
  v4ec += "     4.02           NAVIGATION DATA     M                   RINEX VERSION / TYPE\n";
  v4ec += "                                                            END OF HEADER\n";
  v4ec += "> EPH E11 INAV TEST\n";
  v4ec += "E11 2026 01 24 04 00 00 1.0D-04 2.0D-12 0.0D+00\n";
  v4ec += "    1.0D+00 2.0D+00 3.0D+00 4.0D+00 5.0D+00 6.0D+00 7.0D+00 8.0D+00\n";
  v4ec += "> EPH C03 D1 TEST\n";
  v4ec += "C03 2026 01 24 04 00 00 1.0D-04 2.0D-12 0.0D+00\n";
  v4ec += "    1.0D+00 2.0D+00 3.0D+00 4.0D+00 5.0D+00 6.0D+00 7.0D+00 8.0D+00\n";
  std::istringstream v4ecin(v4ec);
  auto v4ecrf = ceqc::service::rinex::readStream("v4ec.26p", v4ecin);
  assert(v4ecrf.data.navigationRecords.size() == 2);
  assert(v4ecrf.data.navigationRecords[0].system == "E");
  assert(v4ecrf.data.navigationRecords[0].fields.count("IODnav"));
  assert(v4ecrf.data.navigationRecords[0].fields.count("Eccentricity"));
  assert(v4ecrf.data.navigationRecords[1].system == "C");
  assert(v4ecrf.data.navigationRecords[1].fields.count("AODE"));
  assert(v4ecrf.data.navigationRecords[1].fields.count("SqrtA"));

  // Non-EPH RINEX 4 system data records should be grouped instead of being
  // confused with satellite ephemeris records.
  std::string v4sto;
  v4sto += "     4.02           NAVIGATION DATA     M                   RINEX VERSION / TYPE\n";
  v4sto += "                                                            END OF HEADER\n";
  v4sto += "> STO GPUT TEST\n";
  v4sto += "    1.0D-09 2.0D-15 345600 2420\n";
  std::istringstream v4stoin(v4sto);
  auto v4storf = ceqc::service::rinex::readStream("v4sto.26p", v4stoin);
  assert(v4storf.data.navigationRecords.size() == 1);
  assert(v4storf.data.navigationRecords[0].recordType == "STO");
  assert(v4storf.data.navigationRecords[0].messageType == "GPUT");


  // RINEX 4 system data records should round-trip through NAV synthesis as typed
  // STO/EOP/ION records instead of being emitted as raw comment-like lines.
  ceqc::model::RinexFile sysNav; sysNav.header.kind=ceqc::model::RinexKind::Nav; sysNav.header.version=4.02;
  ceqc::model::NavigationRecord sto; sto.recordType="STO"; sto.messageType="GPUT"; sto.messageSubtype="TEST"; sto.values={1e-9,2e-15,345600,2420};
  ceqc::model::NavigationRecord eop; eop.recordType="EOP"; eop.messageType="IERS"; eop.messageSubtype="TEST"; eop.values={0.1,0.2,0.3,0.4,0.5,0.6,345600,2420};
  ceqc::model::NavigationRecord ion; ion.recordType="ION"; ion.messageType="GPSA"; ion.messageSubtype="KLOB"; ion.values={1,2,3,4,5,6,7,8};
  sysNav.data.navigationRecords={sto,eop,ion};
  auto sysMerged=ceqc::service::rinex::merge({sysNav}, ceqc::model::RinexKind::Nav, 4.02);
  bool hasSto=false, hasEop=false, hasIon=false;
  for(auto& l: sysMerged.body){ if(l.find("> STO GPUT TEST")!=std::string::npos) hasSto=true; if(l.find("> EOP IERS TEST")!=std::string::npos) hasEop=true; if(l.find("> ION GPSA KLOB")!=std::string::npos) hasIon=true; }
  assert(hasSto && hasEop && hasIon);

  // RINEX 4.02 GLONASS CDMA and BDS CNAV-family message types should not use
  // legacy GPS field labels.
  std::string v4sub;
  v4sub += "     4.02           NAVIGATION DATA     M                   RINEX VERSION / TYPE\n";
  v4sub += "                                                            END OF HEADER\n";
  v4sub += "> EPH R01 L1OC TEST\n";
  v4sub += "R01 2026 01 24 04 00 00 1.0D-04 2.0D-12 0.0D+00\n";
  v4sub += "    1.0D+00 2.0D+00 3.0D+00 4.0D+00 5.0D+00 6.0D+00 7.0D+00 8.0D+00\n";
  v4sub += "> EPH C03 CNAV TEST\n";
  v4sub += "C03 2026 01 24 04 00 00 1.0D-04 2.0D-12 0.0D+00\n";
  v4sub += "    1.0D+00 2.0D+00 3.0D+00 4.0D+00 5.0D+00 6.0D+00 7.0D+00 8.0D+00\n";
  v4sub += "    9.0D+00 1.0D+01 1.1D+01 1.2D+01\n";
  v4sub += "    1.3D+01 1.4D+01 1.5D+01 1.6D+01\n";
  v4sub += "    1.7D+01 1.8D+01 1.9D+01 2.0D+01\n";
  std::istringstream v4subin(v4sub);
  auto v4subrf=ceqc::service::rinex::readStream("v4sub.26p", v4subin);
  assert(v4subrf.data.navigationRecords.size()==2);
  assert(v4subrf.data.navigationRecords[0].fields.count("FrequencyNumber"));
  assert(v4subrf.data.navigationRecords[1].fields.count("ADot"));
  assert(v4subrf.data.navigationRecords[1].fields.count("ISC_B1Cp"));



  // RINEX 2 OBS must not contain RINEX 3/4-only GLONASS header records.
  // This regressed when Anubis-oriented GLONASS SLOT/BIS records were copied
  // into all OBS versions; teqc rejects such RINEX2 headers as malformed.
  ceqc::model::RinexFile v2obsStrict;
  v2obsStrict.header.kind = ceqc::model::RinexKind::Obs;
  v2obsStrict.header.version = 3.05;
  v2obsStrict.header.lines.push_back({" C1C    0.000 C1P    0.000 C2C    0.000 C2P    0.000        GLONASS COD/PHS/BIS"," C1C    0.000 C1P    0.000 C2C    0.000 C2P    0.000","GLONASS COD/PHS/BIS"});
  ceqc::model::ObservationRecord ro;
  ro.system = "R"; ro.satellite = "R05"; ro.time = ceqc::model::makeUTC(2026,6,5,13,0,30);
  ceqc::model::ObservationValue rov; rov.type = "C1C"; rov.value = 22000000.0; ro.values.push_back(rov);
  v2obsStrict.data.observationRecords.push_back(ro);
  auto v2strictMerged = ceqc::service::rinex::merge({v2obsStrict}, ceqc::model::RinexKind::Obs, 2.11);
  for (const auto& h : v2strictMerged.header.lines) {
    assert(h.label != "GLONASS COD/PHS/BIS");
    assert(h.label != "GLONASS SLOT / FRQ #");
  }

  // RINEX 4 EPH typed NAV markers written from RTCM/UBX translators should not
  // expose internal source tags such as RTCM1042.  gfzrnx treats unknown source
  // tokens as XXXX and rejects subsequent rows.  Use neutral source "0".
  ceqc::model::RinexFile srcNav4; srcNav4.header.kind = ceqc::model::RinexKind::Nav; srcNav4.header.version = 4.02;
  ceqc::model::NavigationRecord cnav; cnav.system="C"; cnav.satellite="C01"; cnav.recordType="EPH"; cnav.messageType="D1"; cnav.messageSubtype="RTCM1042"; cnav.epoch=ceqc::model::makeUTC(2026,6,5,13,0,14);
  cnav.values = {1e-4,2e-12,0,1,2,3,4,5};
  srcNav4.data.navigationRecords.push_back(cnav);
  auto mergedNav4 = ceqc::service::rinex::merge({srcNav4}, ceqc::model::RinexKind::Nav, 4.02);
  bool hasNeutralSource=false, hasTranslatorSource=false;
  for (const auto& l : mergedNav4.body) {
    if (l.find("> EPH C01 D1 0") != std::string::npos) hasNeutralSource=true;
    if (l.find("RTCM1042") != std::string::npos) hasTranslatorSource=true;
  }
  assert(hasNeutralSource && !hasTranslatorSource);


  // RINEX2 OBS/NAV body lines must be padded to at least 80 columns.
  // Anubis 3.11 indexes RINEX2 records by fixed columns and aborts on 79-char
  // LF-only NAV/OBS lines. This test prevents regressions without hardcoding
  // any sample station, satellite, or coordinate.
  ceqc::model::RinexFile v2padObs;
  v2padObs.header.kind = ceqc::model::RinexKind::Obs;
  v2padObs.header.version = 3.05;
  ceqc::model::ObservationRecord po;
  po.system="G"; po.satellite="G01"; po.time=ceqc::model::makeUTC(2026,1,1,0,0,0);
  for (const char* t : {"C1C","L1C","D1C","S1C","C2L","L2L"}) {
    ceqc::model::ObservationValue ov; ov.type=t; ov.value=1.0; po.values.push_back(ov);
  }
  v2padObs.data.observationRecords.push_back(po);
  auto v2padMerged = ceqc::service::rinex::merge({v2padObs}, ceqc::model::RinexKind::Obs, 2.11);
  for (const auto& l : v2padMerged.body) assert(l.size() >= 80);

  ceqc::model::RinexFile v2padNav; v2padNav.header.kind=ceqc::model::RinexKind::Nav; v2padNav.header.version=3.05;
  ceqc::model::NavigationRecord pn; pn.system="G"; pn.satellite="G01"; pn.recordType="EPH"; pn.messageType="LNAV"; pn.epoch=ceqc::model::makeUTC(2026,1,1,0,0,0);
  for (int i=0;i<31;++i) pn.values.push_back(0.001*i);
  v2padNav.data.navigationRecords.push_back(pn);
  auto v2padNavMerged = ceqc::service::rinex::merge({v2padNav}, ceqc::model::RinexKind::Nav, 2.11);
  for (const auto& l : v2padNavMerged.body) assert(l.size() >= 80);


  std::cout << "ceqc tests passed\n";
}

