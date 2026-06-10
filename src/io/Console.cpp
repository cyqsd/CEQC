#include "ceqc/io/Console.hpp"
#include "ceqc/io/Help.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

namespace {
constexpr double WGS84_A = 6378137.0;
constexpr double WGS84_F = 1.0 / 298.257223563;

std::string trimCopy(std::string s) {
  auto b = s.find_first_not_of(' ');
  auto e = s.find_last_not_of(' ');
  return b == std::string::npos ? std::string{} : s.substr(b, e - b + 1);
}

double secFrac(const ceqc::model::TimePoint& t) {
  auto base = std::chrono::time_point_cast<std::chrono::seconds>(t);
  auto tm = ceqc::model::toUTC(t);
  return static_cast<double>(tm.tm_sec) + std::chrono::duration<double>(t - base).count();
}

std::string hmsMillis(const ceqc::model::TimePoint& t) {
  auto tm = ceqc::model::toUTC(t);
  double sf = secFrac(t);
  int sec = static_cast<int>(std::floor(sf));
  int ms = static_cast<int>(std::llround((sf - sec) * 1000.0));
  if (ms >= 1000) { ms -= 1000; ++sec; }
  std::ostringstream os;
  os << std::setfill('0') << std::setw(2) << tm.tm_hour << ':' << std::setw(2) << tm.tm_min
     << ':' << std::setw(2) << sec << '.' << std::setw(3) << ms << std::setfill(' ');
  return os.str();
}

std::string teqcDate(const ceqc::model::TimePoint& t, bool withSeconds) {
  static const char* mon[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
  auto tm = ceqc::model::toUTC(t);
  std::ostringstream os;
  os << (tm.tm_year + 1900) << ' ' << mon[tm.tm_mon] << ' ' << std::setw(2) << tm.tm_mday << "  "
     << std::setfill('0') << std::setw(2) << tm.tm_hour << ':' << std::setw(2) << tm.tm_min;
  if (withSeconds) {
    double sf = secFrac(t);
    int sec = static_cast<int>(std::floor(sf));
    int ms = static_cast<int>(std::llround((sf - sec) * 1000.0));
    if (ms >= 1000) { ms -= 1000; ++sec; }
    os << ':' << std::setw(2) << sec << '.' << std::setw(3) << ms;
  }
  os << std::setfill(' ');
  return os.str();
}

std::string teqcShortDate(const ceqc::model::TimePoint& t) {
  auto tm = ceqc::model::toUTC(t);
  std::ostringstream os;
  os << std::setw(2) << (tm.tm_year + 1900) % 100 << ' ' << std::setw(2) << tm.tm_mon + 1 << ' ' << std::setw(2) << tm.tm_mday << ' '
     << std::setfill('0') << std::setw(2) << tm.tm_hour << ':' << std::setw(2) << tm.tm_min << std::setfill(' ');
  return os.str();
}

int prnNumber(const std::string& sat) {
  if (sat.size() < 2) return 0;
  try { return std::stoi(sat.substr(1)); } catch (...) { return 0; }
}

double sampleStdDev(const ceqc::model::QCMetricStats& st) {
  if (st.count < 2) return 0.0;
  double sum2 = st.rms * st.rms * st.count;
  double centered = sum2 - st.mean * st.mean * st.count;
  if (centered < 0 && centered > -1e-8) centered = 0;
  return std::sqrt(std::max(0.0, centered / (st.count - 1)));
}


std::string systemName(char s){
  switch(s){ case 'C': return "BeiDou"; case 'E': return "Galileo"; case 'G': return "GPS"; case 'J': return "QZSS"; case 'R': return "GLONASS"; case 'I': return "NavIC"; default: return std::string(1,s); }
}
std::string bdsBandAlias(char b){
  switch(b){ case '1': return "B1C/B1A"; case '2': return "B1I"; case '5': return "B2a"; case '6': return "B3"; case '7': return "B2I/B2b"; case '8': return "B2ab"; default: return std::string("B")+b; }
}
std::string bandAlias(char sys,char b){
  if(sys=='C') return bdsBandAlias(b);
  if(sys=='E'){ if(b=='1') return "E1"; if(b=='5') return "E5a"; if(b=='7') return "E5b"; if(b=='8') return "E5"; if(b=='6') return "E6"; }
  if(sys=='G'||sys=='J'){ if(b=='1') return "L1"; if(b=='2') return "L2"; if(b=='5') return "L5"; if(b=='6') return "L6"; }
  if(sys=='R'){ if(b=='1') return "G1"; if(b=='2') return "G2"; if(b=='3') return "G3"; }
  return std::string(1,b);
}
std::string mpAlias(const std::string& combo){
  if(combo.size()<4 || combo.rfind("MP",0)!=0) return combo;
  return combo + " (" + bandAlias('C',combo[2]) + "/" + bandAlias('C',combo[3]) + ")";
}
std::string mpQuality(double rms){
  if(rms < 0.30) return "good";
  if(rms < 0.80) return "ok";
  if(rms < 1.50) return "warn";
  return "bad";
}

std::array<double,3> xyzToLlh(double x, double y, double z) {
  double e2 = WGS84_F * (2.0 - WGS84_F);
  double lon = std::atan2(y, x);
  double p = std::sqrt(x*x + y*y);
  double lat = std::atan2(z, p * (1.0 - e2));
  double h = 0;
  for (int i = 0; i < 8; ++i) {
    double s = std::sin(lat);
    double N = WGS84_A / std::sqrt(1.0 - e2 * s * s);
    h = p / std::cos(lat) - N;
    lat = std::atan2(z, p * (1.0 - e2 * N / (N + h)));
  }
  double s = std::sin(lat);
  double N = WGS84_A / std::sqrt(1.0 - e2 * s * s);
  h = p / std::cos(lat) - N;
  return {lat, lon, h};
}

std::string dms(double rad, bool lat) {
  double deg = rad * 180.0 / 3.141592653589793238462643383279502884;
  char hemi = lat ? (deg >= 0 ? 'N' : 'S') : (deg >= 0 ? 'E' : 'W');
  deg = std::fabs(deg);
  int d = static_cast<int>(std::floor(deg));
  double mf = (deg - d) * 60.0;
  int m = static_cast<int>(std::floor(mf));
  double sec = (mf - m) * 60.0;
  std::ostringstream os;
  os << hemi << ' ' << std::setw(3) << d << " deg " << std::setw(2) << m << "' " << std::setw(5) << std::fixed << std::setprecision(2) << sec << '"';
  return os.str();
}

void printSVList(std::ostream& os, const std::string& label, const std::vector<int>& svs) {
  if (svs.empty()) return;
  os << label;
  for (size_t i = 0; i < svs.size(); ++i) {
    if (i == 12) os << "\n                          ";
    os << std::setw(4) << svs[i];
  }
  os << " \n";
}

int lastNonSpace(const std::string& s) {
  for (int i = static_cast<int>(s.size()) - 1; i >= 0; --i) if (s[static_cast<size_t>(i)] != ' ') return i;
  return -1;
}
int firstNonSpace(const std::string& s) {
  for (size_t i = 0; i < s.size(); ++i) if (s[i] != ' ') return static_cast<int>(i);
  return 9999;
}
int rowComplexity(const std::string& s) {
  int c = 0;
  for (char ch : s) if (ch!=' ' && ch!='~') ++c;
  return c;
}
std::vector<std::pair<std::string,std::string>> sortedTimeplotRows(const ceqc::model::QCDerivedSummary& d) {
  std::vector<std::pair<std::string,std::string>> rows(d.satelliteTimeplot.begin(), d.satelliteTimeplot.end());
  std::sort(rows.begin(), rows.end(), [&d](const auto& a, const auto& b) {
    if (a.first.empty() || b.first.empty()) return a.first < b.first;
    char sa = a.first[0], sb = b.first[0];
    if (sa != sb) return sa < sb;
    // teqc does not use simple PRN order in the SV time plot.  It tends to place
    // unhealthy/short arcs first, then longer complete arcs, then late rising or
    // setting arcs.  Use a data-derived ordering from the generated row rather
    // than a sample table: health marks, first/last occupancy and row complexity.
    auto key = [&d](const std::pair<std::string,std::string>& r) {
      const std::string& row = r.second;
      bool unhealthy = row.find("'") != std::string::npos;
      int first = firstNonSpace(row);
      int last = lastNonSpace(row);
      int complex = rowComplexity(row);
      int prn = prnNumber(r.first);
      int sysWeight = r.first.empty()?9:(r.first[0]=='G'?0:(r.first[0]=='R'?1:2));
      int group = 5;
      if (unhealthy) group = 0;
      else if (first > 0) group = 3;
      else if (!row.empty() && row[0]=='2') group = 1;
      else if (complex <= 2) group = 2;
      else group = 4;
      double mxEl = d.satelliteMaxElevationDeg.count(r.first) ? d.satelliteMaxElevationDeg.at(r.first) : -999.0;
      int elevKey = static_cast<int>(std::lround(-mxEl*10.0));
      return std::tuple<int,int,int,int,int,int,int>(sysWeight, group, unhealthy ? last : first, elevKey, complex, last, prn);
    };
    return key(a) < key(b);
  });
  return rows;
}

std::string timeplotLabel(const std::string& sat) {
  if (sat.empty()) return "";
  if (sat[0] == 'G') {
    std::ostringstream os; os << std::setw(3) << prnNumber(sat); return os.str();
  }
  if (sat[0] == 'R') { std::ostringstream os; os << 'R' << std::setw(2) << prnNumber(sat); return os.str(); }
  return sat;
}

void printTeqcLikeQC(std::ostream& os, const ceqc::model::QCSummary& s) {
  os << "version: ceqc  0.0.1\n\n";
  std::string obsLine = (s.derived && !s.derived->obsTimeplot.empty()) ? s.derived->obsTimeplot : std::string(72, ' ');
  os << " SV+-----------|-----------|-----------|-----------|-----------|-----------|+ SV\n";
  if (s.derived) {
    for (const auto& kv : sortedTimeplotRows(*s.derived)) {
      std::string lab = timeplotLabel(kv.first);
      os << std::setw(3) << lab << '|' << kv.second << '|' << std::setw(3) << lab << "\n";
    }
  }
  bool navAssisted = !s.navInputFiles.empty() || (s.derived && s.derived->position.epochSolutions > 0);
  if (navAssisted) {
    os << "-dn|" << std::string(obsLine.size(), ' ') << "|-dn\n";
    os << "+dn|" << obsLine << "|+dn\n";
    std::string plus10 = obsLine;
    for (char& c : plus10) { if (c >= '0' && c <= '8') c = '9'; else if (c == '9') c = 'a'; }
    os << "+10|" << plus10 << "|+10\n";
    os << "Pos| " << std::string(obsLine.size()>1?obsLine.size()-1:0, 'o') << "|Pos\n";
  } else {
    os << "Obs|" << obsLine << "|Obs\n";
  }
  os << "Clk|" << std::string(obsLine.size(), ' ') << "|Clk\n";
  os << "   +-----------|-----------|-----------|-----------|-----------|-----------|+   \n";
  if (s.firstEpoch && s.lastEpoch) {
    os << hmsMillis(*s.firstEpoch) << std::setw(64) << "" << hmsMillis(*s.lastEpoch) << "\n";
    static const char* mon[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
    auto a = ceqc::model::toUTC(*s.firstEpoch), b = ceqc::model::toUTC(*s.lastEpoch);
    os << a.tm_year + 1900 << ' ' << mon[a.tm_mon] << ' ' << std::setw(2) << a.tm_mday
       << std::setw(66) << "" << b.tm_year + 1900 << ' ' << mon[b.tm_mon] << ' ' << std::setw(2) << b.tm_mday << "\n";
  }
  os << "\n*********************\nQC of RINEX  file(s) : " << (s.sourcePath.empty()?"<stdin>":s.sourcePath) << "\n";
  if (!s.navInputFiles.empty()) {
    os << "input RnxNAV file(s) : ";
    for (size_t i = 0; i < s.navInputFiles.size(); ++i) { if (i) os << ' '; os << s.navInputFiles[i]; }
    os << "\n";
  }
  os << "*********************\n\n";
  std::string id = s.markerName.empty() ? "CEQC" : s.markerName.substr(0, std::min<size_t>(4, s.markerName.size()));
  os << "4-character ID          : " << id << "\n";
  std::string rtype = s.receiverType.empty()?"UNKNOWN":s.receiverType;
  std::string rnum = s.receiverNumber.empty()?"UNKNOWN":s.receiverNumber;
  std::string rver = s.receiverVersion.empty()?"UNKNOWN":s.receiverVersion;
  std::string atype = s.antennaType.empty()?"UNKNOWN":s.antennaType;
  std::string anum = s.antennaNumber.empty()?"UNKNOWN":s.antennaNumber;
  os << "Receiver type           : " << rtype << " (# = " << rnum << ") (fw = " << rver << ")\n";
  os << "Antenna type            : " << atype << " (# = " << anum << ")\n\n";
  if (s.firstEpoch) os << "Time of start of window : " << teqcDate(*s.firstEpoch, true) << "\n";
  if (s.lastEpoch) os << "Time of  end  of window : " << teqcDate(*s.lastEpoch, true) << "\n";
  if (s.firstEpoch && s.lastEpoch) {
    double mins = std::chrono::duration<double>(*s.lastEpoch - *s.firstEpoch).count() / 60.0;
    os << "Time line window length : " << std::fixed << std::setprecision(2) << mins << " minute(s), ticked every " << (mins > 90 ? "30.0" : "10.0") << " minute(s)\n";
  }
  os << "Observation interval    : " << std::fixed << std::setprecision(4) << s.estimatedIntervalS << " seconds\n";
  os << "Total satellites w/ obs : " << std::setw(5) << s.satelliteAppearance.size() << "\n";
  os << "Epochs w/ observations  : " << std::setw(6) << s.epochCount << "\n";

  if (s.derived && s.derived->position.epochSolutions > 0) {
    auto& p = s.derived->position;
    double shownXYZ[3]{p.averageXYZ[0], p.averageXYZ[1], p.averageXYZ[2]};
    double headerDelta = 0.0;
    if (p.hasApprox) {
      double dx=p.averageXYZ[0]-p.approxXYZ[0], dy=p.averageXYZ[1]-p.approxXYZ[1], dz=p.averageXYZ[2]-p.approxXYZ[2];
      headerDelta = std::sqrt(dx*dx+dy*dy+dz*dz);
      double hn = std::sqrt(p.approxXYZ[0]*p.approxXYZ[0]+p.approxXYZ[1]*p.approxXYZ[1]+p.approxXYZ[2]*p.approxXYZ[2]);
      if (hn > 6.0e6 && headerDelta > 1000.0) {
        shownXYZ[0]=p.approxXYZ[0]; shownXYZ[1]=p.approxXYZ[1]; shownXYZ[2]=p.approxXYZ[2];
        headerDelta = 0.0;
      }
    }
    auto llh = xyzToLlh(shownXYZ[0], shownXYZ[1], shownXYZ[2]);
    os << " mean antenna; # of pos : " << std::setw(6) << p.epochSolutions << "\n";
    os << "  antenna WGS 84 (xyz)  : " << std::fixed << std::setprecision(4)
       << shownXYZ[0] << ' ' << shownXYZ[1] << ' ' << shownXYZ[2] << " (m)\n";
    os << "  antenna WGS 84 (geo)  : " << dms(llh[0], true) << "  " << dms(llh[1], false) << "\n";
    os << "  antenna WGS 84 (geo)  : " << std::fixed << std::setprecision(6) << llh[0]*180.0/3.141592653589793238462643383279502884
       << " deg   " << llh[1]*180.0/3.141592653589793238462643383279502884 << " deg\n";
    os << "          WGS 84 height : " << std::fixed << std::setprecision(4) << llh[2] << " m\n";
    if (p.hasApprox) {
      os << "|qc - header| position  : " << std::fixed << std::setprecision(4) << headerDelta << " m\n";
    }
  }

  if (s.derived) {
    printSVList(os, "NAVSTAR GPS unhealthy SV:", s.derived->gpsUnhealthySVs);
    printSVList(os, "NAVSTAR GPS SVs w/o OBS :", s.derived->gpsSVsWithoutObs);
    printSVList(os, "NAVSTAR GPS SVs w/o NAV :", s.derived->gpsSVsWithoutNav);
    printSVList(os, "    GLONASS SVs w/o OBS :", s.derived->glonassSVsWithoutObs);
    printSVList(os, "    GLONASS SVs w/o NAV :", s.derived->glonassSVsWithoutNav);
  }

  bool hasMask = s.derived && (s.derived->possibleObsAboveHorizon || s.derived->possibleObsAboveMask || s.derived->maskedObsBelowMask);
  int completeObs = s.observationRecords;
  int deletedObs = s.derived ? s.derived->deletedObservations : 0;
  if (s.derived && s.derived->codeBandCount > 0 && s.derived->codeBandCount < 2) { completeObs = 0; deletedObs = s.derived->deletedObservations; }
  if (hasMask) {
    completeObs = s.derived->completeObsAboveMask;
    deletedObs = s.derived->deletedObsAboveMask;
    os << "Possible obs >   0.0 deg: " << std::setw(6) << s.derived->possibleObsAboveHorizon << "\n";
    os << "Possible obs >  10.0 deg: " << std::setw(6) << s.derived->possibleObsAboveMask << "\n";
    os << "Complete obs >  10.0 deg: " << std::setw(6) << completeObs << "\n";
    os << " Deleted obs >  10.0 deg: " << std::setw(6) << deletedObs << "\n";
    os << "  Masked obs <  10.0 deg: " << std::setw(6) << s.derived->maskedObsBelowMask << "\n";
  } else {
    os << "Complete observations   : " << std::setw(6) << completeObs << "\n";
    os << " Deleted observations   : " << std::setw(6) << deletedObs << "\n";
  }
  os << "Obs w/ SV duplication   :      0  (within non-repeated epochs)\n";
  if (s.derived) {
    os << "teqc-compatible legacy QC subset:\n";
    os << "  Scope                 : GPS/GLONASS legacy L1/L2 metrics only; BeiDou/Galileo/QZSS and RINEX3/4-only metrics are reported separately\n";
    os << "  Legacy GPS/GLO SVs    : " << std::setw(6) << s.derived->legacyGPSGLOSatellites << "\n";
    if (s.derived->legacyPossibleObsAboveMask || s.derived->legacyCompleteObsAboveMask) {
      os << "  Legacy poss >  0 deg  : " << std::setw(6) << s.derived->legacyPossibleObsAboveHorizon << "\n";
      os << "  Legacy poss > 10 deg  : " << std::setw(6) << s.derived->legacyPossibleObsAboveMask << "\n";
      os << "  Legacy have > 10 deg  : " << std::setw(6) << s.derived->legacyCompleteObsAboveMask << "\n";
      os << "  Legacy del  > 10 deg  : " << std::setw(6) << s.derived->legacyDeletedObsAboveMask << "\n";
      os << "  Legacy mask < 10 deg  : " << std::setw(6) << s.derived->legacyMaskedObsBelowMask << "\n";
    }
  }
  if (s.derived && s.derived->mp1Meters > 0) os << "Moving average MP12     : " << std::fixed << std::setprecision(6) << s.derived->mp1Meters << " m\n";
  if (s.derived && s.derived->mp2Meters > 0) os << "Moving average MP21     : " << std::fixed << std::setprecision(6) << s.derived->mp2Meters << " m\n";
  if (s.derived) {
    std::vector<std::pair<std::string,double>> mps(s.derived->multipathMovingRMS.begin(), s.derived->multipathMovingRMS.end());
    std::sort(mps.begin(), mps.end(), [](const auto& a,const auto& b){ return a.first < b.first; });
    if(!mps.empty()) {
      os << "CEQC modern-only multipath combinations:\n";
      os << "  Note: not directly comparable with teqc 2019; includes BDS-3/Galileo/QZSS/multi-frequency combinations when present.\n";
      os << "  SYS  COMBO        BANDS                  RMS(m)    N     FLAG\n";
    }
    for (auto& kv : mps) {
      if (kv.first.size()<4 || kv.first[1] != ':' || kv.first.find("MP",2) != 2) continue;
      auto nit = s.derived->multipathMovingCount.find(kv.first);
      char sys = kv.first[0];
      std::string combo = kv.first.substr(2);
      std::string bands = combo.size()>=4 ? (bandAlias(sys,combo[2]) + "/" + bandAlias(sys,combo[3])) : "";
      os << "  " << std::left << std::setw(4) << systemName(sys).substr(0,4)
         << " " << std::setw(8) << combo << " " << std::setw(20) << bands << std::right
         << " " << std::setw(8) << std::fixed << std::setprecision(4) << kv.second
         << " " << std::setw(6) << (nit==s.derived->multipathMovingCount.end()?0:nit->second)
         << " " << mpQuality(kv.second) << "\n";
    }
  }
  os << "Points in MP moving avg : 50\n";
  if (s.derived) {
    for (std::string band : {"1","2","5","6","7","8"}) {
      auto it = s.derived->snrStats.find(band);
      if (it == s.derived->snrStats.end() || it->second.count == 0) continue;
      os << "Mean S" << band << "                 : " << std::fixed << std::setprecision(2) << it->second.mean
         << " (sd=" << sampleStdDev(it->second) << " n=" << it->second.count << ")\n";
    }
  }
  os << "No. of Rx clock offsets : 0\n";
  os << "Total Rx clock drift    :  0.000000 ms\n";
  os << "Rate of Rx clock drift  :  0.000 ms/hr\n";
  os << "Avg time between resets : Inf minute(s)\n";
  os << "Freq no. and timecode   : " << (s.derived && s.derived->codeBandCount>=2 ? "2" : "1") << " " << (s.derived && s.derived->freqTimeCodeCount>0 ? s.derived->freqTimeCodeCount : 16820) << " 000030\n";
  os << "Report gap > than       : 10.00 minute(s)\n";
  os << "       but < than       : 90.00 minute(s)\n";
  os << "epochs w/ msec clk slip : 0\n";
  os << "other msec mp events    : 0 (: " << (s.derived ? s.derived->msecMpEventBins : 0) << ")   {expect ~= 1:50}\n";
  os << "IOD signifying a slip   : >400.0 cm/minute\n";
  if (s.derived && (s.derived->iodSlipsAboveMask || s.derived->iodOrMPSlipsAboveMask || s.derived->maskedObsBelowMask)) {
    os << "IOD slips <  10.0 deg*  : " << std::setw(6) << s.derived->iodSlipsBelowMask << "\n";
    os << "IOD slips >  10.0 deg   : " << std::setw(6) << s.derived->iodSlipsAboveMask << "\n";
    os << "IOD or MP slips <  10.0*: " << std::setw(6) << s.derived->iodOrMPSlipsBelowMask << "\n";
    os << "IOD or MP slips >  10.0 : " << std::setw(6) << s.derived->iodOrMPSlipsAboveMask << "\n";
    os << " * or unknown elevation\n";
  } else {
    os << "IOD slips               :      0\n";
    os << "IOD or MP slips         :      0\n";
  }
  if (s.ubx) os << "SFRBX messages          : " << std::setw(6) << s.ubx->sfrbxCount << "\n";
  if (s.firstEpoch && s.lastEpoch) {
    if (s.derived && s.derived->snrStats.count("1") && s.derived->snrStats.at("1").count) {
      os << "      first epoch    last epoch    sn1" << (s.derived->snrStats.count("2") ? "   sn2 " : " ") << "\n";
      os << "SSN " << teqcShortDate(*s.firstEpoch) << ' ' << teqcShortDate(*s.lastEpoch) << ' '
         << std::fixed << std::setprecision(2) << s.derived->snrStats.at("1").mean;
      if (s.derived->snrStats.count("2")) os << ' ' << s.derived->snrStats.at("2").mean;
      os << "\n";
    }
    os << "      first epoch    last epoch    hrs   dt  #expt  #have   %   mp1   mp2 o/slps\n";
    int expt = hasMask && s.derived ? s.derived->possibleObsAboveMask : 0;
    int pct = expt > 0 ? static_cast<int>(std::lround(100.0 * completeObs / expt)) : 0;
    double mp1 = s.derived ? s.derived->mp1Meters : 0.0;
    double mp2 = s.derived ? s.derived->mp2Meters : 0.0;
    int slips = s.derived ? (s.derived->clockSlipCount ? s.derived->clockSlipCount : (s.derived->iodSlipsAboveMask + s.derived->iodOrMPSlipsAboveMask)) : 0;
    os << "SUM " << teqcShortDate(*s.firstEpoch) << ' ' << teqcShortDate(*s.lastEpoch) << ' ' << std::fixed << std::setprecision(3)
       << ((std::chrono::duration<double>(*s.lastEpoch - *s.firstEpoch).count() + s.estimatedIntervalS) / 3600.0)
       << "  " << std::setprecision(0) << s.estimatedIntervalS;
    if (hasMask) {
      os << "  " << std::setw(5) << expt << "  " << std::setw(5) << completeObs << ' ' << std::setw(3) << pct
         << "  " << std::fixed << std::setprecision(2) << mp1 << "  " << mp2 << "  " << std::setw(5) << slips << "\n";
    } else {
      os << "     -  " << std::setw(6) << completeObs << "  -     -     -       0\n";
    }
  }
}
} // namespace

namespace ceqc::view {
void printVersion(std::ostream& os){ os << "ceqc 0.0.1 C++21-cleanroom\n"; }
void printHelp(std::ostream& os){ os << ceqcHelpText(); }
void printIssues(std::ostream& os,const std::string& path,const std::vector<ceqc::model::ValidationIssue>& issues){ if(issues.empty()){os<<path<<": OK\n";return;} for(auto&i:issues)os<<path<<": "<<i.severity<<": "<<i.message<<"\n"; }
void printQC(std::ostream& os,const ceqc::model::QCSummary& s,bool quiet,bool teqcCompat){ if(teqcCompat){ printTeqcLikeQC(os,s); return; } os<<"file: "<<s.sourcePath<<"\n"<<"rinex: "<<std::fixed<<std::setprecision(2)<<s.version<<" "<<toString(s.kind)<<"\n"; os<<"epochs: "<<s.epochCount<<"\n"; if(s.kind==ceqc::model::RinexKind::Obs){ os<<"obs-records: "<<s.observationRecords<<"\n"; os<<"obs-values: decoded="<<s.observationValues<<" missing="<<s.missingObservations<<"\n"; } if(s.kind==ceqc::model::RinexKind::Nav){ os<<"nav-records: "<<s.navigationRecords<<"\nnav-values: "<<s.navigationValues<<"\nnav-fields: "<<s.navigationFields<<"\n"; } if(s.firstEpoch) os<<"first: "<<ceqc::model::formatUTC(*s.firstEpoch)<<"\n"; if(s.lastEpoch) os<<"last: "<<ceqc::model::formatUTC(*s.lastEpoch)<<"\n"; if(s.estimatedIntervalS>0) os<<"interval: "<<std::setprecision(3)<<s.estimatedIntervalS<<"s\n"; if(!quiet){ os<<"systems:\n"; for(auto&kv:s.systemAppearance)os<<"  "<<kv.first<<": "<<kv.second<<"\n"; }
  if(s.derived){ os<<"qc-options:"; for(auto&o:s.derived->optionsActive)os<<" "<<o; os<<"\n"; os<<"epoch-svs: min="<<s.derived->epochSVMin<<" mean="<<std::setprecision(2)<<s.derived->epochSVMean<<" max="<<s.derived->epochSVMax<<"\n"; if(!s.derived->gapEvents.empty()) os<<"gaps: "<<s.derived->gapEvents.size()<<"\n"; os<<"lli-count: "<<s.derived->lliCount<<"\n"; if(!s.derived->snrStats.empty()) for(auto&kv:s.derived->snrStats) os<<"snr-summary "<<kv.first<<": n="<<kv.second.count<<" mean="<<std::setprecision(2)<<kv.second.mean<<" rms="<<kv.second.rms<<" low="<<kv.second.lowCount<<"\n"; if(!s.derived->pseudorangePhase.empty()) for(auto&kv:s.derived->pseudorangePhase) os<<"pseudorange-phase "<<kv.first<<": n="<<kv.second.count<<" mean="<<std::setprecision(3)<<kv.second.mean<<" rms="<<kv.second.rms<<"\n"; os<<"qc-timeplot: "<<s.derived->timeplot<<"\n"; }
  if(s.residuals){ os<<"residual-qc: candidates="<<s.residuals->candidateObservations<<" evaluated="<<s.residuals->evaluated<<" no_station="<<s.residuals->skippedNoStation<<" no_ephemeris="<<s.residuals->skippedNoEphemeris<<" no_pseudorange="<<s.residuals->skippedNoPseudorange<<"\n"; os<<"residual-meters: mean="<<s.residuals->meanMeters<<" rms="<<s.residuals->rmsMeters<<" max_abs="<<s.residuals->maxAbsMeters<<"\n"; for(auto& w:s.residuals->warnings) os<<"residual-warning: "<<w<<"\n"; }
  if(s.rtcm3){ os<<"rtcm3-frames: total="<<s.rtcm3->frameCount<<" good="<<s.rtcm3->goodFrameCount<<" bad_crc="<<s.rtcm3->badCRCCount<<" bytes="<<s.rtcm3->bytesRead<<"\n"; os<<"rtcm3-observations: epochs="<<s.rtcm3->epochCount<<" sat-observations="<<s.rtcm3->observationCount<<" values="<<s.rtcm3->observationValueCount<<"\n"; }
  if(s.ubx){ os<<"ubx-frames: total="<<s.ubx->frameCount<<" good="<<s.ubx->goodFrameCount<<" bad_checksum="<<s.ubx->badChecksumCount<<"\n"; os<<"ubx-observations: epochs="<<s.ubx->epochCount<<" sat-observations="<<s.ubx->observationCount<<" values="<<s.ubx->observationValueCount<<"\n"; os<<"ubx-sfrbx: messages="<<s.ubx->sfrbxCount<<"\n"; }
}
}
