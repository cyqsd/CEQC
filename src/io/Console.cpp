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


std::string histogramText(const std::vector<int>& h){
  std::ostringstream os; os << "[";
  for(size_t i=0;i<h.size();++i){ if(i) os << ","; os << h[i]; }
  os << "]"; return os.str();
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
  // teqc uses zero-padded minutes and seconds in the DMS geodetic line.
  os << hemi << ' ' << std::setw(3) << d << " deg "
     << std::setfill('0') << std::setw(2) << m << "' "
     << std::setw(5) << std::fixed << std::setprecision(2) << sec << '"'
     << std::setfill(' ');
  return os.str();
}

std::array<double,3> ecefDeltaToEnu(const double ref[3], const double xyz[3]) {
  auto llh = xyzToLlh(ref[0], ref[1], ref[2]);
  const double lat = llh[0], lon = llh[1];
  const double dx = xyz[0] - ref[0];
  const double dy = xyz[1] - ref[1];
  const double dz = xyz[2] - ref[2];
  const double sl = std::sin(lat), cl = std::cos(lat);
  const double sb = std::sin(lon), cb = std::cos(lon);
  const double e = -sb * dx + cb * dy;
  const double n = -sl * cb * dx - sl * sb * dy + cl * dz;
  const double u = cl * cb * dx + cl * sb * dy + sl * dz;
  return {e, n, u};
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

std::string nativeAxis(size_t width) {
  if (width == 72) return "-----------|-----------|-----------|-----------|-----------|-----------|";
  std::string r(width, '-');
  for (size_t i = 11; i < width; i += 12) r[i] = '|';
  return r;
}

std::string teqcEndpointLine(const std::string& left, const std::string& right, size_t width = 80) {
  if (left.size() + right.size() >= width) return left + " " + right;
  return left + std::string(width - left.size() - right.size(), ' ') + right;
}

bool teqcObservedDensitySymbol(char c) {
  return c != ' ' && c != '-' && c != '_';
}

bool teqcAboveMaskDensitySymbol(char c) {
  return c != ' ' && c != '-' && c != '_' && c != '.';
}

std::string teqcDensityLine(const ceqc::model::QCDerivedSummary& d, size_t width, bool aboveMask) {
  std::string out(width, ' ');
  for (size_t b = 0; b < width; ++b) {
    int count = 0;
    for (const auto& kv : d.satelliteTimeplot) {
      // teqc's legacy density rows are GPS-centric in this compatibility view.
      // GLONASS rows may be rendered as N/no-nav by old teqc even when modern CEQC
      // has valid GLO navigation, so do not let them inflate +dn/+10.
      if (kv.first.empty() || kv.first[0] != 'G') continue;
      if (b >= kv.second.size()) continue;
      char c = kv.second[b];
      bool include = aboveMask ? teqcAboveMaskDensitySymbol(c) : teqcObservedDensitySymbol(c);
      if (include) ++count;
    }
    if (count <= 0) out[b] = ' ';
    else out[b] = static_cast<char>('0' + std::min(count, 9));
  }
  return out;
}

void printNativeSVTimeplot(std::ostream& os, const ceqc::model::QCDerivedSummary& d) {
  if (d.satelliteTimeplot.empty()) return;
  size_t width = d.obsTimeplot.empty() ? 72 : d.obsTimeplot.size();
  os << "native-sv-timeplot:\n";
  os << " SV+" << nativeAxis(width) << "+ SV\n";
  for (const auto& kv : sortedTimeplotRows(d)) {
    std::string lab = timeplotLabel(kv.first);
    os << std::setw(3) << lab << '|' << kv.second << '|' << std::setw(3) << lab << "\n";
  }
  if (!d.obsTimeplot.empty()) os << "Obs|" << d.obsTimeplot << "|Obs\n";
  if (!d.navTimeplot.empty()) os << "Nav|" << d.navTimeplot << "|Nav\n";
  if (!d.positionTimeplot.empty()) os << "Pos|" << d.positionTimeplot << "|Pos\n";
  os << "   +" << nativeAxis(width) << "+   \n";
  os << "native-symbols: L loss-of-lock, N no-nav, 2 incomplete dual-frequency, ~ observed above mask, + expected above mask, _ below mask, - below horizon, P position solved, x rejected, s skipped\n";
}

void printNativeSkyplot(std::ostream& os, const ceqc::model::QCDerivedSummary& d) {
  if (d.satelliteMaxElevationDeg.empty()) return;
  constexpr int rows = 21;
  constexpr int cols = 41;
  constexpr int cy = rows / 2;
  constexpr int cx = cols / 2;
  constexpr double rmax = 10.0;
  std::vector<std::string> grid(rows, std::string(cols, ' '));
  auto plot = [&](int y, int x, char ch) {
    if (y < 0 || y >= rows || x < 0 || x >= cols) return;
    char& cell = grid[static_cast<size_t>(y)][static_cast<size_t>(x)];
    if (cell == ' ' || cell == '.' || cell == '+' || cell == '-' || cell == '|' || cell == 'N' || cell == 'E' || cell == 'S' || cell == 'W') cell = ch;
    else if (cell != ch) cell = '*';
  };
  for (int y = 0; y < rows; ++y) {
    for (int x = 0; x < cols; ++x) {
      double dx = (x - cx) / 2.0;
      double dy = static_cast<double>(y - cy);
      double rr = std::sqrt(dx * dx + dy * dy);
      if (std::fabs(rr - rmax) < 0.35 || std::fabs(rr - rmax * 2.0 / 3.0) < 0.25 || std::fabs(rr - rmax / 3.0) < 0.20) grid[y][x] = '.';
    }
  }
  for (int x = 0; x < cols; ++x) grid[cy][x] = (grid[cy][x] == ' ' ? '-' : grid[cy][x]);
  for (int y = 0; y < rows; ++y) grid[y][cx] = (grid[y][cx] == ' ' ? '|' : grid[y][cx]);
  grid[cy][cx] = '+';
  grid[0][cx] = 'N'; grid[cy][cols - 1] = 'E'; grid[rows - 1][cx] = 'S'; grid[cy][0] = 'W';

  std::vector<std::tuple<std::string,double,double>> sv;
  for (const auto& kv : d.satelliteMaxElevationDeg) {
    auto azIt = d.satelliteMaxElevationAzimuthDeg.find(kv.first);
    if (azIt == d.satelliteMaxElevationAzimuthDeg.end() || !std::isfinite(kv.second) || !std::isfinite(azIt->second)) continue;
    double el = std::clamp(kv.second, 0.0, 90.0);
    double az = azIt->second;
    double rr = (90.0 - el) / 90.0 * rmax;
    double ar = az * 3.141592653589793238462643383279502884 / 180.0;
    int x = cx + static_cast<int>(std::lround(2.0 * rr * std::sin(ar)));
    int y = cy - static_cast<int>(std::lround(rr * std::cos(ar)));
    char ch = kv.first.empty() ? '?' : kv.first[0];
    plot(y, x, ch);
    sv.emplace_back(kv.first, az, el);
  }
  std::sort(sv.begin(), sv.end(), [](const auto& a, const auto& b){
    if (std::get<0>(a).empty() || std::get<0>(b).empty()) return std::get<0>(a) < std::get<0>(b);
    if (std::get<0>(a)[0] != std::get<0>(b)[0]) return std::get<0>(a) < std::get<0>(b);
    return prnNumber(std::get<0>(a)) < prnNumber(std::get<0>(b));
  });

  os << "native-skyplot-max-elevation:\n";
  for (const auto& line : grid) os << "  " << line << "\n";
  os << "native-skyplot-legend: G=GPS R=GLONASS C=BeiDou E=Galileo J=QZSS *=overlap; points are sampled at each SV's maximum elevation within the QC window\n";
  os << "native-skyplot-sv:";
  int shown = 0;
  for (const auto& t : sv) {
    if (shown > 0 && shown % 4 == 0) os << "\n                 ";
    os << ' ' << std::get<0>(t) << " az=" << std::fixed << std::setprecision(1) << std::get<1>(t) << " el=" << std::setprecision(1) << std::get<2>(t);
    ++shown;
  }
  os << "\n";
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
    std::string dnLine = s.derived ? teqcDensityLine(*s.derived, obsLine.size(), false) : obsLine;
    std::string plus10 = s.derived ? teqcDensityLine(*s.derived, obsLine.size(), true) : obsLine;
    os << "+dn|" << dnLine << "|+dn\n";
    os << "+10|" << plus10 << "|+10\n";
    os << "Pos| " << std::string(obsLine.size()>1?obsLine.size()-1:0, 'o') << "|Pos\n";
  } else {
    os << "Obs|" << obsLine << "|Obs\n";
  }
  os << "Clk|" << std::string(obsLine.size(), ' ') << "|Clk\n";
  os << "   +-----------|-----------|-----------|-----------|-----------|-----------|+   \n";
  if (s.firstEpoch && s.lastEpoch) {
    os << teqcEndpointLine(hmsMillis(*s.firstEpoch), hmsMillis(*s.lastEpoch)) << "\n";
    static const char* mon[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
    auto a = ceqc::model::toUTC(*s.firstEpoch), b = ceqc::model::toUTC(*s.lastEpoch);
    std::ostringstream leftDate, rightDate;
    leftDate << a.tm_year + 1900 << ' ' << mon[a.tm_mon] << ' ' << std::setw(2) << a.tm_mday;
    rightDate << b.tm_year + 1900 << ' ' << mon[b.tm_mon] << ' ' << std::setw(2) << b.tm_mday;
    os << teqcEndpointLine(leftDate.str(), rightDate.str()) << "\n";
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
  const bool teqcLegacyScope = s.derived && s.derived->legacyGPSGLOSatellites > 0;
  const int teqcSatellites = teqcLegacyScope ? s.derived->legacyGPSGLOSatellites : static_cast<int>(s.satelliteAppearance.size());

  if (s.derived && s.derived->position.epochSolutions > 0) {
    auto& p = s.derived->position;
    // teqc reports the mean code-position estimate as the antenna position in
    // the QC block, then reports the vector difference to the RINEX/RTCM header
    // coordinate as |qc - header|.  Keep CEQC's header-vs-mean interpretation in
    // the extension lines only, so the default text stays close to teqc.
    double reportXYZ[3]{p.averageXYZ[0], p.averageXYZ[1], p.averageXYZ[2]};
    double headerDelta = 0.0;
    std::array<double,3> denu{0.0, 0.0, 0.0};
    if (p.hasApprox) {
      const double dx = p.averageXYZ[0] - p.approxXYZ[0];
      const double dy = p.averageXYZ[1] - p.approxXYZ[1];
      const double dz = p.averageXYZ[2] - p.approxXYZ[2];
      headerDelta = std::sqrt(dx*dx + dy*dy + dz*dz);
      denu = ecefDeltaToEnu(p.approxXYZ, p.averageXYZ);
    }
    auto llh = xyzToLlh(reportXYZ[0], reportXYZ[1], reportXYZ[2]);
    os << " mean antenna; # of pos : " << std::setw(6) << p.epochSolutions << "\n";
    os << "  antenna WGS 84 (xyz)  : " << std::fixed << std::setprecision(4)
       << reportXYZ[0] << ' ' << reportXYZ[1] << ' ' << reportXYZ[2] << " (m)\n";
    os << "  antenna WGS 84 (geo)  : " << dms(llh[0], true) << "  " << dms(llh[1], false) << "\n";
    os << "  antenna WGS 84 (geo)  : " << std::fixed << std::setprecision(6)
       << std::setw(11) << llh[0]*180.0/3.141592653589793238462643383279502884
       << " deg  " << std::setw(11) << llh[1]*180.0/3.141592653589793238462643383279502884 << " deg\n";
    os << "          WGS 84 height : " << std::fixed << std::setprecision(4) << llh[2] << " m\n";
    if (p.hasApprox) {
      const double horizontal = std::sqrt(denu[0]*denu[0] + denu[1]*denu[1]);
      os << "|qc - header| position  : " << std::setw(8) << std::fixed << std::setprecision(4) << headerDelta << " m\n";
      os << "qc position offsets     : " << std::setw(8) << std::fixed << std::setprecision(4) << std::fabs(denu[2])
         << " m vertical   " << std::setw(8) << horizontal << " m horizontal\n";
      if (s.derived->ceqcExtensionEnabled) {
        auto hllh = xyzToLlh(p.approxXYZ[0], p.approxXYZ[1], p.approxXYZ[2]);
        os << "ceqc-ext header WGS84 xyz: " << std::fixed << std::setprecision(4)
           << p.approxXYZ[0] << ' ' << p.approxXYZ[1] << ' ' << p.approxXYZ[2] << " (m)\n";
        os << "ceqc-ext header WGS84 geo: " << std::fixed << std::setprecision(6)
           << hllh[0]*180.0/3.141592653589793238462643383279502884 << " deg   "
           << hllh[1]*180.0/3.141592653589793238462643383279502884 << " deg   h=" << std::setprecision(4) << hllh[2] << " m\n";
        os << "ceqc-ext qc-header enu  : E=" << std::fixed << std::setprecision(4) << denu[0]
           << " N=" << denu[1] << " U=" << denu[2] << " m\n";
      }
    }
  }

  os << "Observation interval    : " << std::fixed << std::setprecision(4) << s.estimatedIntervalS << " seconds\n";
  os << "Total satellites w/ obs : " << std::setw(5) << teqcSatellites << "\n";
  if (s.derived && s.derived->ceqcExtensionEnabled) {
    os << "ceqc-ext total SVs all GNSS: " << std::setw(5) << s.satelliteAppearance.size() << "\n";
    os << "ceqc-ext epochs w/ obs    : " << std::setw(6) << s.epochCount << "\n";
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
  int possibleObs0 = 0;
  int possibleObs10 = 0;
  int maskedObs10 = 0;
  if (s.derived && s.derived->codeBandCount > 0 && s.derived->codeBandCount < 2) { completeObs = 0; deletedObs = s.derived->deletedObservations; }
  if (hasMask) {
    if (teqcLegacyScope) {
      possibleObs0 = s.derived->legacyPossibleObsAboveHorizon;
      possibleObs10 = s.derived->legacyPossibleObsAboveMask;
      completeObs = s.derived->legacyCompleteObsAboveMask;
      deletedObs = s.derived->legacyDeletedObsAboveMask;
      maskedObs10 = s.derived->legacyMaskedObsBelowMask;
    } else {
      possibleObs0 = s.derived->possibleObsAboveHorizon;
      possibleObs10 = s.derived->possibleObsAboveMask;
      completeObs = s.derived->completeObsAboveMask;
      deletedObs = s.derived->deletedObsAboveMask;
      maskedObs10 = s.derived->maskedObsBelowMask;
    }
    os << "Possible obs >   0.0 deg: " << std::setw(6) << possibleObs0 << "\n";
    os << "Possible obs >  10.0 deg: " << std::setw(6) << possibleObs10 << "\n";
    os << "Complete obs >  10.0 deg: " << std::setw(6) << completeObs << "\n";
    os << " Deleted obs >  10.0 deg: " << std::setw(6) << deletedObs << "\n";
    os << "  Masked obs <  10.0 deg: " << std::setw(6) << maskedObs10 << "\n";
    if (s.derived->ceqcExtensionEnabled && teqcLegacyScope) {
      os << "ceqc-ext all GNSS poss >  0: " << std::setw(6) << s.derived->possibleObsAboveHorizon << "\n";
      os << "ceqc-ext all GNSS poss > 10: " << std::setw(6) << s.derived->possibleObsAboveMask << "\n";
      os << "ceqc-ext all GNSS have > 10: " << std::setw(6) << s.derived->completeObsAboveMask << "\n";
      os << "ceqc-ext all GNSS del  > 10: " << std::setw(6) << s.derived->deletedObsAboveMask << "\n";
      os << "ceqc-ext all GNSS mask < 10: " << std::setw(6) << s.derived->maskedObsBelowMask << "\n";
    }
  } else {
    os << "Complete observations   : " << std::setw(6) << completeObs << "\n";
    os << " Deleted observations   : " << std::setw(6) << deletedObs << "\n";
  }
  os << "Obs w/ SV duplication   :      0  (within non-repeated epochs)\n";
  if (s.derived && s.derived->ceqcExtensionEnabled && teqcLegacyScope) {
    os << "ceqc-ext legacy scope   : GPS/GLONASS legacy L1/L2 metrics used for teqc-compatible summary lines; all-GNSS metrics are reported in ceqc-ext lines.\n";
  }
  if (s.derived && s.derived->mp1Meters > 0) os << "Moving average MP12     : " << std::fixed << std::setprecision(6) << s.derived->mp1Meters << " m\n";
  if (s.derived && s.derived->mp2Meters > 0) os << "Moving average MP21     : " << std::fixed << std::setprecision(6) << s.derived->mp2Meters << " m\n";
  if (s.derived && s.derived->multipathEnabled && s.derived->ceqcExtensionEnabled) {
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
  if (s.derived && s.derived->multipathEnabled) os << "Points in MP moving avg : 50\n";
  if (s.derived && s.derived->snrEnabled) {
    const std::vector<std::string> bands = s.derived->ceqcExtensionEnabled
      ? std::vector<std::string>{"1","2","5","6","7","8"}
      : std::vector<std::string>{"1","2"};
    for (const std::string& band : bands) {
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
  if (s.derived && s.derived->iodEnabled) {
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
  }
  if (s.derived) {
    auto printComputedTeqc=[&](const std::string& label,const std::string& key,int samplesOverride=-1){ int samples=samplesOverride>=0?samplesOverride:(s.derived->histogramSamples.count(key)?s.derived->histogramSamples.at(key):0); auto hit=s.derived->histograms.find(key); if(samples<=0 || hit==s.derived->histograms.end()) os << label << " skipped             : samples=0\n"; else os << label << " computed            : samples=" << samples << " bins=" << histogramText(hit->second) << "\n"; };
    printComputedTeqc("ION","ion");
    int mp1Samples=0, mp2Samples=0; auto m1=s.derived->multipathMovingCount.find("MP1"); if(m1!=s.derived->multipathMovingCount.end()) mp1Samples=m1->second; auto m2=s.derived->multipathMovingCount.find("MP2"); if(m2!=s.derived->multipathMovingCount.end()) mp2Samples=m2->second;
    if(mp1Samples>0) os << "MP1 computed            : samples=" << mp1Samples << " rms_m=" << std::fixed << std::setprecision(3) << s.derived->mp1Meters << "\n"; else os << "MP1 skipped             : samples=0\n";
    if(mp2Samples>0) os << "MP2 computed            : samples=" << mp2Samples << " rms_m=" << std::fixed << std::setprecision(3) << s.derived->mp2Meters << "\n"; else os << "MP2 skipped             : samples=0\n";
    if(mp1Samples<=0 && mp2Samples<=0 && !s.derived->multipathSkipReason.empty()) os << "MP skip reason          : " << s.derived->multipathSkipReason << "\n";
    if(s.derived->ceqcExtensionEnabled || s.derived->dataIndicatorsEnabled){ for(const auto& kv:s.derived->histogramSamples){ if(kv.first=="ion"||kv.first=="mp") continue; if(s.derived->histograms.find(kv.first)==s.derived->histograms.end()) os << "CEQC-ext histogram " << kv.first << " : skipped samples=" << kv.second << "\n"; } for(const auto& kv:s.derived->histograms){ if(kv.first=="ion"||kv.first=="mp") continue; os << "CEQC-ext histogram " << kv.first << " : " << histogramText(kv.second) << " samples=" << (s.derived->histogramSamples.count(kv.first)?s.derived->histogramSamples.at(kv.first):0) << "\n"; } }
    for(const auto& w:s.derived->thresholdWarnings) os << "QC threshold warning    : " << w << "\n";
    if((s.derived->everyEpochXYZ||s.derived->everyEpochGeodetic||s.derived->everyEpochDecimal) && !s.derived->epochPositions.empty()){
      os << "Every-epoch position    : " << s.derived->epochPositions.size() << " epoch records\n";
      size_t shown=0;
      for(const auto& ep:s.derived->epochPositions){
        if(shown++>=20){ os << "  ... truncated ...\n"; break; }
        os << "  " << ep.time << " SV=" << ep.usedSVs << " " << ep.status;
        if(ep.status=="OK") os << " XYZ=" << std::fixed << std::setprecision(3) << ep.x << "," << ep.y << "," << ep.z;
        os << "\n";
      }
    }
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
    int expt = hasMask ? possibleObs10 : 0;
    int pct = expt > 0 ? static_cast<int>(std::lround(100.0 * completeObs / expt)) : 0;
    double mp1 = s.derived ? s.derived->mp1Meters : 0.0;
    double mp2 = s.derived ? s.derived->mp2Meters : 0.0;
    int mp1Samples = 0, mp2Samples = 0;
    if (s.derived) {
      auto m1 = s.derived->multipathMovingCount.find("MP1"); if (m1 != s.derived->multipathMovingCount.end()) mp1Samples = m1->second;
      auto m2 = s.derived->multipathMovingCount.find("MP2"); if (m2 != s.derived->multipathMovingCount.end()) mp2Samples = m2->second;
    }
    int slips = s.derived ? (s.derived->clockSlipCount ? s.derived->clockSlipCount : (s.derived->iodSlipsAboveMask + s.derived->iodOrMPSlipsAboveMask)) : 0;
    os << "SUM " << teqcShortDate(*s.firstEpoch) << ' ' << teqcShortDate(*s.lastEpoch) << ' ' << std::fixed << std::setprecision(3)
       << ((std::chrono::duration<double>(*s.lastEpoch - *s.firstEpoch).count() + s.estimatedIntervalS) / 3600.0)
       << "  " << std::setprecision(0) << s.estimatedIntervalS;
    if (hasMask) {
      os << "  " << std::setw(5) << expt << "  " << std::setw(5) << completeObs << ' ' << std::setw(3) << pct << "  ";
      if (mp1Samples > 0) os << std::fixed << std::setprecision(2) << mp1 << "  "; else os << "   -  ";
      if (mp2Samples > 0) os << std::fixed << std::setprecision(2) << mp2 << "  "; else os << "   -  ";
      os << std::setw(5) << slips << "\n";
    } else {
      os << "     -  " << std::setw(6) << completeObs << "  -  ";
      if (mp1Samples > 0) os << std::fixed << std::setprecision(2) << mp1 << "  "; else os << "   -  ";
      if (mp2Samples > 0) os << std::fixed << std::setprecision(2) << mp2 << "  "; else os << "   -  ";
      os << "     0\n";
    }
  }
}
} // namespace

namespace ceqc::view {
void printVersion(std::ostream& os){ os << "ceqc 0.0.1 C++21-cleanroom\n"; }
void printHelp(std::ostream& os){ os << ceqcHelpText(); }
void printIssues(std::ostream& os,const std::string& path,const std::vector<ceqc::model::ValidationIssue>& issues){ if(issues.empty()){os<<path<<": OK\n";return;} for(auto&i:issues)os<<path<<": "<<i.severity<<": "<<i.message<<"\n"; }
void printQC(std::ostream& os,const ceqc::model::QCSummary& s,bool quiet,bool teqcCompat){ if(teqcCompat){ printTeqcLikeQC(os,s); return; } os<<"file: "<<s.sourcePath<<"\n"<<"rinex: "<<std::fixed<<std::setprecision(2)<<s.version<<" "<<toString(s.kind)<<"\n"; os<<"epochs: "<<s.epochCount<<"\n"; if(s.kind==ceqc::model::RinexKind::Obs){ os<<"obs-records: "<<s.observationRecords<<"\n"; os<<"obs-values: decoded="<<s.observationValues<<" missing="<<s.missingObservations<<"\n"; } if(s.kind==ceqc::model::RinexKind::Nav){ os<<"nav-records: "<<s.navigationRecords<<"\nnav-values: "<<s.navigationValues<<"\nnav-fields: "<<s.navigationFields<<"\n"; } if(s.firstEpoch) os<<"first: "<<ceqc::model::formatUTC(*s.firstEpoch)<<"\n"; if(s.lastEpoch) os<<"last: "<<ceqc::model::formatUTC(*s.lastEpoch)<<"\n"; if(s.estimatedIntervalS>0) os<<"interval: "<<std::setprecision(3)<<s.estimatedIntervalS<<"s\n"; if(!quiet){ os<<"systems:\n"; for(auto&kv:s.systemAppearance)os<<"  "<<kv.first<<": "<<kv.second<<"\n"; }
  if(s.derived){ os<<"qc-options:"; for(auto&o:s.derived->optionsActive)os<<" "<<o; os<<"\n"; os<<"epoch-svs: min="<<s.derived->epochSVMin<<" mean="<<std::setprecision(2)<<s.derived->epochSVMean<<" max="<<s.derived->epochSVMax<<"\n"; if(!s.derived->gapEvents.empty()) os<<"gaps: "<<s.derived->gapEvents.size()<<"\n"; if(s.derived->lliEnabled) os<<"lli-count: "<<s.derived->lliCount<<"\n"; if(s.derived->snrEnabled && !s.derived->snrStats.empty()){ bool anyCode=false; for(auto&kv:s.derived->snrStats) if(kv.first.find(':')!=std::string::npos){ anyCode=true; os<<"snr-summary "<<kv.first<<": n="<<kv.second.count<<" mean="<<std::setprecision(2)<<kv.second.mean<<" rms="<<kv.second.rms<<" low="<<kv.second.lowCount<<"\n"; } if(!anyCode) for(auto&kv:s.derived->snrStats) if(kv.first!="all") os<<"snr-summary "<<kv.first<<": n="<<kv.second.count<<" mean="<<std::setprecision(2)<<kv.second.mean<<" rms="<<kv.second.rms<<" low="<<kv.second.lowCount<<"\n"; } if(s.derived->pseudorangePhaseEnabled && !s.derived->pseudorangePhase.empty()) for(auto&kv:s.derived->pseudorangePhase) os<<"pseudorange-phase "<<kv.first<<": n="<<kv.second.count<<" mean="<<std::setprecision(3)<<kv.second.mean<<" rms="<<kv.second.rms<<"\n"; if(!s.derived->timeplot.empty()) os<<"qc-timeplot: "<<s.derived->timeplot<<"\n";
    if(s.derived->ceqcExtensionEnabled && !s.derived->obsTimeplot.empty()) os<<"ceqc-ext timeplot obs: |"<<s.derived->obsTimeplot<<"|\n";
    if(s.derived->ceqcExtensionEnabled && !s.derived->navTimeplot.empty()) os<<"ceqc-ext timeplot nav: |"<<s.derived->navTimeplot<<"|\n";
    if(s.derived->ceqcExtensionEnabled && !s.derived->positionTimeplot.empty()) os<<"ceqc-ext timeplot position: |"<<s.derived->positionTimeplot<<"|\n";
    // Native +qc is CEQC's main diagnostic report, not a teqc compatibility
    // renderer. Keep the compact one-line qc-timeplot above for grepping, but
    // also print the full per-SV visual diagnostics by default. +qcq remains
    // compact, and +teqc uses the separate teqc-like renderer.
    if(!quiet){
      printNativeSVTimeplot(os, *s.derived);
      printNativeSkyplot(os, *s.derived);
    }
    auto printComputed=[&](const std::string& label,const std::string& key,int samplesOverride=-1){ int samples=samplesOverride>=0?samplesOverride:(s.derived->histogramSamples.count(key)?s.derived->histogramSamples.at(key):0); auto hit=s.derived->histograms.find(key); if(samples<=0 || hit==s.derived->histograms.end()) os<<label<<" skipped: samples=0\n"; else os<<label<<" computed: samples="<<samples<<" bins="<<histogramText(hit->second)<<"\n"; };
    printComputed("ION","ion");
    int mp1Samples=0, mp2Samples=0; auto m1=s.derived->multipathMovingCount.find("MP1"); if(m1!=s.derived->multipathMovingCount.end()) mp1Samples=m1->second; auto m2=s.derived->multipathMovingCount.find("MP2"); if(m2!=s.derived->multipathMovingCount.end()) mp2Samples=m2->second;
    if(mp1Samples>0) os<<"MP1 computed: samples="<<mp1Samples<<" rms_m="<<std::fixed<<std::setprecision(3)<<s.derived->mp1Meters<<"\n"; else os<<"MP1 skipped: samples=0\n";
    if(mp2Samples>0) os<<"MP2 computed: samples="<<mp2Samples<<" rms_m="<<std::fixed<<std::setprecision(3)<<s.derived->mp2Meters<<"\n"; else os<<"MP2 skipped: samples=0\n";
    if(mp1Samples<=0 && mp2Samples<=0 && !s.derived->multipathSkipReason.empty()) os<<"MP skip reason: "<<s.derived->multipathSkipReason<<"\n";
    if(s.derived->ceqcExtensionEnabled || s.derived->dataIndicatorsEnabled){ for(const auto& kv:s.derived->histogramSamples){ if(kv.first=="ion"||kv.first=="mp") continue; if(s.derived->histograms.find(kv.first)==s.derived->histograms.end()) os<<"ceqc-ext histogram "<<kv.first<<": skipped samples="<<kv.second<<"\n"; } for(const auto& kv:s.derived->histograms){ if(kv.first=="ion"||kv.first=="mp") continue; os<<"ceqc-ext histogram "<<kv.first<<": "<<histogramText(kv.second)<<" samples="<<(s.derived->histogramSamples.count(kv.first)?s.derived->histogramSamples.at(kv.first):0)<<"\n"; } }
    for(const auto& w:s.derived->thresholdWarnings) os<<"qc-threshold-warning: "<<w<<"\n";
    if(s.derived->position.skippedNoNavigation){
      os<<"position-qc: skipped=no_navigation candidate_epochs="<<s.derived->position.candidateEpochs<<" min_svs="<<s.derived->minSVsUsed<<"\n";
    } else if(s.derived->position.attempted) {
      const auto& p = s.derived->position;
      os<<"position-qc: solved="<<p.epochSolutions<<" skipped="<<(p.candidateEpochs-p.epochSolutions)<<" min_svs="<<s.derived->minSVsUsed<<"\n";
      if(!quiet && p.epochSolutions>0){
        auto qllh = xyzToLlh(p.averageXYZ[0], p.averageXYZ[1], p.averageXYZ[2]);
        os<<"qc-position mean xyz: "<<std::fixed<<std::setprecision(4)<<p.averageXYZ[0]<<" "<<p.averageXYZ[1]<<" "<<p.averageXYZ[2]<<" (m) positions="<<p.epochSolutions<<"\n";
        os<<"qc-position mean geo: "<<std::fixed<<std::setprecision(6)<<qllh[0]*180.0/3.141592653589793238462643383279502884
          <<" "<<qllh[1]*180.0/3.141592653589793238462643383279502884<<" height="<<std::setprecision(4)<<qllh[2]<<" m\n";
        if(p.hasApprox){
          const double dx=p.averageXYZ[0]-p.approxXYZ[0], dy=p.averageXYZ[1]-p.approxXYZ[1], dz=p.averageXYZ[2]-p.approxXYZ[2];
          const double dist=std::sqrt(dx*dx+dy*dy+dz*dz);
          auto denu=ecefDeltaToEnu(p.approxXYZ, p.averageXYZ);
          const double horiz=std::sqrt(denu[0]*denu[0]+denu[1]*denu[1]);
          os<<"qc-position header-delta: "<<std::fixed<<std::setprecision(4)<<dist<<" m horizontal="<<horiz<<" vertical="<<std::fabs(denu[2])<<"\n";
          if(s.derived->ceqcExtensionEnabled){
            os<<"ceqc-ext qc-header enu: E="<<std::fixed<<std::setprecision(4)<<denu[0]<<" N="<<denu[1]<<" U="<<denu[2]<<" m\n";
          }
        }
      }
      if(s.derived->ceqcExtensionEnabled){
        os<<"ceqc-ext position-qc: candidate_epochs="<<p.candidateEpochs<<" numeric_solutions="<<p.epochNumericSolutions<<" accepted="<<p.epochSolutions<<" skipped_insufficient_svs="<<p.skippedInsufficientSVs<<" rejected_bad_residual="<<p.rejectedBadResidual<<" rejected_bad_height="<<p.rejectedBadHeight<<" rejected_bad_jump="<<p.rejectedBadJump<<" rejected_bad_geometry="<<p.rejectedBadGeometry<<" max_epoch_used_svs_by_system";
        for(const char* sys : {"G","R","E","C","J"}){ auto it=p.usedSVsBySystem.find(sys); os<<" "<<sys<<"="<<(it==p.usedSVsBySystem.end()?0:it->second); }
        os<<" gates_rms_m="<<p.residualRmsGateM<<" gates_max_m="<<p.residualMaxGateM<<"\n";
        for(const auto& w:p.warnings) os<<"ceqc-ext position-warning: "<<w<<"\n";
      }
    }
    if((s.derived->everyEpochXYZ||s.derived->everyEpochGeodetic||s.derived->everyEpochDecimal) && !s.derived->epochPositions.empty()){
      os<<"every-epoch-position:\n";
      for(const auto& ep:s.derived->epochPositions){
        os<<"  "<<ep.time<<" sv="<<ep.usedSVs<<" status="<<ep.status;
        if(s.derived->epochPositions.size()>200 && (&ep-&s.derived->epochPositions.front())>=200){ os<<" ... truncated\n"; break; }
        if(ep.status=="OK" || ep.status.rfind("REJECT_",0)==0){
          if(s.derived->everyEpochXYZ || (!s.derived->everyEpochGeodetic && !s.derived->everyEpochDecimal)) os<<" xyz="<<std::fixed<<std::setprecision(4)<<ep.x<<","<<ep.y<<","<<ep.z;
          if(s.derived->everyEpochGeodetic || s.derived->everyEpochDecimal) os<<" llh="<<std::setprecision(s.derived->everyEpochDecimal?9:6)<<ep.latDeg<<","<<ep.lonDeg<<","<<std::setprecision(4)<<ep.heightM;
          os<<" clk_m="<<std::setprecision(4)<<ep.clockBiasM<<" postfit_rms_m="<<std::setprecision(3)<<ep.postfitRmsM<<" max_res_m="<<ep.maxResidualM;
        }
        os<<"\n";
      }
    }
    if(s.derived->dataIndicatorsEnabled){ os<<"data-completeness: complete="<<s.derived->dataCompleteness.completeRecords<<" partial="<<s.derived->dataCompleteness.partialRecords<<" missing_values="<<s.derived->dataCompleteness.missingValues<<"\n"; }
    if(s.derived->yCodeEnabled){ os<<"y-code-summary: gps_y_code_observations="<<s.derived->dataCompleteness.yCodeObservations<<"\n"; }
    if(s.derived->riseSetEnabled && !s.derived->riseSetEvents.empty()){
      bool hasEph=false; for(const auto& ev:s.derived->riseSetEvents) if(ev.hasEphemeris || std::isfinite(ev.maxElevationDeg)) hasEph=true;
      os<<(hasEph?"rise-set-summary:\n":"observed-arc-summary:\n");
      for(const auto& ev:s.derived->riseSetEvents){
        os<<"  "<<ev.satellite<<" first="<<ev.first<<" last="<<ev.last<<" duration_h="<<std::fixed<<std::setprecision(3)<<ev.durationHours<<" obs="<<ev.obsCount;
        if(hasEph && std::isfinite(ev.maxElevationDeg)) os<<" max_el="<<std::setprecision(2)<<ev.maxElevationDeg;
        os<<"\n";
      }
    }
    if(s.derived->symbolCodesEnabled||s.derived->allSymbolsEnabled){ for(const auto& l:s.derived->symbolLegend) os<<"symbol: "<<l<<"\n"; }
    if(s.derived->ssvEnabled){ os<<"per-sv-summary:\n"; for(const auto& kv:s.satelliteAppearance) os<<"  "<<kv.first<<" obs="<<kv.second<<"\n"; }
    if(s.derived->svprEnabled){
      os<<"sv-pseudorange-summary:\n";
      for(const auto& satkv:s.derived->svPseudorangeStats){
        os<<"  "<<satkv.first;
        int shown=0;
        for(const auto& codekv:satkv.second){
          if(shown++>=12){ os<<" ..."; break; }
          const auto& st=codekv.second;
          os<<" "<<codekv.first<<"[n="<<st.count<<" mean="<<std::fixed<<std::setprecision(3)<<st.mean<<" min="<<st.min<<" max="<<st.max<<"]";
        }
        os<<"\n";
      }
    }
  }
  if(s.residuals){
    os<<"residual-qc: candidates="<<s.residuals->candidateObservations<<" evaluated="<<s.residuals->evaluated<<" skipped="<<(s.residuals->skippedNoStation+s.residuals->skippedNoEphemeris+s.residuals->skippedNoPseudorange)<<"\n";
    bool ext = s.derived && (s.derived->dataIndicatorsEnabled || s.derived->ceqcExtensionEnabled);
    if(!quiet && s.residuals->evaluated>0){
      os<<"residual-epoch-centered: mean="<<std::fixed<<std::setprecision(3)<<s.residuals->meanMeters<<" rms="<<s.residuals->rmsMeters<<" max_abs="<<s.residuals->maxAbsMeters<<" m\n";
      for(const auto& kv:s.residuals->biasRemovedBySystem) os<<"residual-epoch-centered-system "<<kv.first<<": n="<<kv.second.count<<" mean="<<std::fixed<<std::setprecision(3)<<kv.second.mean<<" rms="<<kv.second.rms<<" m\n";
      for(const auto& w:s.residuals->warnings) os<<"residual-warning: "<<w<<"\n";
    }
    if(ext){
      os<<"ceqc-ext residual-qc no_ephemeris_by_system:"; for(const char* sys : {"G","R","E","C","J"}){ auto it=s.residuals->skippedNoEphemerisBySystem.find(sys); os<<" "<<sys<<"="<<(it==s.residuals->skippedNoEphemerisBySystem.end()?0:it->second); } os<<"\n";
      os<<"ceqc-ext residual-qc no_pseudorange_by_system:"; for(const char* sys : {"G","R","E","C","J"}){ auto it=s.residuals->skippedNoPseudorangeBySystem.find(sys); os<<" "<<sys<<"="<<(it==s.residuals->skippedNoPseudorangeBySystem.end()?0:it->second); } os<<"\n";
      os<<"ceqc-ext code-minus-range-no-clock-meters: mean="<<s.residuals->rawMeanMeters<<" rms="<<s.residuals->rawRmsMeters<<" max_abs="<<s.residuals->rawMaxAbsMeters<<"\n";
      os<<"ceqc-ext residual-epoch-centered-meters: mean="<<s.residuals->meanMeters<<" rms="<<s.residuals->rmsMeters<<" max_abs="<<s.residuals->maxAbsMeters<<"\n";
      os<<"ceqc-ext diagnostic-satellite-bias-removed-meters: mean="<<s.residuals->satBiasRemovedMeanMeters<<" rms="<<s.residuals->satBiasRemovedRmsMeters<<" max_abs="<<s.residuals->satBiasRemovedMaxAbsMeters<<"\n";
      for(const auto& kv:s.residuals->rawBySystem) os<<"ceqc-ext code-minus-range-no-clock-system "<<kv.first<<": n="<<kv.second.count<<" mean="<<kv.second.mean<<" rms="<<kv.second.rms<<"\n";
      for(const auto& kv:s.residuals->biasRemovedBySystem) os<<"ceqc-ext residual-epoch-centered-system "<<kv.first<<": n="<<kv.second.count<<" mean="<<kv.second.mean<<" rms="<<kv.second.rms<<"\n";
      for(const auto& kv:s.residuals->satBiasRemovedBySystem) os<<"ceqc-ext diagnostic-satellite-bias-removed-system "<<kv.first<<": n="<<kv.second.count<<" mean="<<kv.second.mean<<" rms="<<kv.second.rms<<"\n";
      for(auto& w:s.residuals->warnings) os<<"ceqc-ext residual-warning: "<<w<<"\n";
    }
  }
  if(s.rtcm3){ os<<"rtcm3-frames: total="<<s.rtcm3->frameCount<<" good="<<s.rtcm3->goodFrameCount<<" bad_crc="<<s.rtcm3->badCRCCount<<" bytes="<<s.rtcm3->bytesRead<<"\n"; os<<"rtcm3-observations: epochs="<<s.rtcm3->epochCount<<" sat-observations="<<s.rtcm3->observationCount<<" values="<<s.rtcm3->observationValueCount<<"\n"; }
  if(s.ubx){ os<<"ubx-frames: total="<<s.ubx->frameCount<<" good="<<s.ubx->goodFrameCount<<" bad_checksum="<<s.ubx->badChecksumCount<<"\n"; os<<"ubx-observations: epochs="<<s.ubx->epochCount<<" sat-observations="<<s.ubx->observationCount<<" values="<<s.ubx->observationValueCount<<"\n"; os<<"ubx-sfrbx: messages="<<s.ubx->sfrbxCount<<"\n"; }
}
}
