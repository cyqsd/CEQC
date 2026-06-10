#include "ceqc/rinex/RinexService.hpp"
#include "ceqc/qc/QCService.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <ctime>
#include <cctype>

namespace ceqc::service::rinex {
namespace {

WriterOptions g_writerOptions{};

bool rtklibCompat() { return g_writerOptions.rtklibCompat; }

int systemOrderKey(const std::string& sys) {
  if (!rtklibCompat()) return 100;
  if (sys == "G") return 0;
  if (sys == "R") return 1;
  if (sys == "C") return 2;
  if (sys == "E") return 3;
  if (sys == "J") return 4;
  if (sys == "S") return 5;
  if (sys == "I") return 6;
  return 50;
}

double normVer(double v, RinexKind k, double current = 0.0) {
  auto canonical = [&](double major) {
    if (major < 3) return 2.11;
    if (major < 4) return 3.05;
    return 4.02;
  };
  if (v >= 2 && v < 5) {
    double major = std::floor(v + 1e-9);
    if (std::fabs(v - major) < 1e-9) return canonical(v);
    return v; // explicit subversion, e.g. +v2.10, +v3.03, +v4.00
  }
  if (current >= 2 && current < 5) return current;
  return k == RinexKind::Nav ? 4.02 : 3.05;
}

std::string kindFirst(RinexKind k, double v) {
  v = normVer(v, k, 0.0);
  std::ostringstream os;
  os << std::fixed << std::setw(9) << std::setprecision(2) << v << "           ";
  if (k == RinexKind::Obs) {
    if (rtklibCompat() && v >= 3.0) os << std::left << std::setw(20) << "OBSERVATION DATA" << "M: Mixed            ";
    else os << std::left << std::setw(20) << "OBSERVATION DATA" << "M                   ";
  } else if (k == RinexKind::Nav) {
    if (v < 3.0) {
      // CEQC's RINEX 2 NAV writer emits GPS LNAV records only.  RINEX 2
      // does not have the RINEX 3/4 mixed-navigation system flag, so strict
      // readers such as Anubis expect the legacy GPS navigation file type to
      // be declared explicitly in columns 21--40.  Leaving this as the generic
      // "NAVIGATION DATA" makes Anubis warn "RINEXN system not defined,
      // used GPS" even though the content is valid GPS NAV.
      os << std::left << std::setw(20) << "N: GPS NAV DATA" << "                    ";
    } else if (rtklibCompat() && v >= 3.0) {
      // RTKLIB_EX/RTKPLOT-compatible NAV uses a RINEX 3.05 mixed-GNSS header
      // and the traditional untyped NAV body.  Do not emit RINEX-4 typed
      // markers in this branch.
      os << std::left << std::setw(20) << "N: GNSS NAV DATA" << "M: Mixed            ";
    } else {
      os << std::left << std::setw(20) << "NAVIGATION DATA" << "M                   ";
    }
  } else {
    os << std::left << std::setw(20) << "METEOROLOGICAL DATA" << "M                   ";
  }
  return os.str();
}


HeaderLine line(std::string v, std::string l);

std::string utcNowCompact() {
  std::time_t now = std::time(nullptr);
  std::tm tm{};
#if defined(_WIN32)
  gmtime_s(&tm, &now);
#else
  gmtime_r(&now, &tm);
#endif
  std::ostringstream os;
  os << std::setfill('0') << std::setw(4) << tm.tm_year + 1900
     << std::setw(2) << tm.tm_mon + 1 << std::setw(2) << tm.tm_mday
     << " " << std::setw(2) << tm.tm_hour << std::setw(2) << tm.tm_min
     << std::setw(2) << tm.tm_sec << " UTC";
  return os.str();
}

std::string rtklibPgmRunByDate() {
  std::ostringstream os;
  // RTKPlot compatibility mode intentionally mirrors RTKCONV-EX's header
  // fingerprint.  Several RTKLIB_EX GUI code paths are more permissive when
  // files look exactly like its own converter output.  This string is emitted
  // only for explicit +rtklib/+rtkplot projections; strict CEQC output keeps
  // the normal program label.
  os << std::left << std::setw(20) << "RTKCONV-EX 2.5.0"
     << std::setw(20) << ""
     << std::setw(20) << utcNowCompact();
  return os.str();
}

std::string lowerCopy(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  return s;
}

std::string rtklibInputFormat(const std::vector<RinexFile>& files) {
  for (const auto& f : files) {
    if (f.rtcm3) return "RTCM 3";
    if (f.ubx) return "UBX";
  }
  for (const auto& f : files) {
    std::string p = lowerCopy(f.path);
    if (p.find("rtcm") != std::string::npos || p.find("rtcm3") != std::string::npos) return "RTCM 3";
    if (p.find("ubx") != std::string::npos || p.rfind(".ubx") != std::string::npos) return "UBX";
  }
  return "RINEX";
}

std::string rtklibLogPath(const std::vector<RinexFile>& files) {
  for (const auto& f : files) {
    if (f.rtcm3 && !f.rtcm3->sourcePath.empty()) return f.rtcm3->sourcePath;
    if (f.ubx && !f.ubx->sourcePath.empty()) return f.ubx->sourcePath;
  }
  for (const auto& f : files) if (!f.path.empty()) return f.path;
  return "";
}

void addRtklibComments(RinexFile& out, const std::vector<RinexFile>& files) {
  if (!rtklibCompat() || out.header.version < 3.0) return;
  // RTKCONV-EX always emits COMMENT records immediately after the PGM line.
  // RTKPlot is surprisingly sensitive to this Windows/RTKCONV style in some
  // import paths, so the explicit +rtklib/+rtkplot projection emits these
  // records even when there is no meaningful external log text available.
  std::string fmt = "format: " + rtklibInputFormat(files);
  std::string log = "log: " + rtklibLogPath(files);
  out.header.lines.push_back(line(fmt, "COMMENT"));
  out.header.lines.push_back(line(log, "COMMENT"));
}

HeaderLine line(std::string v, std::string l) {
  if (v.size() > 60) v = v.substr(0, 60);
  std::ostringstream os;
  os << std::left << std::setw(60) << v << std::setw(20) << l;
  std::string raw = os.str();
  return {raw, raw.substr(0, 60), l};
}

std::string padMin80(std::string s) {
  // RINEX records are traditionally transported as 80-column records.  RINEX
  // 3/4 OBS data records may legitimately exceed 80 columns, but NAV records
  // and RINEX2 OBS records should still be at least 80 characters for strict
  // tools such as Anubis that index fixed columns. Padding never changes the
  // numeric content; it only adds trailing blanks.
  if (s.size() < 80) s.append(80 - s.size(), ' ');
  return s;
}

double secOf(const TimePoint& t) {
  auto base = std::chrono::time_point_cast<std::chrono::seconds>(t);
  auto frac = std::chrono::duration<double>(t - base).count();
  auto tm = toUTC(t);
  return static_cast<double>(tm.tm_sec) + frac;
}

std::string slot(const ObservationValue* v) {
  // RTKLIB_EX/RTKPLOT is much more tolerant of RTKCONV-style blank LLI/SSI
  // columns than of a carrier arc where every initial MSM lock state is marked
  // with LLI=1.  Keep strict/mainline RINEX LLI/SSI intact, but in the explicit
  // +rtklib/+rtkplot branch emit blank flag columns to match RTKCONV-EX golden
  // plotting files.  This is a compatibility projection only; the source
  // ObservationValue remains unchanged and normal output is unaffected.
  auto lliCol = [&](const ObservationValue* vv) -> std::string {
    if (rtklibCompat()) return " ";
    return (!vv || vv->lli.empty()) ? " " : vv->lli.substr(0, 1);
  };
  auto ssiCol = [&](const ObservationValue* vv) -> std::string {
    if (rtklibCompat()) return " ";
    return (!vv || vv->ssi.empty()) ? " " : vv->ssi.substr(0, 1);
  };
  if (!v || !v->value || std::isnan(*v->value)) {
    if (!v) return "                ";
    std::ostringstream os;
    os << std::setw(14) << "" << std::setw(1) << lliCol(v) << std::setw(1) << ssiCol(v);
    return os.str();
  }
  std::ostringstream os;
  double vv = std::round((*v->value) * 1000.0) / 1000.0;
  os << std::setw(14) << std::fixed << std::setprecision(3) << vv
     << std::setw(1) << lliCol(v)
     << std::setw(1) << ssiCol(v);
  return os.str();
}

bool keepMeta(const std::string& l) {
  static const std::vector<std::string> labels = {
      "MARKER NAME", "MARKER NUMBER", "MARKER TYPE", "OBSERVER / AGENCY", "REC # / TYPE / VERS",
      "ANT # / TYPE", "APPROX POSITION XYZ", "ANTENNA: DELTA H/E/N", "WAVELENGTH FACT L1/2",
      "TIME OF FIRST OBS", "TIME OF LAST OBS", "INTERVAL", "LEAP SECONDS", "IONOSPHERIC CORR",
      "TIME SYSTEM CORR", "GLONASS SLOT / FRQ #", "GLONASS COD/PHS/BIS", "COMMENT"};
  return std::find(labels.begin(), labels.end(), l) != labels.end();
}

std::vector<HeaderLine> copiedMeta(const std::vector<RinexFile>& files, RinexKind k) {
  std::vector<HeaderLine> out;
  std::set<std::string> seen;
  for (auto& f : files) {
    if (f.header.kind != k) continue;
    for (auto& h : f.header.lines) {
      if (keepMeta(h.label)) {
        std::string key = h.label + ":" + h.value;
        if (!seen.count(key)) {
          out.push_back(h);
          seen.insert(key);
        }
      }
    }
    if (!out.empty()) break;
  }
  return out;
}

std::string navSci(double v) {
  if (rtklibCompat()) {
    // RTKLIB/RTKCONV writes RINEX NAV numeric fields in the traditional
    // D19.12 form: leading sign/blank, decimal mantissa in [-1,1), D exponent
    // with two digits (e.g. " -.335284043103D-03"), not normalized C++
    // scientific notation ("-3.352840431030D-04").  RTKPLOT's parser is more
    // reliable with this exact style, so use it only in the explicit
    // +rtklib/+rtkplot compatibility branch.
    if (!std::isfinite(v) || v == 0.0) return "  .000000000000D+00";
    int exp10 = static_cast<int>(std::floor(std::log10(std::fabs(v)))) + 1;
    double mant = v / std::pow(10.0, exp10);
    // Round to 12 digits after the decimal point. If rounding reaches 1.0,
    // renormalize so the mantissa still fits the D19.12 field.
    double scale = 1e12;
    mant = std::round(mant * scale) / scale;
    if (std::fabs(mant) >= 1.0) { mant /= 10.0; ++exp10; }
    std::string sign = mant < 0.0 ? " -" : "  ";
    mant = std::fabs(mant);
    long long frac = static_cast<long long>(std::llround(mant * scale));
    if (frac >= static_cast<long long>(scale)) { frac /= 10; ++exp10; }
    std::ostringstream os;
    os << sign << "." << std::setw(12) << std::setfill('0') << frac
       << "D" << (exp10 < 0 ? '-' : '+') << std::setw(2) << std::setfill('0') << std::abs(exp10);
    return os.str();
  }
  std::ostringstream os;
  os << std::uppercase << std::scientific << std::setprecision(12) << std::setw(19) << v;
  return os.str();
}

double secondsOfWeekFromEpoch(const ceqc::model::TimePoint& tp, int y, int m, int d) {
  auto epoch0 = ceqc::model::makeUTC(y, m, d, 0, 0, 0.0);
  double seconds = std::chrono::duration<double>(tp - epoch0).count();
  double sow = std::fmod(seconds, 604800.0);
  if (sow < 0.0) sow += 604800.0;
  return sow;
}

int gpsWeekFromEpoch(const ceqc::model::TimePoint& tp) {
  auto gps0 = ceqc::model::makeUTC(1980, 1, 6, 0, 0, 0.0);
  double seconds = std::chrono::duration<double>(tp - gps0).count();
  return static_cast<int>(std::floor(seconds / 604800.0 + 1e-9));
}

double expandModuloWeek(double rawWeek, int modulo, const ceqc::model::TimePoint& epoch) {
  if (!std::isfinite(rawWeek) || rawWeek < 0 || modulo <= 0) return rawWeek;
  int ref = gpsWeekFromEpoch(epoch);
  int raw = static_cast<int>(std::llround(rawWeek));
  if (raw >= ref - modulo / 2 && raw <= ref + modulo / 2) return rawWeek;
  int k = static_cast<int>(std::llround(static_cast<double>(ref - raw) / static_cast<double>(modulo)));
  int full = raw + k * modulo;
  if (full < 0) full = raw;
  return static_cast<double>(full);
}

void normalizeNavValuesForOutput(const NavigationRecord& r, std::vector<double>& navVals) {
  if (!r.epoch) return;
  auto patchField = [&](const std::string& name, int modulo) {
    auto it = r.fields.find(name);
    if (it == r.fields.end()) return;
    size_t idx = it->second.index;
    if (idx >= navVals.size()) return;
    navVals[idx] = expandModuloWeek(navVals[idx], modulo, *r.epoch);
  };
  auto patchTransmissionTime = [&](double fallbackTow) {
    auto it = r.fields.find("TransmissionTime");
    if (it == r.fields.end()) return;
    size_t idx = it->second.index;
    if (idx >= navVals.size()) return;
    double cur = navVals[idx];
    // RTKLIB_EX/RTKPLOT compatibility must preserve the converter-style raw
    // transmission time even when it is outside the strict 0..604800 week
    // interval.  Example: some GPS LNAV records exported by RTKCONV-EX carry
    // TtoM=696759; replacing it with Toe makes RTKPlot fail to compute the same
    // satellite positions/sky plot.  Strict gfzrnx/Anubis output still repairs
    // absent/invalid TtoM below.
    if (rtklibCompat() && std::isfinite(cur) && cur > 0.0) return;
    bool bad = !std::isfinite(cur) || cur <= 0.0 || cur >= 604800.0;
    if (!bad) return;
    double candidate = std::numeric_limits<double>::quiet_NaN();
    auto preferField = [&](const std::string& name) {
      if (std::isfinite(candidate) && candidate > 0.0 && candidate < 604800.0) return;
      auto jt = r.fields.find(name);
      if (jt == r.fields.end()) return;
      double v = jt->second.value;
      if (std::isfinite(v) && v > 0.0 && v < 604800.0) candidate = v;
    };
    // RTCM GPS/QZSS 1019/1044 ephemeris does not reliably provide a RINEX
    // TtoM value in every stream.  When the decoded value is absent/zero,
    // derive it from the same record's time tags rather than writing a literal
    // zero and making strict validators patch it afterwards.
    preferField("TransmissionTimeRaw");
    preferField("Toe");
    preferField("Toc");
    if (!(std::isfinite(candidate) && candidate > 0.0 && candidate < 604800.0)) candidate = fallbackTow;
    if (std::isfinite(candidate) && candidate > 0.0 && candidate < 604800.0) navVals[idx] = candidate;
  };
  if (r.system == "G" || r.system == "J") {
    // RTCM 1019/1044 carry 10-bit GPS/QZSS week numbers.  RINEX 3/4 NAV
    // expects the continuous full GPS week.  gfzrnx rejects the raw modulo
    // value (e.g. 374 when the record epoch is GPS week 2422).
    patchField("Week", 1024);
    patchField("GPSWeek", 1024);
    patchTransmissionTime(secondsOfWeekFromEpoch(*r.epoch, 1980, 1, 6));
  } else if (r.system == "E") {
    // Galileo GST week in RTCM is 12-bit; keep it continuous when an epoch is known.
    patchField("Week", 4096);
    patchTransmissionTime(secondsOfWeekFromEpoch(*r.epoch, 1999, 8, 22));
  } else if (r.system == "C") {
    patchTransmissionTime(secondsOfWeekFromEpoch(*r.epoch, 2006, 1, 1));
  }
}

std::vector<std::string> formatNavRecord(const NavigationRecord& r, double version) {
  std::vector<std::string> out;
  std::string rt = r.recordType.empty() ? "EPH" : r.recordType;
  std::string mt = r.messageType.empty() ? (r.system == "R" ? "FDMA" : "LNAV") : r.messageType;
  std::string st = r.messageSubtype.empty() ? "BNK" : r.messageSubtype;
  // RINEX 4 system data records (STO/EOP/ION) are not satellite ephemerides.
  // Preserve the typed marker and write numeric payload in stable 4-value rows.
  if (rt != "EPH") {
    if (version < 4.0) return out;
    std::string marker = "> " + rt + " " + mt;
    if (!st.empty() && st != "BNK") marker += " " + st;
    out.push_back(marker);
    for (size_t i = 0; i < r.values.size(); i += 4) {
      std::ostringstream ln;
      ln << "    ";
      for (size_t j = 0; j < 4 && i + j < r.values.size(); ++j) ln << navSci(r.values[i + j]);
      out.push_back(ln.str());
    }
    if (r.values.empty() && r.rawLines.size() > 1) out.insert(out.end(), r.rawLines.begin()+1, r.rawLines.end());
    for (auto& navLine : out) navLine = padMin80(navLine);
    return out;
  }
  if (r.satellite.empty() || !r.epoch) {
    if (version < 4.0) return out;
    out.insert(out.end(), r.rawLines.begin(), r.rawLines.end());
    for (auto& navLine : out) navLine = padMin80(navLine);
    return out;
  }
  if (version < 3.0 && r.system != "G") return out;
  std::vector<double> navVals = r.values;
  normalizeNavValuesForOutput(r, navVals);
  // Anubis 3.11 evaluates GLONASS ephemeris usability from the extended
  // FDMA numeric record.  If RINEX3 writes only the legacy 15-value set while
  // RINEX4 writes 19 values, the same observations produce different GLO
  // expected counts.  Keep RINEX3/4 same-source QC comparable by writing the
  // extended 19-value GLONASS FDMA payload for RINEX3+; the RTCM1020-unavailable
  // tail terms are neutral zeros and are not sample-specific.
  if (version >= 3.0 && r.system == "R" && (mt == "FDMA" || mt.empty()) && navVals.size() < 19) {
    navVals.resize(19, 0.0);
  }
  if (version >= 4.0 && !rtklibCompat()) {
    // RINEX 4 NAV message subtype/source must be a value understood by strict
    // readers. Diagnostic source tags such as RTCM1042, UBX-SFRBX, RAW or
    // SFRBX-ASSEMBLED are useful internally but make gfzrnx map the message to
    // XXXX and then reject the following numeric rows.  Keep real RINEX source
    // identifiers when they are already standards-like (e.g. "0" or tests
    // that explicitly round-trip a token), but collapse CEQC translator source
    // annotations to the neutral source "0".
    std::string outSubtype = st;
    auto startsWith = [](const std::string& x, const std::string& p) {
      return x.rfind(p, 0) == 0;
    };
    if (outSubtype.empty() || outSubtype == "BNK" || outSubtype == "SFRBX" ||
        outSubtype == "SFRBX-ASSEMBLED" || outSubtype == "RAW" ||
        startsWith(outSubtype, "RTCM") || startsWith(outSubtype, "UBX") ||
        startsWith(outSubtype, "CEQC")) {
      outSubtype = "0";
    }
    out.push_back("> EPH " + r.satellite + " " + mt + " " + outSubtype);
  }
  auto tm = toUTC(*r.epoch);
  auto val = [&](size_t i) { return i < navVals.size() ? navVals[i] : 0.0; };
  std::ostringstream first;
  if (version < 3.0) {
    std::string prn = r.satellite.size() > 1 ? r.satellite.substr(1) : r.satellite;
    first << std::setw(2) << prn << " " << std::setw(2) << (tm.tm_year + 1900) % 100
          << " " << std::setw(2) << tm.tm_mon + 1 << " " << std::setw(2) << tm.tm_mday
          << " " << std::setw(2) << tm.tm_hour << " " << std::setw(2) << tm.tm_min << " " << std::setw(4) << tm.tm_sec;
  } else {
    first << std::left << std::setw(3) << r.satellite << std::right << " " << std::setw(4) << tm.tm_year + 1900
          << " " << std::setfill('0') << std::setw(2) << tm.tm_mon + 1 << " " << std::setw(2) << tm.tm_mday
          << " " << std::setw(2) << tm.tm_hour << " " << std::setw(2) << tm.tm_min << " " << std::setw(2) << tm.tm_sec
          << std::setfill(' ');
  }
  first << navSci(val(0)) << navSci(val(1)) << navSci(val(2));
  out.push_back(first.str());
  for (size_t i = 3; i < navVals.size(); i += 4) {
    std::ostringstream ln;
    // RINEX 2 NAV continuation lines are 3 blanks + 4*D19.12.
    // RINEX 3/4 NAV continuation lines are 4 blanks + 4*D19.12.
    // Using 4 blanks for RINEX 2 shifts every field by one column; gfzrnx then
    // reads the GPS week field as e.g. "0 2.422..." and rejects it.
    ln << (version < 3.0 ? "   " : "    ");
    for (size_t j = 0; j < 4 && i + j < navVals.size(); ++j) ln << navSci(navVals[i + j]);
    out.push_back(ln.str());
  }
  if (!rtklibCompat()) {
    for (auto& navLine : out) navLine = padMin80(navLine);
  }
  return out;
}

std::map<std::string, std::vector<ObservationRecord>> buckets(const std::vector<ObservationRecord>& rs) {
  std::map<std::string, std::vector<ObservationRecord>> b;
  for (auto& r : rs) b[formatUTC(r.time)].push_back(r);
  return b;
}

std::string rinex2Code(std::string c) {
  if (c.size() >= 2) return c.substr(0, 2);
  return c;
}

std::string outputObsCode(const std::string& sys, std::string code) {
  if (!rtklibCompat() || code.size() != 3) return code;
  // RTKLIB_EX 2.5.0 uses its own legacy plotting-friendly signal names for
  // several RTCM MSM signals.  Keep this mapping only in +rtklib/+rtkplot so
  // the strict gfzrnx/Anubis branch keeps the original CEQC/RINEX mapping.
  if (sys == "G" && code[1] == '2' && code[2] == 'L') code[2] = 'W';
  if (sys == "J" && code[1] == '2' && code[2] == 'L') code[2] = 'W';
  if (sys == "C" && code[1] == '7' && code[2] == 'X') code[2] = 'I';
  return code;
}

int rtklibSignalRank(const std::string& sys, const std::string& code) {
  if (code.size() != 3) return 1000;
  const std::string sig = code.substr(1, 2);
  auto rankIn = [&](const std::vector<std::string>& order) {
    auto it = std::find(order.begin(), order.end(), sig);
    return it == order.end() ? 900 : static_cast<int>(it - order.begin());
  };
  if (sys == "G" || sys == "J") return rankIn({"1C", "2W", "5Q", "1W", "2L", "5X"});
  if (sys == "R") return rankIn({"1C", "2C", "1P", "2P"});
  if (sys == "C") return rankIn({"2I", "7I", "6I", "1P", "5P", "7X"});
  if (sys == "E") return rankIn({"1X", "5X", "7X", "8X", "6X", "1C"});
  return rankIn({code.substr(1, 2)});
}

int obsKindRank(char c) {
  if (c == 'C' || c == 'P') return 0;
  if (c == 'L') return 1;
  if (c == 'D') return 2;
  if (c == 'S') return 3;
  return 9;
}

void sortObsTypes(std::vector<std::string>& v, const std::string& sys) {
  std::sort(v.begin(), v.end(), [&](const auto& a, const auto& b) {
    if (rtklibCompat() && a.size() == 3 && b.size() == 3) {
      int fa = rtklibSignalRank(sys, a), fb = rtklibSignalRank(sys, b);
      if (fa != fb) return fa < fb;
      int ka = obsKindRank(a[0]), kb = obsKindRank(b[0]);
      if (ka != kb) return ka < kb;
      return a < b;
    }
    int oa = obsKindRank(a.empty() ? '\0' : a[0]);
    int ob = obsKindRank(b.empty() ? '\0' : b[0]);
    return oa == ob ? a < b : oa < ob;
  });
}

bool hasFiniteObservationValue(const ObservationValue& val) {
  return val.value && std::isfinite(*val.value);
}

std::map<std::string, std::vector<std::string>> observationTypes(const std::vector<ObservationRecord>& recs, bool v2) {
  std::map<std::string, std::vector<std::string>> types;
  for (auto& r : recs) {
    auto sys = r.system.empty() ? "G" : r.system;
    auto& tv = types[v2 ? "G" : sys];
    for (auto& val : r.values) {
      // Declare only observation types for which at least one numeric value is
      // actually written.  Some translators preserve LLI/SSI state on an empty
      // carrier slot; that is useful inside an existing declared type, but it
      // must not create a new SYS / # / OBS TYPES entry by itself.  This keeps
      // RINEX4 OBS TYPES aligned with RINEX3 for the same source observations
      // and prevents validators from expecting data that CEQC never outputs.
      if (!hasFiniteObservationValue(val)) continue;
      auto mapped = outputObsCode(sys, val.type);
      auto code = v2 ? rinex2Code(mapped) : mapped;
      if (std::find(tv.begin(), tv.end(), code) == tv.end()) tv.push_back(code);
    }
  }
  for (auto& [sys, tv] : types) sortObsTypes(tv, sys);
  return types;
}

std::map<std::string, const ObservationValue*> valueMap(const ObservationRecord& r, bool v2) {
  std::map<std::string, const ObservationValue*> m;
  for (auto& v : r.values) {
    auto mapped = outputObsCode(r.system.empty() ? "G" : r.system, v.type);
    auto code = v2 ? rinex2Code(mapped) : mapped;
    if (!m.count(code)) m[code] = &v;
  }
  return m;
}

std::vector<std::string> orderedSystems(const std::map<std::string, std::vector<std::string>>& types) {
  std::vector<std::string> keys;
  for (auto& kv : types) keys.push_back(kv.first);
  if (rtklibCompat()) {
    std::sort(keys.begin(), keys.end(), [](const auto& a, const auto& b) {
      int oa = systemOrderKey(a), ob = systemOrderKey(b);
      return oa == ob ? a < b : oa < ob;
    });
  }
  return keys;
}

void addV3ObsTypes(RinexFile& out, const std::map<std::string, std::vector<std::string>>& types) {
  for (auto& sys : orderedSystems(types)) {
    const auto& tv = types.at(sys);
    for (size_t i = 0; i < tv.size(); i += 13) {
      std::ostringstream v;
      if (i == 0) {
        if (rtklibCompat()) v << sys << std::setw(5) << tv.size(); // RTKCONV: "G    9"
        else v << sys << " " << std::setw(3) << tv.size();
      } else {
        v << "      ";
      }
      for (size_t j = i; j < tv.size() && j < i + 13; ++j) v << " " << tv[j];
      out.header.lines.push_back(line(v.str(), "SYS / # / OBS TYPES"));
    }
  }
}

void addV3PhaseShifts(RinexFile& out, const std::map<std::string, std::vector<std::string>>& types) {
  for (auto& sys : orderedSystems(types)) {
    const auto& tv = types.at(sys);
    for (auto& code : tv) {
      if (code.size() == 3 && code[0] == 'L') {
        std::ostringstream v;
        // Strict RINEX uses an explicit zero phase shift.  RTKCONV-EX leaves
        // most phase-shift records blank and writes 0.00000 only for GPS L2W/L5Q;
        // match that style in +rtklib/+rtkplot because RTKPlot imports these
        // files through its own converter-oriented path.
        v << sys << " " << code;
        if (rtklibCompat()) {
          if (sys == "G" && (code == "L2W" || code == "L5Q")) {
            v << "  " << std::setw(7) << std::fixed << std::setprecision(5) << 0.0;
          }
        } else {
          v << " " << std::setw(8) << std::fixed << std::setprecision(5) << 0.0 << " " << std::setw(2) << 0;
        }
        out.header.lines.push_back(line(v.str(), "SYS / PHASE SHIFT"));
      }
    }
  }
}

void addV2ObsTypes(RinexFile& out, const std::vector<std::string>& tv) {
  for (size_t i = 0; i < tv.size(); i += 9) {
    std::ostringstream v;
    if (i == 0) v << std::setw(6) << tv.size();
    else v << "      ";
    for (size_t j = i; j < tv.size() && j < i + 9; ++j) v << std::setw(6) << tv[j];
    out.header.lines.push_back(line(v.str(), "# / TYPES OF OBSERV"));
  }
}

std::string rinexTimeHeader(const TimePoint& t, const std::string& scale) {
  auto tm = toUTC(t);
  std::ostringstream os;
  if (rtklibCompat()) {
    // RTKCONV-EX pads calendar fields with zero-filled two-digit month/day/hour/min.
    os << "  " << std::setw(4) << tm.tm_year + 1900
       << "    " << std::setfill('0') << std::setw(2) << tm.tm_mon + 1
       << "    " << std::setw(2) << tm.tm_mday
       << "    " << std::setw(2) << tm.tm_hour
       << "    " << std::setw(2) << tm.tm_min
       << std::setfill(' ') << std::setw(13) << std::fixed << std::setprecision(7) << secOf(t)
       << "     " << scale;
  } else {
    os << "  " << std::setw(4) << tm.tm_year + 1900 << std::setw(6) << tm.tm_mon + 1
       << std::setw(6) << tm.tm_mday << std::setw(6) << tm.tm_hour << std::setw(6) << tm.tm_min
       << std::setw(13) << std::fixed << std::setprecision(7) << secOf(t) << "     " << scale;
  }
  return os.str();
}

bool hasHeaderLabel(const RinexFile& out, const std::string& label);

double estimateObsIntervalSeconds(const std::vector<ObservationRecord>& recs) {
  // teqc reports a nominal observation interval derived from the actual epoch
  // sequence.  Do the same for generated/merged OBS output: collect unique
  // epoch times, build a histogram of positive adjacent deltas, and use the
  // dominant delta.  This deliberately avoids fixed values such as 15 s and
  // also avoids blindly preserving a stale input INTERVAL after windowing or
  // decimation.
  if (recs.size() < 2) return 0.0;
  std::vector<TimePoint> ts;
  ts.reserve(recs.size());
  for (const auto& r : recs) ts.push_back(r.time);
  std::sort(ts.begin(), ts.end());
  ts.erase(std::unique(ts.begin(), ts.end()), ts.end());
  if (ts.size() < 2) return 0.0;

  struct Bin { int count = 0; double total = 0.0; };
  std::map<long long, Bin> hist;
  for (size_t i = 1; i < ts.size(); ++i) {
    double dt = std::chrono::duration<double>(ts[i] - ts[i-1]).count();
    if (!(dt > 0.0) || dt >= 86400.0) continue;
    // Millisecond bins are fine for RINEX header INTERVAL precision and keep
    // epochs such as 4.994/9.994/14.994 grouped as exactly 5.000 s.
    long long key = static_cast<long long>(std::llround(dt * 1000.0));
    if (key <= 0) continue;
    hist[key].count += 1;
    hist[key].total += dt;
  }
  if (hist.empty()) return 0.0;

  auto best = hist.begin();
  for (auto it = hist.begin(); it != hist.end(); ++it) {
    if (it->second.count > best->second.count ||
        (it->second.count == best->second.count && it->first < best->first)) {
      best = it;
    }
  }
  return best->second.total / static_cast<double>(best->second.count);
}

void addObsTimeHeaders(RinexFile& out, const std::vector<ObservationRecord>& recs) {
  if (recs.empty()) return;
  auto first = recs.front().time;
  auto last = recs.front().time;
  for (auto& r : recs) { if (r.time < first) first = r.time; if (r.time > last) last = r.time; }
  out.header.lines.push_back(line(rinexTimeHeader(first, "GPS"), "TIME OF FIRST OBS"));
  out.header.lines.push_back(line(rinexTimeHeader(last, "GPS"), "TIME OF LAST OBS"));
  if (!rtklibCompat() && !hasHeaderLabel(out, "INTERVAL")) {
    double interval = estimateObsIntervalSeconds(recs);
    if (interval > 0.0) {
      std::ostringstream v;
      v << std::fixed << std::setprecision(3) << std::setw(10) << interval;
      out.header.lines.push_back(line(v.str(), "INTERVAL"));
    }
  }
}

bool hasHeaderLabel(const RinexFile& out, const std::string& label) {
  for (auto& h : out.header.lines) if (h.label == label) return true;
  return false;
}

void setHeaderLine(RinexFile& out, const std::string& label, const std::string& value) {
  HeaderLine h = line(value, label);
  for (auto& e : out.header.lines) {
    if (e.label == label) { e = h; return; }
  }
  out.header.lines.push_back(h);
}

void removeHeaderLabel(RinexFile& out, const std::string& label) {
  out.header.lines.erase(std::remove_if(out.header.lines.begin(), out.header.lines.end(), [&](const HeaderLine& h) {
    return h.label == label;
  }), out.header.lines.end());
}

void applyRtklibObsMetadata(RinexFile& out) {
  if (!rtklibCompat() || out.header.kind != RinexKind::Obs || out.header.version < 3.0) return;
  // RTKCONV-EX writes a station-id style MARKER NAME and blank MARKER NUMBER /
  // MARKER TYPE / OBSERVER fields when the RTCM stream has no human station
  // metadata.  Do this only in +rtklib/+rtkplot mode so strict RINEX output can
  // keep CEQC's normal metadata defaults or user-supplied header edits.
  setHeaderLine(out, "MARKER NAME", "0000");
  setHeaderLine(out, "MARKER NUMBER", "");
  setHeaderLine(out, "MARKER TYPE", "");
  setHeaderLine(out, "OBSERVER / AGENCY", "");
}

void reorderRtklibObsHeader(RinexFile& out) {
  if (!rtklibCompat() || out.header.kind != RinexKind::Obs || out.header.version < 3.0) return;
  static const std::vector<std::string> order = {
    "RINEX VERSION / TYPE", "PGM / RUN BY / DATE", "COMMENT",
    "MARKER NAME", "MARKER NUMBER", "MARKER TYPE", "OBSERVER / AGENCY",
    "REC # / TYPE / VERS", "ANT # / TYPE", "APPROX POSITION XYZ", "ANTENNA: DELTA H/E/N",
    "SYS / # / OBS TYPES", "TIME OF FIRST OBS", "TIME OF LAST OBS",
    "SYS / PHASE SHIFT", "GLONASS SLOT / FRQ #", "GLONASS COD/PHS/BIS", "INTERVAL"
  };
  std::vector<HeaderLine> old = out.header.lines;
  std::vector<HeaderLine> nw;
  std::vector<char> used(old.size(), 0);
  for (const auto& label : order) {
    for (size_t i = 0; i < old.size(); ++i) {
      if (!used[i] && old[i].label == label) { nw.push_back(old[i]); used[i] = 1; }
    }
  }
  for (size_t i = 0; i < old.size(); ++i) if (!used[i] && old[i].label != "END OF HEADER") nw.push_back(old[i]);
  out.header.lines.swap(nw);
}

bool isZeroApproxHeader(const RinexFile& out) {
  for (auto& h : out.header.lines) {
    if (h.label != "APPROX POSITION XYZ") continue;
    std::istringstream is(h.value); double x=0,y=0,z=0;
    if(is>>x>>y>>z) return std::sqrt(x*x+y*y+z*z) < 1.0e6;
  }
  return true;
}

void maybeFillApproxFromNav(RinexFile& out, const std::vector<RinexFile>& files) {
  if (out.header.kind != RinexKind::Obs) return;
  if (!isZeroApproxHeader(out)) return;
  std::vector<NavigationRecord> navs;
  for (auto& f : files) if (f.header.kind == RinexKind::Nav) navs.insert(navs.end(), f.data.navigationRecords.begin(), f.data.navigationRecords.end());
  if (navs.empty()) return;
  if (auto xyz = ceqc::service::qc::estimateApproxPosition(out, navs)) {
    std::ostringstream v;
    v << std::fixed << std::setprecision(4) << std::setw(14) << (*xyz)[0] << std::setw(14) << (*xyz)[1] << std::setw(14) << (*xyz)[2];
    setHeaderLine(out, "APPROX POSITION XYZ", v.str());
  }
}

void addRequiredObsDefaults(RinexFile& out) {
  if (!hasHeaderLabel(out, "MARKER NAME")) out.header.lines.push_back(line(rtklibCompat() ? "0000" : "CEQC", "MARKER NAME"));
  if (!hasHeaderLabel(out, "MARKER NUMBER")) out.header.lines.push_back(line(rtklibCompat() ? "" : "0000", "MARKER NUMBER"));
  if (rtklibCompat() && !hasHeaderLabel(out, "MARKER TYPE")) out.header.lines.push_back(line("", "MARKER TYPE"));
  if (!hasHeaderLabel(out, "OBSERVER / AGENCY")) out.header.lines.push_back(line(rtklibCompat() ? "" : "UNKNOWN             CEQC", "OBSERVER / AGENCY"));
  if (!hasHeaderLabel(out, "REC # / TYPE / VERS")) out.header.lines.push_back(line("UNKNOWN             UNKNOWN             UNKNOWN", "REC # / TYPE / VERS"));
  if (!hasHeaderLabel(out, "ANT # / TYPE")) out.header.lines.push_back(line("UNKNOWN             UNKNOWN", "ANT # / TYPE"));
  if (!hasHeaderLabel(out, "APPROX POSITION XYZ")) out.header.lines.push_back(line("        0.0000        0.0000        0.0000", "APPROX POSITION XYZ"));
  if (!hasHeaderLabel(out, "ANTENNA: DELTA H/E/N")) out.header.lines.push_back(line("        0.0000        0.0000        0.0000", "ANTENNA: DELTA H/E/N"));
}

std::optional<double> navField(const NavigationRecord& n, const std::string& name) {
  auto it = n.fields.find(name);
  if (it == n.fields.end()) return std::nullopt;
  return it->second.value;
}

void addGlonassSlotAndBiasHeaders(RinexFile& out, const std::vector<RinexFile>& files) {
  if (out.header.kind != RinexKind::Obs || out.header.version < 3.0) return;
  std::map<std::string,int> slots;
  for (const auto& f : files) {
    if (f.header.kind != RinexKind::Nav) continue;
    for (const auto& n : f.data.navigationRecords) {
      if (n.satellite.size() < 2 || n.satellite[0] != 'R') continue;
      auto k = navField(n, "FrequencyNumber");
      if (!k || !std::isfinite(*k)) continue;
      if (!slots.count(n.satellite)) slots[n.satellite] = static_cast<int>(std::lround(*k));
    }
  }
  if (!slots.empty() && !hasHeaderLabel(out, "GLONASS SLOT / FRQ #")) {
    std::vector<std::pair<std::string,int>> pairs(slots.begin(), slots.end());
    for (size_t start = 0; start < pairs.size(); start += 8) {
      std::ostringstream v;
      if (start == 0) v << std::setw(3) << pairs.size(); else v << "   ";
      for (size_t i = start; i < pairs.size() && i < start + 8; ++i) {
        v << ' ' << std::setw(3) << pairs[i].first << std::setw(3) << pairs[i].second;
      }
      out.header.lines.push_back(line(v.str(), "GLONASS SLOT / FRQ #"));
    }
  }
  bool hasGloObs = out.header.obsTypes.count("R") || std::any_of(out.data.observationRecords.begin(), out.data.observationRecords.end(), [](const ObservationRecord& r){ return r.system == "R"; });
  if (hasGloObs && !hasHeaderLabel(out, "GLONASS COD/PHS/BIS")) {
    // No receiver-specific GLONASS code/phase bias is known from RTCM/UBX here.
    // Writing explicit zero biases keeps RINEX3/4 headers complete for Anubis
    // without inventing sample-specific corrections.
    out.header.lines.push_back(line(" C1C    0.000 C1P    0.000 C2C    0.000 C2P    0.000", "GLONASS COD/PHS/BIS"));
  }
}

void indexHeader(RinexHeader& h) {
  h.byLabel.clear();
  for (size_t i = 0; i < h.lines.size(); ++i) h.byLabel[h.lines[i].label].push_back(i);
}

bool hasMappedValues(const ObservationRecord& r, const std::vector<std::string>& types, bool v2);

void makeV3Body(RinexFile& out, const std::vector<ObservationRecord>& recs, const std::map<std::string, std::vector<std::string>>& types) {
  auto b = buckets(recs);
  for (auto& [_, rs] : b) {
    std::vector<const ObservationRecord*> writable;
    writable.reserve(rs.size());
    for (auto& r : rs) {
      auto it = types.find(r.system);
      if (it == types.end() || it->second.empty()) continue;
      // The epoch satellite count must match the records actually written.
      // If an edit/format projection drops all observables for a satellite,
      // do not declare it in the epoch header; gfzrnx otherwise reports
      // mismatch of number of satellites(N) and records got(M).
      if (!hasMappedValues(r, it->second, false)) continue;
      writable.push_back(&r);
    }
    if (writable.empty()) continue;
    if (rtklibCompat()) {
      std::stable_sort(writable.begin(), writable.end(), [](const ObservationRecord* a, const ObservationRecord* b) {
        int oa = systemOrderKey(a->system), ob = systemOrderKey(b->system);
        if (oa != ob) return oa < ob;
        return a->satellite < b->satellite;
      });
    }
    auto tm = toUTC(writable.front()->time);
    std::ostringstream e;
    e << "> " << std::setw(4) << tm.tm_year + 1900 << " " << std::setfill('0') << std::setw(2) << tm.tm_mon + 1
      << " " << std::setw(2) << tm.tm_mday << " " << std::setw(2) << tm.tm_hour << " " << std::setw(2) << tm.tm_min
      << std::setfill(' ') << " " << std::setw(10) << std::fixed << std::setprecision(7) << secOf(writable.front()->time)
      // RINEX 3/4 epoch flag is one digit and number-of-satellites is the
      // following I3 field.  Do not insert an extra blank between them:
      // gfzrnx reads columns 33-35 for the I3 value and otherwise reads only
      // the tens digit (e.g. 31 -> 3), causing satellite-count mismatch.
      << "  0" << std::setw(3) << writable.size();
    out.body.push_back(e.str());
    for (auto* rp : writable) {
      const auto& r = *rp;
      auto it = types.find(r.system);
      std::ostringstream l;
      l << std::left << std::setw(3) << r.satellite;
      auto vm = valueMap(r, false);
      for (auto& t : it->second) l << slot(vm.count(t) ? vm[t] : nullptr);
      // RINEX 3/4 observation records are one logical record per satellite.
      // Do not wrap at 80 characters: that splits the 16-character observable
      // fields and makes strict validators treat continuation fragments such as
      // "   0" as new satellite records with system '0'.  Long body records
      // are accepted by RINEX 3/4 tools including gfzrnx.
      out.body.push_back(l.str());
    }
  }
}


bool hasMappedValues(const ObservationRecord& r, const std::vector<std::string>& types, bool v2) {
  auto vm = valueMap(r, v2);
  for (auto& t : types) {
    auto it = vm.find(t);
    if (it != vm.end() && it->second && it->second->value && std::isfinite(*(it->second->value))) return true;
  }
  return false;
}

void makeV2Body(RinexFile& out, const std::vector<ObservationRecord>& recs, const std::vector<std::string>& types) {
  auto b = buckets(recs);
  for (auto& [_, rs0] : b) {
    std::vector<ObservationRecord> rs;
    rs.reserve(rs0.size());
    for (auto& r : rs0) if (hasMappedValues(r, types, true)) rs.push_back(r);
    if (rs.empty()) continue;
    auto tm = toUTC(rs[0].time);
    for (size_t start = 0; start < rs.size(); start += 999) {
      size_t end = std::min(rs.size(), start + static_cast<size_t>(999));
      std::ostringstream e;
      e << " " << std::setw(2) << std::setfill('0') << (tm.tm_year + 1900) % 100 << std::setfill(' ')
        << std::setw(3) << tm.tm_mon + 1 << std::setw(3) << tm.tm_mday << std::setw(3) << tm.tm_hour
        << std::setw(3) << tm.tm_min << std::setw(11) << std::fixed << std::setprecision(7) << secOf(rs[0].time)
        << std::setw(3) << 0 << std::setw(3) << (end - start);
      std::string eline = e.str();
      size_t countInLine = 0;
      for (size_t i = start; i < end; ++i) {
        if (countInLine == 12) {
          out.body.push_back(padMin80(eline));
          eline = "                                ";
          countInLine = 0;
        }
        eline += rs[i].satellite.substr(0, 3);
        ++countInLine;
      }
      out.body.push_back(padMin80(eline));
      int linesPerSat = static_cast<int>((types.size() + 4) / 5);
      (void)linesPerSat;
      for (size_t i = start; i < end; ++i) {
        auto vm = valueMap(rs[i], true);
        std::string linePart;
        int n = 0;
        for (auto& t : types) {
          linePart += slot(vm.count(t) ? vm[t] : nullptr);
          if (++n == 5) {
            out.body.push_back(padMin80(linePart));
            linePart.clear();
            n = 0;
          }
        }
        if (!linePart.empty()) out.body.push_back(padMin80(linePart));
      }
    }
  }
}

} // namespace

RinexFile windowObservation(const RinexFile& rf, const std::optional<TimePoint>& start, const std::optional<TimePoint>& end) {
  if (rf.header.kind != RinexKind::Obs || (!start && !end)) return rf;
  RinexFile cp = rf;
  cp.data.observationRecords.clear();
  for (auto& r : rf.data.observationRecords) {
    if (start && r.time < *start) continue;
    if (end && r.time > *end) continue;
    cp.data.observationRecords.push_back(r);
  }
  return merge({cp}, RinexKind::Obs, rf.header.version);
}

RinexFile decimate(const RinexFile& rf, const DecimationSpec& spec) {
  if (!spec.enabled) return rf;
  auto iv = spec.interval.count();
  auto off = spec.offset.count();
  if (iv <= 0) return rf;
  RinexFile cp = rf;
  auto keepTime = [&](const TimePoint& t) {
    auto s = std::chrono::duration_cast<std::chrono::seconds>(t.time_since_epoch()).count();
    return ((s - off) % iv + iv) % iv == 0;
  };
  if (rf.header.kind == RinexKind::Obs) {
    cp.data.observationRecords.clear();
    for (auto& r : rf.data.observationRecords) if (keepTime(r.time)) cp.data.observationRecords.push_back(r);
    return merge({cp}, RinexKind::Obs, rf.header.version);
  }
  if (rf.header.kind == RinexKind::Nav) {
    cp.data.navigationRecords.clear();
    for (auto& r : rf.data.navigationRecords) if (!r.epoch || keepTime(*r.epoch)) cp.data.navigationRecords.push_back(r);
    return merge({cp}, RinexKind::Nav, rf.header.version);
  }
  if (rf.header.kind == RinexKind::Met) {
    cp.data.meteorologicalRecords.clear();
    cp.body.clear();
    for (auto& r : rf.data.meteorologicalRecords) if (keepTime(r.time)) { cp.data.meteorologicalRecords.push_back(r); cp.body.push_back(r.rawLine); }
    return merge({cp}, RinexKind::Met, rf.header.version);
  }
  return rf;
}

RinexFile merge(const std::vector<RinexFile>& files, RinexKind kind, double targetVersion) {
  RinexFile out;
  for (const auto& f : files) { if (f.header.kind == kind && !f.path.empty()) { out.path = f.path; break; } }
  out.header.kind = kind;
  double currentVersion = 0.0; for (const auto& f : files) if (f.header.kind == kind && f.header.version >= 2.0) { currentVersion = f.header.version; break; }
  out.header.version = normVer(targetVersion, kind, currentVersion);
  const bool v2obs = kind == RinexKind::Obs && out.header.version < 3.0;
  out.header.lines.push_back(line(kindFirst(kind, out.header.version), "RINEX VERSION / TYPE"));
  out.header.lines.push_back(line(rtklibCompat() ? rtklibPgmRunByDate() : "ceqc            ceqc            ", "PGM / RUN BY / DATE"));
  addRtklibComments(out, files);
  for (auto& meta : copiedMeta(files, kind)) {
    if (kind == RinexKind::Obs && meta.label == "TIME OF FIRST OBS") continue;
    if (kind == RinexKind::Obs && meta.label == "TIME OF LAST OBS") continue;
    if (kind == RinexKind::Obs && meta.label == "INTERVAL") continue;
    if (kind == RinexKind::Obs && meta.label == "SYS / PHASE SHIFT") continue;
    // RINEX 3/4-only GLONASS header records must not leak into RINEX 2 OBS.
    // teqc and strict RINEX 2 readers treat them as malformed header records.
    if (kind == RinexKind::Obs && v2obs &&
        (meta.label == "GLONASS SLOT / FRQ #" || meta.label == "GLONASS COD/PHS/BIS")) continue;
    out.header.lines.push_back(meta);
  }
  if (kind == RinexKind::Obs) {
    std::vector<ObservationRecord> recs;
    for (auto& f : files) if (f.header.kind == RinexKind::Obs) recs.insert(recs.end(), f.data.observationRecords.begin(), f.data.observationRecords.end());
    if (v2obs) {
      // RINEX 2 is a legacy observation format.  Modern constellations such as
      // Galileo/BDS/QZSS can be represented only awkwardly in RINEX2 mixed files,
      // and Anubis 3.11's RINEX2 reader treats long modern mixed epochs as
      // incomplete/buffer-limited.  Keep RINEX2 output to the legacy systems that
      // the RINEX2 QC path is intended to validate (GPS/GLONASS).  RINEX3/4
      // remain the canonical output for BDS/GAL/QZS and keep all systems.
      recs.erase(std::remove_if(recs.begin(), recs.end(), [](const ObservationRecord& r) {
        return !(r.system == "G" || r.system == "R");
      }), recs.end());
    }
    std::stable_sort(recs.begin(), recs.end(), [](auto& a, auto& b) { return a.time < b.time; });
    out.data.observationRecords = recs;
    auto types = observationTypes(recs, v2obs);
    maybeFillApproxFromNav(out, files);
    addRequiredObsDefaults(out);
    applyRtklibObsMetadata(out);
    addGlonassSlotAndBiasHeaders(out, files);
    if (v2obs) {
      auto tv = types.empty() ? std::vector<std::string>{} : types.begin()->second;
      addV2ObsTypes(out, tv);
      out.header.obsTypes["G"] = tv;
    } else {
      addV3ObsTypes(out, types);
      addV3PhaseShifts(out, types);
      out.header.obsTypes = types;
    }
    addObsTimeHeaders(out, recs);
    reorderRtklibObsHeader(out);
    out.header.lines.push_back(line("", "END OF HEADER"));
    indexHeader(out.header);
    if (v2obs) makeV2Body(out, recs, out.header.obsTypes["G"]);
    else makeV3Body(out, recs, types);
    parseContent(out);
    return out;
  }

  for (auto& f : files) {
    if (f.header.kind != kind) continue;
    if (kind == RinexKind::Nav) out.data.navigationRecords.insert(out.data.navigationRecords.end(), f.data.navigationRecords.begin(), f.data.navigationRecords.end());
    if (kind == RinexKind::Met) out.data.meteorologicalRecords.insert(out.data.meteorologicalRecords.end(), f.data.meteorologicalRecords.begin(), f.data.meteorologicalRecords.end());
    if (kind != RinexKind::Nav) out.body.insert(out.body.end(), f.body.begin(), f.body.end());
  }
  if (kind == RinexKind::Met && !out.data.meteorologicalRecords.empty()) {
    std::vector<std::string> mt;
    for (const auto& r : out.data.meteorologicalRecords) {
      for (const auto& kv : r.values) if (std::find(mt.begin(), mt.end(), kv.first) == mt.end()) mt.push_back(kv.first);
    }
    if (mt.empty()) {
      for (const auto& f : files) if (f.header.kind == RinexKind::Met) { mt = f.header.metTypes; if (!mt.empty()) break; }
    }
    if (!mt.empty()) {
      std::ostringstream v;
      v << std::setw(6) << mt.size();
      for (auto& t : mt) v << std::setw(6) << t;
      out.header.lines.push_back(line(v.str(), "# / TYPES OF OBSERV"));
      out.header.metTypes = mt;
    }
    out.body.clear();
    for (const auto& r : out.data.meteorologicalRecords) {
      auto tm = toUTC(r.time);
      std::ostringstream row;
      row << std::setw(3) << (tm.tm_year % 100) << std::setw(3) << (tm.tm_mon + 1) << std::setw(3) << tm.tm_mday
          << std::setw(3) << tm.tm_hour << std::setw(3) << tm.tm_min << std::setw(3) << tm.tm_sec;
      for (auto& t : out.header.metTypes) {
        auto it = r.values.find(t);
        if (it != r.values.end()) row << std::setw(8) << std::fixed << std::setprecision(3) << it->second;
        else row << std::setw(8) << "";
      }
      out.body.push_back(row.str());
    }
  }
  if (kind == RinexKind::Nav && !out.data.navigationRecords.empty()) {
    if (!rtklibCompat()) {
      std::sort(out.data.navigationRecords.begin(), out.data.navigationRecords.end(), [](const auto& a, const auto& b) {
        if (a.satellite != b.satellite) return a.satellite < b.satellite;
        if (a.epoch && b.epoch && *a.epoch != *b.epoch) return *a.epoch < *b.epoch;
        return a.messageType < b.messageType;
      });
    }
    // In +rtklib/+rtkplot mode preserve translator/message order.  RTKCONV-EX
    // outputs mixed G/R/C ephemerides in stream order, and RTKPlot has fewer
    // surprises when CEQC does not regroup the whole NAV by satellite system.
    std::set<std::string> seen;
    for (auto& r : out.data.navigationRecords) {
      std::string key = r.satellite + ":" + r.messageType + ":" + r.messageSubtype + ":" + (r.epoch ? formatUTC(*r.epoch) : "");
      if (seen.count(key)) continue;
      seen.insert(key);
      auto lines = formatNavRecord(r, out.header.version);
      out.body.insert(out.body.end(), lines.begin(), lines.end());
    }
  }
  out.header.lines.push_back(line("", "END OF HEADER"));
  indexHeader(out.header);
  return out;
}

void setWriterOptions(WriterOptions options){ g_writerOptions = options; }
WriterOptions writerOptions(){ return g_writerOptions; }

} // namespace ceqc::service::rinex
