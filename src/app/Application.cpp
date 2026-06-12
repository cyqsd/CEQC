#include "ceqc/app/Application.hpp"
#include "ceqc/cli/CommandLine.hpp"
#include "ceqc/qc/QCService.hpp"
#include "ceqc/rinex/RinexService.hpp"
#include "ceqc/translate/Translator.hpp"
#include "ceqc/io/Console.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <filesystem>
#include <set>


namespace {
uint32_t binexCRC32(const std::string& data) {
  uint32_t crc = 0xFFFFFFFFu;
  for(unsigned char c: data){ crc ^= c; for(int i=0;i<8;++i) crc = (crc&1) ? (crc>>1)^0xEDB88320u : (crc>>1); }
  return ~crc;
}
void writeBinexVL(std::ostream& os, uint32_t v) {
  // BINEX-style variable length: high bit continues, 7 payload bits per byte.
  unsigned char buf[5]; int n=0; do { buf[n++] = static_cast<unsigned char>(v & 0x7F); v >>= 7; } while(v);
  for(int i=n-1;i>=0;--i){ unsigned char b=buf[i]; if(i) b|=0x80; os.put(static_cast<char>(b)); }
}

std::string lowerPathExt(std::string s) {
  for(char& c:s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}
std::vector<std::string> candidateNavFilesForObs(const std::string& obsPath) {
  namespace fs = std::filesystem;
  std::vector<std::string> out;
  fs::path p(obsPath);
  fs::path dir = p.parent_path();
  std::string stem = p.stem().string();
  if (stem.empty()) return out;
  std::vector<std::string> exts = {".nav", ".NAV", ".n", ".N", ".rnx", ".RNX"};
  for (auto& e : exts) {
    fs::path c = dir / (stem + e);
    if (fs::exists(c)) out.push_back(c.string());
  }
  // RINEX 2 compressed/legacy naming can be station/DOY.??n; search same directory
  // for files with the same stem and a navigation-looking extension.
  if (fs::exists(dir.empty()?fs::path("."):dir)) {
    std::error_code ec;
    for (auto& ent : fs::directory_iterator(dir.empty()?fs::path("."):dir, ec)) {
      if (ec || !ent.is_regular_file()) continue;
      auto q = ent.path();
      if (q.stem().string() != stem) continue;
      auto ext = lowerPathExt(q.extension().string());
      bool navExt = (ext == ".nav" || ext == ".rnx" || ext == ".gnav" || ext == ".glo" || ext == ".hnav" || ext == ".qnav" || ext == ".lnav");
      if (!navExt && ext.size()==4) {
        char last = ext.back();
        navExt = (last=='n' || last=='g' || last=='l' || last=='h' || last=='q' || last=='p');
      }
      if (navExt) out.push_back(q.string());
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}
bool alreadyLoadedPath(const std::vector<ceqc::model::RinexFile>& files, const std::string& path) {
  namespace fs = std::filesystem;
  std::error_code ec;
  auto a = fs::weakly_canonical(fs::path(path), ec).string();
  for (auto& rf : files) {
    std::error_code ec2;
    auto b = fs::weakly_canonical(fs::path(rf.path), ec2).string();
    if (!a.empty() && a == b) return true;
  }
  return false;
}

void writeBinex00Metadata(const std::string& path, const ceqc::model::Operation& op, const std::vector<ceqc::model::RinexFile>& files) {
  std::ostringstream payload;
  payload << "ceqc BINEX 0x00 site metadata\n";
  for (auto& kv : op.binexMetadata) payload << kv.first << "=" << kv.second << "\n";
  for (auto& rf : files) {
    if (rf.header.kind == ceqc::model::RinexKind::Obs || rf.header.kind == ceqc::model::RinexKind::Nav || rf.header.kind == ceqc::model::RinexKind::Met) {
      for (auto& h : rf.header.lines) {
        if (h.label == "MARKER NAME" || h.label == "MARKER NUMBER" || h.label == "MARKER TYPE" || h.label == "OBSERVER / AGENCY" || h.label == "REC # / TYPE / VERS" || h.label == "ANT # / TYPE" || h.label == "APPROX POSITION XYZ" || h.label == "ANTENNA: DELTA H/E/N")
          payload << h.label << "=" << h.value << "\n";
      }
      break;
    }
  }
  std::string pl = payload.str();
  std::ofstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot create " + path);
  std::ostringstream record;
  record.put(static_cast<char>(0x00));
  writeBinexVL(record, static_cast<uint32_t>(pl.size()));
  record.write(pl.data(), static_cast<std::streamsize>(pl.size()));
  std::string rec = record.str();
  uint32_t crc = binexCRC32(rec);
  // sync + record 0x00 + variable length + payload + CRC32 little-endian.
  f.put(static_cast<char>(0xE2));
  f.write(rec.data(), static_cast<std::streamsize>(rec.size()));
  for(int i=0;i<4;++i) f.put(static_cast<char>((crc>>(8*i))&0xFF));
}

bool hasKey(const std::map<std::string,std::string>& m, const std::string& token, std::string* value=nullptr) {
  auto up = [](std::string x){ for(char& c:x) c=(char)std::toupper((unsigned char)c); return x; };
  std::string t = up(token);
  for (auto& kv : m) {
    std::string k = up(kv.first);
    if (k.find(t) != std::string::npos) { if(value) *value = kv.second; return true; }
  }
  return false;
}
std::map<std::string,std::string> headerOnly(const std::map<std::string,std::string>& m) {
  std::map<std::string,std::string> out;
  auto up = [](std::string x){ for(char& c:x) c=(char)std::toupper((unsigned char)c); return x; };
  for (auto& kv : m) {
    std::string k = up(kv.first);
    if (k.find(".OBS")!=std::string::npos || k.find(".-OBS")!=std::string::npos || k.find(".RENAME_OBS")!=std::string::npos || k.find(".DEC")!=std::string::npos || k.find(".SUMMARY")!=std::string::npos) continue;
    out[kv.first] = kv.second;
  }
  return out;
}
void writeObsSummary(const std::string& target, const ceqc::model::RinexFile& rf, bool append) {
  if (target.empty()) return;
  std::ostream* out = nullptr;
  std::ofstream file;
  if (target != "e" && target != "stdout" && target != "-") { file.open(target, append ? std::ios::app : std::ios::out); out=&file; }
  else out=&std::cerr;
  if(!out || !(*out)) return;
  std::map<std::string,int> counts;
  for(auto& r: rf.data.observationRecords) for(auto& v:r.values) counts[v.type]++;
  *out << "# ceqc observable summary\n";
  *out << "epochs " << rf.data.observationEpochs.size() << "\n";
  *out << "records " << rf.data.observationRecords.size() << "\n";
  for(auto& kv:counts) *out << kv.first << " " << kv.second << "\n";
}

std::string readAllText(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot read golden file: " + path);
  std::ostringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

std::string normalizeEOL(std::string s, const std::string& eol) {
  std::string tmp;
  tmp.reserve(s.size());
  for (size_t i = 0; i < s.size(); ++i) {
    if (s[i] == '\r') {
      if (i + 1 < s.size() && s[i + 1] == '\n') ++i;
      tmp.push_back('\n');
    } else {
      tmp.push_back(s[i]);
    }
  }
  if (eol == "crlf") {
    std::string out;
    out.reserve(tmp.size() + 64);
    for (char c : tmp) {
      if (c == '\n') out += "\r\n";
      else out.push_back(c);
    }
    return out;
  }
  return tmp;
}

void compareGolden(const std::string& actual, const ceqc::model::Operation& op, std::ostream& err) {
  if (op.teqcGolden.empty()) return;
  auto golden = normalizeEOL(readAllText(op.teqcGolden), op.teqcEOL);
  if (actual == golden) {
    err << "teqc-compat: byte-identical to " << op.teqcGolden << "\n";
    return;
  }
  size_t n = std::min(actual.size(), golden.size());
  size_t off = 0;
  while (off < n && actual[off] == golden[off]) ++off;
  err << "teqc-compat: byte mismatch at offset " << off;
  if (off < actual.size()) err << ": ceqc=0x" << std::hex << std::setw(2) << std::setfill('0') << (static_cast<unsigned>(static_cast<unsigned char>(actual[off]))) << std::dec << std::setfill(' ');
  else err << ": ceqc=<EOF>";
  if (off < golden.size()) err << " teqc=0x" << std::hex << std::setw(2) << std::setfill('0') << (static_cast<unsigned>(static_cast<unsigned char>(golden[off]))) << std::dec << std::setfill(' ');
  else err << " teqc=<EOF>";
  err << "\n";
  if (!op.teqcDiff.empty()) {
    std::ofstream d(op.teqcDiff);
    size_t a = off > 64 ? off - 64 : 0;
    size_t b = std::min(std::max(actual.size(), golden.size()), off + 64);
    d << "offset " << off << "\n";
    d << "--- ceqc context ---\n" << actual.substr(a, std::min(actual.size(), b) - a) << "\n";
    d << "--- teqc context ---\n" << golden.substr(a, std::min(golden.size(), b) - a) << "\n";
  }
}


std::string systemLongToShort(const std::string& x) {
  if (x == "GPS") return "G"; if (x == "GLONASS") return "R"; if (x == "SBAS") return "S";
  if (x == "Galileo" || x == "GALILEO") return "E"; if (x == "Beidou" || x == "BEIDOU" || x == "BDS") return "C";
  if (x == "QZSS") return "J"; if (x == "IRNSS" || x == "NavIC") return "I";
  return x.empty()?"":std::string(1, x[0]);
}
int satPRN(const std::string& sat) {
  if (sat.size() < 2) return 0;
  try { return std::stoi(sat.substr(1)); } catch(...) { return 0; }
}
bool keepSatellite(const std::string& sat, const ceqc::model::Operation& op) {
  if (sat.empty()) return true;
  std::string sys(1, sat[0]);
  int prn = satPRN(sat);
  for (auto& kv : op.maxExpectedSVs) if (systemLongToShort(kv.first) == sys && prn > kv.second) return false;
  bool anyInclusive = false;
  for (const auto& kv : op.svSelections) if (kv.second.specified && kv.second.includeMode) anyInclusive = true;
  auto it = op.svSelections.find(sys);
  if (it == op.svSelections.end() || !it->second.specified) return !anyInclusive;
  const auto& sel = it->second;
  if (sel.all) return sel.includeMode;
  bool listed = std::find(sel.prns.begin(), sel.prns.end(), prn) != sel.prns.end();
  return sel.includeMode ? listed : !listed;
}
void applySVFilters(ceqc::model::RinexFile& rf, const ceqc::model::Operation& op) {
  using namespace ceqc::model;
  if (rf.header.kind == RinexKind::Obs) {
    std::vector<ObservationRecord> kept;
    for (const auto& r : rf.data.observationRecords) if (keepSatellite(r.satellite, op)) kept.push_back(r);
    if (!op.svDuplicates) {
      std::map<std::string,ObservationRecord> last;
      for (const auto& r : kept) last[ceqc::model::formatUTC(r.time)+":"+r.satellite] = r;
      kept.clear(); for (auto& kv : last) kept.push_back(kv.second);
    }
    if (op.maxRxSVsSpecified && op.maxRxSVs > 0) {
      std::map<std::string,int> perEpoch;
      std::vector<ObservationRecord> limited;
      for (const auto& r : kept) { auto key = ceqc::model::formatUTC(r.time); if (++perEpoch[key] <= op.maxRxSVs) limited.push_back(r); }
      kept.swap(limited);
    }
    if (op.orderSVByPRN) std::sort(kept.begin(), kept.end(), [](const auto& a,const auto& b){ return a.time==b.time ? a.satellite < b.satellite : a.time < b.time; });
    rf.data.observationRecords.swap(kept);
    rf = ceqc::service::rinex::merge({rf}, RinexKind::Obs, rf.header.version);
  } else if (rf.header.kind == RinexKind::Nav) {
    std::vector<NavigationRecord> kept;
    for (const auto& r : rf.data.navigationRecords) if (keepSatellite(r.satellite, op)) kept.push_back(r);
    rf.data.navigationRecords.swap(kept);
    rf = ceqc::service::rinex::merge({rf}, RinexKind::Nav, rf.header.version);
  }
}

bool messageMatch(const ceqc::model::NavigationRecord& r, const std::string& token) {
  auto up=[](std::string x){ for(char& c:x)c=(char)std::toupper((unsigned char)c); return x; };
  std::string t=up(token), m=up(r.messageType), sub=up(r.messageSubtype), sys=up(r.system), sat=up(r.satellite);
  return t==m || t==sub || t==sys || t==sat || t==(m+":"+sub) || t==(sys+":"+m);
}
void applyEphFilters(ceqc::model::RinexFile& rf, const ceqc::model::Operation& op) {
  using namespace ceqc::model;
  if(rf.header.kind != RinexKind::Nav) return;
  if(op.ephIncludeTypes.empty() && op.ephExcludeTypes.empty()) return;
  std::vector<NavigationRecord> kept;
  for(const auto& r: rf.data.navigationRecords){
    bool keep = op.ephIncludeTypes.empty();
    for(const auto& t: op.ephIncludeTypes) if(messageMatch(r,t)) keep=true;
    for(const auto& t: op.ephExcludeTypes) if(messageMatch(r,t)) keep=false;
    if(keep) kept.push_back(r);
  }
  rf.data.navigationRecords.swap(kept);
  rf = ceqc::service::rinex::merge({rf}, RinexKind::Nav, rf.header.version);
}

void printConfig(std::ostream& os, const ceqc::model::Operation& op, bool all, bool bcf) {
  os << (bcf ? "# ceqc BINEX configuration export\n" : "# ceqc configuration\n");
  os << (op.qc ? "+qc" : "-qc") << "\n";
  if (!op.translatorName.empty()) os << "-tr " << op.translatorName << "\n";
  if (op.targetVersion >= 2) os << "+v" << static_cast<int>(op.targetVersion) << "\n";
  for (auto& kv : op.obsHeaderEdits) os << kv.first << " \"" << kv.second << "\"\n";
  for (auto& kv : op.navHeaderEdits) os << kv.first << " \"" << kv.second << "\"\n";
  for (auto& kv : op.metHeaderEdits) os << kv.first << " \"" << kv.second << "\"\n";
  for (auto& kv : op.binexMetadata) os << kv.first << " \"" << kv.second << "\"\n";
  for (auto& kv : op.svSelections) {
    os << (kv.second.includeMode?"+":"-") << kv.first;
    if (!kv.second.all) for (size_t i=0;i<kv.second.prns.size();++i) os << (i?",":"") << kv.second.prns[i];
    os << "\n";
  }
  if (op.maxRxChannels != 12 || all) os << "-max_rx_ch " << op.maxRxChannels << "\n";
  if (op.maxRxSVs != 12 || all) os << "-max_rx_SVs " << op.maxRxSVs << "\n";
  if (all) os << "# all known major switches\n+obs <file>\n+nav <file>\n+met <file>\n+qc\n+plot\n+slips <file>\n+teqc\n";
}

void printFormats(std::ostream& os, const ceqc::service::translator::Registry& reg) {
  os << "ceqc supported input translators:\n";
  for (const auto& f : reg.formats()) {
    os << "  " << (f.implemented ? "* " : "- ") << f.name;
    if (!f.aliases.empty()) {
      os << " (";
      for (size_t i = 0; i < f.aliases.size(); ++i) { if (i) os << ","; os << f.aliases[i]; }
      os << ")";
    }
    os << " : " << f.description << "\n";
  }
}


} // namespace

namespace ceqc::app {
using namespace ceqc::model;

int Application::run(const std::vector<std::string>& args) {
  try {
    auto op = ceqc::cli::parseArgs(args);
    if (op.showHelp || args.empty()) { view::printHelp(err_); return 0; }
    if (op.showVersion) { view::printVersion(err_); return 0; }
    if (op.showID) { err_ << "ceqc 0.0.1 C++21-cleanroom\n"; return 0; }
    if (op.showConfig || op.showBCF) { printConfig(out_, op, op.showAllConfig, op.showBCF); return 0; }

    service::translator::Registry reg;
    if (op.showFormats) { printFormats(out_, reg); return 0; }
    if (op.inputFiles.empty()) { err_ << "ceqc: no input files\n"; return 2; }
    service::rinex::setParserOptions({op.relaxRinex, op.extendRinex2, op.reformatRinex});
    std::vector<RinexFile> files;
    for (auto& p : op.inputFiles) {
      const service::translator::Translator* tr = nullptr;
      if (!op.translatorName.empty()) tr = reg.find(op.translatorName);
      else tr = reg.detect(p);
      if (!tr) throw std::runtime_error("cannot detect input format: " + p);
      auto decoded = tr->decode(p);
      for (auto& d : decoded) if (d.path.empty()) d.path = p;
      files.insert(files.end(), decoded.begin(), decoded.end());
    }
    for (auto& n : op.qcOptions.navFiles) files.push_back(service::rinex::readFile(n));
    if (op.qc) {
      std::vector<std::string> autoNavs;
      for (auto& rf : files) {
        if (rf.header.kind != RinexKind::Obs) continue;
        for (auto& n : candidateNavFilesForObs(rf.path)) {
          if (!alreadyLoadedPath(files, n)) autoNavs.push_back(n);
        }
      }
      std::sort(autoNavs.begin(), autoNavs.end());
      std::vector<std::string> dedupNavs;
      std::set<std::string> canonSeen;
      for (auto& n : autoNavs) {
        std::error_code ec;
        std::string c = std::filesystem::weakly_canonical(std::filesystem::path(n), ec).string();
        if (c.empty()) c = n;
        if (canonSeen.insert(c).second) dedupNavs.push_back(n.rfind("./",0)==0 ? n.substr(2) : n);
      }
      for (auto& n : dedupNavs) {
        try { files.push_back(service::rinex::readFile(n)); }
        catch(const std::exception& e) { err_ << "ceqc: warning: cannot read auto NAV file " << n << ": " << e.what() << "\n"; }
      }
    }

    for (auto& rf : files) {
      if (rf.header.kind == RinexKind::Obs) {
        service::rinex::applyHeaderEdits(rf, headerOnly(op.obsHeaderEdits));
        std::string v;
        if (hasKey(op.obsHeaderEdits, ".RENAME_OBS", &v)) service::rinex::applyObsTypeFilter(rf, v, false, true);
        if (hasKey(op.obsHeaderEdits, ".OBS", &v) && !hasKey(op.obsHeaderEdits, ".-OBS", nullptr)) service::rinex::applyObsTypeFilter(rf, v, false, false);
        if (hasKey(op.obsHeaderEdits, ".-OBS", &v)) service::rinex::applyObsTypeFilter(rf, v, true, false);
        rf = service::rinex::windowObservation(rf, op.windowStart, op.windowEnd);
        rf = service::rinex::decimate(rf, op.obsDecimation);
        if (op.obsSummaryTarget != "e") writeObsSummary(op.obsSummaryTarget, rf, op.obsSummaryAppend);
      } else if (rf.header.kind == RinexKind::Nav) {
        service::rinex::applyHeaderEdits(rf, headerOnly(op.navHeaderEdits));
        rf = service::rinex::decimate(rf, op.navDecimation);
      } else if (rf.header.kind == RinexKind::Met) {
        service::rinex::applyHeaderEdits(rf, headerOnly(op.metHeaderEdits));
        std::string v;
        if (hasKey(op.metHeaderEdits, ".RENAME_OBS", &v)) service::rinex::applyMetTypeFilter(rf, v, false, true);
        if (hasKey(op.metHeaderEdits, ".OBS", &v) && !hasKey(op.metHeaderEdits, ".-OBS", nullptr)) service::rinex::applyMetTypeFilter(rf, v, false, false);
        if (hasKey(op.metHeaderEdits, ".-OBS", &v)) service::rinex::applyMetTypeFilter(rf, v, true, false);
        rf = service::rinex::decimate(rf, op.metDecimation);
      }
      applySVFilters(rf, op);
      applyEphFilters(rf, op);
      if (op.targetVersion > 0.0 && rf.header.kind != RinexKind::Unknown) {
        rf = service::rinex::merge({rf}, rf.header.kind, op.targetVersion);
      }
    }

    if (op.verifyOnly) {
      for (auto& rf : files) view::printIssues(out_, rf.path, service::rinex::validate(rf));
    }

    std::string renderedQC;
    std::vector<QCSummary> qcSummariesForAux;
    if (op.qc) {
      std::vector<NavigationRecord> navs;
      std::vector<std::string> navPaths;
      std::map<std::string,int> navSatAppear;
      for (auto& rf : files) if (rf.header.kind == RinexKind::Nav) {
        navPaths.push_back(rf.path);
        navs.insert(navs.end(), rf.data.navigationRecords.begin(), rf.data.navigationRecords.end());
        for (auto& nr : rf.data.navigationRecords) if(!nr.satellite.empty()) navSatAppear[nr.satellite]++;
      }
      bool hasObsForQC = false;
      for (auto& rf : files) if (rf.header.kind == RinexKind::Obs) hasObsForQC = true;
      std::ostringstream qout;
      for (auto& rf : files) {
        if (hasObsForQC && rf.header.kind == RinexKind::Nav) continue;
        auto s = navs.empty() ? service::qc::analyze(rf, op.qcOptions) : service::qc::analyzeWithNavigation(rf, navs, op.qcOptions);
        s.navInputFiles = navPaths;
        s.navigationSatelliteAppearance = navSatAppear;
        if (rf.header.kind == RinexKind::Obs) qcSummariesForAux.push_back(s);
        view::printQC(qout, s, op.quietQC, op.teqcCompat || !op.quietQC);
      }
      renderedQC = normalizeEOL(qout.str(), op.teqcEOL);
      out_ << renderedQC;
      compareGolden(renderedQC, op, err_);
    }

    auto rtklibTargetVersion = [&](RinexKind kind) {
      double requested = op.targetVersion > 0.0 ? op.targetVersion : (kind == RinexKind::Nav ? 4.02 : 3.05);
      if (!op.rtklibCompat) return requested;
      // +rtklib/+rtkplot is a compatibility projection for RTKLIB_EX/RTKPLOT,
      // not a replacement for strict mainline RINEX.  RTKLIB_EX 2.5.0 can plot
      // both its own RINEX 3.05 and RINEX 4.02 exports, provided the NAV body is
      // the legacy mixed body (no RINEX-4 > EPH typed records), signal aliases
      // and blank LLI/SSI follow RTKCONV-EX style.  Therefore preserve the
      // explicit +v requested by the user.  If no +v is provided, default the
      // RTKLIB compatibility export to 3.05, which is RTKCONV-EX's safest common
      // plotting format.
      if (op.targetVersion > 0.0) return requested;
      (void)kind;
      return 3.05;
    };

    service::rinex::setWriterOptions({op.rtklibCompat});
    if (!op.outputObs.empty()) {
      auto m = service::rinex::merge(files, RinexKind::Obs, rtklibTargetVersion(RinexKind::Obs));
      service::rinex::applyHeaderEdits(m, headerOnly(op.obsHeaderEdits));
      service::rinex::writeFile(op.outputObs, m);
      if (op.obsSummaryTarget == "e") writeObsSummary(op.outputObs, m, true);
    }
    if (!op.outputNav.empty()) {
      // Always remove a previous target before deciding whether a fresh NAV can
      // be written.  Several validation scripts reuse the same filenames; if a
      // later decode has no ephemeris records, leaving an older header-only NAV
      // behind makes gfzrnx/Anubis report a false CEQC format error.
      std::error_code rmEc;
      std::filesystem::remove(op.outputNav, rmEc);

      double navTarget = rtklibTargetVersion(RinexKind::Nav);
      bool haveNavRecords = false;
      for (const auto& rf : files) {
        if (rf.header.kind == RinexKind::Nav && !rf.data.navigationRecords.empty()) { haveNavRecords = true; break; }
      }
      if (haveNavRecords) {
        auto m = service::rinex::merge(files, RinexKind::Nav, navTarget);
        service::rinex::applyHeaderEdits(m, headerOnly(op.navHeaderEdits));
        if (m.data.navigationRecords.empty()) {
          err_ << "ceqc: warning: navigation records were decoded but merge produced no output records; skip writing " << op.outputNav << "\n";
        } else {
          service::rinex::writeFile(op.outputNav, m);
        }
      } else {
        // Do not create a syntactically valid but data-empty NAV file.  gfzrnx and
        // Anubis correctly treat header-only NAV as an error ("no input data got").
        // For RAWX-only UBX logs, OBS is still valid and NAV validation must be skipped.
        err_ << "ceqc: warning: no navigation ephemeris records decoded; not writing " << op.outputNav << "\n";
      }
    }
    service::rinex::setWriterOptions({false});
    if (!op.outputMet.empty()) { auto m = service::rinex::merge(files, RinexKind::Met, op.targetVersion); service::rinex::applyHeaderEdits(m, headerOnly(op.metHeaderEdits)); service::rinex::writeFile(op.outputMet, m); }
    if (!op.outputBinex.empty()) { writeBinex00Metadata(op.outputBinex, op, files); }

    if (op.qcOptions.plot || op.qcOptions.slipsEnabled) {
      auto emitAux = [&](const QCSummary& s) {
        if (op.qcOptions.plot) {
          std::string root = op.qcOptions.root.empty() ? "ceqc" : op.qcOptions.root;
          std::ofstream f(root + ".plot");
          f << service::qc::makePlot(s);
        }
        if (op.qcOptions.slipsEnabled && !op.qcOptions.slipsTarget.empty()) {
          std::ofstream slips(op.qcOptions.slipsTarget, op.qcOptions.slipsAppend ? std::ios::app : std::ios::out);
          if (s.derived) {
            for (auto& ev : s.derived->slipEvents) slips << ev.time << " " << ev.satellite << " " << ev.type << " " << ev.detail << "\n";
          }
        }
      };
      if (!qcSummariesForAux.empty()) {
        for (const auto& s : qcSummariesForAux) emitAux(s);
      } else {
        std::vector<NavigationRecord> navs;
        for (auto& rf : files) if (rf.header.kind == RinexKind::Nav) navs.insert(navs.end(), rf.data.navigationRecords.begin(), rf.data.navigationRecords.end());
        for (auto& rf : files) {
          if (rf.header.kind != RinexKind::Obs) continue;
          auto s = navs.empty() ? service::qc::analyze(rf, op.qcOptions) : service::qc::analyzeWithNavigation(rf, navs, op.qcOptions);
          emitAux(s);
        }
      }
    }

    return 0;
  } catch (const std::exception& e) {
    err_ << "ceqc: " << e.what() << "\n";
    return 1;
  }
}
} // namespace ceqc::app
