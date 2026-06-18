#include "ceqc/cli/CommandLine.hpp"
#include "ceqc/rinex/RinexService.hpp"
#include <cctype>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace ceqc::cli {
using namespace ceqc::model;
namespace {

bool has(const std::string& s, const std::string& p) { return s.rfind(p, 0) == 0; }
double num(const std::string& s) { return std::stod(s); }

bool nearVersion(double a, double b) { return std::fabs(a - b) < 5e-4; }

double exactRinexVersionValue(double v) {
  // normVer() treats mathematically integral values (3.00, 4.00) as a major
  // request and maps them to CEQC defaults.  Add a tiny, non-printing epsilon
  // for explicit .00 requests so the header still prints 3.00/4.00.
  double major = std::floor(v + 1e-9);
  if (std::fabs(v - major) < 1e-9) return v + 1e-6;
  return v;
}

double parseRinexVersionOption(const std::string& opt) {
  std::string vs = opt.substr(2);
  if (vs.empty()) throw std::runtime_error("invalid RINEX version option " + opt);
  if (vs == "2") return 2.11;
  if (vs == "3") return 3.05;
  if (vs == "4") return 4.02;
  double v = 0.0;
  try { v = std::stod(vs); }
  catch (...) { throw std::runtime_error("invalid RINEX version option " + opt); }
  // Keep the accepted list intentionally explicit.  CEQC must not silently emit
  // a header version whose body/header profile has not been accepted for validation.
  static const double supported[] = {
    2.10, 2.11,
    3.00, 3.01, 3.02, 3.03, 3.04, 3.05,
    4.00, 4.01, 4.02
  };
  for (double ok : supported) if (nearVersion(v, ok)) return exactRinexVersionValue(ok);
  throw std::runtime_error(
    "unsupported RINEX version option " + opt +
    "; supported versioned profiles: +v2/+v2.10/+v2.11, +v3/+v3.00..+v3.05, +v4/+v4.00..+v4.02");
}

std::vector<std::string> split(const std::string& s, char d=',') {
  std::vector<std::string> r;
  std::string x;
  std::istringstream is(s);
  while (std::getline(is, x, d)) {
    if (!x.empty()) r.push_back(x);
  }
  return r;
}

std::vector<std::string> tokenizeConfigLine(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool quote = false;
  char qchar = 0;
  for (size_t i = 0; i < line.size(); ++i) {
    char c = line[i];
    if (!quote && (c == '#' || c == ';')) break;
    if (quote) {
      if (c == qchar) { quote = false; continue; }
      if (c == '\\' && i + 1 < line.size()) { cur.push_back(line[++i]); continue; }
      cur.push_back(c);
      continue;
    }
    if (c == '\'' || c == '"') { quote = true; qchar = c; continue; }
    if (std::isspace(static_cast<unsigned char>(c))) {
      if (!cur.empty()) { out.push_back(cur); cur.clear(); }
    } else {
      cur.push_back(c);
    }
  }
  if (!cur.empty()) out.push_back(cur);
  return out;
}

std::vector<std::string> readConfigTokens(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot read config file: " + path);
  std::vector<std::string> out;
  std::string line;
  while (std::getline(f, line)) {
    auto v = tokenizeConfigLine(line);
    out.insert(out.end(), v.begin(), v.end());
  }
  return out;
}

std::vector<std::string> expandConfigArgs(const std::vector<std::string>& args, std::vector<std::string>& configs) {
  std::vector<std::string> out;
  for (size_t i = 0; i < args.size(); ++i) {
    const auto& a = args[i];
    if ((a == "-config" || a == "--config") && i + 1 < args.size()) {
      std::string cfg = args[++i];
      configs.push_back(cfg);
      auto toks = readConfigTokens(cfg);
      auto nested = expandConfigArgs(toks, configs);
      out.insert(out.end(), nested.begin(), nested.end());
    } else {
      out.push_back(a);
    }
  }
  return out;
}

DecimationSpec dec(std::string s) {
  DecimationSpec d;
  d.enabled = true;
  d.raw = s;
  auto p = s.find(':');
  std::string a = p == std::string::npos ? s : s.substr(0, p);
  std::string b = p == std::string::npos ? "0" : s.substr(p + 1);
  auto parse = [](std::string x) {
    if (x.ends_with("h")) return std::chrono::seconds(static_cast<int>(std::stod(x.substr(0, x.size() - 1)) * 3600));
    if (x.ends_with("m")) return std::chrono::seconds(static_cast<int>(std::stod(x.substr(0, x.size() - 1)) * 60));
    if (x.ends_with("s")) x.pop_back();
    return std::chrono::seconds(static_cast<int>(std::stod(x)));
  };
  d.interval = parse(a);
  d.offset = parse(b);
  return d;
}

TimePoint parseTime(const std::string& s) {
  int y=1970, m=1, d=1, hh=0, mi=0;
  double sec=0;
  char c=0;
  std::istringstream is(s);
  is >> y >> c >> m >> c >> d >> c >> hh >> c >> mi >> c >> sec;
  if (is) return makeUTC(y, m, d, hh, mi, sec);

  // teqc window strings: [[[[[[YY]YY]MM]DD]hh]mm]ss[.sssss]
  std::string t = s;
  std::string frac;
  auto dot = t.find('.');
  if (dot != std::string::npos) { frac = t.substr(dot); t = t.substr(0, dot); }
  bool digits = !t.empty() && std::all_of(t.begin(), t.end(), [](char ch){ return std::isdigit(static_cast<unsigned char>(ch)); });
  if (!digits) return makeUTC(1970, 1, 1, 0, 0, 0);
  auto toInt = [&](size_t pos, size_t n, int def)->int { return pos+n<=t.size()?std::stoi(t.substr(pos,n)):def; };
  // Stable defaults; relative short forms are resolved to Unix epoch date rather than local date for reproducible tests.
  y=1970; m=1; d=1; hh=0; mi=0; sec=0;
  if (t.size()==2) { sec = toInt(0,2,0); }
  else if (t.size()==4) { mi=toInt(0,2,0); sec=toInt(2,2,0); }
  else if (t.size()==6) { hh=toInt(0,2,0); mi=toInt(2,2,0); sec=toInt(4,2,0); }
  else if (t.size()==8) { d=toInt(0,2,1); hh=toInt(2,2,0); mi=toInt(4,2,0); sec=toInt(6,2,0); }
  else if (t.size()==10) { m=toInt(0,2,1); d=toInt(2,2,1); hh=toInt(4,2,0); mi=toInt(6,2,0); sec=toInt(8,2,0); }
  else if (t.size()==12) { y=toInt(0,2,70); y += y<80?2000:1900; m=toInt(2,2,1); d=toInt(4,2,1); hh=toInt(6,2,0); mi=toInt(8,2,0); sec=toInt(10,2,0); }
  else if (t.size()==14) { y=toInt(0,4,1970); m=toInt(4,2,1); d=toInt(6,2,1); hh=toInt(8,2,0); mi=toInt(10,2,0); sec=toInt(12,2,0); }
  else return makeUTC(1970,1,1,0,0,0);
  if (!frac.empty()) sec += std::stod("0"+frac);
  return makeUTC(y,m,d,hh,mi,sec);
}

std::string collectXYZ(const std::vector<std::string>& args, size_t& i) {
  if (i + 3 >= args.size()) throw std::runtime_error("not enough values");
  std::string a = args[++i];
  std::string b = args[++i];
  std::string c = args[++i];
  return a + " " + b + " " + c;
}

std::string collectN(const std::vector<std::string>& args, size_t& i, int n) {
  if (i + static_cast<size_t>(n) >= args.size()) throw std::runtime_error("not enough values");
  std::string out;
  for (int k = 0; k < n; ++k) {
    if (k) out += " ";
    out += args[++i];
  }
  return out;
}

std::string collectUntilOption(const std::vector<std::string>& args, size_t& i) {
  std::string out;
  while (i + 1 < args.size()) {
    const auto& nx = args[i + 1];
    if (!nx.empty() && (nx[0] == '+' || nx[0] == '-')) break;
    if (!out.empty()) out += " ";
    out += args[++i];
  }
  return out;
}

std::string collectModWF(const std::vector<std::string>& args, size_t& i) {
  if (i + 3 >= args.size()) throw std::runtime_error("-O.mod_wf requires i j n {SV...}");
  std::string out;
  std::string i1 = args[++i];
  std::string j1 = args[++i];
  int n = std::stoi(args[++i]);
  out += i1 + std::string(" ") + j1 + std::string(" ") + std::to_string(n);
  for (int k = 0; k < n && i + 1 < args.size(); ++k) out += std::string(" ") + args[++i];
  return out;
}

bool nextIsValue(const std::vector<std::string>& args, size_t i) {
  return i + 1 < args.size() && !args[i + 1].empty() && args[i + 1][0] != '+' && args[i + 1][0] != '-';
}


std::vector<int> parsePRNList(std::string spec) {
  std::vector<int> out;
  std::replace(spec.begin(), spec.end(), ';', ',');
  std::istringstream is(spec);
  std::string tok;
  while (std::getline(is, tok, ',')) {
    if (tok.empty()) continue;
    auto dash = tok.find('-');
    if (dash == std::string::npos) {
      out.push_back(std::stoi(tok));
    } else {
      int a = std::stoi(tok.substr(0, dash));
      int b = std::stoi(tok.substr(dash + 1));
      if (b < a) std::swap(a, b);
      for (int v = a; v <= b; ++v) out.push_back(v);
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

bool isSystemFilter(const std::string& a) {
  if (a.size() < 2) return false;
  if (a[0] != '+' && a[0] != '-') return false;
  char c = a[1];
  return std::string("GRSECJI").find(c) != std::string::npos;
}

std::string sysName(char c) {
  switch(c) { case 'G': return "G"; case 'R': return "R"; case 'S': return "S"; case 'E': return "E"; case 'C': return "C"; case 'J': return "J"; case 'I': return "I"; default: return std::string(1,c); }
}

std::string lower(std::string s) { for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c))); return s; }

std::string normalizeSystemToken(std::string s) {
  s = lower(s);
  if (s == "g" || s == "gps") return "G";
  if (s == "r" || s == "glo" || s == "glonass") return "R";
  if (s == "e" || s == "gal" || s == "galileo") return "E";
  if (s == "c" || s == "bds" || s == "bei" || s == "beidou") return "C";
  if (s == "j" || s == "qzs" || s == "qzss") return "J";
  if (s == "s" || s == "sbas") return "S";
  if (s == "i" || s == "irn" || s == "navic") return "I";
  if (!s.empty()) {
    char c = static_cast<char>(std::toupper(static_cast<unsigned char>(s[0])));
    if (std::string("GRECJSI").find(c) != std::string::npos) return std::string(1, c);
  }
  return s;
}


} // namespace

Operation parseArgs(const std::vector<std::string>& rawArgs) {
  std::vector<std::string> configFiles;
  auto args = expandConfigArgs(rawArgs, configFiles);
  Operation op;
  op.configFiles = configFiles;
  bool qcRequested = false;
  bool translationRequested = false;
  for (size_t k = 0; k < args.size(); ++k) {
    const auto& x = args[k];
    if (x == "+qc" || x == "+qcq" || x == "+plot" || x == "+slips" || x == "++slips") qcRequested = true;
    if (x == "-tr") translationRequested = true;
  }
  for (size_t i = 0; i < args.size(); ++i) {
    auto a = args[i];
    if (a == "+help" || a == "-help" || a == "--help" || a == "-h") op.showHelp = true;
    else if (a == "+version" || a == "-version" || a == "--version") op.showVersion = true;
    else if (a == "+id" || a == "-id") op.showID = true;
    else if (a == "+config") op.showConfig = true;
    else if (a == "++config") { op.showConfig = true; op.showAllConfig = true; }
    else if (a == "+bcf") op.showBCF = true;
    else if (a == "+formats" || a == "++formats" || a == "--formats") op.showFormats = true;
    else if (a == "+rtklib" || a == "+rtkplot") op.rtklibCompat = true;
    else if (a == "-rtklib" || a == "-rtkplot") { /* teqc-style compatibility is enabled only with +rtklib/+rtkplot; ignore minus-form aliases */ }
    else if (a == "+verify") op.verifyOnly = true;
    else if (a == "+qc") op.qc = true;
    else if (a == "+qcq") { op.qc = true; op.quietQC = true; }
    else if (a == "-qc" || a == "-qcq") op.qc = false;
    else if (a == "+teqc" || a == "--teqc-compat") op.teqcCompat = true;
    else if (a == "-teqc_golden" && nextIsValue(args, i)) op.teqcGolden = args[++i];
    else if (a == "-teqc_diff" && nextIsValue(args, i)) op.teqcDiff = args[++i];
    else if (a == "-teqc_eol" && nextIsValue(args, i)) op.teqcEOL = lower(args[++i]);
    else if ((a == "+qc_json" || a == "+qc-json" || a == "+json_qc" || a == "+qcjson") && nextIsValue(args, i)) { op.qcJsonPath = args[++i]; op.qc = true; }
    else if ((a == "-qc_json" || a == "-qc-json") && nextIsValue(args, i)) { op.qcJsonPath = args[++i]; op.qc = true; }
    else if (a == "+obs" && nextIsValue(args, i)) op.outputObs = args[++i];
    else if (a == "+nav" && nextIsValue(args, i)) {
      std::string v = args[++i];
      // teqc compatibility: in QC/plot mode, +nav names the auxiliary NAV file(s).
      // In translation/export mode, keep CEQC/RTKLIB meaning: +nav is the output
      // NAV path.  This lets both of these work without ambiguity:
      //   ceqc -tr ubx +obs out.o +nav out.p raw.ubx
      //   ceqc +qc +plot +nav brdc.nav obs.o
      if (qcRequested && !translationRequested) {
        auto vals = split(v);
        op.qcOptions.navFiles.insert(op.qcOptions.navFiles.end(), vals.begin(), vals.end());
      } else {
        op.outputNav = v;
      }
    }
    else if (a == "+met" && nextIsValue(args, i)) op.outputMet = args[++i];
    else if (a == "+binex" && nextIsValue(args, i)) op.outputBinex = args[++i];
    else if (a == "-tr" && nextIsValue(args, i)) op.translatorName = args[++i];
    else if (a.rfind("+v",0)==0 && a.size()>2 && std::isdigit((unsigned char)a[2])) {
      op.targetVersion = parseRinexVersionOption(a);
    }
    else if (a == "-max_rx_ch" || a == "-max_rx_channels") { if (nextIsValue(args,i)) op.maxRxChannels = std::stoi(args[++i]); }
    else if (a == "-max_rx_SVs") { if (nextIsValue(args,i)) { op.maxRxSVs = std::stoi(args[++i]); op.maxRxSVsSpecified = true; } }
    else if (has(a,"-n_")) { if (nextIsValue(args,i)) op.maxExpectedSVs[a.substr(3)] = std::stoi(args[++i]); }
    else if (a == "+ch") op.useAllChannels = true;
    else if (has(a,"-ch")) op.useAllChannels = false;
    else if (a == "+NaN_obs") op.allowNaNObs = true;
    else if (a == "-NaN_obs") { op.allowNaNObs = false; if (nextIsValue(args,i)) ++i; }
    else if (a == "+sv_duplicates") op.svDuplicates = true;
    else if (a == "-sv_duplicates") op.svDuplicates = false;
    else if (a == "+svo") op.orderSVByPRN = true;
    else if (a == "-svo") op.orderSVByPRN = false;
    else if (a == "+relax") op.relaxRinex = true;
    else if (a == "-relax") op.relaxRinex = false;
    else if (a == "+extend") op.extendRinex2 = true;
    else if (a == "-extend") op.extendRinex2 = false;
    else if (a == "+reformat") op.reformatRinex = true;
    else if (a == "-reformat") op.reformatRinex = false;
    else if ((a == "-st" || a == "+st") && nextIsValue(args, i)) op.windowStart = parseTime(args[++i]);
    else if ((a == "-e" || a == "+e") && nextIsValue(args, i)) op.windowEnd = parseTime(args[++i]);
    else if (a == "-nav" && nextIsValue(args, i)) op.qcOptions.navFiles = split(args[++i]);
    else if (a == "-no_orbit" && nextIsValue(args, i)) { for (auto& s : split(args[++i], '+')) op.qcOptions.noOrbitSystems[normalizeSystemToken(s)] = true; }
    else if ((a == "-no_position" || a == "-no_pos") && nextIsValue(args, i)) { for (auto& s : split(args[++i], '+')) op.qcOptions.noPositionSystems[normalizeSystemToken(s)] = true; }
    else if (a == "+ap") op.qcOptions.averagePosition = true;
    else if (a == "-ap") op.qcOptions.averagePosition = false;
    else if (a == "+pos" || a == "+position") op.qcOptions.positionOnly = true;
    else if (a == "-pos" || a == "-position") op.qcOptions.positionOnly = false;
    else if (a == "+eep") op.qcOptions.everyEpochPosition = true;
    else if (a == "+eepx") { op.qcOptions.everyEpochPosition = true; op.qcOptions.everyEpochXYZ = true; }
    else if (a == "+eepg") { op.qcOptions.everyEpochPosition = true; op.qcOptions.everyEpochGeodetic = true; }
    else if (a == "+eepd") { op.qcOptions.everyEpochPosition = true; op.qcOptions.everyEpochDecimal = true; }
    else if (a == "+cl" || a == "+clock_slips") op.qcOptions.clockSlips = true;
    else if (a == "-cl" || a == "-clock_slips") op.qcOptions.clockSlips = false;
    else if (a == "+data") op.qcOptions.dataIndicators = true;
    else if (a == "-data") op.qcOptions.dataIndicators = false;
    else if (a == "+ceqc_ext" || a == "+ceqc-ext") op.qcOptions.ceqcExtension = true;
    else if (a == "-ceqc_ext" || a == "-ceqc-ext") op.qcOptions.ceqcExtension = false;
    else if (a == "+ion") op.qcOptions.ion = true;
    else if (a == "-ion") op.qcOptions.ion = false;
    else if (a == "+iod") op.qcOptions.iod = true;
    else if (a == "-iod") op.qcOptions.iod = false;
    else if (a == "+mp") op.qcOptions.multipath = true;
    else if (a == "-mp") op.qcOptions.multipath = false;
    else if (a == "+mp_raw") op.qcOptions.mpRaw = true;
    else if (a == "-mp_raw") op.qcOptions.mpRaw = false;
    else if (a == "+sn" || a == "+snr") op.qcOptions.snr = true;
    else if (a == "-sn" || a == "-snr") op.qcOptions.snr = false;
    else if (a == "+lli") op.qcOptions.lli = true;
    else if (a == "-lli") op.qcOptions.lli = false;
    else if (a == "+pl") op.qcOptions.pseudorangePhase = true;
    else if (a == "-pl") op.qcOptions.pseudorangePhase = false;
    else if (a == "+plot") op.qcOptions.plot = true;
    else if (a == "-plot") op.qcOptions.plot = false;
    else if (a == "+mask") op.qcOptions.mask = true;
    else if (a == "-mask") op.qcOptions.mask = false;
    else if (a == "+rs") op.qcOptions.riseSet = true;
    else if (a == "-rs") op.qcOptions.riseSet = false;
    else if (a == "+ssv") op.qcOptions.ssv = true;
    else if (a == "+svpr") op.qcOptions.svpr = true;
    else if (a == "+Y-code") op.qcOptions.yCode = true;
    else if (a == "-Y-code") op.qcOptions.yCode = false;
    else if (a == "+sym") op.qcOptions.symbolCodes = true;
    else if (a == "++sym") op.qcOptions.allSymbols = true;
    else if (a == "+slips") { op.qcOptions.slipsEnabled = true; if (nextIsValue(args, i)) op.qcOptions.slipsTarget = args[++i]; }
    else if (a == "++slips") { op.qcOptions.slipsEnabled = true; op.qcOptions.slipsAppend = true; if (nextIsValue(args, i)) op.qcOptions.slipsTarget = args[++i]; }
    else if (a == "-slips") { op.qcOptions.slipsEnabled = false; op.qcOptions.slipsAppend = false; op.qcOptions.slipsTarget.clear(); }
    else if ((a == "+w" || a == "-w") && nextIsValue(args, i)) { int w = std::stoi(args[++i]); if (w > 0) op.qcOptions.width = w; }
    else if ((a == "-bins" || a == "-ion_bins" || a == "-mp_bins" || a == "-sn_bins" || a == "-min_SVs") && nextIsValue(args, i)) {
      int v = std::stoi(args[++i]);
      if (a == "-bins") op.qcOptions.bins = v;
      else if (a == "-ion_bins") op.qcOptions.ionBins = v;
      else if (a == "-mp_bins") op.qcOptions.mpBins = v;
      else if (a == "-sn_bins") op.qcOptions.snBins = v;
      else op.qcOptions.minSVs = v;
    }
    else if ((a == "-set_mask" || a == "-set_maskdeg") && nextIsValue(args, i)) op.qcOptions.setMaskDeg = num(args[++i]);
    else if ((a == "-set_horizon" || a == "-set_hor") && nextIsValue(args, i)) op.qcOptions.setHorizonDeg = num(args[++i]);
    else if ((a == "-set_comp" || a == "-set_comparison") && nextIsValue(args, i)) op.qcOptions.setComparisonDeg = num(args[++i]);
    else if (a == "-root" && nextIsValue(args, i)) op.qcOptions.root = args[++i];
    else if (a == "-mp_sigmas" && nextIsValue(args, i)) op.qcOptions.mpSigmas = num(args[++i]);
    else if (a == "-mp_win" && nextIsValue(args, i)) op.qcOptions.mpWindow = std::stoi(args[++i]);
    else if (a == "-msec_tol" && nextIsValue(args, i)) op.qcOptions.msecTol = num(args[++i]);
    else if (a == "-pos_conv" && nextIsValue(args, i)) op.qcOptions.positionConvM = num(args[++i]);
    else if (a == "-pos_h_min" && nextIsValue(args, i)) op.qcOptions.positionHMinM = num(args[++i]);
    else if (a == "-pos_h_max" && nextIsValue(args, i)) op.qcOptions.positionHMaxM = num(args[++i]);
    else if (a == "-pos_jump" && nextIsValue(args, i)) op.qcOptions.positionJumpM = num(args[++i]);
    else if ((a == "-ion_jump" || a == "-iod_jump" || a == "-gap_mn" || a == "-gap_mx" || a == "-code_sigmas" || a == "-mpca_as") && nextIsValue(args, i)) {
      double v = num(args[++i]);
      if (a == "-ion_jump") op.qcOptions.ionJumpCM = v;
      else if (a == "-iod_jump") op.qcOptions.iodJumpCMPerMin = v;
      else if (a == "-gap_mn") op.qcOptions.gapMinMinutes = v;
      else if (a == "-gap_mx") op.qcOptions.gapMaxNoNavMinutes = v;
      else if (a == "-code_sigmas") op.qcOptions.codeSigmas = v;
      else if (a == "-mpca_as") op.qcOptions.mpCAASPercent = v;
    }
    else if (has(a, "-min_L") && nextIsValue(args, i)) op.qcOptions.minSNR[a.substr(6, 1)] = std::stoi(args[++i]);
    else if (has(a, "-mp") && a.find("_rms") != std::string::npos && nextIsValue(args, i)) op.qcOptions.mpRMSCM[a.substr(3,2)] = num(args[++i]);
    else if ((a == "+eph" || a == "-eph") && nextIsValue(args, i)) {
      auto vals = split(args[++i]);
      auto& dst = (a[0] == '+') ? op.ephIncludeTypes : op.ephExcludeTypes;
      for (auto v : vals) { std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return std::toupper(c); }); dst.push_back(v); }
    }
    else if (isSystemFilter(a)) {
      auto sys = sysName(a[1]);
      auto& sel = op.svSelections[sys];
      sel.specified = true;
      sel.includeMode = a[0] == '+';
      if (a.size() == 2) {
        sel.all = true;
        sel.prns.clear();
      } else {
        sel.all = false;
        sel.prns = parsePRNList(a.substr(2));
      }
    }
    else if (has(a, "-O.")) {
      auto au = lower(a);
      if (au.find("mod_wf") != std::string::npos) op.obsHeaderEdits[a] = collectModWF(args, i);
      else if (au.find("def_wf") != std::string::npos) op.obsHeaderEdits[a] = collectN(args, i, 2);
      else if (au.find("px") != std::string::npos || au.find("pg") != std::string::npos || au.find("pe") != std::string::npos || au.find("slant") != std::string::npos || au == "+o.sl" || au == "-o.sl") op.obsHeaderEdits[a] = collectXYZ(args, i);
      else if (au.find("start") != std::string::npos || au == "-o.st") op.obsHeaderEdits[a] = collectN(args, i, 6);
      else if (au.find("dec") != std::string::npos && nextIsValue(args, i)) op.obsDecimation = dec(args[++i]);
      else if (nextIsValue(args, i)) op.obsHeaderEdits[a] = args[++i];
    }
    else if (has(a, "+O.")) {
      auto au = lower(a);
      if (au.find("mod_wf") != std::string::npos) op.obsHeaderEdits[a] = collectModWF(args, i);
      else if (au.find("def_wf") != std::string::npos) op.obsHeaderEdits[a] = collectN(args, i, 2);
      else if (au.find("px") != std::string::npos || au.find("pg") != std::string::npos || au.find("pe") != std::string::npos || au.find("slant") != std::string::npos || au == "+o.sl") op.obsHeaderEdits[a] = collectXYZ(args, i);
      else if (au.find("start") != std::string::npos || au == "+o.st") op.obsHeaderEdits[a] = collectN(args, i, 6);
      else if (au.find("summary") != std::string::npos) { op.obsSummaryAppend = true; if (nextIsValue(args, i)) op.obsSummaryTarget = args[++i]; }
      else if (nextIsValue(args, i)) op.obsHeaderEdits[a] = args[++i];
    }
    else if (has(a, "-N.") || has(a, "+N.")) {
      auto au = lower(a);
      if (au.find("dec") != std::string::npos && nextIsValue(args, i)) op.navDecimation = dec(args[++i]);
      else if (au.find("dutc") != std::string::npos || au.find("iona") != std::string::npos || au.find("ionb") != std::string::npos) op.navHeaderEdits[a] = collectUntilOption(args, i);
      else if (nextIsValue(args, i)) op.navHeaderEdits[a] = args[++i];
    }
    else if (has(a, "-M.") || has(a, "+M.")) {
      auto au = lower(a);
      if (au.find("dec") != std::string::npos && nextIsValue(args, i)) op.metDecimation = dec(args[++i]);
      else if ((au.find("model") != std::string::npos || au.find("position") != std::string::npos || au.find("comment") != std::string::npos) && nextIsValue(args, i)) op.metHeaderEdits[a] = args[++i];
      else if (nextIsValue(args, i)) op.metHeaderEdits[a] = args[++i];
    }
    else if (has(a, "-B.") || has(a, "+B.")) {
      auto au = lower(a);
      if (au.find("px") != std::string::npos || au.find("pg") != std::string::npos || au.find("pe") != std::string::npos || au.find("slant") != std::string::npos || au=="-b.sl") op.binexMetadata[a] = collectUntilOption(args, i);
      else if (nextIsValue(args, i)) op.binexMetadata[a] = args[++i];
    }
    else if (a == "+relax") op.relaxRinex = true;
    else if (a == "-relax") op.relaxRinex = false;
    else if (a == "+extend") op.extendRinex2 = true;
    else if (a == "-extend") op.extendRinex2 = false;
    else if (a == "+reformat") op.reformatRinex = true;
    else if (a == "-reformat") op.reformatRinex = false;
    else if (!a.empty() && (a[0] == '-' || a[0] == '+')) { /* accept unimplemented teqc switch for compatibility */ }
    else op.inputFiles.push_back(a);
  }
  for (auto& kv : op.obsHeaderEdits) op.headerEdits[kv.first] = kv.second;
  return op;
}

} // namespace ceqc::cli
