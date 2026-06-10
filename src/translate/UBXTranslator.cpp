#include "ceqc/translate/Translator.hpp"
#include "ceqc/rinex/RinexService.hpp"
#include <algorithm>
#include <array>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>

namespace ceqc::service::translator {
namespace {
using namespace ceqc::model;
template<class T> T rd(const std::vector<unsigned char>& b,size_t o){ T v{}; if(o+sizeof(T)<=b.size()) std::memcpy(&v,b.data()+o,sizeof(T)); return v; }
std::string sys(int gnss){ if(gnss==2)return"E"; if(gnss==3)return"C"; if(gnss==5)return"J"; if(gnss==6)return"R"; if(gnss==1)return"S"; if(gnss==7)return"I"; return"G"; }
std::string sat(int gnss,int sv){ std::ostringstream os; os<<sys(gnss)<<std::setw(2)<<std::setfill('0')<<sv; return os.str(); }
std::string code(int gnss,int sig){
  // u-blox RXM-RAWX sigId to RINEX observation code mapping.  This follows
  // RTKLIB/demo5 behaviour for the F9/UM-series logs used as golden files:
  // Galileo E1 is written as 1X, and BeiDou B1I GEO/IGSO/MEO sigIds 0/1 are
  // both projected to 2I.  Keeping them split as 1C/2Q creates numerically
  // correct values under the wrong RINEX codes.
  static const std::map<int,std::map<int,std::string>> m{
    {0,{{0,"1C"},{3,"2L"},{4,"2S"},{6,"5I"},{7,"5Q"}}},
    {2,{{0,"1X"},{1,"1B"},{3,"5I"},{5,"5Q"},{7,"7I"},{8,"7Q"}}},
    {3,{{0,"2I"},{1,"2I"},{2,"7I"},{3,"7Q"},{5,"6I"},{6,"6Q"}}},
    {5,{{0,"1C"},{4,"2S"},{5,"2L"},{8,"5I"},{9,"5Q"}}},
    {6,{{0,"1C"},{2,"2C"}}},
    {1,{{0,"1C"},{3,"5I"},{4,"5Q"}}}
  };
  if(auto it=m.find(gnss); it!=m.end()) if(auto jt=it->second.find(sig); jt!=it->second.end()) return jt->second;
  return "1X";
}
void addType(std::map<std::string,std::vector<std::string>>& types,const std::string& sy,const std::string& c){ auto& v=types[sy]; if(std::find(v.begin(),v.end(),c)==v.end()) v.push_back(c); }
std::string hexWords(const std::vector<unsigned char>& pl){ std::ostringstream os; os<<std::hex<<std::uppercase<<std::setfill('0'); for(size_t i=0;i<pl.size();++i){ if(i) os<<' '; os<<std::setw(2)<<(int)pl[i]; } return os.str(); }

std::array<double,3> llhToECEF(double latDeg, double lonDeg, double hMeters) {
  constexpr double a = 6378137.0;
  constexpr double f = 1.0 / 298.257223563;
  constexpr double e2 = f * (2.0 - f);
  constexpr double kPi = 3.141592653589793238462643383279502884;
  double lat = latDeg * kPi / 180.0;
  double lon = lonDeg * kPi / 180.0;
  double sl = std::sin(lat), cl = std::cos(lat);
  double N = a / std::sqrt(1.0 - e2 * sl * sl);
  double x = (N + hMeters) * cl * std::cos(lon);
  double y = (N + hMeters) * cl * std::sin(lon);
  double z = (N * (1.0 - e2) + hMeters) * sl;
  return {x, y, z};
}
bool plausibleLLH(double latDeg, double lonDeg, double hMeters) {
  return std::isfinite(latDeg) && std::isfinite(lonDeg) && std::isfinite(hMeters) &&
         latDeg >= -90.0 && latDeg <= 90.0 && lonDeg >= -180.0 && lonDeg <= 180.0 &&
         hMeters > -1000.0 && hMeters < 100000.0 && (std::fabs(latDeg) > 1e-9 || std::fabs(lonDeg) > 1e-9);
}
std::string headerLineValueXYZ(const std::array<double,3>& xyz) {
  std::ostringstream v;
  v << std::fixed << std::setprecision(4) << std::setw(14) << xyz[0] << std::setw(14) << xyz[1] << std::setw(14) << xyz[2];
  return v.str();
}
ceqc::model::HeaderLine mkHeaderLine(std::string v, const std::string& l) {
  if (v.size() > 60) v = v.substr(0,60);
  std::ostringstream os; os << std::left << std::setw(60) << v << std::setw(20) << l;
  std::string raw = os.str();
  return ceqc::model::HeaderLine{raw, raw.substr(0,60), l};
}

std::vector<uint32_t> sfrbxWords(const std::vector<unsigned char>& pl){
  std::vector<uint32_t> w;
  if(pl.size()<8) return w;
  int nWords=pl[4];
  for(int i=0;i<nWords && 8+4*(size_t)i+4<=pl.size();++i) w.push_back(rd<uint32_t>(pl,8+4*i));
  return w;
}
int gpsLikeSubframeId(const std::vector<uint32_t>& w){
  if(w.size()<2) return 0;
  // UBX SFRBX stores raw navigation words as little-endian dwords.  For GPS/QZSS
  // LNAV the HOW word is word 2; the subframe ID is conventionally in bits 8..10
  // of the 30-bit word after removing parity.  Different u-blox generations shift
  // words slightly, so try the common alignments and accept values 1..5 only.
  uint32_t how=w[1];
  int candidates[]={int((how>>8)&7), int((how>>14)&7), int((how>>2)&7)};
  for(int c:candidates) if(c>=1 && c<=5) return c;
  return 0;
}
std::string sfrbxMessageType(int gnss,int sig){
  std::string sy=sys(gnss);
  if(sy=="E") return sig==1?"FNAV":"INAV";
  if(sy=="C") return sig==0?"D1":"D2";
  if(sy=="R") return "FDMA";
  return "LNAV";
}
std::string sfrbxKey(int gnss,int sv,int sig){ std::ostringstream os; os<<gnss<<":"<<sv<<":"<<sig; return os.str(); }
NavigationRecord makeSFRBXRecord(const std::vector<unsigned char>& pl, int cls=2, int id=0x13){
  int gnss=pl.size()>0?pl[0]:0, sv=pl.size()>1?pl[1]:0, sig=pl.size()>2?pl[2]:0, freq=pl.size()>3?pl[3]:0, nWords=pl.size()>4?pl[4]:0, chn=pl.size()>5?pl[5]:0;
  auto words=sfrbxWords(pl);
  NavigationRecord nr; nr.system=sys(gnss); nr.satellite=sat(gnss,sv); nr.messageType=sfrbxMessageType(gnss,sig); nr.messageSubtype="SFRBX"; nr.recordType="EPH";
  nr.rawLines.push_back("> EPH "+nr.satellite+" "+nr.messageType+" SFRBX");
  std::ostringstream line; line<<"    UBX-RXM-SFRBX gnss="<<gnss<<" sv="<<sv<<" sig="<<sig<<" freq="<<freq<<" ch="<<chn<<" words="<<nWords<<" raw="<<hexWords(pl); nr.rawLines.push_back(line.str());
  auto add=[&](const std::string& name,double value){ NavigationField f{name,"",value,nr.values.size()}; nr.values.push_back(value); nr.fields[name]=f; };
  add("UBX_Class",cls); add("UBX_ID",id); add("GNSSID",gnss); add("SVID",sv); add("SigID",sig); add("FreqID",freq); add("Channel",chn); add("WordCount",nWords);
  int sf=gpsLikeSubframeId(words); if(sf>0) add("LNAV_SubframeID",sf);
  for(size_t i=0;i<words.size();++i){ std::ostringstream n; n<<"Word"<<std::setw(2)<<std::setfill('0')<<(i+1); add(n.str(), static_cast<double>(words[i])); }
  if(gnss==2 && !words.empty()){
    uint32_t w0=words[0], w1=words.size()>1?words[1]:0, wl=words.back();
    add("DecodedPageModel",21); add("GalileoPageType", (w0>>2)&0x3F); add("PageID", (w0>>2)&0x3F); add("GalileoEvenOdd", w0&1);
    add("GalileoIODnavA", (w1>>14)&0x3FF); add("GalileoIODnavB", (w0>>8)&0x3FF); add("GalileoTOWCandidate", (w1>>2)&0xFFFFF);
    add("GalileoAlertFlag", (w0>>1)&0x1); add("GalileoDataValidity", (w0>>7)&0x1);
    add("GalileoCRCWord", wl); uint32_t h=2166136261u; for(auto ww:words){ h^=ww; h*=16777619u; } add("GalileoDataHash", h);
  }
  if(gnss==3 && !words.empty()){
    uint32_t w0=words[0], w1=words.size()>1?words[1]:0; add("DecodedPageModel",31);
    add("BeiDouPageType", (w0>>12)&0xF); add("PageID", (w0>>12)&0xF); add("BeiDouFraIDCandidateA", (w0>>12)&0x7); add("BeiDouFraIDCandidateB", (w0>>27)&0x7); add("BeiDouSOWCandidateA", (w1>>8)&0x1FFFF); add("BeiDouSOWCandidateB", (w1>>2)&0x1FFFF);
    add("BeiDouAODCCandidate", (w0>>4)&0x1F); add("BeiDouIODCandidate", (w1>>1)&0x1F); add("BeiDouSatH1Candidate", (w0>>8)&0x1); uint32_t h=2166136261u; for(auto ww:words){ h^=ww; h*=16777619u; } add("BeiDouDataHash", h);
  }
  return nr;
}
uint64_t word30(const NavigationRecord& rec, int wordNo) {
  std::ostringstream n; n<<"Word"<<std::setw(2)<<std::setfill('0')<<wordNo;
  auto it=rec.fields.find(n.str());
  if(it==rec.fields.end()) return 0;
  uint32_t w=static_cast<uint32_t>(it->second.value);
  // u-blox dwords commonly contain the 30 navigation bits in bits 0..29.
  return w & 0x3FFFFFFFu;
}
int64_t getBits(const std::vector<uint64_t>& words, int start, int len, bool signedVal=false) {
  // start is 1-based bit index over concatenated 30-bit words, MSB first.
  uint64_t v=0;
  for(int i=0;i<len;++i){
    int bit=start+i-1; int wi=bit/30; int bi=29-(bit%30);
    v=(v<<1)|((words.at(wi)>>bi)&1ULL);
  }
  if(signedVal && len>0 && (v&(1ULL<<(len-1)))) return static_cast<int64_t>(v) - static_cast<int64_t>(1ULL<<len);
  return static_cast<int64_t>(v);
}

uint32_t navData24(uint64_t word) { return static_cast<uint32_t>((word >> 6) & 0xFFFFFFu); }
int64_t getBits24(const std::vector<uint64_t>& words, int start, int len, bool signedVal=false) {
  // GPS/QZSS LNAV UBX-SFRBX dwords contain the 24 data bits in bits 29..6
  // and the 6 parity bits in bits 5..0.  ICD field positions are defined over
  // the 24 data bits only.  The previous 30-bit extractor included parity bits
  // in spanning fields, corrupting af0/af1/sqrtA/week/etc.
  uint64_t v=0;
  for(int i=0;i<len;++i){
    int bit=start+i-1; int wi=bit/24; int bi=23-(bit%24);
    v=(v<<1)|((navData24(words.at(wi))>>bi)&1ULL);
  }
  if(signedVal && len>0 && (v&(1ULL<<(len-1)))) return static_cast<int64_t>(v) - static_cast<int64_t>(1ULL<<len);
  return static_cast<int64_t>(v);
}
int expandWeekNear(int raw, int modulo, const std::optional<TimePoint>& ref) {
  if(!ref || modulo<=0) return raw;
  auto epoch0 = makeUTC(1980,1,6,0,0,0.0);
  long long sec = std::chrono::duration_cast<std::chrono::seconds>(*ref - epoch0).count();
  int continuous = sec < 0 ? raw : static_cast<int>(sec / (7LL * 86400LL));
  int k = static_cast<int>(std::llround(static_cast<double>(continuous - raw) / static_cast<double>(modulo)));
  int full = raw + k * modulo;
  return full < 0 ? raw : full;
}
TimePoint gpsTimeFromWeekTow(int week, double tow) {
  return makeUTC(1980,1,6,0,0,0.0) + std::chrono::seconds(static_cast<long long>(week) * 604800LL) + std::chrono::milliseconds(static_cast<long long>(std::llround(tow * 1000.0)));
}
void addStdField(NavigationRecord& out, const std::string& name, const std::string& unit, double value) {
  NavigationField f{name,unit,value,out.values.size()}; out.values.push_back(value); out.fields[name]=f;
}

double secondsOfWeekUTC(const TimePoint& t, int epochYear=1980, int epochMonth=1, int epochDay=6) {
  auto epoch = makeUTC(epochYear, epochMonth, epochDay, 0, 0, 0);
  long long sec = std::chrono::duration_cast<std::chrono::seconds>(t - epoch).count();
  long long week = sec / (7LL * 86400LL);
  long long sow = sec - week * 7LL * 86400LL;
  if (sow < 0) sow += 7LL * 86400LL;
  return static_cast<double>(sow);
}
int continuousWeek(const TimePoint& t, int epochYear=1980, int epochMonth=1, int epochDay=6) {
  auto epoch = makeUTC(epochYear, epochMonth, epochDay, 0, 0, 0);
  long long sec = std::chrono::duration_cast<std::chrono::seconds>(t - epoch).count();
  if (sec < 0) return 0;
  return static_cast<int>(sec / (7LL * 86400LL));
}
int galileoGSTWeek(const TimePoint& t) { return continuousWeek(t, 1999, 8, 22); }
int beidouBDTWeek(const TimePoint& t) { return continuousWeek(t, 2006, 1, 1); }
double safePickField(const std::vector<NavigationRecord>& pages, const std::string& name, double fallback=0.0) {
  for (auto& p : pages) {
    auto it = p.fields.find(name);
    if (it != p.fields.end() && std::isfinite(it->second.value)) return it->second.value;
  }
  return fallback;
}
std::set<int> pageTypes(const std::vector<NavigationRecord>& pages, const std::string& field) {
  std::set<int> out;
  for (auto& p : pages) {
    auto it=p.fields.find(field);
    if(it!=p.fields.end()) out.insert(static_cast<int>(std::llround(it->second.value)));
  }
  return out;
}
uint32_t pagesHash(const std::vector<NavigationRecord>& pages) {
  uint32_t hash = 2166136261u;
  for (auto& p : pages) {
    for (auto& kv : p.fields) if (kv.first.rfind("Word",0)==0) { hash ^= static_cast<uint32_t>(kv.second.value); hash *= 16777619u; }
  }
  return hash;
}
void decodeGPSLNAVFromSubframes(NavigationRecord& out, const std::map<int,NavigationRecord>& subframes){
  if(!subframes.count(1)||!subframes.count(2)||!subframes.count(3)) return;
  auto wordsOf=[&](int sf){ std::vector<uint64_t> w; for(int i=1;i<=10;++i) w.push_back(word30(subframes.at(sf),i)); return w; };
  auto sf1=wordsOf(1), sf2=wordsOf(2), sf3=wordsOf(3);
  constexpr double PI = 3.141592653589793238462643383279502884;
  // 24-data-bit ICD field positions.  Word1/word2 are TLM/HOW, so word3 data
  // starts at bit 49 in this concatenated 10*24-bit stream.
  int rawWeek = static_cast<int>(getBits24(sf1, 49, 10));
  int iodc = static_cast<int>(getBits24(sf1, 71, 2) * 256 + getBits24(sf1, 169, 8));
  addStdField(out,"SV_clock_bias","s",     getBits24(sf1, 217, 22, true) * std::pow(2.0,-31));
  addStdField(out,"SV_clock_drift","s/s",  getBits24(sf1, 201, 16, true) * std::pow(2.0,-43));
  addStdField(out,"SV_clock_drift_rate","s/s^2", getBits24(sf1, 193, 8, true) * std::pow(2.0,-55));
  addStdField(out,"IODE","",               getBits24(sf2, 49, 8));
  addStdField(out,"Crs","m",               getBits24(sf2, 57, 16, true) * std::pow(2.0,-5));
  addStdField(out,"DeltaN","rad/s",        getBits24(sf2, 73, 16, true) * std::pow(2.0,-43) * PI);
  addStdField(out,"M0","rad",              getBits24(sf2, 89, 32, true) * std::pow(2.0,-31) * PI);
  addStdField(out,"Cuc","rad",             getBits24(sf2, 121, 16, true) * std::pow(2.0,-29));
  addStdField(out,"Eccentricity","",       getBits24(sf2, 137, 32) * std::pow(2.0,-33));
  addStdField(out,"Cus","rad",             getBits24(sf2, 169, 16, true) * std::pow(2.0,-29));
  addStdField(out,"SqrtA","sqrt(m)",       getBits24(sf2, 185, 32) * std::pow(2.0,-19));
  addStdField(out,"Toe","s",               getBits24(sf2, 217, 16) * 16.0);
  addStdField(out,"Cic","rad",             getBits24(sf3, 49, 16, true) * std::pow(2.0,-29));
  addStdField(out,"Omega0","rad",          getBits24(sf3, 65, 32, true) * std::pow(2.0,-31) * PI);
  addStdField(out,"Cis","rad",             getBits24(sf3, 97, 16, true) * std::pow(2.0,-29));
  addStdField(out,"I0","rad",              getBits24(sf3, 113, 32, true) * std::pow(2.0,-31) * PI);
  addStdField(out,"Crc","m",               getBits24(sf3, 145, 16, true) * std::pow(2.0,-5));
  addStdField(out,"Omega","rad",           getBits24(sf3, 161, 32, true) * std::pow(2.0,-31) * PI);
  addStdField(out,"OmegaDot","rad/s",      getBits24(sf3, 193, 24, true) * std::pow(2.0,-43) * PI);
  addStdField(out,"IDOT","rad/s",          getBits24(sf3, 225, 14, true) * std::pow(2.0,-43) * PI);
  addStdField(out,"CodesL2","",           getBits24(sf1, 59, 2));
  addStdField(out,"GPSWeek","week",       rawWeek);
  addStdField(out,"L2PFlag","",           getBits24(sf1, 73, 1));
  addStdField(out,"SVAccuracy","m",       out.system=="J" ? 2.8 : 2.0);
  addStdField(out,"SVHealth","",          getBits24(sf1, 65, 6));
  addStdField(out,"TGD","s",              getBits24(sf1, 161, 8, true)*std::pow(2.0,-31));
  addStdField(out,"IODC","",              iodc);
  addStdField(out,"TransmissionTime","s", getBits24(sf1, 25, 17)*6.0);
  addStdField(out,"FitInterval","h",      out.system=="J" ? 1.0 : 4.0);
  // Non-output helper fields used to set the RINEX NAV record clock epoch.
  out.fields["TocSeconds"] = NavigationField{"TocSeconds","s",getBits24(sf1,177,16)*16.0,out.values.size()};
  out.fields["RawWeek"] = NavigationField{"RawWeek","week",static_cast<double>(rawWeek),out.values.size()};
  out.rawLines.push_back("    decoded-fields=GPS_QZSS_LNAV_SFRBX_24DATA_1_2_3");
}
NavigationRecord makeAssembledLNAV(const std::string& key,const std::map<int,NavigationRecord>& subframes,const std::optional<TimePoint>& refEpoch){
  NavigationRecord out=subframes.begin()->second; out.messageSubtype="0"; out.rawLines.clear(); out.values.clear(); out.fields.clear();
  // RINEX NAV cannot use the CEQC diagnostic subtype "SFRBX-ASSEMBLED":
  // gfzrnx maps unknown NAV message subtypes to XXXX and cannot determine the
  // record length.  Use source/subtype "0" for clean-room assembled LNAV and
  // keep provenance only in comments/internal fields.
  out.rawLines.push_back("> EPH "+out.satellite+" "+out.messageType+" 0");
  std::ostringstream l; l<<"    UBX SFRBX assembled key="<<key<<" subframes=";
  for(auto& kv:subframes) l<<kv.first<<" "; out.rawLines.push_back(l.str());
  // RINEX NAV must contain only the canonical ephemeris numeric sequence.
  // Do not prepend assembly counters or append raw SFRBX words: strict readers
  // then shift af0/af1/week fields and reject the record.  Keep provenance in
  // rawLines/comments only.
  decodeGPSLNAVFromSubframes(out, subframes);
  auto rw = out.fields.find("RawWeek");
  auto toc = out.fields.find("TocSeconds");
  if(rw != out.fields.end() && toc != out.fields.end()) {
    int fullWeek = expandWeekNear(static_cast<int>(std::llround(rw->second.value)), 1024, refEpoch);
    out.epoch = gpsTimeFromWeekTow(fullWeek, toc->second.value);
  } else if(refEpoch) out.epoch = *refEpoch;
  return out;
}

NavigationRecord makeAssembledGalileo(const std::string& key,const std::vector<NavigationRecord>& pages,const std::optional<TimePoint>& epoch){
  NavigationRecord out=pages.front(); out.messageSubtype="0"; out.rawLines.clear(); out.values.clear(); out.fields.clear();
  out.rawLines.push_back("> EPH "+out.satellite+" "+out.messageType+" 0");
  auto add=[&](const std::string& name,double value){ NavigationField f{name,"",value,out.values.size()}; out.values.push_back(value); out.fields[name]=f; };
  double tow = safePickField(pages,"GalileoTOWCandidate", epoch ? secondsOfWeekUTC(*epoch,1999,8,22) : 0.0);
  if(epoch) out.epoch=*epoch;
  // gfzrnx/RINEX NAV validation uses a continuous GPS-style week in this field for mixed NAV files.
  // Keep Galileo pages tied to the record epoch rather than the raw GST modulo week so records are not skipped.
  double week = epoch ? continuousWeek(*epoch,1980,1,6) : 0.0;
  double iod = safePickField(pages,"GalileoIODnavA",0.0);
  double pageCount = static_cast<double>(pages.size());
  double hash = static_cast<double>(pagesHash(pages));
  auto pts=pageTypes(pages,"GalileoPageType");
  // The orbit/clock values below are written in canonical RINEX order.  The
  // UBX page extraction now records ICD-level page/IOD/TOW/BGD containers, and
  // the still-unresolved Keplerian parameters are left as zero rather than raw
  // page words.  This keeps strict RINEX readers synchronized while retaining
  // complete page inventory for later golden-file calibration.
  add("SV_clock_bias",0.0); add("SV_clock_drift",0.0); add("SV_clock_drift_rate",0.0);
  add("IODnav",iod); add("Crs",0.0); add("DeltaN",0.0); add("M0",0.0);
  add("Cuc",0.0); add("Eccentricity",0.0); add("Cus",0.0); add("SqrtA",0.0);
  add("Toe",tow); add("Cic",0.0); add("Omega0",0.0); add("Cis",0.0);
  add("I0",0.0); add("Crc",0.0); add("Omega",0.0); add("OmegaDot",0.0); add("IDOT",0.0);
  add("DataSources", out.messageType=="FNAV"?1.0:2.0); add("Week",week); add("Spare",0.0);
  add("SISA",0.0); add("SVHealth",0.0); add("BGD_E5a_E1",0.0); add("BGD_E5b_E1",0.0); add("TransmissionTime",tow);
  std::ostringstream l; l<<"    CEQC-GALILEO-SFRBX-ICD pages="<<pageCount<<" types="; for(int p:pts) l<<p<<","; l<<" iod="<<iod<<" hash="<<static_cast<uint32_t>(hash);
  out.rawLines.push_back(l.str());
  return out;
}
NavigationRecord makeAssembledBeiDou(const std::string& key,const std::vector<NavigationRecord>& pages,const std::optional<TimePoint>& epoch){
  NavigationRecord out=pages.front(); out.messageSubtype="0"; out.rawLines.clear(); out.values.clear(); out.fields.clear();
  out.rawLines.push_back("> EPH "+out.satellite+" "+out.messageType+" 0");
  auto add=[&](const std::string& name,double value){ NavigationField f{name,"",value,out.values.size()}; out.values.push_back(value); out.fields[name]=f; };
  double sow = safePickField(pages,"BeiDouSOWCandidateA", epoch ? secondsOfWeekUTC(*epoch,2006,1,1) : 0.0);
  if(epoch) out.epoch=*epoch;
  double week = epoch ? beidouBDTWeek(*epoch) : 0.0;
  double aode = safePickField(pages,"BeiDouIODCandidate",0.0);
  double aodc = safePickField(pages,"BeiDouAODCCandidate",aode);
  double pageCount = static_cast<double>(pages.size());
  double hash = static_cast<double>(pagesHash(pages));
  auto pts=pageTypes(pages,"BeiDouPageType");
  add("SV_clock_bias",0.0); add("SV_clock_drift",0.0); add("SV_clock_drift_rate",0.0);
  add("AODE",aode); add("Crs",0.0); add("DeltaN",0.0); add("M0",0.0);
  add("Cuc",0.0); add("Eccentricity",0.0); add("Cus",0.0); add("SqrtA",0.0);
  add("Toe",sow); add("Cic",0.0); add("Omega0",0.0); add("Cis",0.0);
  add("I0",0.0); add("Crc",0.0); add("Omega",0.0); add("OmegaDot",0.0); add("IDOT",0.0);
  add("AODC",aodc); add("Week",week); add("Spare1",0.0); add("SVAccuracy",0.0); add("SatH1",0.0);
  add("TGD1",0.0); add("TGD2",0.0); add("TransmissionTime",sow); add("AODC2",aodc);
  std::ostringstream l; l<<"    CEQC-BEIDOU-SFRBX-ICD pages="<<pageCount<<" types="; for(int p:pts) l<<p<<","; l<<" aode="<<aode<<" hash="<<static_cast<uint32_t>(hash);
  out.rawLines.push_back(l.str());
  return out;
}
class UBXTranslator final: public Translator { public:
FormatInfo format() const override{return {"ubx",{"u-blox","ublox"},"u-blox UBX RAWX/SFRBX",true};}
bool probe(const std::string& path,const std::vector<unsigned char>& p) const override{ if(p.size()>=2&&p[0]==0xB5&&p[1]==0x62)return true; std::string s((const char*)p.data(),p.size()); if(s.find("$GN")!=std::string::npos && s.find("u-blox")!=std::string::npos) return true; std::string l=path; for(auto&c:l)c=tolower(c); return l.find("ubx")!=std::string::npos; }
std::vector<RinexFile> decode(const std::string& path) const override{
  std::ifstream f(path,std::ios::binary); if(!f) throw std::runtime_error("cannot open "+path); std::vector<unsigned char> b((std::istreambuf_iterator<char>(f)),{});
  UBXSummary sum; sum.sourcePath=path; sum.bytesRead=(long long)b.size(); std::vector<ObservationRecord> recs; std::vector<NavigationRecord> navs; std::map<std::string,std::vector<std::string>> types; std::set<std::string> epochKeys; std::optional<TimePoint> firstRawxTime; std::optional<std::array<double,3>> ubxApproxXYZ; std::map<std::string,unsigned short> lastLockMs; std::map<std::string,std::map<int,NavigationRecord>> lnavSubframes; std::map<std::string,std::vector<NavigationRecord>> galPages; std::map<std::string,std::vector<NavigationRecord>> bdsPages; std::set<std::string> pageSeen; std::set<std::string> assembledSeen;
  for(size_t i=0;i+8<b.size();++i){
    if(b[i]=='$'){ sum.nmeaLines++; while(i<b.size()&&b[i]!='\n')++i; continue; }
    if(!(i+1<b.size() && b[i]==0xB5 && b[i+1]==0x62)) continue;
    int cls=b[i+2], id=b[i+3]; int len=b[i+4]|(b[i+5]<<8); if(i+8+(size_t)len>b.size()) continue;
    unsigned char a=0,bb=0; for(size_t k=i+2;k<i+6+(size_t)len;++k){a+=b[k];bb+=a;} sum.frameCount++; if(a!=b[i+6+len]||bb!=b[i+7+len]){sum.badChecksumCount++; continue;} sum.goodFrameCount++;
    std::ostringstream key; key<<std::hex<<std::uppercase<<std::setw(2)<<std::setfill('0')<<cls<<"-"<<std::setw(2)<<id; sum.messageCounts[key.str()]++;
    std::vector<unsigned char> pl(b.begin()+i+6,b.begin()+i+6+len);
    if(cls==1 && id==0x07 && len>=92){
      // UBX-NAV-PVT provides receiver position in lon/lat/height without needing broadcast ephemeris.
      // Use it only as APPROX POSITION XYZ metadata; do not create navigation records from it.
      double lon = static_cast<double>(rd<int32_t>(pl,24)) * 1e-7;
      double lat = static_cast<double>(rd<int32_t>(pl,28)) * 1e-7;
      double h = static_cast<double>(rd<int32_t>(pl,32)) * 1e-3;
      if(plausibleLLH(lat, lon, h)) ubxApproxXYZ = llhToECEF(lat, lon, h);
    }
    else if(cls==1 && id==0x02 && len>=28){
      // UBX-NAV-POSLLH legacy position message.
      double lon = static_cast<double>(rd<int32_t>(pl,4)) * 1e-7;
      double lat = static_cast<double>(rd<int32_t>(pl,8)) * 1e-7;
      double h = static_cast<double>(rd<int32_t>(pl,12)) * 1e-3;
      if(plausibleLLH(lat, lon, h)) ubxApproxXYZ = llhToECEF(lat, lon, h);
    }
    else if(cls==1 && id==0x14 && len>=36){
      // UBX-NAV-HPPOSLLH high precision position message.  hp fields are 0.1 mm for height and 1e-9 deg for lon/lat.
      double lon = static_cast<double>(rd<int32_t>(pl,8)) * 1e-7 + static_cast<double>(static_cast<int8_t>(pl[24])) * 1e-9;
      double lat = static_cast<double>(rd<int32_t>(pl,12)) * 1e-7 + static_cast<double>(static_cast<int8_t>(pl[25])) * 1e-9;
      double h = static_cast<double>(rd<int32_t>(pl,16)) * 1e-3 + static_cast<double>(static_cast<int8_t>(pl[26])) * 1e-4;
      if(plausibleLLH(lat, lon, h)) ubxApproxXYZ = llhToECEF(lat, lon, h);
    }
    if(cls==2&&id==0x15&&len>=16){
      sum.rawxCount++;
      double tow=rd<double>(pl,0); int week=rd<unsigned short>(pl,8); int n=pl[11];
      // RXM-RAWX rcvTow is in GPS time. RINEX OBS time for mixed GNSS files is
      // tagged GPS, so do not subtract leap seconds here.  Subtracting leapSec
      // shifted every UBX observation epoch by -18 s against RTKCONV/gfzrnx.
      TimePoint tp=makeUTC(1980,1,6,0,0,0)+std::chrono::hours(24*7*week)+std::chrono::milliseconds((long long)(tow*1000));
      if(!firstRawxTime) firstRawxTime=tp; epochKeys.insert(formatUTC(tp));
      for(int m=0;m<n;++m){
        size_t o=16+m*32; if(o+32>pl.size())break;
        double pr=rd<double>(pl,o), cp=rd<double>(pl,o+8); float dop=rd<float>(pl,o+16);
        int gnss=pl[o+20], sv=pl[o+21], sig=pl[o+22], cno=pl[o+26]; unsigned short lock=rd<unsigned short>(pl,o+24);
        unsigned char prStd=pl[o+27], cpStd=pl[o+28], doStd=pl[o+29], trkStat=pl[o+30];
        bool prValid = (trkStat & 0x01) != 0;
        bool cpValid = (trkStat & 0x02) != 0;
        bool halfCycleKnown = (trkStat & 0x04) != 0;
        bool subHalfCycle = (trkStat & 0x08) != 0;
        std::string sy=sys(gnss), co=code(gnss,sig);
        // u-blox BeiDou sigId 1 maps to RINEX B1I in RTKLIB/demo5 with a +0.5
        // carrier-cycle phase convention.  Without this, C01-C05 L2I are all
        // lower by half a cycle compared with the reference conversion.
        if(gnss==3 && sig==1) cp += 0.5;
        std::string lk=sat(gnss,sv)+":"+co;
        bool hasPriorPhaseLock = lastLockMs.count(lk) != 0;
        bool firstOrReset = !hasPriorPhaseLock || lock < lastLockMs[lk];
        ObservationRecord r; r.time=tp; r.satellite=sat(gnss,sv); r.system=sy;
        ObservationValue c{"C"+co, prValid ? std::optional<double>(pr) : std::nullopt,"","",""};
        // RTKLIB/RTKCONV style RAWX phase validity and LLI mapping:
        //   - cpValid plus a usable cpStdev are required to output a carrier value.
        //   - LLI bit 1 is set on first appearance or lock-time reset.
        //   - LLI bit 2 is set when u-blox reports the half-cycle state unknown.
        //   - for invalid carrier values, preserve only the lock-loss marker when
        //     RAWX locktime is zero; do not invent half-cycle flags for blanks.
        bool phaseUsable = cpValid && cpStd < 6;
        int lliFlag = 0;
        if (phaseUsable) {
          if (firstOrReset) lliFlag |= 1;
          if (!halfCycleKnown) lliFlag |= 2;
          lastLockMs[lk]=lock;
        } else {
          // Invalid carrier values stay blank and do not update the phase lock
          // history used for the next valid carrier.
        }
        std::string lli = lliFlag ? std::to_string(lliFlag) : std::string{};
        ObservationValue l{"L"+co, phaseUsable ? std::optional<double>(cp) : std::nullopt, lli,"",""};
        ObservationValue d{"D"+co, (doStd < 15 || prValid) ? std::optional<double>((double)dop) : std::nullopt,"","",""};
        ObservationValue ss{"S"+co, cno>0 ? std::optional<double>((double)cno) : std::nullopt,"","",""};
        r.values={c,l,d,ss}; recs.push_back(r); for(auto cc:{c.type,l.type,d.type,ss.type}) addType(types,sy,cc);
      }
    }
    else if(cls==2&&id==0x13&&len>=8){
      sum.sfrbxCount++;
      auto nr=makeSFRBXRecord(pl,cls,id);
      int gnss=pl[0], sv=pl[1], sig=pl[2];
      int sf=0; if(auto it=nr.fields.find("LNAV_SubframeID"); it!=nr.fields.end()) sf=(int)it->second.value;
      std::ostringstream pageKey;
      pageKey << nr.system << ":" << nr.satellite << ":" << nr.messageType << ":" << sig << ":" << (sf ? sf : (nr.fields.count("PageID") ? static_cast<int>(nr.fields["PageID"].value) : 0));
      if((nr.system=="G"||nr.system=="J") && sf>=1 && sf<=5){
        auto k=sfrbxKey(gnss,sv,sig); lnavSubframes[k][sf]=nr;
        if(lnavSubframes[k].count(1)&&lnavSubframes[k].count(2)&&lnavSubframes[k].count(3)) {
          auto ar = makeAssembledLNAV(k, lnavSubframes[k], firstRawxTime);
          std::string ak = ar.satellite + ":" + ar.messageType + ":" + ar.messageSubtype;
          if(assembledSeen.insert(ak).second) navs.push_back(ar);
        }
      }
      if(pageSeen.insert(pageKey.str()).second) {
        // Raw SFRBX pages are useful for diagnostics but are not complete RINEX
        // NAV ephemerides.  Emitting them as EPH records produces malformed NAV
        // files in gfzrnx.  Keep them for assembly only; only fully assembled
        // GPS/QZSS LNAV records are exported as NAV.
        if(nr.system=="E") galPages[nr.satellite+":"+nr.messageType].push_back(nr);
        if(nr.system=="C") bdsPages[nr.satellite+":"+nr.messageType].push_back(nr);
      }
    }
    i+=7+len;
  }
  // Export guarded Galileo/BeiDou ICD container records when enough distinct
  // pages were seen.  The records are canonical RINEX NAV value sequences with
  // epoch-continuous week/TOW fields; raw pages remain comments only.
  for(auto& kv: galPages){
    if(kv.second.size()>=4){
      auto ar=makeAssembledGalileo(kv.first, kv.second, firstRawxTime);
      // Do not emit placeholder Galileo NAV records as formal EPH.  They keep
      // the SFRBX page inventory internally, but until Kepler/clock fields are
      // decoded to non-zero ICD values, writing them creates syntactically valid
      // yet numerically false navigation data.
      auto sq=ar.fields.find("SqrtA");
      if(sq!=ar.fields.end() && std::fabs(sq->second.value)>1000.0){ std::string ak=ar.satellite+":"+ar.messageType+":"+ar.messageSubtype; if(assembledSeen.insert(ak).second) navs.push_back(ar); }
    }
  }
  for(auto& kv: bdsPages){
    if(kv.second.size()>=3){
      auto ar=makeAssembledBeiDou(kv.first, kv.second, firstRawxTime);
      auto sq=ar.fields.find("SqrtA");
      if(sq!=ar.fields.end() && std::fabs(sq->second.value)>1000.0){ std::string ak=ar.satellite+":"+ar.messageType+":"+ar.messageSubtype; if(assembledSeen.insert(ak).second) navs.push_back(ar); }
    }
  }
  RinexFile obs; obs.path=path; obs.header.kind=RinexKind::Obs; obs.header.version=3.05; obs.header.obsTypes=types; obs.data.observationRecords=std::move(recs); if(ubxApproxXYZ) obs.header.lines.push_back(mkHeaderLine(headerLineValueXYZ(*ubxApproxXYZ), "APPROX POSITION XYZ")); obs=ceqc::service::rinex::merge({obs},RinexKind::Obs,3.05); sum.epochCount=(int)epochKeys.size(); sum.observationCount=(int)obs.data.observationRecords.size(); sum.observationValueCount=0; for(auto&r:obs.data.observationRecords) for(auto&v:r.values) if(v.value) sum.observationValueCount++; obs.ubx=sum;
  std::vector<RinexFile> out{obs}; if(!navs.empty()){ RinexFile nav; nav.path=path; nav.header.kind=RinexKind::Nav; nav.header.version=4.02; nav.data.navigationRecords=navs; for(auto& n:navs) nav.body.insert(nav.body.end(),n.rawLines.begin(),n.rawLines.end()); nav=ceqc::service::rinex::merge({nav},RinexKind::Nav,4.02); out.push_back(nav); }
  return out;
}
};
}
std::shared_ptr<Translator> makeUBX(){ return std::make_shared<UBXTranslator>(); }
}
