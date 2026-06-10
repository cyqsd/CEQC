#include "ceqc/translate/Translator.hpp"
#include "ceqc/rinex/RinexService.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <optional>
#include <regex>
#include <set>
#include <sstream>

namespace ceqc::service::translator {
namespace {
using namespace ceqc::model;

uint32_t crc24q(const std::vector<unsigned char>& data){
  uint32_t crc=0;
  for(auto b:data){ crc^=(uint32_t)b<<16; for(int i=0;i<8;++i){ crc<<=1; if(crc&0x1000000) crc^=0x1864CFB; } }
  return crc&0xFFFFFF;
}
struct Bits{
  const std::vector<unsigned char>& d; size_t p=0;
  int rem() const { return static_cast<int>(d.size()*8-p); }
  uint64_t u(int n){ uint64_t v=0; for(int i=0;i<n&&p<d.size()*8;++i){ v=(v<<1)|((d[p/8]>>(7-p%8))&1); ++p; } return v; }
  int64_t s(int n){ auto v=u(n); if(n && (v&(1ULL<<(n-1)))) return (int64_t)v-(int64_t)(1ULL<<n); return (int64_t)v; }
  int64_t sm(int n){ auto sign=u(1); auto mag=u(n-1); return sign?-static_cast<int64_t>(mag):static_cast<int64_t>(mag); }
  void skip(int n){ p=std::min(p+static_cast<size_t>(n), d.size()*8); }
};
std::string sat(const std::string& sys,int prn){ std::ostringstream os; os<<sys<<std::setw(2)<<std::setfill('0')<<prn; return os.str(); }
std::string hexWords(const std::vector<unsigned char>& pl){ std::ostringstream os; os<<std::hex<<std::uppercase<<std::setfill('0'); for(size_t i=0;i<pl.size();++i){ if(i) os<<' '; os<<std::setw(2)<<(int)pl[i]; } return os.str(); }
TimePoint inferAnchor(const std::string& path){
  std::regex yd(R"((20[0-9]{2})([0-9]{3}))"); std::smatch m; int hour=0;
  std::string base=path.substr(path.find_last_of("/\\")==std::string::npos?0:path.find_last_of("/\\")+1);
  std::regex pre(R"(^([0-9]{2})([0-9]{2})([0-9]{2})([0-9]{2}))"); std::smatch hm; if(std::regex_search(base,hm,pre)) hour=std::stoi(hm[4]);
  if(std::regex_search(base,m,yd)){ int y=std::stoi(m[1]); int doy=std::stoi(m[2]); return makeUTC(y,1,1,hour,0,0)+std::chrono::hours(24*(doy-1)); }
  return makeUTC(2026,1,1,0,0,0);
}
TimePoint rtcmTime(TimePoint anchor,const std::string& sys,long long rawMs){
  constexpr long long dayMs=86400000LL; long long ms=rawMs;
  if(sys=="R"){ ms=(rawMs & ((1LL<<27)-1)) - 3LL*3600*1000 + 18*1000; }
  else if(sys=="C") ms=rawMs+14*1000;
  ms%=dayMs; if(ms<0) ms+=dayMs;
  auto day=std::chrono::floor<std::chrono::hours>(anchor.time_since_epoch());
  auto days=std::chrono::duration_cast<std::chrono::hours>(day).count()/24;
  auto midnight=TimePoint{}+std::chrono::hours(days*24);
  return midnight+std::chrono::milliseconds(ms);
}
std::vector<int> maskIdx(uint64_t mask,int width){ std::vector<int> out; for(int i=0;i<width;++i) if(mask&(1ULL<<(width-1-i))) out.push_back(i+1); return out; }
std::string msmCode(const std::string& sys,int sig){
  static const std::map<std::string,std::map<int,std::string>> m{
    {"G",{{2,"1C"},{10,"2L"},{15,"2S"},{16,"2L"},{17,"2X"},{22,"5I"},{23,"5Q"},{24,"5X"}}},
    {"R",{{2,"1C"},{3,"1P"},{8,"2C"},{9,"2P"}}},
    {"E",{{2,"1C"},{3,"1A"},{4,"1B"},{5,"1X"},{12,"5I"},{13,"5Q"},{14,"5X"},{18,"7I"},{19,"7Q"},{20,"7X"},{22,"8I"},{23,"8Q"},{24,"8X"}}},
    {"C",{{2,"2I"},{3,"2Q"},{4,"2X"},{8,"6I"},{9,"6Q"},{10,"6X"},{12,"7I"},{13,"7Q"},{14,"7X"}}},
    {"J",{{2,"1C"},{8,"2S"},{9,"2L"},{10,"2X"},{15,"5I"},{16,"5Q"},{17,"5X"}}}
  };
  if(auto it=m.find(sys); it!=m.end()) if(auto jt=it->second.find(sig); jt!=it->second.end()) return jt->second;
  return "1X";
}
double frequency(const std::string& sys,const std::string& code){ if(code.empty()) return 0; char b=code[0]; if(sys=="G"||sys=="J"){ if(b=='1')return 1575.42e6; if(b=='2')return 1227.60e6; if(b=='5')return 1176.45e6; } if(sys=="R"){ if(b=='1')return 1602e6; if(b=='2')return 1246e6; } if(sys=="E"){ if(b=='1')return 1575.42e6; if(b=='5')return 1176.45e6; if(b=='7')return 1207.14e6; if(b=='8')return 1191.795e6; } if(sys=="C"){ if(b=='2')return 1561.098e6; if(b=='6')return 1268.52e6; if(b=='7')return 1207.14e6; } return 0; }
std::string obsKey(const TimePoint& t,const std::string& sat){ return formatUTC(t)+":"+sat; }
void addType(std::map<std::string,std::vector<std::string>>& types,const std::string& sys,const std::string& code){ auto& v=types[sys]; if(std::find(v.begin(),v.end(),code)==v.end()) v.push_back(code); }
void addVal(ObservationRecord& r,const std::string& type,double val){ ObservationValue v; v.type=type; v.value=val; r.values.push_back(v); }
void addVal(ObservationRecord& r,const std::string& type,double val,const std::string& lli,const std::string& ssi=""){
  ObservationValue v; v.type=type; v.value=val; v.lli=lli; v.ssi=ssi; r.values.push_back(v);
}
struct RtcmLockState { int lock=-1; int half=-1; };
std::string rtcmLliFor(const std::string& key,int lock,int half,std::map<std::string,RtcmLockState>& states){
  int flag=0;
  auto it=states.find(key);
  if(it==states.end()) flag |= 1;
  else {
    // MSM lock indicators are logarithmic/saturated counters.  A decrease means
    // loss of lock; a half-cycle ambiguity state change also starts a new arc.
    if(lock < it->second.lock) flag |= 1;
    if(half != it->second.half) flag |= 1;
  }
  if(half) flag |= 2;
  states[key] = {lock,half};
  return flag ? std::to_string(flag) : std::string{};
}


void addField(NavigationRecord& nr,const std::string& name,const std::string& unit,double value){
  NavigationField f{name,unit,value,nr.values.size()}; nr.values.push_back(value); nr.fields[name]=f;
}

std::optional<double> navFieldValue(const NavigationRecord& nr, const std::string& name){
  auto it=nr.fields.find(name);
  if(it==nr.fields.end()) return std::nullopt;
  return it->second.value;
}
double rtcmObsNominalFrequency(const std::string& sys,const std::string& code){ return frequency(sys,code); }
double glonassFrequencyFromSlot(char band, int k){
  if(band=='1') return (1602.0 + 0.5625*k)*1e6;
  if(band=='2') return (1246.0 + 0.4375*k)*1e6;
  return 0.0;
}
void normalizeGlonassFdmaObservationCycles(std::map<std::string,ObservationRecord>& recMap, const std::vector<NavigationRecord>& navs){
  std::map<std::string,int> slotBySat;
  for(const auto& n:navs){
    if(n.system!="R" || n.satellite.empty()) continue;
    auto fk=navFieldValue(n,"FrequencyNumber");
    if(fk) slotBySat[n.satellite]=(int)std::lround(*fk);
  }
  if(slotBySat.empty()) return;
  for(auto& kv:recMap){
    auto& r=kv.second;
    if(r.system!="R") continue;
    auto sit=slotBySat.find(r.satellite);
    if(sit==slotBySat.end()) continue;
    int k=sit->second;
    for(auto& v:r.values){
      if(!v.value || v.type.size()<2) continue;
      if(v.type[0]!='L' && v.type[0]!='D') continue;
      char band=v.type[1];
      double nominal=rtcmObsNominalFrequency("R",std::string(1,band));
      double actual=glonassFrequencyFromSlot(band,k);
      if(nominal>0.0 && actual>0.0) *v.value *= actual/nominal;
    }
  }
}
double sod(double sec){
  while(sec < 0.0) sec += 86400.0;
  while(sec >= 86400.0) sec -= 86400.0;
  return sec;
}
TimePoint dayTime(TimePoint anchor,double sec){
  auto hrs=std::chrono::floor<std::chrono::hours>(anchor.time_since_epoch());
  auto days=std::chrono::duration_cast<std::chrono::hours>(hrs).count()/24;
  auto midnight=TimePoint{}+std::chrono::hours(days*24);
  sec = sod(sec);
  return midnight+std::chrono::milliseconds((long long)std::llround(sec*1000.0));
}
NavigationRecord parseGPS1019(Bits& br, TimePoint anchor){
  constexpr double PI=3.14159265358979323846;
  int prn=(int)br.u(6); int week=(int)br.u(10); int ura=(int)br.u(4); int codeL2=(int)br.u(2); (void)codeL2;
  double idot=br.s(14)*std::pow(2.0,-43)*PI; int iode=(int)br.u(8); double toc=br.u(16)*16.0;
  double af2=br.s(8)*std::pow(2.0,-55), af1=br.s(16)*std::pow(2.0,-43), af0=br.s(22)*std::pow(2.0,-31); int iodc=(int)br.u(10);
  double crs=br.s(16)*std::pow(2.0,-5); double dn=br.s(16)*std::pow(2.0,-43)*PI; double m0=br.s(32)*std::pow(2.0,-31)*PI;
  double cuc=br.s(16)*std::pow(2.0,-29); double ecc=br.u(32)*std::pow(2.0,-33); double cus=br.s(16)*std::pow(2.0,-29); double sqrtA=br.u(32)*std::pow(2.0,-19); double toe=br.u(16)*16.0;
  double cic=br.s(16)*std::pow(2.0,-29); double om0=br.s(32)*std::pow(2.0,-31)*PI; double cis=br.s(16)*std::pow(2.0,-29); double i0=br.s(32)*std::pow(2.0,-31)*PI;
  double crc=br.s(16)*std::pow(2.0,-5); double omega=br.s(32)*std::pow(2.0,-31)*PI; double omdot=br.s(24)*std::pow(2.0,-43)*PI; double tgd=br.s(8)*std::pow(2.0,-31); int health=(int)br.u(6); br.skip(1); br.skip(1); double tt=br.u(17); br.skip(2); br.skip(1);
  NavigationRecord nr; nr.system="G"; nr.satellite=sat("G",prn); nr.recordType="EPH"; nr.messageType="LNAV"; nr.messageSubtype="RTCM1019"; nr.epoch=dayTime(anchor,toc);
  addField(nr,"SV_clock_bias","s",af0); addField(nr,"SV_clock_drift","s/s",af1); addField(nr,"SV_clock_drift_rate","s/s^2",af2);
  addField(nr,"IODE","",iode); addField(nr,"Crs","m",crs); addField(nr,"DeltaN","rad/s",dn); addField(nr,"M0","rad",m0);
  addField(nr,"Cuc","rad",cuc); addField(nr,"Eccentricity","",ecc); addField(nr,"Cus","rad",cus); addField(nr,"SqrtA","sqrt(m)",sqrtA);
  addField(nr,"Toe","s",toe); addField(nr,"Cic","rad",cic); addField(nr,"Omega0","rad",om0); addField(nr,"Cis","rad",cis);
  addField(nr,"I0","rad",i0); addField(nr,"Crc","m",crc); addField(nr,"Omega","rad",omega); addField(nr,"OmegaDot","rad/s",omdot);
  addField(nr,"IDOT","rad/s",idot); addField(nr,"CodesL2","",codeL2); addField(nr,"Week","week",week); addField(nr,"L2PFlag","",0);
  addField(nr,"SVAccuracy","",ura); addField(nr,"SVHealth","",health); addField(nr,"TGD","s",tgd); addField(nr,"IODC","",iodc); addField(nr,"TransmissionTime","s",tt); addField(nr,"FitInterval","",0);
  nr.rawLines.push_back("> EPH "+nr.satellite+" LNAV RTCM1019");
  return nr;
}

NavigationRecord parseBDS1042(Bits& br, TimePoint anchor){
  constexpr double PI=3.14159265358979323846;
  int prn=(int)br.u(6); int week=(int)br.u(13); int ura=(int)br.u(4); double idot=br.s(14)*std::pow(2.0,-43)*PI; int aode=(int)br.u(5); double toc=br.u(17)*8.0;
  double af2=br.s(11)*std::pow(2.0,-66), af1=br.s(22)*std::pow(2.0,-50), af0=br.s(24)*std::pow(2.0,-33); int aodc=(int)br.u(5);
  double crs=br.s(18)*std::pow(2.0,-6); double dn=br.s(16)*std::pow(2.0,-43)*PI; double m0=br.s(32)*std::pow(2.0,-31)*PI;
  double cuc=br.s(18)*std::pow(2.0,-31); double ecc=br.u(32)*std::pow(2.0,-33); double cus=br.s(18)*std::pow(2.0,-31); double sqrtA=br.u(32)*std::pow(2.0,-19); double toe=br.u(17)*8.0;
  double cic=br.s(18)*std::pow(2.0,-31); double om0=br.s(32)*std::pow(2.0,-31)*PI; double cis=br.s(18)*std::pow(2.0,-31); double i0=br.s(32)*std::pow(2.0,-31)*PI;
  double crc=br.s(18)*std::pow(2.0,-6); double omega=br.s(32)*std::pow(2.0,-31)*PI; double omdot=br.s(24)*std::pow(2.0,-43)*PI;
  double tgd1=br.s(10)*1e-10, tgd2=br.s(10)*1e-10; int health=(int)br.u(1);
  NavigationRecord nr; nr.system="C"; nr.satellite=sat("C",prn); nr.recordType="EPH"; nr.messageType="D1"; nr.messageSubtype="RTCM1042"; nr.epoch=dayTime(anchor,toc);
  addField(nr,"SV_clock_bias","s",af0); addField(nr,"SV_clock_drift","s/s",af1); addField(nr,"SV_clock_drift_rate","s/s^2",af2);
  addField(nr,"AODE","",aode); addField(nr,"Crs","m",crs); addField(nr,"DeltaN","rad/s",dn); addField(nr,"M0","rad",m0);
  addField(nr,"Cuc","rad",cuc); addField(nr,"Eccentricity","",ecc); addField(nr,"Cus","rad",cus); addField(nr,"SqrtA","sqrt(m)",sqrtA);
  addField(nr,"Toe","s",toe); addField(nr,"Cic","rad",cic); addField(nr,"Omega0","rad",om0); addField(nr,"Cis","rad",cis);
  addField(nr,"I0","rad",i0); addField(nr,"Crc","m",crc); addField(nr,"Omega","rad",omega); addField(nr,"OmegaDot","rad/s",omdot);
  addField(nr,"IDOT","rad/s",idot); addField(nr,"AODC","",aodc); addField(nr,"Week","week",week); addField(nr,"Spare1","",0);
  addField(nr,"SVAccuracy","",ura); addField(nr,"SatH1","",health); addField(nr,"TGD1","s",tgd1); addField(nr,"TGD2","s",tgd2); addField(nr,"TransmissionTime","s",toc);
  nr.rawLines.push_back("> EPH "+nr.satellite+" D1 RTCM1042"); return nr;
}

NavigationRecord parseGLO1020(Bits& br, TimePoint anchor){
  int prn=(int)br.u(6); int freqRaw=(int)br.u(5); br.skip(2+2); int tkH=(int)br.u(5), tkM=(int)br.u(6), tkS=(int)br.u(1); int bn=(int)br.u(1); br.skip(1); int tb=(int)br.u(7);
  double vx=br.sm(24)*std::pow(2.0,-20); double x=br.sm(27)*std::pow(2.0,-11); double ax=br.sm(5)*std::pow(2.0,-30);
  double vy=br.sm(24)*std::pow(2.0,-20); double y=br.sm(27)*std::pow(2.0,-11); double ay=br.sm(5)*std::pow(2.0,-30);
  double vz=br.sm(24)*std::pow(2.0,-20); double z=br.sm(27)*std::pow(2.0,-11); double az=br.sm(5)*std::pow(2.0,-30); br.skip(1);
  double gamma=br.sm(11)*std::pow(2.0,-40); br.skip(3); double tau=br.sm(22)*std::pow(2.0,-30);
  double local=tb*900.0, toeUtc=sod(local-3*3600.0);
  double frameUtc=sod(tkH*3600.0+tkM*60.0+tkS*30.0-3*3600.0);
  NavigationRecord nr; nr.system="R"; nr.satellite=sat("R",prn); nr.recordType="EPH"; nr.messageType="FDMA"; nr.messageSubtype="RTCM1020"; nr.epoch=dayTime(anchor,toeUtc);
  addField(nr,"TauN","s",tau); addField(nr,"GammaN","",gamma); addField(nr,"MessageFrameTime","s",frameUtc);
  addField(nr,"X","km",x); addField(nr,"VX","km/s",vx); addField(nr,"AX","km/s^2",ax); addField(nr,"Health","",bn);
  addField(nr,"Y","km",y); addField(nr,"VY","km/s",vy); addField(nr,"AY","km/s^2",ay); addField(nr,"FrequencyNumber","",freqRaw-7);
  addField(nr,"Z","km",z); addField(nr,"VZ","km/s",vz); addField(nr,"AZ","km/s^2",az); addField(nr,"AgeOfOperation","d",0);
  nr.rawLines.push_back("> EPH "+nr.satellite+" FDMA RTCM1020");
  return nr;
}


NavigationRecord parseQZSS1044(Bits& br, TimePoint anchor){
  constexpr double PI=3.14159265358979323846;
  int prn=(int)br.u(4); int tocRaw=(int)br.u(16); double af2=br.s(8)*std::pow(2.0,-55), af1=br.s(16)*std::pow(2.0,-43), af0=br.s(22)*std::pow(2.0,-31); int iode=(int)br.u(8); double crs=br.s(16)*std::pow(2.0,-5); double dn=br.s(16)*std::pow(2.0,-43)*PI; double m0=br.s(32)*std::pow(2.0,-31)*PI;
  double cuc=br.s(16)*std::pow(2.0,-29); double ecc=br.u(32)*std::pow(2.0,-33); double cus=br.s(16)*std::pow(2.0,-29); double sqrtA=br.u(32)*std::pow(2.0,-19); double toe=br.u(16)*16.0;
  double cic=br.s(16)*std::pow(2.0,-29); double om0=br.s(32)*std::pow(2.0,-31)*PI; double cis=br.s(16)*std::pow(2.0,-29); double i0=br.s(32)*std::pow(2.0,-31)*PI;
  double crc=br.s(16)*std::pow(2.0,-5); double omega=br.s(32)*std::pow(2.0,-31)*PI; double omdot=br.s(24)*std::pow(2.0,-43)*PI; double idot=br.s(14)*std::pow(2.0,-43)*PI; int codeL2=(int)br.u(2); int week=(int)br.u(10); int ura=(int)br.u(4); int health=(int)br.u(6); double tgd=br.s(8)*std::pow(2.0,-31); int iodc=(int)br.u(10); double tt=br.u(17);
  NavigationRecord nr; nr.system="J"; nr.satellite=sat("J",prn); nr.recordType="EPH"; nr.messageType="LNAV"; nr.messageSubtype="RTCM1044"; nr.epoch=dayTime(anchor,tocRaw*16.0);
  addField(nr,"SV_clock_bias","s",af0); addField(nr,"SV_clock_drift","s/s",af1); addField(nr,"SV_clock_drift_rate","s/s^2",af2); addField(nr,"IODE","",iode); addField(nr,"Crs","m",crs); addField(nr,"DeltaN","rad/s",dn); addField(nr,"M0","rad",m0); addField(nr,"Cuc","rad",cuc); addField(nr,"Eccentricity","",ecc); addField(nr,"Cus","rad",cus); addField(nr,"SqrtA","sqrt(m)",sqrtA); addField(nr,"Toe","s",toe); addField(nr,"Cic","rad",cic); addField(nr,"Omega0","rad",om0); addField(nr,"Cis","rad",cis); addField(nr,"I0","rad",i0); addField(nr,"Crc","m",crc); addField(nr,"Omega","rad",omega); addField(nr,"OmegaDot","rad/s",omdot); addField(nr,"IDOT","rad/s",idot); addField(nr,"CodesL2","",codeL2); addField(nr,"Week","week",week); addField(nr,"L2PFlag","",0); addField(nr,"SVAccuracy","",ura); addField(nr,"SVHealth","",health); addField(nr,"TGD","s",tgd); addField(nr,"IODC","",iodc); addField(nr,"TransmissionTime","s",tt); addField(nr,"Toc","s",tocRaw*16.0);
  nr.rawLines.push_back("> EPH "+nr.satellite+" LNAV RTCM1044"); return nr;
}

NavigationRecord parseGAL104x(Bits& br, TimePoint anchor, int msg){
  constexpr double PI=3.14159265358979323846;
  int prn=(int)br.u(6); int week=(int)br.u(12); int iod=(int)br.u(10); int sisa=(int)br.u(8); double idot=br.s(14)*std::pow(2.0,-43)*PI; double toc=br.u(14)*60.0; double af2=br.s(6)*std::pow(2.0,-59); double af1=br.s(21)*std::pow(2.0,-46); double af0=br.s(31)*std::pow(2.0,-34);
  double crs=br.s(16)*std::pow(2.0,-5); double dn=br.s(16)*std::pow(2.0,-43)*PI; double m0=br.s(32)*std::pow(2.0,-31)*PI; double cuc=br.s(16)*std::pow(2.0,-29); double ecc=br.u(32)*std::pow(2.0,-33); double cus=br.s(16)*std::pow(2.0,-29); double sqrtA=br.u(32)*std::pow(2.0,-19); double toe=br.u(14)*60.0; double cic=br.s(16)*std::pow(2.0,-29); double om0=br.s(32)*std::pow(2.0,-31)*PI; double cis=br.s(16)*std::pow(2.0,-29); double i0=br.s(32)*std::pow(2.0,-31)*PI; double crc=br.s(16)*std::pow(2.0,-5); double omega=br.s(32)*std::pow(2.0,-31)*PI; double omdot=br.s(24)*std::pow(2.0,-43)*PI; double bgdA=br.s(10)*std::pow(2.0,-32); double bgdB=br.s(10)*std::pow(2.0,-32); int health=(int)br.u(6);
  NavigationRecord nr; nr.system="E"; nr.satellite=sat("E",prn); nr.recordType="EPH"; nr.messageType=(msg==1045?"FNAV":"INAV"); nr.messageSubtype=(msg==1045?"RTCM1045":"RTCM1046"); nr.epoch=dayTime(anchor,toc);
  addField(nr,"SV_clock_bias","s",af0); addField(nr,"SV_clock_drift","s/s",af1); addField(nr,"SV_clock_drift_rate","s/s^2",af2); addField(nr,"IODnav","",iod); addField(nr,"Crs","m",crs); addField(nr,"DeltaN","rad/s",dn); addField(nr,"M0","rad",m0); addField(nr,"Cuc","rad",cuc); addField(nr,"Eccentricity","",ecc); addField(nr,"Cus","rad",cus); addField(nr,"SqrtA","sqrt(m)",sqrtA); addField(nr,"Toe","s",toe); addField(nr,"Cic","rad",cic); addField(nr,"Omega0","rad",om0); addField(nr,"Cis","rad",cis); addField(nr,"I0","rad",i0); addField(nr,"Crc","m",crc); addField(nr,"Omega","rad",omega); addField(nr,"OmegaDot","rad/s",omdot); addField(nr,"IDOT","rad/s",idot); addField(nr,"DataSources","",msg==1045?1:2); addField(nr,"Week","week",week); addField(nr,"Spare","",0); addField(nr,"SISA","m",sisa); addField(nr,"SVHealth","",health); addField(nr,"BGD_E5a_E1","s",bgdA); addField(nr,"BGD_E5b_E1","s",bgdB); addField(nr,"TransmissionTime","s",toc);
  nr.rawLines.push_back("> EPH "+nr.satellite+" "+nr.messageType+" "+nr.messageSubtype); return nr;
}

class RTCM3Translator final: public Translator { public:
FormatInfo format() const override{return {"rtcm3",{"rtcm","rtcm-3","rtcmv3"},"RTCM version 3",true};}
bool probe(const std::string& path,const std::vector<unsigned char>& p) const override{ if(p.size()>=3&&p[0]==0xD3&&(p[1]&0xFC)==0)return true; std::string l=path; for(auto&c:l)c=tolower(c); return l.find("rtcm")!=std::string::npos; }
std::vector<RinexFile> decode(const std::string& path) const override{
  std::ifstream f(path,std::ios::binary); if(!f) throw std::runtime_error("cannot open "+path); std::vector<unsigned char> b((std::istreambuf_iterator<char>(f)),{});
  RTCM3Summary sum; sum.sourcePath=path; sum.bytesRead=(long long)b.size(); TimePoint anchor=inferAnchor(path);
  std::map<std::string,ObservationRecord> recMap; std::map<std::string,std::vector<std::string>> types; std::set<std::string> epochs;
  std::vector<NavigationRecord> navs; RTCM3Station station; bool haveStation=false; std::string antenna="UNKNOWN", antSerial="", receiver="UNKNOWN", receiverVer="", receiverSerial="";
  std::map<std::string,RtcmLockState> msmLockStates;
  for(size_t i=0;i+6<b.size();++i){ if(b[i]!=0xD3) continue; int len=((b[i+1]&3)<<8)|b[i+2]; if(len<=0||i+6+(size_t)len>b.size()){sum.badLengthCount++; continue;} std::vector<unsigned char> msg(b.begin()+i,b.begin()+i+3+len); uint32_t want=((uint32_t)b[i+3+len]<<16)|((uint32_t)b[i+4+len]<<8)|b[i+5+len]; int mt=(b[i+3]<<4)|(b[i+4]>>4); sum.frameCount++; if(crc24q(msg)!=want){sum.badCRCCount++; continue;} sum.goodFrameCount++; sum.messageCounts[mt]++;
    std::vector<unsigned char> pl(b.begin()+i+3,b.begin()+i+3+len); Bits br{pl}; int type=(int)br.u(12);
    if(type==1005||type==1006){ int st=(int)br.u(12); int itrf=(int)br.u(6); br.skip(4); double X=br.s(38)*0.0001; br.skip(2); double Y=br.s(38)*0.0001; br.skip(2); double Z=br.s(38)*0.0001; double h=0; if(type==1006) h=br.u(16)*0.0001; station={st,itrf,X,Y,Z,h,true}; haveStation=true; sum.stations[st]=station; sum.stationCounts[st]++; }
    else if(type==1033){ int st=(int)br.u(12); (void)st; int n=(int)br.u(8); antenna.clear(); for(int k=0;k<n && br.rem()>=8;++k) antenna.push_back((char)br.u(8)); br.skip(8); int ns=(int)br.u(8); antSerial.clear(); for(int k=0;k<ns && br.rem()>=8;++k) antSerial.push_back((char)br.u(8)); int nr=(int)br.u(8); receiver.clear(); for(int k=0;k<nr && br.rem()>=8;++k) receiver.push_back((char)br.u(8)); int nv=(int)br.u(8); receiverVer.clear(); for(int k=0;k<nv && br.rem()>=8;++k) receiverVer.push_back((char)br.u(8)); int nrs=(int)br.u(8); receiverSerial.clear(); for(int k=0;k<nrs && br.rem()>=8;++k) receiverSerial.push_back((char)br.u(8)); }
    else if(type>=1071&&type<=1137){ int group=type/10, var=type%10; std::string sys = group==107?"G":(group==108?"R":(group==109?"E":(group==111?"J":(group==112?"C":"S")))); br.u(12); long long rawMs=(long long)br.u(30); br.skip(1+3+7+2+2+1+3); auto sats=maskIdx(br.u(64),64); auto sigs=maskIdx(br.u(32),32); std::vector<std::pair<int,int>> cells; for(size_t si=0; si<sats.size(); ++si) for(size_t gi=0; gi<sigs.size(); ++gi) if(br.u(1)) cells.push_back({(int)si,(int)gi});
      std::vector<int> rough(sats.size()), mod(sats.size()); std::vector<int64_t> rate(sats.size()); for(auto&x:rough)x=(int)br.u(8); if(var==5||var==7) for(size_t k=0;k<sats.size();++k) br.u(4); for(auto&x:mod)x=(int)br.u(10); if(var==5||var==7) for(auto&x:rate)x=br.s(14);
      TimePoint tp=rtcmTime(anchor,sys,rawMs); epochs.insert(formatUTC(tp));
      int finePrBits=(var==6||var==7)?20:15, finePhBits=(var==6||var==7)?24:22, lockBits=(var==6||var==7)?10:4, cnrBits=(var==6||var==7)?10:6;
      double finePrDen=(var==6||var==7)?std::pow(2.0,29):std::pow(2.0,24);
      double finePhDen=(var==6||var==7)?std::pow(2.0,31):std::pow(2.0,29);
      // RTCM MSM cell data are grouped by field, not interleaved by cell:
      //   all fine pseudoranges, all fine phase-ranges, all locks, all half-cycle
      //   indicators, all C/N0 values, then optional all fine phase-rate values.
      // Reading them per-cell shifts every field after the first cell and produces
      // apparently plausible code values but completely broken carrier arcs, which
      // shows up as kilometre-level MP residuals.
      const size_t ncells = cells.size();
      std::vector<int64_t> finePrs(ncells), finePhs(ncells), fineRates(ncells,0);
      std::vector<int> locks(ncells), halves(ncells), cnrs(ncells);
      for(size_t ci=0; ci<ncells; ++ci) finePrs[ci]=br.s(finePrBits);
      for(size_t ci=0; ci<ncells; ++ci) finePhs[ci]=br.s(finePhBits);
      for(size_t ci=0; ci<ncells; ++ci) locks[ci]=(int)br.u(lockBits);
      for(size_t ci=0; ci<ncells; ++ci) halves[ci]=(int)br.u(1);
      for(size_t ci=0; ci<ncells; ++ci) cnrs[ci]=(int)br.u(cnrBits);
      if(var==5||var==7) for(size_t ci=0; ci<ncells; ++ci) fineRates[ci]=br.s(15);
      for(size_t ci=0; ci<ncells; ++ci){
        auto [si,gi]=cells[ci];
        auto finePr=finePrs[ci]; auto finePh=finePhs[ci];
        int lock=locks[ci]; int half=halves[ci]; int cnr=cnrs[ci]; int64_t fineRate=fineRates[ci];
        std::string sid=sat(sys,sats[si]); auto& r=recMap[obsKey(tp,sid)]; r.time=tp; r.satellite=sid; r.system=sys; std::string code=msmCode(sys,sigs[gi]);
        bool validRange = rough[si]!=255 && mod[si]!=1023;
        bool validFinePr = finePr != -(1LL<<(finePrBits-1));
        bool validFinePh = finePh != -(1LL<<(finePhBits-1));
        if(validRange && validFinePr){
          double ms=rough[si]+mod[si]/1024.0; double c_ms=299792.458; double pr=(ms+finePr/finePrDen)*c_ms;
          addVal(r,"C"+code,pr); addType(types,sys,"C"+code);
          double freq=frequency(sys,code);
          if(freq>0 && validFinePh){
            double lambda=299792458.0/freq; double ph_ms=ms+finePh/finePhDen;
            std::string lli=rtcmLliFor(sid+":"+code,lock,half,msmLockStates);
            addVal(r,"L"+code,(ph_ms*c_ms)/lambda,lli); addType(types,sys,"L"+code);
            if(var==5||var==7){ addVal(r,"D"+code,-(rate[si]+fineRate*0.0001)/lambda); addType(types,sys,"D"+code); }
          }
        }
        double snr=(var==6||var==7)?cnr/16.0:cnr; addVal(r,"S"+code,snr); addType(types,sys,"S"+code);
      }
    }
    else if(type==1019||type==1020||type==1042||type==1044||type==1045||type==1046){
      NavigationRecord nr;
      if(type==1019) nr=parseGPS1019(br,anchor);
      else if(type==1020) nr=parseGLO1020(br,anchor);
      else if(type==1042) nr=parseBDS1042(br,anchor);
      else if(type==1044) nr=parseQZSS1044(br,anchor);
      else if(type==1045||type==1046) nr=parseGAL104x(br,anchor,type);
      else { nr.system=type==1042?"C":(type==1044?"J":(type>=1045?"E":"G")); nr.satellite=sat(nr.system,(int)br.u(type==1044?4:6)); nr.messageType=type==1020?"FDMA":(type>=1045?"INAV":"LNAV"); nr.messageSubtype="RAW"; nr.recordType="EPH"; nr.rawLines.push_back("> EPH "+nr.satellite+" "+nr.messageType+" RAW"); }
      nr.rawLines.push_back("    RTCM3 "+std::to_string(type)+" "+hexWords(pl)); navs.push_back(nr); sum.ephemeris[std::to_string(type)+":"+nr.satellite]={type,nr.system,nr.satellite,++sum.ephemeris[std::to_string(type)+":"+nr.satellite].count}; }
    i+=5+len;
  }
  normalizeGlonassFdmaObservationCycles(recMap, navs);
  RinexFile obs; obs.path=path; obs.header.kind=RinexKind::Obs; obs.header.version=3.05; obs.header.obsTypes=types; obs.data.observationRecords.reserve(recMap.size()); for(auto&kv:recMap) obs.data.observationRecords.push_back(kv.second); if(haveStation){ auto mk=[](std::string v,std::string l){ if(v.size()>60)v=v.substr(0,60); std::ostringstream os; os<<std::left<<std::setw(60)<<v<<std::setw(20)<<l; std::string raw=os.str(); return HeaderLine{raw,raw.substr(0,60),l};}; std::ostringstream pos; pos<<std::fixed<<std::setprecision(4)<<std::setw(14)<<station.x<<std::setw(14)<<station.y<<std::setw(14)<<station.z; obs.header.lines.push_back(mk(pos.str(),"APPROX POSITION XYZ")); std::ostringstream recLine; recLine<<std::left<<std::setw(20)<<(receiverSerial.empty()?"UNKNOWN":receiverSerial.substr(0,20))<<std::setw(20)<<receiver.substr(0,20)<<std::setw(20)<<receiverVer.substr(0,20); obs.header.lines.push_back(mk(recLine.str(),"REC # / TYPE / VERS")); std::ostringstream antLine; antLine<<std::left<<std::setw(20)<<(antSerial.empty()?"UNKNOWN":antSerial.substr(0,20))<<std::setw(20)<<antenna.substr(0,20); obs.header.lines.push_back(mk(antLine.str(),"ANT # / TYPE")); }
  obs=ceqc::service::rinex::merge({obs},RinexKind::Obs,3.05); sum.epochCount=(int)epochs.size(); sum.observationCount=(int)obs.data.observationRecords.size(); sum.observationValueCount=0; for(auto&r:obs.data.observationRecords) sum.observationValueCount+=(int)r.values.size(); obs.rtcm3=sum;
  std::vector<RinexFile> out{obs}; if(!navs.empty()){ RinexFile nav; nav.path=path; nav.header.kind=RinexKind::Nav; nav.header.version=4.02; nav.data.navigationRecords=navs; for(auto& n:navs) nav.body.insert(nav.body.end(),n.rawLines.begin(),n.rawLines.end()); nav=ceqc::service::rinex::merge({nav},RinexKind::Nav,4.02); out.push_back(nav); } return out;
}
};
}
std::shared_ptr<Translator> makeRTCM3(){ return std::make_shared<RTCM3Translator>(); }
}
