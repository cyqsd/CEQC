#include "ceqc/rinex/RinexService.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace ceqc::service::rinex {
namespace {
ParserOptions g_parserOptions{};
std::string trim(const std::string& s){ size_t a=0,b=s.size(); while(a<b && std::isspace((unsigned char)s[a])) ++a; while(b>a && std::isspace((unsigned char)s[b-1])) --b; return s.substr(a,b-a); }
std::string upper(std::string s){ for(auto& c:s)c=std::toupper((unsigned char)c); return s; }
std::vector<std::string> fields(const std::string& s){ std::istringstream is(s); std::vector<std::string> v; std::string x; while(is>>x)v.push_back(x); return v; }
std::string substrSafe(const std::string& s,size_t a,size_t b){ if(a>=s.size()) return {}; if(b>s.size()) b=s.size(); if(b<a) b=a; return s.substr(a,b-a); }
double rinexFloat(std::string s){ for(auto& c:s) if(c=='D'||c=='d') c='E'; auto t=trim(s); if(t=="NaN"||t=="NAN"||t=="nan") return std::numeric_limits<double>::quiet_NaN(); return std::stod(t); }
std::string guessLabel(const std::string& line){
  static const std::vector<std::string> labels={"RINEX VERSION / TYPE","PGM / RUN BY / DATE","MARKER NAME","MARKER NUMBER","MARKER TYPE","OBSERVER / AGENCY","REC # / TYPE / VERS","ANT # / TYPE","APPROX POSITION XYZ","ANTENNA: DELTA H/E/N","WAVELENGTH FACT L1/2","# / TYPES OF OBSERV","SYS / # / OBS TYPES","TIME OF FIRST OBS","TIME OF LAST OBS","INTERVAL","LEAP SECONDS","IONOSPHERIC CORR","TIME SYSTEM CORR","COMMENT","END OF HEADER"};
  auto up=upper(line);
  for(auto& l:labels){ auto p=up.find(l); if(p!=std::string::npos) return l; }
  return "";
}
HeaderLine parseHeaderLine(const std::string& line){
  HeaderLine h; h.raw=line;
  if(line.size()>=60){ h.value=line.substr(0,60); h.label=upper(trim(line.substr(60))); }
  else { h.value=trim(line); h.label=g_parserOptions.relax?guessLabel(line):""; }
  auto guessed = guessLabel(line);
  // Many real-world RINEX files and hand-built tests are slightly short, which
  // can split the 20-character label field (e.g. "SYS / # / OBS TYPES" becomes
  // "YS / # / OBS TYPES" at column 61).  Recognize the known label in the whole
  // line and cut the value before the label while still preserving strict parsing
  // for unknown labels.
  if(!guessed.empty() && (h.label.empty() || h.label.find(guessed)==std::string::npos)) {
    h.label = guessed;
    auto p = upper(line).find(guessed);
    if(p != std::string::npos) h.value = line.substr(0, p);
  }
  if(h.label.empty() && g_parserOptions.relax) h.label=guessed;
  if(!h.label.empty()){
    auto p=upper(h.value).find(h.label);
    if(p!=std::string::npos) h.value=trim(h.value.substr(0,p));
  }
  return h;
}
std::string normalizeSat(const std::string& id, const std::string& defaultSys="G"){
  auto t=upper(trim(id));
  if(t.empty()) return {};
  auto pad2=[](std::string n){ if(n.size()==1) return std::string("0")+n; return n; };
  if(std::isdigit((unsigned char)t[0])) return (defaultSys.empty()?"G":defaultSys.substr(0,1))+pad2(t);
  if(t.size()==2 && std::isalpha((unsigned char)t[0]) && std::isdigit((unsigned char)t[1])) return t.substr(0,1)+"0"+t.substr(1,1);
  return t;
}
std::string fixedSatID(const std::string& line){
  if(line.size()<3) return {};
  auto id=upper(trim(line.substr(0,3)));
  if(id.size()<2) return {};
  if(std::string("GRECSIJL").find(id[0])!=std::string::npos) return normalizeSat(id, std::string(1,id[0]));
  return {};
}
TimePoint parseEpoch(int y,int m,int d,int hh,int mm,double ss){ return makeUTC(y,m,d,hh,mm,ss); }
int fullYear(int yy){ if(yy<80) return 2000+yy; if(yy<100) return 1900+yy; return yy; }
std::optional<double> rinexFloatOpt(const std::string& s){ auto t=trim(s); if(t.empty()) return {}; try{ return rinexFloat(t); }catch(...){ return {}; } }
std::vector<double> fixedNavFloats(const std::string& line,size_t start){
  std::vector<double> vals;
  if(start>=line.size()) return vals;
  for(size_t pos=start; pos<line.size(); pos+=19){
    auto chunk=substrSafe(line,pos,pos+19);
    if(trim(chunk).empty()) continue;
    auto v=rinexFloatOpt(chunk);
    if(v) vals.push_back(*v);
  }
  return vals;
}
std::optional<TimePoint> navTimeFixed(const std::string& line,double ver){
  try{
    if(ver>=3.0){
      int y=std::stoi(trim(substrSafe(line,4,8)));
      int mo=std::stoi(trim(substrSafe(line,9,11)));
      int d=std::stoi(trim(substrSafe(line,12,14)));
      int hh=std::stoi(trim(substrSafe(line,15,17)));
      int mi=std::stoi(trim(substrSafe(line,18,20)));
      double ss=std::stod(trim(substrSafe(line,21,23)));
      return parseEpoch(y,mo,d,hh,mi,ss);
    }
    int y=fullYear(std::stoi(trim(substrSafe(line,3,5))));
    int mo=std::stoi(trim(substrSafe(line,6,8)));
    int d=std::stoi(trim(substrSafe(line,9,11)));
    int hh=std::stoi(trim(substrSafe(line,12,14)));
    int mi=std::stoi(trim(substrSafe(line,15,17)));
    double ss=std::stod(trim(substrSafe(line,18,22)));
    return parseEpoch(y,mo,d,hh,mi,ss);
  }catch(...){ return {}; }
}
void buildHeaderIndex(RinexHeader& h){ h.byLabel.clear(); for(size_t i=0;i<h.lines.size();++i) if(!h.lines[i].label.empty()) h.byLabel[h.lines[i].label].push_back(i); }
void parseObsTypes(RinexHeader& h){
  h.obsTypes.clear();
  if(h.kind!=RinexKind::Obs) return;
  if(h.version>=3.0){
    std::string lastSys;
    std::map<std::string,int> expected;
    for(auto idx: h.byLabel["SYS / # / OBS TYPES"]){
      auto rawValue = h.lines[idx].value;
      auto f=fields(rawValue);
      if(f.empty()) continue;
      std::string sys;
      size_t start=0;
      if(f[0].size()==1 && std::string("GRECSIJL").find((char)std::toupper((unsigned char)f[0][0]))!=std::string::npos){
        sys=upper(f[0]);
        start=1;
        if(f.size()>1){ try{ expected[sys]=std::stoi(f[1]); start=2; }catch(...){ } }
        lastSys=sys;
      } else if(!lastSys.empty()) {
        sys=lastSys;
        start=0;
        // Fixed-column continuation lines can have the system in column 1 blank
        // and observation codes starting around column 7.
      } else continue;
      for(size_t i=start;i<f.size();++i) if(f[i].size()==3) h.obsTypes[sys].push_back(upper(f[i]));
    }
    for(auto& kv: expected){
      auto& v=h.obsTypes[kv.first];
      if(kv.second>0 && (int)v.size()>kv.second) v.resize(kv.second);
    }
  } else {
    std::vector<std::string> types;
    int expected=0;
    for(auto idx: h.byLabel["# / TYPES OF OBSERV"]){
      auto f=fields(h.lines[idx].value);
      size_t start=0;
      if(!f.empty()){
        try{ expected=std::stoi(f[0]); start=1; }catch(...){ }
      }
      for(size_t i=start;i<f.size();++i) {
        auto code=upper(f[i]);
        if(!code.empty()) types.push_back(code);
      }
    }
    if(expected>0 && (int)types.size()>expected) types.resize(expected);
    if(!types.empty()){
      std::string sys=h.satelliteSystem.empty()?"G":upper(h.satelliteSystem.substr(0,1));
      if(sys==" " || sys=="M") sys="G";
      h.obsTypes[sys]=types;
      h.obsTypes["G"]=types;
    }
  }
}
void parseMetTypes(RinexHeader& h){ h.metTypes.clear(); if(h.kind!=RinexKind::Met) return; for(auto idx: h.byLabel["# / TYPES OF OBSERV"]){ auto f=fields(h.lines[idx].value); size_t start=0; if(!f.empty()){ try{ std::stoi(f[0]); start=1; }catch(...){ } } for(size_t i=start;i<f.size();++i) h.metTypes.push_back(upper(f[i])); } }
RinexHeader buildHeader(std::vector<HeaderLine> lines){
  RinexHeader h; h.lines=std::move(lines); buildHeaderIndex(h);
  if(!h.lines.empty()){
    auto f=fields(h.lines[0].raw);
    if(!f.empty()){ try{ h.version=std::stod(f[0]); }catch(...){ } }
    auto first=h.lines[0].raw; auto up=upper(first);
    if(first.size()>20){
      std::string t=trim(first.substr(20,1));
      if(t=="O" || (t=="M" && up.find("OBS")!=std::string::npos)) h.kind=RinexKind::Obs;
      else if(t=="N"||t=="G"||t=="R"||t=="E"||t=="C"||t=="J"||t=="I"||t=="S"||t=="L") h.kind=RinexKind::Nav;
      else if(t=="M" && up.find("MET")!=std::string::npos) h.kind=RinexKind::Met;
    }
    if(up.find("OBSERVATION")!=std::string::npos) h.kind=RinexKind::Obs;
    else if(up.find("NAVIGATION")!=std::string::npos || up.find("NAV DATA")!=std::string::npos) h.kind=RinexKind::Nav;
    else if(up.find("METEOROLOGICAL")!=std::string::npos) h.kind=RinexKind::Met;
    if(first.size()>40) h.satelliteSystem=trim(first.substr(40,1));
    if(up.find("GLONASS")!=std::string::npos) h.satelliteSystem="R";
    else if(up.find("GALILEO")!=std::string::npos) h.satelliteSystem="E";
    else if(up.find("BEIDOU")!=std::string::npos || up.find("BDS")!=std::string::npos) h.satelliteSystem="C";
    else if(up.find("QZSS")!=std::string::npos) h.satelliteSystem="J";
    else if(up.find("SBAS")!=std::string::npos || up.find("GEOSTATIONARY")!=std::string::npos) h.satelliteSystem="S";
    else if(up.find("IRNSS")!=std::string::npos || up.find("NAVIC")!=std::string::npos) h.satelliteSystem="I";
    else if(h.satelliteSystem.empty() && h.kind==RinexKind::Nav && h.version<3.0) h.satelliteSystem="G";
  }
  parseObsTypes(h); parseMetTypes(h); return h;
}
std::vector<ObservationValue> parseValues(const std::string& text,const std::vector<std::string>& types){ std::vector<ObservationValue> vals; for(size_t i=0;i<types.size();++i){ std::string chunk=substrSafe(text,i*16,i*16+16); ObservationValue v; v.type=types[i]; v.raw=chunk; auto vs=trim(substrSafe(chunk,0,14)); if(!vs.empty()){ try{ v.value=rinexFloat(vs); }catch(...){ } } v.lli=trim(substrSafe(chunk,14,15)); v.ssi=trim(substrSafe(chunk,15,16)); vals.push_back(v); } return vals; }
std::vector<std::string> v2SatList(const std::vector<std::string>& lines,size_t i,int nsat,const std::string& defaultSys){
  std::string joined;
  for(int k=0;k<1+(nsat-1)/12 && i+(size_t)k<lines.size();++k) if(lines[i+k].size()>32) joined+=lines[i+k].substr(32);
  std::vector<std::string> sats;
  for(size_t off=0;off+3<=joined.size() && (int)sats.size()<nsat;off+=3){ auto id=trim(joined.substr(off,3)); if(!id.empty()) sats.push_back(normalizeSat(id, defaultSys)); }
  return sats;
}
std::vector<std::string> firstObsTypes(const RinexHeader& h,const std::string& sys="G"){ auto it=h.obsTypes.find(sys); if(it!=h.obsTypes.end()) return it->second; auto g=h.obsTypes.find("G"); if(g!=h.obsTypes.end()) return g->second; if(!h.obsTypes.empty()) return h.obsTypes.begin()->second; return {}; }
void parseObs(RinexFile& rf){ if(rf.header.version>=3.0){ for(size_t i=0;i<rf.body.size();++i){ auto line=rf.body[i]; if(line.empty()||line[0]!='>') continue; auto f=fields(line.substr(1)); if(f.size()<7) continue; int y=std::stoi(f[0]),mo=std::stoi(f[1]),d=std::stoi(f[2]),hh=std::stoi(f[3]),mi=std::stoi(f[4]); double sec=std::stod(f[5]); int flag=std::stoi(f[6]); int nsat=f.size()>7?std::stoi(f[7]):0; ObservationEpoch ep; ep.time=parseEpoch(y,mo,d,hh,mi,sec); ep.flag=flag; ep.lineIndex=i; ep.rawLine=line; size_t epochIndex=rf.data.observationEpochs.size(); size_t j=i+1; while(j<rf.body.size() && (int)ep.satellites.size()<nsat){ if(!rf.body[j].empty()&&rf.body[j][0]=='>') break; auto sat=fixedSatID(rf.body[j]); if(sat.empty()){ ++j; continue; } auto sys=sat.substr(0,1); auto types=firstObsTypes(rf.header,sys); std::string obsText=substrSafe(rf.body[j],3,rf.body[j].size()); std::vector<std::string> raw{rf.body[j]}; ++j; while(!types.empty() && obsText.size() < types.size()*16 && j<rf.body.size() && (rf.body[j].empty()||rf.body[j][0]!='>') && fixedSatID(rf.body[j]).empty()){ auto cont = rf.body[j]; obsText += cont.size() > 3 ? cont.substr(3) : std::string{}; raw.push_back(rf.body[j]); ++j; } ObservationRecord rec; rec.epochIndex=epochIndex; rec.time=ep.time; rec.flag=flag; rec.satellite=sat; rec.system=sys; rec.rawLines=raw; rec.values=parseValues(obsText,types); rf.data.observationRecords.push_back(rec); ep.satellites.push_back(sat); }
      rf.data.observationEpochs.push_back(ep); i=j-1; } } else { auto types=firstObsTypes(rf.header); for(size_t i=0;i<rf.body.size();++i){ auto f=fields(rf.body[i]); if(f.size()<8) continue; int yy; try{ yy=std::stoi(f[0]); }catch(...){ continue; } if(yy<0||yy>99) continue; int y=fullYear(yy),mo=std::stoi(f[1]),d=std::stoi(f[2]),hh=std::stoi(f[3]),mi=std::stoi(f[4]); double sec=std::stod(f[5]); int flag=std::stoi(f[6]), nsat=std::stoi(f[7]); ObservationEpoch ep; ep.time=parseEpoch(y,mo,d,hh,mi,sec); ep.flag=flag; ep.lineIndex=i; ep.rawLine=rf.body[i]; ep.satellites=v2SatList(rf.body,i,nsat, (rf.header.satelliteSystem.empty()?"G":rf.header.satelliteSystem.substr(0,1))); size_t epochIndex=rf.data.observationEpochs.size(); rf.data.observationEpochs.push_back(ep); size_t j=i+1+(nsat>12?(nsat-1)/12:0); int linesPerSat=(int)(types.size()+4)/5; for(int s=0;s<nsat && j<rf.body.size();++s){ std::string obsText; std::vector<std::string> raw; for(int k=0;k<linesPerSat && j<rf.body.size();++k,++j){ raw.push_back(rf.body[j]); obsText+=rf.body[j]; } std::string sat=s<(int)ep.satellites.size()?ep.satellites[s]:"G"+std::to_string(s+1); ObservationRecord rec; rec.epochIndex=epochIndex; rec.time=ep.time; rec.flag=flag; rec.satellite=sat; rec.system=sat.substr(0,1); rec.rawLines=raw; rec.values=parseValues(obsText,types); rf.data.observationRecords.push_back(rec); } i=j-1; } } }
std::optional<TimePoint> navTime(const std::vector<std::string>& f,size_t off,bool longYear){ if(f.size()<off+5) return {}; try{ int y=std::stoi(f[off]),mo=std::stoi(f[off+1]),d=std::stoi(f[off+2]),hh=std::stoi(f[off+3]),mi=std::stoi(f[off+4]); double ss=f.size()>off+5?std::stod(f[off+5]):0; if(!longYear) y=fullYear(y); return parseEpoch(y,mo,d,hh,mi,ss); }catch(...){ return {}; } }
bool isNavStart(const std::string& line,double ver){ auto t=trim(line); if(t.rfind(">",0)==0) return true; auto f=fields(line); if(ver>=3.0) return f.size()>=7 && f[0].size()>=2 && std::string("GRECIJSL").find(f[0][0])!=std::string::npos; return f.size()>=6 && std::isdigit((unsigned char)f[0][0]); }
std::vector<std::string> navNamesFor(const NavigationRecord& r){
  static const std::vector<std::string> gps={"SV_clock_bias","SV_clock_drift","SV_clock_drift_rate","IODE","Crs","DeltaN","M0","Cuc","Eccentricity","Cus","SqrtA","Toe","Cic","Omega0","Cis","I0","Crc","Omega","OmegaDot","IDOT","CodesL2","Week","L2PFlag","SVAccuracy","SVHealth","TGD","IODC","TransmissionTime","FitInterval"};
  static const std::vector<std::string> gpsCnav={"SV_clock_bias","SV_clock_drift","SV_clock_drift_rate","ADot","Crs","DeltaN","M0","Cuc","Eccentricity","Cus","SqrtA","Top","Cic","Omega0","Cis","I0","Crc","Omega","OmegaDot","IDOT","DeltaNDot","URAI_NED0","URAI_NED1","URAI_ED","SVHealth","TGD","URAI_NED2","ISC_L1CA","ISC_L2C","ISC_L5I5","ISC_L5Q5","TransmissionTime","Week","Flags"};
  static const std::vector<std::string> glo={"TauN","GammaN","MessageFrameTime","X","VX","AX","Health","Y","VY","AY","FrequencyNumber","Z","VZ","AZ","AgeOfOperation"};
  static const std::vector<std::string> gloCdma={"SV_clock_bias","SV_clock_drift","SV_clock_drift_rate","Health","FrequencyNumber","Tk","Tb","X","VX","AX","Y","VY","AY","Z","VZ","AZ","DeltaTau","En","P1","P2","P3","P4","NT","M"};
  static const std::vector<std::string> gal={"SV_clock_bias","SV_clock_drift","SV_clock_drift_rate","IODnav","Crs","DeltaN","M0","Cuc","Eccentricity","Cus","SqrtA","Toe","Cic","Omega0","Cis","I0","Crc","Omega","OmegaDot","IDOT","DataSources","Week","Spare","SISA","SVHealth","BGD_E5a_E1","BGD_E5b_E1","TransmissionTime"};
  static const std::vector<std::string> bds={"SV_clock_bias","SV_clock_drift","SV_clock_drift_rate","AODE","Crs","DeltaN","M0","Cuc","Eccentricity","Cus","SqrtA","Toe","Cic","Omega0","Cis","I0","Crc","Omega","OmegaDot","IDOT","AODC","Week","Spare1","SVAccuracy","SatH1","TGD1","TGD2","TransmissionTime","AODC2"};
  static const std::vector<std::string> bdsCnav={"SV_clock_bias","SV_clock_drift","SV_clock_drift_rate","ADot","Crs","DeltaN","M0","Cuc","Eccentricity","Cus","SqrtA","Toe","Cic","Omega0","Cis","I0","Crc","Omega","OmegaDot","IDOT","DeltaNDot","ISC_B1Cp","ISC_B2ap","ISC_B2bp","TGD_B1Cp","TGD_B2ap","TGD_B2bp","SISAI_oe","SISAI_ocb","SISAI_oc1","SISAI_oc2","SatType","IODC"};
  static const std::vector<std::string> navic={"SV_clock_bias","SV_clock_drift","SV_clock_drift_rate","IODEC","Crs","DeltaN","M0","Cuc","Eccentricity","Cus","SqrtA","Toe","Cic","Omega0","Cis","I0","Crc","Omega","OmegaDot","IDOT","Week","URA","Health","TGD","TransmissionTime"};
  static const std::vector<std::string> sto={"A0","A1","ReferenceTime","ReferenceWeek","Provider","UTC_ID","SISAccuracy"};
  static const std::vector<std::string> eop={"Xp","Xpdot","Yp","Ypdot","UT1_UTC","LOD","ReferenceTime","ReferenceWeek"};
  static const std::vector<std::string> ionKlob={"Alpha0","Alpha1","Alpha2","Alpha3","Beta0","Beta1","Beta2","Beta3","ReferenceTime","ReferenceWeek"};
  static const std::vector<std::string> ionNequick={"Ai0","Ai1","Ai2","Region1","Region2","Region3","Region4","Region5","ReferenceTime","ReferenceWeek"};
  static const std::vector<std::string> ionBds={"Alpha0","Alpha1","Alpha2","Alpha3","Beta0","Beta1","Beta2","Beta3","BDGIM0","BDGIM1","BDGIM2","BDGIM3","BDGIM4","BDGIM5","BDGIM6","BDGIM7","BDGIM8","BDGIM9"};
  auto rt=upper(r.recordType), mt=upper(r.messageType);
  if(rt=="STO") return sto;
  if(rt=="EOP") return eop;
  if(rt=="ION"){
    if(mt.find("NEQ")!=std::string::npos || mt.find("GAL")!=std::string::npos) return ionNequick;
    if(r.system=="C" || mt.find("BDS")!=std::string::npos || mt.find("B1")!=std::string::npos) return ionBds;
    return ionKlob;
  }
  if(r.system=="R") return (mt=="L1OC" || mt=="L3OC" || mt=="CDMA") ? gloCdma : glo;
  if(r.system=="E") return gal;
  if(r.system=="C") return (mt=="CNAV" || mt=="CNV1" || mt=="CNV2" || mt=="CNV3" || mt=="B1C" || mt=="B2A" || mt=="B2B") ? bdsCnav : bds;
  if(r.system=="I") return navic;
  if((r.system=="G" || r.system=="J") && (r.messageType=="CNAV" || r.messageType=="CNAV2" || r.messageType=="CNV2" || r.messageType=="CNAV-2")) return gpsCnav;
  return gps;
}
void fieldModel(NavigationRecord& r){
  auto names=navNamesFor(r);
  r.fields.clear();
  for(size_t i=0;i<r.values.size();++i){
    std::string n=i<names.size()?names[i]:"value_"+std::to_string(i);
    r.fields[n]={n,"",r.values[i],i};
  }
}

void parseNav(RinexFile& rf){
  for(size_t i=0;i<rf.body.size();){
    if(trim(rf.body[i]).empty()){ ++i; continue; }
    if(!isNavStart(rf.body[i],rf.header.version)){ ++i; continue; }
    size_t start=i++;
    bool typed = trim(rf.body[start]).rfind(">",0)==0;
    if(typed){
      while(i<rf.body.size() && trim(rf.body[i]).rfind(">",0)!=0) ++i;
    } else {
      while(i<rf.body.size() && !isNavStart(rf.body[i],rf.header.version)) ++i;
    }
    std::vector<std::string> lines(rf.body.begin()+start,rf.body.begin()+i);
    NavigationRecord rec; rec.rawLines=lines; rec.lineIndex=start; rec.recordType="EPH"; rec.messageSubtype="BNK"; size_t dataLine=0;
    if(!lines.empty() && trim(lines[0]).rfind(">",0)==0){
      auto f=fields(trim(lines[0]).substr(1));
      if(!f.empty()) rec.recordType=upper(f[0]);
      auto looksSat=[](const std::string& x){ return x.size()>=2 && std::string("GRECSIJL").find((char)std::toupper((unsigned char)x[0]))!=std::string::npos && std::isdigit((unsigned char)x[1]); };
      if(f.size()>1){
        std::string tok=upper(f[1]);
        if(looksSat(tok)){
          rec.satellite=normalizeSat(tok,std::string(1,tok[0])); rec.system=rec.satellite.substr(0,1);
          if(f.size()>2) rec.messageType=upper(f[2]);
          if(f.size()>3) rec.messageSubtype=upper(f[3]);
        } else {
          rec.messageType=tok;
          rec.system=tok.empty()?std::string{}:tok.substr(0,1);
          if(f.size()>2) rec.messageSubtype=upper(f[2]);
        }
      }
      dataLine=1;
    }
    if(dataLine<lines.size()){
      auto f=fields(lines[dataLine]);
      if(rf.header.version>=3.0 && rec.recordType=="EPH"){
        if(rec.satellite.empty()){
          auto sat=fixedSatID(lines[dataLine]);
          if(!sat.empty()){ rec.satellite=sat; rec.system=sat.substr(0,1); }
          else if(!f.empty()){ rec.satellite=upper(f[0]); rec.system=rec.satellite.empty()?std::string{}:rec.satellite.substr(0,1); }
        }
        rec.epoch=navTimeFixed(lines[dataLine], rf.header.version);
        auto nums=fixedNavFloats(lines[dataLine], 23);
        rec.values.insert(rec.values.end(), nums.begin(), nums.end());
      } else if(rf.header.version<3.0 && f.size()>=1){
        rec.system=(rf.header.satelliteSystem.empty()?"G":rf.header.satelliteSystem.substr(0,1));
        rec.satellite=normalizeSat(substrSafe(lines[dataLine],0,2), rec.system);
        rec.epoch=navTimeFixed(lines[dataLine], rf.header.version);
        auto nums=fixedNavFloats(lines[dataLine], 22);
        rec.values.insert(rec.values.end(), nums.begin(), nums.end());
      } else {
        // RINEX 4 non-EPH system data records (STO/EOP/ION/etc.) are numeric
        // payload records without a satellite epoch line.  Preserve them as
        // NavigationRecord values for QC inventory and round-trip output.
        for(auto& x: f) try{ rec.values.push_back(rinexFloat(x)); }catch(...){ }
      }
    }
    for(size_t l=dataLine+1;l<lines.size();++l){
      if(rec.recordType=="EPH"){
        size_t startCol = rf.header.version<3.0 ? 3 : 4;
        auto nums=fixedNavFloats(lines[l], startCol);
        rec.values.insert(rec.values.end(), nums.begin(), nums.end());
      } else {
        for(auto& x: fields(lines[l])) try{ rec.values.push_back(rinexFloat(x)); }catch(...){ }
      }
    }
    if(rec.messageType.empty()) rec.messageType=(rec.system=="R"?"FDMA":"LNAV");
    fieldModel(rec);
    rf.data.navigationRecords.push_back(rec);
  }
}
void parseMet(RinexFile& rf){ for(size_t i=0;i<rf.body.size();++i){ auto f=fields(rf.body[i]); if(f.size()<6) continue; try{ int y=fullYear(std::stoi(f[0])),mo=std::stoi(f[1]),d=std::stoi(f[2]),hh=std::stoi(f[3]),mi=std::stoi(f[4]); double ss=std::stod(f[5]); MeteorologicalRecord r; r.time=parseEpoch(y,mo,d,hh,mi,ss); r.rawLine=rf.body[i]; r.lineIndex=i; for(size_t k=6;k<f.size() && k-6<rf.header.metTypes.size();++k) r.values[rf.header.metTypes[k-6]]=rinexFloat(f[k]); rf.data.meteorologicalRecords.push_back(r); }catch(...){ } } }
}
void setParserOptions(ParserOptions options){ g_parserOptions = options; }
ParserOptions parserOptions(){ return g_parserOptions; }
RinexFile readFile(const std::string& path){ std::ifstream f(path); if(!f) throw std::runtime_error("cannot open "+path); return readStream(path,f); }
RinexFile readStream(const std::string& path,std::istream& is){ RinexFile rf; rf.path=path; std::string line; bool inHeader=true; std::vector<HeaderLine> h; while(std::getline(is,line)){ if(!line.empty()&&line.back()=='\r') line.pop_back(); if(inHeader){ auto hl=parseHeaderLine(line); h.push_back(hl); if(hl.label=="END OF HEADER") inHeader=false; } else rf.body.push_back(line); } rf.header=buildHeader(h); parseContent(rf); return rf; }
void parseContent(RinexFile& rf){ rf.data={}; if(rf.header.kind==RinexKind::Obs) parseObs(rf); else if(rf.header.kind==RinexKind::Nav) parseNav(rf); else if(rf.header.kind==RinexKind::Met) parseMet(rf); }
std::vector<ValidationIssue> validate(const RinexFile& rf){ std::vector<ValidationIssue> issues; auto miss=[&](const std::string& l){ if(!rf.header.byLabel.count(l)) issues.push_back({"ERROR","missing required RINEX header label \""+l+"\"",l});}; miss("RINEX VERSION / TYPE"); if(rf.header.kind==RinexKind::Obs){ miss("END OF HEADER"); if(rf.header.version>=3.0) miss("SYS / # / OBS TYPES"); else miss("# / TYPES OF OBSERV"); } else miss("END OF HEADER"); if(rf.header.version==0) issues.push_back({"ERROR","cannot parse RINEX version","RINEX VERSION / TYPE"}); return issues; }
}
