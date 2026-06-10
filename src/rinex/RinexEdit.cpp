#include "ceqc/rinex/RinexService.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <sstream>
#include <set>
#include <vector>
namespace ceqc::service::rinex {
namespace {
std::string trim(const std::string& s){size_t a=0,b=s.size();while(a<b&&isspace((unsigned char)s[a]))++a;while(b>a&&isspace((unsigned char)s[b-1]))--b;return s.substr(a,b-a);} 
std::string upper(std::string s){for(auto& c:s)c=toupper((unsigned char)c);return s;} 
HeaderLine makeLine(std::string value,std::string label){ if(value.size()>60)value=value.substr(0,60); std::ostringstream os; os<<std::left<<std::setw(60)<<value<<std::setw(20)<<upper(label); std::string raw=os.str(); return {raw,raw.substr(0,60),upper(label)}; } 
void rebuild(RinexHeader& h){ h.byLabel.clear(); for(size_t i=0;i<h.lines.size();++i) h.byLabel[h.lines[i].label].push_back(i); } 
std::vector<std::string> splitCodes(const std::string& s){ std::vector<std::string> out; std::string token; for(char c:s){ if(std::isspace((unsigned char)c)||c==','||c=='+'){ if(!token.empty()){out.push_back(upper(token));token.clear();} } else token.push_back(c); } if(!token.empty()) out.push_back(upper(token)); return out;}
std::vector<std::string> splitWords(const std::string& s){ std::istringstream is(s); std::vector<std::string> v; std::string x; while(is>>x)v.push_back(x); return v; }
std::string keyOf(const std::string& opt){ auto x=opt; auto pos=x.find('.'); if(pos!=std::string::npos) x=x.substr(pos+1); return upper(x); }
bool isSpecial(const std::string& k){ static const std::set<std::string> s={"OBS","-OBS","RENAME_OBS","DECIMATE","SUMMARY"}; return s.count(k)>0; }
std::string existingValue(const RinexFile& rf,const std::string& label){ auto it=rf.header.byLabel.find(label); if(it==rf.header.byLabel.end()||it->second.empty()) return {}; return rf.header.lines[it->second.front()].value; }
std::string col(const std::string& s,size_t a,size_t n){ if(a>=s.size()) return {}; return trim(s.substr(a,std::min(n,s.size()-a))); }
std::string fixed3(const std::string& value){ auto f=splitWords(value); std::ostringstream os; for(int i=0;i<3;++i){ double v=0; if(i<(int)f.size()) { try{v=std::stod(f[i]);}catch(...){}} os<<std::setw(14)<<std::fixed<<std::setprecision(4)<<v; } return os.str(); }

std::string fixedXYZ(double x, double y, double z){ std::ostringstream os; os<<std::setw(14)<<std::fixed<<std::setprecision(4)<<x<<std::setw(14)<<y<<std::setw(14)<<z; return os.str(); }
std::string geodeticToXYZ(const std::string& value){
  auto f=splitWords(value); if(f.size()<3) return fixed3(value);
  double lat=0,lon=0,h=0; try{ lat=std::stod(f[0]); lon=std::stod(f[1]); h=std::stod(f[2]); }catch(...){ return fixed3(value); }
  constexpr double a=6378137.0; constexpr double invf=298.257223563; double ff=1.0/invf; double e2=ff*(2-ff);
  constexpr double PI=3.141592653589793238462643383279502884; double phi=lat*PI/180.0, lam=lon*PI/180.0; double sp=std::sin(phi), cp=std::cos(phi); double N=a/std::sqrt(1-e2*sp*sp);
  return fixedXYZ((N+h)*cp*std::cos(lam),(N+h)*cp*std::sin(lam),(N*(1-e2)+h)*sp);
}
std::string slantToHEN(const std::string& value){
  // teqc-style -O.sl[ant] takes: dh slant_distance diameter.
  // Convert to antenna H/E/N with E=N=0 and H = dh + sqrt(s^2 - (d/2)^2).
  auto f=splitWords(value); if(f.size()<3) return fixed3(value); double dh=0,s=0,d=0; try{dh=std::stod(f[0]); s=std::stod(f[1]); d=std::stod(f[2]);}catch(...){return fixed3(value);} double h=dh+std::sqrt(std::max(0.0,s*s-(d/2.0)*(d/2.0))); return fixedXYZ(h,0,0);
}
std::string formatRinexObsTime(const std::string& value){ auto f=splitWords(value); if(f.size()<6) return value; int y=0,mo=0,da=0,hh=0,mi=0; double sec=0; try{y=std::stoi(f[0]);mo=std::stoi(f[1]);da=std::stoi(f[2]);hh=std::stoi(f[3]);mi=std::stoi(f[4]);sec=std::stod(f[5]);}catch(...){return value;} std::ostringstream os; os<<std::setw(6)<<y<<std::setw(6)<<mo<<std::setw(6)<<da<<std::setw(6)<<hh<<std::setw(6)<<mi<<std::setw(13)<<std::fixed<<std::setprecision(7)<<sec<<"     GPS"; return os.str(); }
std::string formatWF(const std::string& value){ auto f=splitWords(value); int i=1,j=1; try{ if(f.size()>0)i=std::stoi(f[0]); if(f.size()>1)j=std::stoi(f[1]); }catch(...){ } std::ostringstream os; os<<std::setw(6)<<i<<std::setw(6)<<j; return os.str(); }
std::string formatModWF(const std::string& value){ auto f=splitWords(value); std::ostringstream os; int i=1,j=1,n=0; try{ if(f.size()>0)i=std::stoi(f[0]); if(f.size()>1)j=std::stoi(f[1]); if(f.size()>2)n=std::stoi(f[2]); }catch(...){ } os<<std::setw(6)<<i<<std::setw(6)<<j<<std::setw(6)<<n; for(size_t k=3;k<f.size();++k) os<<std::setw(6)<<f[k]; return os.str(); }
std::string formatIonCorr(const std::string& value, bool alpha){ auto f=splitWords(value); if(f.empty()) return value; std::string sys=upper(f[0]); std::string typ; if(sys=="G") typ=alpha?"GPSA":"GPSB"; else if(sys=="C") typ=alpha?"BDSA":"BDSB"; else if(sys=="J") typ=alpha?"QZSA":"QZSB"; else if(sys=="I") typ=alpha?"IRNA":"IRNB"; else if(sys=="E") typ="GAL"; else typ=sys+(alpha?"A":"B"); std::ostringstream os; os<<std::left<<std::setw(5)<<typ; for(size_t k=1;k<f.size() && k<5;++k){ double v=0; try{v=std::stod(f[k]);}catch(...){ } os<<std::right<<std::setw(12)<<std::scientific<<std::uppercase<<std::setprecision(4)<<v; } return os.str(); }
std::string formatTimeSystemCorr(const std::string& value){ auto f=splitWords(value); if(f.empty()) return value; std::string sys=upper(f[0]); std::string typ=(sys=="G"?"GPUT":sys=="R"?"GLUT":sys=="E"?"GAUT":sys=="C"?"BDUT":sys=="J"?"QZUT":sys=="I"?"IRUT":sys+"UT"); std::ostringstream os; os<<std::left<<std::setw(5)<<typ; for(size_t k=1;k<f.size() && k<7;++k){ double v=0; try{v=std::stod(f[k]);}catch(...){ } if(k<=2) os<<std::right<<std::setw(17)<<std::scientific<<std::uppercase<<std::setprecision(9)<<v; else os<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(0)<<v; } return os.str(); }
std::string fixedPGM(const std::string& runBy){ std::ostringstream os; os<<std::left<<std::setw(20)<<"ceqc"<<std::setw(20)<<runBy.substr(0,20)<<std::setw(20)<<""; return os.str(); }
void setLine(RinexFile& rf,const std::string& label,const std::string& value,bool append=false){ auto line=makeLine(value,label); size_t end=rf.header.lines.size(); auto eh=rf.header.byLabel.find("END OF HEADER"); if(eh!=rf.header.byLabel.end()&&!eh->second.empty()) end=eh->second[0]; auto it=rf.header.byLabel.find(label); if(append || it==rf.header.byLabel.end()||it->second.empty()){ rf.header.lines.insert(rf.header.lines.begin()+static_cast<long>(end),line); } else { for(auto idx:it->second) rf.header.lines[idx]=line; } rebuild(rf.header); }
}
std::string normalizeHeaderLabel(const std::string& opt){ auto x=keyOf(opt); static std::map<std::string,std::string> m{{"MO","MARKER NAME"},{"MONUMENT","MARKER NAME"},{"MN","MARKER NUMBER"},{"MT","MARKER TYPE"},{"AG","OBSERVER / AGENCY"},{"AGENCY","OBSERVER / AGENCY"},{"O","OBSERVER / AGENCY"},{"OI","OBSERVER / AGENCY"},{"OPERATOR","OBSERVER / AGENCY"},{"RN","REC # / TYPE / VERS"},{"RT","REC # / TYPE / VERS"},{"RV","REC # / TYPE / VERS"},{"AN","ANT # / TYPE"},{"AT","ANT # / TYPE"},{"PX","APPROX POSITION XYZ"},{"PG","APPROX POSITION XYZ"},{"PE","ANTENNA: DELTA H/E/N"},{"DH","ANTENNA: DELTA H/E/N"},{"SLANT","ANTENNA: DELTA H/E/N"},{"START","TIME OF FIRST OBS"},{"ST","TIME OF FIRST OBS"},{"DEF_WF","WAVELENGTH FACT L1/2"},{"MOD_WF","WAVELENGTH FACT L1/2"},{"INT","INTERVAL"},{"INTERVAL","INTERVAL"},{"LEAP","LEAP SECONDS"},{"RUN_BY","PGM / RUN BY / DATE"},{"R","PGM / RUN BY / DATE"},{"SYSTEM","RINEX VERSION / TYPE"},{"COMMENT","COMMENT"},{"C","COMMENT"},{"OBS","SYS / # / OBS TYPES"},{"OBS_TYPES","SYS / # / OBS TYPES"},{"IONA","IONOSPHERIC CORR"},{"IONB","IONOSPHERIC CORR"},{"DUTC","TIME SYSTEM CORR"},{"MODEL","MET SENSOR MOD/TYPE/ACC"},{"POSITION","MET SENSOR POS XYZ/H"}}; auto it=m.find(x); if(it!=m.end()) return it->second; std::replace(x.begin(),x.end(),'_',' '); return x; }
void applyHeaderEdits(RinexFile& rf,const std::map<std::string,std::string>& edits){ 
  if(edits.empty()) return; 
  std::map<std::string,std::string> k2v; for(auto& [raw,val]:edits) k2v[keyOf(raw)]=val;
  // Structured receiver fields: teqc -O.rn/-O.rt/-O.rv modify columns of REC # / TYPE / VERS.
  if(k2v.count("RN")||k2v.count("RT")||k2v.count("RV")){
    auto old=existingValue(rf,"REC # / TYPE / VERS"); std::string rn=col(old,0,20), rt=col(old,20,20), rv=col(old,40,20);
    if(k2v.count("RN")) rn=k2v["RN"]; if(k2v.count("RT")) rt=k2v["RT"]; if(k2v.count("RV")) rv=k2v["RV"];
    std::ostringstream v; v<<std::left<<std::setw(20)<<rn.substr(0,20)<<std::setw(20)<<rt.substr(0,20)<<std::setw(20)<<rv.substr(0,20); setLine(rf,"REC # / TYPE / VERS",v.str());
  }
  if(k2v.count("AN")||k2v.count("AT")){
    auto old=existingValue(rf,"ANT # / TYPE"); std::string an=col(old,0,20), at=col(old,20,40);
    if(k2v.count("AN")) an=k2v["AN"]; if(k2v.count("AT")) at=k2v["AT"];
    std::ostringstream v; v<<std::left<<std::setw(20)<<an.substr(0,20)<<std::setw(40)<<at.substr(0,40); setLine(rf,"ANT # / TYPE",v.str());
  }
  for(auto& [rawLabel,val]:edits){ 
    auto k=keyOf(rawLabel); if(isSpecial(k)) continue;
    if(k=="RN"||k=="RT"||k=="RV"||k=="AN"||k=="AT") continue;
    auto label=normalizeHeaderLabel(rawLabel); 
    if(label=="COMMENT") { setLine(rf,label,val,true); continue; }
    if(rf.header.kind==RinexKind::Met && (k=="MODEL"||k=="POSITION")) {
      setLine(rf,"COMMENT",std::string(k=="MODEL"?"MET model: ":"MET position: ")+val,true);
      continue;
    }
    if(k=="PG") { setLine(rf,"APPROX POSITION XYZ",geodeticToXYZ(val)); continue; }
    if(k=="SLANT"||k=="SL") { setLine(rf,"ANTENNA: DELTA H/E/N",slantToHEN(val)); continue; }
    if(label=="APPROX POSITION XYZ"||label=="ANTENNA: DELTA H/E/N") { setLine(rf,label,fixed3(val)); continue; }
    if(label=="TIME OF FIRST OBS") { setLine(rf,label,formatRinexObsTime(val)); continue; }
    if(label=="WAVELENGTH FACT L1/2") { setLine(rf,label,(k=="MOD_WF"?formatModWF(val):formatWF(val)), k=="MOD_WF"); continue; }
    if(label=="IONOSPHERIC CORR") { setLine(rf,label,formatIonCorr(val,k=="IONA"),true); continue; }
    if(label=="TIME SYSTEM CORR") { setLine(rf,label,formatTimeSystemCorr(val),true); continue; }
    if(label=="PGM / RUN BY / DATE") { setLine(rf,label,fixedPGM(val)); continue; }
    setLine(rf,label,val); 
  } 
  rebuild(rf.header); parseContent(rf); 
}
void applyObsTypeFilter(RinexFile& rf,const std::string& obsList,bool exclude,bool renameOnly){ auto codes=splitCodes(obsList); if(codes.empty()) return; for(auto& rec:rf.data.observationRecords){ std::vector<ObservationValue> nv; for(size_t i=0;i<rec.values.size();++i){ bool in=std::find(codes.begin(),codes.end(),rec.values[i].type)!=codes.end(); if(renameOnly && i<codes.size()){ auto v=rec.values[i]; v.type=codes[i]; nv.push_back(v); } else if(exclude ? !in : in) nv.push_back(rec.values[i]); } rec.values=nv; }
  rf.header.obsTypes.clear(); for(auto& r:rf.data.observationRecords) for(auto& v:r.values) rf.header.obsTypes[r.system].push_back(v.type); for(auto& [s,v]:rf.header.obsTypes){ std::sort(v.begin(),v.end()); v.erase(std::unique(v.begin(),v.end()),v.end()); } }
void applyMetTypeFilter(RinexFile& rf,const std::string& obsList,bool exclude,bool renameOnly){ auto codes=splitCodes(obsList); if(codes.empty()) return; for(auto& r:rf.data.meteorologicalRecords){ if(renameOnly){ std::map<std::string,double> nv; size_t i=0; for(auto& kv:r.values){ nv[i<codes.size()?codes[i++]:kv.first]=kv.second; } r.values=std::move(nv); } else { for(auto it=r.values.begin(); it!=r.values.end();){ bool in=std::find(codes.begin(),codes.end(),it->first)!=codes.end(); if(exclude ? in : !in) it=r.values.erase(it); else ++it; } } } }
}
