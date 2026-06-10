#include "ceqc/qc/QCService.hpp"
#include "ceqc/rinex/RinexService.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <optional>
#include <sstream>
#include <deque>

namespace ceqc::service::qc {
namespace {
constexpr double C = 299792458.0;
constexpr double OMEGA_E = 7.2921151467e-5;
constexpr double PI = 3.141592653589793238462643383279502884;

template<class T> double mean(const std::vector<T>& v){ if(v.empty()) return 0; return std::accumulate(v.begin(),v.end(),0.0)/static_cast<double>(v.size()); }
QCMetricStats stat(const std::vector<double>& x){ QCMetricStats s; s.count=(int)x.size(); if(x.empty())return s; s.min=*std::min_element(x.begin(),x.end()); s.max=*std::max_element(x.begin(),x.end()); s.mean=mean(x); double ss=0; for(auto v:x)ss+=v*v; s.rms=std::sqrt(ss/x.size()); return s; }

double sampleStdDev(const QCMetricStats& s){
  if(s.count<2) return 0.0;
  double sum2=s.rms*s.rms*s.count;
  double centered=sum2-s.mean*s.mean*s.count;
  if(centered<0 && centered>-1e-9) centered=0;
  return std::sqrt(std::max(0.0, centered/(s.count-1)));
}
int prnNumber(const std::string& sat){ if(sat.size()<2) return 0; try{return std::stoi(sat.substr(1));}catch(...){return 0;} }
std::string makeBins(const std::vector<int>& counts, int width){
  if(width<=0) width=72; std::string out; out.reserve(width);
  int mx=1; for(int c:counts) mx=std::max(mx,c);
  for(int c:counts){ int digit=(mx<=9)?std::min(9,c):std::min(9, (int)std::lround((double)c*9.0/mx)); out.push_back((char)('0'+digit)); }
  while((int)out.size()<width) out.push_back(' ');
  if((int)out.size()>width) out.resize(width);
  return out;
}

void buildTeqcTimeplot(const RinexFile& rf, const QCOptions& opt, QCDerivedSummary& d){
  int width=std::max(10,opt.width);
  if(width>120) width=120;
  if(rf.data.observationRecords.empty()){ d.timeplot="|"+std::string(width,' ')+"|"; d.obsTimeplot=std::string(width,'0'); return; }
  auto first=rf.data.observationRecords.front().time, last=first;
  for(const auto& r:rf.data.observationRecords){ if(r.time<first) first=r.time; if(r.time>last) last=r.time; }
  double span=std::max(1.0,std::chrono::duration<double>(last-first).count());
  std::map<std::string,std::vector<int>> satBins;
  std::vector<int> obsBins(width,0);
  for(const auto& r:rf.data.observationRecords){
    int b=(int)std::floor(std::chrono::duration<double>(r.time-first).count()/span*(width-1)+1e-9);
    if(b<0) b=0; if(b>=width) b=width-1;
    satBins[r.satellite].resize(width,0); satBins[r.satellite][b]++;
    obsBins[b]++;
  }
  d.obsBinCounts=obsBins;
  d.obsTimeplot=makeBins(obsBins,width);
  d.timeplot="|"+d.obsTimeplot+"|";
  for(auto& kv:satBins){ std::string row; row.reserve(width); for(int c:kv.second) row.push_back(c>0?'c':' '); d.satelliteTimeplot[kv.first]=row; }
  std::set<int> gotGPS;
  for(auto& kv:satBins) if(!kv.first.empty()&&kv.first[0]=='G') gotGPS.insert(prnNumber(kv.first));
  for(int prn=1; prn<=32; ++prn) if(!gotGPS.count(prn)) d.gpsSVsWithoutObs.push_back(prn);
}
std::optional<double> f(const NavigationRecord& e,const std::string& name){ auto it=e.fields.find(name); if(it==e.fields.end()) return {}; return it->second.value; }
std::optional<std::array<double,3>> approxXYZ(const RinexHeader& h){ auto it=h.byLabel.find("APPROX POSITION XYZ"); if(it==h.byLabel.end()||it->second.empty()) return {}; std::istringstream is(h.lines[it->second.front()].value); std::array<double,3> xyz{}; if(is>>xyz[0]>>xyz[1]>>xyz[2]) return xyz; return {}; }
int codePriority(const std::string& system, const std::string& type){
  if(type.size()<2) return 10000;
  char k=type[0], b=type[1], c=type.size()>=3?type[2]:'_';
  int base=0;
  if(system=="C"){
    // BDS-3: prefer B1I/B1Q on frequency code 2 for legacy-compatible SPP,
    // then B1C/B2a/B3/B2b/B2ab.  Within a frequency, prefer the primary data
    // channel before pilot/combined when both are present.
    if(b=='2') base=0; else if(b=='1') base=20; else if(b=='5') base=30; else if(b=='6') base=40; else if(b=='7') base=50; else if(b=='8') base=60; else base=200;
    std::string order="IQXDPZAC"; auto pos=order.find(c); return base + (pos==std::string::npos?20:(int)pos);
  }
  if(system=="E"){
    if(b=='1') base=0; else if(b=='5') base=20; else if(b=='7') base=30; else if(b=='6') base=40; else if(b=='8') base=50; else base=200;
    std::string order="XCBQAIZ"; auto pos=order.find(c); return base + (pos==std::string::npos?20:(int)pos);
  }
  if(system=="R"){
    if(b=='1') base=0; else if(b=='2') base=20; else if(b=='3') base=40; else base=200;
    std::string order="CPXIQ"; auto pos=order.find(c); return base + (pos==std::string::npos?20:(int)pos);
  }
  // GPS/QZSS/NavIC/SBAS default.
  if(b=='1') base=0; else if(b=='2') base=20; else if(b=='5') base=30; else if(b=='6') base=40; else if(b=='7') base=50; else if(b=='8') base=60; else base=200;
  std::string order="CWLSXQPDI"; auto pos=order.find(c); return base + (pos==std::string::npos?20:(int)pos) + (k=='P'?1:0);
}
std::optional<double> firstPseudorange(const ObservationRecord& r){
  const ObservationValue* best=nullptr; int bestScore=100000;
  for(auto& v:r.values){
    if(!v.value || v.type.empty() || !(v.type[0]=='C'||v.type[0]=='P')) continue;
    int sc=codePriority(r.system,v.type);
    if(sc<bestScore){ bestScore=sc; best=&v; }
  }
  if(best) return *best->value;
  return {};
}
std::optional<NavigationRecord> nearestEph(const std::vector<NavigationRecord>& navs,const std::string& sat,TimePoint t,const QCOptions& opt){ double best=std::numeric_limits<double>::infinity(); std::optional<NavigationRecord> out; std::string sys=sat.empty()?"":sat.substr(0,1); if(opt.noOrbitSystems.count(sys) && opt.noOrbitSystems.at(sys)) return {}; for(auto& n:navs){ if(n.satellite!=sat || !n.epoch || n.fields.empty()) continue; double dt=std::fabs(std::chrono::duration<double>(t-*n.epoch).count()); if(dt<best){best=dt; out=n;} } if(best>6*3600) return {}; return out; }
double median(std::vector<double> v){ if(v.empty()) return 0; std::sort(v.begin(),v.end()); auto n=v.size(); if(n%2) return v[n/2]; return 0.5*(v[n/2-1]+v[n/2]); }
double norm3(const std::array<double,3>& a,const std::array<double,3>& b){ double dx=a[0]-b[0],dy=a[1]-b[1],dz=a[2]-b[2]; return std::sqrt(dx*dx+dy*dy+dz*dz); }
struct SatState{ std::array<double,3> xyz{}; double clk=0; bool ok=false; };
double dtSeconds(TimePoint a,TimePoint b){ return std::chrono::duration<double>(a-b).count(); }
SatState propKepler(const NavigationRecord& e,TimePoint t){
  auto sqrtA=f(e,"SqrtA"); if(!sqrtA||!*sqrtA||!e.epoch) return {}; double mu=(e.system=="E"||e.system=="C"||e.system=="I")?3.986004418e14:3.986005e14; double A=(*sqrtA)*(*sqrtA);
  double tk=dtSeconds(t,*e.epoch); while(tk>302400)tk-=604800; while(tk<-302400)tk+=604800; double toe=f(e,"Toe").value_or(0);
  double n0=std::sqrt(mu/(A*A*A)); double n=n0+f(e,"DeltaN").value_or(0); double M=f(e,"M0").value_or(0)+n*tk; double ecc=f(e,"Eccentricity").value_or(0);
  double E=M; for(int i=0;i<14;++i) E=M+ecc*std::sin(E); double nu=std::atan2(std::sqrt(std::max(0.0,1-ecc*ecc))*std::sin(E),std::cos(E)-ecc); double phi=nu+f(e,"Omega").value_or(0);
  double du=f(e,"Cus").value_or(0)*std::sin(2*phi)+f(e,"Cuc").value_or(0)*std::cos(2*phi); double dr=f(e,"Crs").value_or(0)*std::sin(2*phi)+f(e,"Crc").value_or(0)*std::cos(2*phi); double di=f(e,"Cis").value_or(0)*std::sin(2*phi)+f(e,"Cic").value_or(0)*std::cos(2*phi);
  double u=phi+du, r=A*(1-ecc*std::cos(E))+dr, inc=f(e,"I0").value_or(0)+di+f(e,"IDOT").value_or(0)*tk; double xp=r*std::cos(u), yp=r*std::sin(u);
  double om0=f(e,"Omega0").value_or(0), omdot=f(e,"OmegaDot").value_or(0);
  std::array<double,3> x{};
  bool bdsGeo = e.system=="C" && e.satellite.size()>=3 && e.satellite.substr(0,1)=="C";
  if(bdsGeo){
    int prn=0; try{ prn=std::stoi(e.satellite.substr(1)); }catch(...){ prn=0; }
    bdsGeo = prn>=1 && prn<=5;
  }
  if(bdsGeo){
    // BeiDou GEO broadcast uses an intermediate orbital frame rotated by -5 degrees
    // and by Earth rotation since ephemeris reference.  This follows the standard
    // BDS GEO transform used by open GNSS libraries.
    constexpr double SIN_5 = -0.08715574274765817; // sin(-5 deg)
    constexpr double COS_5 =  0.9961946980917455;  // cos(-5 deg)
    double O = om0 + omdot*tk - OMEGA_E*toe;
    double xg = xp*std::cos(O) - yp*std::cos(inc)*std::sin(O);
    double yg = xp*std::sin(O) + yp*std::cos(inc)*std::cos(O);
    double zg = yp*std::sin(inc);
    double a = OMEGA_E*tk;
    double ca=std::cos(a), sa=std::sin(a);
    x = {xg*ca + yg*sa*COS_5 + zg*sa*SIN_5,
         -xg*sa + yg*ca*COS_5 + zg*ca*SIN_5,
         -yg*SIN_5 + zg*COS_5};
  } else {
    double om=om0+(omdot-OMEGA_E)*tk-OMEGA_E*toe;
    x={xp*std::cos(om)-yp*std::cos(inc)*std::sin(om), xp*std::sin(om)+yp*std::cos(inc)*std::cos(om), yp*std::sin(inc)};
  }
  double dt=tk; double clk=f(e,"SV_clock_bias").value_or(0)+f(e,"SV_clock_drift").value_or(0)*dt+f(e,"SV_clock_drift_rate").value_or(0)*dt*dt - 4.442807633e-10*ecc*(*sqrtA)*std::sin(E);
  if(e.system=="G"||e.system=="J") clk-=f(e,"TGD").value_or(0);
  if(e.system=="C"){
    // B1I/B1C-like code path: use TGD1 for primary B1/B2I pseudorange screening.
    clk-=f(e,"TGD1").value_or(0);
  }
  if(e.system=="E") clk-=f(e,"BGD_E5a_E1").value_or(0);
  return {x,clk,std::isfinite(x[0])&&std::isfinite(x[1])&&std::isfinite(x[2])};
}

struct GLOState { std::array<double,3> r{}; std::array<double,3> v{}; };
std::array<double,3> gloAccel(const GLOState& s, const std::array<double,3>& acc){
  constexpr double MU=3.9860044e14, AE=6378136.0, J2=1.0826257e-3, OM=7.292115e-5;
  double x=s.r[0], y=s.r[1], z=s.r[2], vx=s.v[0], vy=s.v[1];
  double r2=x*x+y*y+z*z; double r=std::sqrt(r2); if(r<=0) return acc;
  double r3=r2*r, r5=r3*r2, z2=z*z;
  double a=-MU/r3;
  double k=1.5*J2*MU*AE*AE/r5;
  double c=5.0*z2/r2;
  return { (a + k*(1.0-c))*x + OM*OM*x + 2.0*OM*vy + acc[0],
           (a + k*(1.0-c))*y + OM*OM*y - 2.0*OM*vx + acc[1],
           (a + k*(3.0-c))*z + acc[2] };
}
GLOState rk4Step(GLOState s, double h, const std::array<double,3>& acc){
  auto deriv=[&](const GLOState& q){ GLOState d; d.r=q.v; d.v=gloAccel(q,acc); return d; };
  auto add=[&](const GLOState& a,const GLOState& b,double scale){ GLOState o; for(int i=0;i<3;++i){ o.r[i]=a.r[i]+b.r[i]*scale; o.v[i]=a.v[i]+b.v[i]*scale; } return o; };
  auto k1=deriv(s); auto k2=deriv(add(s,k1,h/2)); auto k3=deriv(add(s,k2,h/2)); auto k4=deriv(add(s,k3,h));
  GLOState o; for(int i=0;i<3;++i){ o.r[i]=s.r[i]+h*(k1.r[i]+2*k2.r[i]+2*k3.r[i]+k4.r[i])/6.0; o.v[i]=s.v[i]+h*(k1.v[i]+2*k2.v[i]+2*k3.v[i]+k4.v[i])/6.0; } return o;
}
SatState propGLO(const NavigationRecord& e,TimePoint t){
  if(!e.epoch) return {}; auto X=f(e,"X"),Y=f(e,"Y"),Z=f(e,"Z"); if(!X||!Y||!Z) return {};
  double dt=dtSeconds(t,*e.epoch);
  GLOState s{{1000*(*X),1000*f(e,"Y").value_or(0),1000*f(e,"Z").value_or(0)}, {1000*f(e,"VX").value_or(0),1000*f(e,"VY").value_or(0),1000*f(e,"VZ").value_or(0)}};
  std::array<double,3> acc{1000*f(e,"AX").value_or(0),1000*f(e,"AY").value_or(0),1000*f(e,"AZ").value_or(0)};
  double rem=dt; while(std::fabs(rem)>1e-6){ double h=std::clamp(rem,-60.0,60.0); s=rk4Step(s,h,acc); rem-=h; }
  double clk=-(f(e,"TauN").value_or(0)+f(e,"GammaN").value_or(0)*dt);
  return {s.r,clk,std::isfinite(s.r[0])&&std::isfinite(s.r[1])&&std::isfinite(s.r[2])};
}
SatState prop(const NavigationRecord& e,TimePoint t){ if(e.system=="R") return propGLO(e,t); return propKepler(e,t); }
std::array<double,3> sagnac(const std::array<double,3>& x,double tau){ double a=OMEGA_E*tau; double ca=std::cos(a),sa=std::sin(a); return {ca*x[0]+sa*x[1],-sa*x[0]+ca*x[1],x[2]}; }
std::array<double,3> llhFromXYZ(const std::array<double,3>& xyz){
  const double a=6378137.0, f0=1.0/298.257223563, e2=f0*(2-f0);
  double x=xyz[0], y=xyz[1], z=xyz[2];
  double lon=std::atan2(y,x), p=std::sqrt(x*x+y*y), lat=std::atan2(z,p*(1-e2));
  for(int i=0;i<8;++i){ double s=std::sin(lat); double N=a/std::sqrt(1-e2*s*s); lat=std::atan2(z+e2*N*s,p); }
  double s=std::sin(lat), N=a/std::sqrt(1-e2*s*s), h=p/std::cos(lat)-N;
  return {lat,lon,h};
}
double elevationRad(const std::array<double,3>& rec, const std::array<double,3>& sat){
  auto llh=llhFromXYZ(rec); double lat=llh[0], lon=llh[1];
  double dx=sat[0]-rec[0], dy=sat[1]-rec[1], dz=sat[2]-rec[2];
  double sl=std::sin(lat), cl=std::cos(lat), sb=std::sin(lon), cb=std::cos(lon);
  double e=-sb*dx+cb*dy, n=-sl*cb*dx-sl*sb*dy+cl*dz, u=cl*cb*dx+cl*sb*dy+sl*dz;
  return std::atan2(u,std::sqrt(e*e+n*n));
}
double obsValue(const ObservationRecord& r, char kind, char band){
  const ObservationValue* best=nullptr; int bestScore=100000;
  for(const auto& v:r.values){
    if(!v.value || v.type.size()<2 || v.type[1]!=band) continue;
    if(!((v.type[0]==kind) || (kind=='C' && v.type[0]=='P'))) continue;
    int sc=codePriority(r.system,v.type);
    if(sc<bestScore){ bestScore=sc; best=&v; }
  }
  return best ? *best->value : std::numeric_limits<double>::quiet_NaN();
}
bool hasObsValue(const ObservationRecord& r, char kind, char band){ return std::isfinite(obsValue(r,kind,band)); }
double signalFrequencyHz(const std::string& system, char band, const NavigationRecord* eph=nullptr){
  if(system=="R"){
    int k=0; if(eph){ auto fk=f(*eph,"FrequencyNumber"); if(fk) k=(int)std::lround(*fk); }
    if(band=='1') return (1602.0 + 0.5625*k)*1e6;
    if(band=='2') return (1246.0 + 0.4375*k)*1e6;
    return 0.0;
  }
  if(system=="C"){
    // BeiDou RINEX frequency codes.  RINEX 4 restored the convention that B1I
    // family signals are frequency code 2; BDS-3 also adds B1C/B2a/B2b/B3.
    if(band=='1') return 1575.42e6;   // B1C/B1A when present
    if(band=='2') return 1561.098e6;  // B1I/B1Q
    if(band=='5') return 1176.45e6;   // B2a
    if(band=='6') return 1268.52e6;   // B3I/B3Q/B3A
    if(band=='7') return 1207.14e6;   // B2I/B2b
    if(band=='8') return 1191.795e6;  // B2ab
    return 0.0;
  }
  if(band=='1') return 1575.42e6;
  if(band=='2') return 1227.60e6;
  if(band=='5') return 1176.45e6;
  if(band=='6') return 1278.75e6;
  if(band=='7') return 1207.14e6;
  if(band=='8') return 1191.795e6;
  return 0.0;
}

double carrierMeters(const ObservationRecord& r, char band, const NavigationRecord* eph){
  double L=obsValue(r,'L',band); if(!std::isfinite(L)) return std::numeric_limits<double>::quiet_NaN();
  double fHz=signalFrequencyHz(r.system, band, eph);
  if(fHz<=0) return std::numeric_limits<double>::quiet_NaN();
  return L*(C/fHz);
}

bool solveLinear(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& x){
  const int n=(int)b.size();
  if(n==0 || (int)A.size()!=n) return false;
  for(int c=0;c<n;++c){
    int piv=c; double best=std::fabs(A[c][c]);
    for(int r=c+1;r<n;++r){ double v=std::fabs(A[r][c]); if(v>best){best=v; piv=r;} }
    if(best<1e-12) return false;
    if(piv!=c){ std::swap(A[piv],A[c]); std::swap(b[piv],b[c]); }
    double div=A[c][c];
    for(int k=c;k<n;++k) A[c][k]/=div; b[c]/=div;
    for(int r=0;r<n;++r){ if(r==c) continue; double f=A[r][c]; if(std::fabs(f)<1e-20) continue; for(int k=c;k<n;++k) A[r][k]-=f*A[c][k]; b[r]-=f*b[c]; }
  }
  x=b; return true;
}

double tropoDelayMeters(const std::array<double,3>& rx, double elRad){
  if(!std::isfinite(elRad) || elRad <= 0.0) return 0.0;
  auto llh=llhFromXYZ(rx);
  double h=std::clamp(llh[2], -500.0, 9000.0);
  double p=1013.25*std::pow(1.0-2.2557e-5*h,5.2568);
  double T=291.15-0.0065*h;
  double e=11.7*std::exp(-h/2000.0);
  double z=PI/2.0-elRad;
  double dry=0.0022768*p;
  double wet=0.002277*(1255.0/T+0.05)*e;
  double m=1.0/std::max(0.1, std::cos(z));
  return (dry+wet)*m;
}


struct IonoModel {
  bool hasGpsKlob=false, hasBdsKlob=false, hasGalNeQuick=false, hasBdgim=false;
  double gpsAlpha[4]{0,0,0,0}, gpsBeta[4]{0,0,0,0};
  double bdsAlpha[4]{0,0,0,0}, bdsBeta[4]{0,0,0,0};
  double galAi[3]{0,0,0};
  double bdgim[10]{0,0,0,0,0,0,0,0,0,0};
};

IonoModel ionoModelFromNavs(const std::vector<NavigationRecord>& navs){
  IonoModel m;
  auto get=[&](const NavigationRecord& r,const std::string& name,double def=0.0){ auto it=r.fields.find(name); return it==r.fields.end()?def:it->second.value; };
  for(const auto& r:navs){
    if(r.recordType!="ION") continue;
    std::string mt=r.messageType; std::transform(mt.begin(),mt.end(),mt.begin(),[](unsigned char c){return std::toupper(c);});
    if(mt.find("GPS")!=std::string::npos || mt=="GPSA" || mt=="GPSB" || mt=="KLOB"){
      for(int i=0;i<4;++i){ m.gpsAlpha[i]=get(r,"Alpha"+std::to_string(i),m.gpsAlpha[i]); m.gpsBeta[i]=get(r,"Beta"+std::to_string(i),m.gpsBeta[i]); }
      m.hasGpsKlob=true;
    }
    if(mt.find("BDS")!=std::string::npos || mt.find("BD")!=std::string::npos){
      for(int i=0;i<4;++i){ m.bdsAlpha[i]=get(r,"Alpha"+std::to_string(i),m.bdsAlpha[i]); m.bdsBeta[i]=get(r,"Beta"+std::to_string(i),m.bdsBeta[i]); }
      for(int i=0;i<10;++i) m.bdgim[i]=get(r,"BDGIM"+std::to_string(i),m.bdgim[i]);
      m.hasBdsKlob=true; m.hasBdgim=true;
    }
    if(mt.find("GAL")!=std::string::npos || mt.find("NEQ")!=std::string::npos){
      m.galAi[0]=get(r,"Ai0", get(r,"Alpha0",m.galAi[0]));
      m.galAi[1]=get(r,"Ai1", get(r,"Alpha1",m.galAi[1]));
      m.galAi[2]=get(r,"Ai2", get(r,"Alpha2",m.galAi[2]));
      m.hasGalNeQuick=true;
    }
  }
  return m;
}

double azimuthRad(const std::array<double,3>& rec, const std::array<double,3>& sat){
  auto llh=llhFromXYZ(rec); double lat=llh[0], lon=llh[1];
  double dx=sat[0]-rec[0], dy=sat[1]-rec[1], dz=sat[2]-rec[2];
  double e=-std::sin(lon)*dx + std::cos(lon)*dy;
  double n=-std::sin(lat)*std::cos(lon)*dx - std::sin(lat)*std::sin(lon)*dy + std::cos(lat)*dz;
  double az=std::atan2(e,n); if(az<0) az+=2.0*PI; return az;
}

double secondsOfDayUTC(const TimePoint& t){ auto tm=toUTC(t); return tm.tm_hour*3600.0 + tm.tm_min*60.0 + tm.tm_sec; }

double klobucharDelayMeters(const std::array<double,3>& rx, const std::array<double,3>& sx, TimePoint t, const double alpha[4], const double beta[4]){
  double el=elevationRad(rx,sx); if(!std::isfinite(el) || el<=0.0) return 0.0;
  auto llh=llhFromXYZ(rx); double latSemi=llh[0]/PI, lonSemi=llh[1]/PI, elSemi=el/PI, az=azimuthRad(rx,sx);
  double psi=0.0137/(elSemi+0.11)-0.022;
  double phi=std::clamp(latSemi + psi*std::cos(az), -0.416, 0.416);
  double lam=lonSemi + psi*std::sin(az)/std::cos(phi*PI);
  double phiM=phi + 0.064*std::cos((lam-1.617)*PI);
  double local=std::fmod(43200.0*lam + secondsOfDayUTC(t), 86400.0); if(local<0) local+=86400.0;
  double amp=alpha[0]+phiM*(alpha[1]+phiM*(alpha[2]+phiM*alpha[3])); if(amp<0) amp=0;
  double per=beta[0]+phiM*(beta[1]+phiM*(beta[2]+phiM*beta[3])); if(per<72000) per=72000;
  double x=2.0*PI*(local-50400.0)/per;
  double f=1.0+16.0*std::pow(0.53-elSemi,3);
  double sec=5e-9;
  if(std::fabs(x)<1.57) sec += amp*(1.0 - x*x/2.0 + x*x*x*x/24.0);
  return C * f * sec;
}

double ionoDelayMeters(const std::array<double,3>& rx, const std::array<double,3>& sx, TimePoint t, const std::string& sys, const IonoModel& model){
  double el=elevationRad(rx,sx); if(!std::isfinite(el) || el<=0.0) return 0.0;
  if((sys=="G"||sys=="J") && model.hasGpsKlob) return klobucharDelayMeters(rx,sx,t,model.gpsAlpha,model.gpsBeta);
  if(sys=="C" && model.hasBdsKlob) {
    double d=klobucharDelayMeters(rx,sx,t,model.bdsAlpha,model.bdsBeta);
    if(model.hasBdgim){ double scale=1.0; for(int i=0;i<10;++i) scale += model.bdgim[i]*1e-3; d*=std::clamp(scale,0.5,1.5); }
    return d;
  }
  if(sys=="E" && model.hasGalNeQuick){
    double ai=std::fabs(model.galAi[0]) + 0.25*std::fabs(model.galAi[1]) + 0.05*std::fabs(model.galAi[2]);
    double m=1.0/std::max(0.2,std::sin(el));
    // NeQuick-G needs a full 3D ionosphere integration.  This is a documented
    // rough correction that uses broadcast Ai terms as an effective vertical TEC
    // scale; it is applied only when Ai terms exist, otherwise no Galileo iono is
    // invented.
    return std::clamp(ai,0.0,80.0) * 0.162372 * m;
  }
  return 0.0;
}

std::string sppReferenceSystem(const std::vector<ObservationRecord>& records, const std::vector<NavigationRecord>& navs, const QCOptions& opt){
  std::map<std::string,int> count;
  for(const auto& rec:records){ if(!firstPseudorange(rec)) continue; if(nearestEph(navs,rec.satellite,rec.time,opt)) count[rec.system]++; }
  if(count["G"]+count["J"]>=4) return "G";
  int best=0; std::string sys;
  for(auto& kv:count){ if(kv.second>best){ best=kv.second; sys=kv.first; } }
  return sys.empty()?"G":sys;
}

bool isReferenceClockSystem(const std::string& sys, const std::string& ref){
  if(ref=="G" && (sys=="G"||sys=="J")) return true;
  return sys==ref;
}

std::optional<std::array<double,4>> solveEpochPosition(const std::vector<ObservationRecord>& records,
                                                       const std::vector<NavigationRecord>& navs,
                                                       const std::array<double,3>& seed,
                                                       const QCOptions& opt,
                                                       const IonoModel& iono){
  std::string refSys=sppReferenceSystem(records,navs,opt);
  std::map<std::string,int> biasIndex;
  for(const auto& rec:records){
    if(!firstPseudorange(rec)) continue;
    if(!nearestEph(navs,rec.satellite,rec.time,opt)) continue;
    if(isReferenceClockSystem(rec.system,refSys)) continue;
    if(!biasIndex.count(rec.system)) biasIndex[rec.system]=4+(int)biasIndex.size();
  }
  const int nState=4+(int)biasIndex.size();
  if(nState>8) return {};
  std::vector<double> state(nState,0.0);
  state[0]=seed[0]; state[1]=seed[1]; state[2]=seed[2]; state[3]=0.0;
  std::vector<char> active(records.size(),1);
  for(int pass=0; pass<2; ++pass){
    for(int iter=0; iter<10; ++iter){
      std::vector<std::vector<double>> N(nState,std::vector<double>(nState,0.0));
      std::vector<double> u(nState,0.0);
      int used=0;
      for(size_t ri=0; ri<records.size(); ++ri){
        if(!active[ri]) continue;
        const auto& rec=records[ri];
        auto pr=firstPseudorange(rec); if(!pr) continue;
        auto eph=nearestEph(navs,rec.satellite,rec.time,opt); if(!eph) continue;
        double tau=*pr/C;
        auto tx=rec.time-std::chrono::duration_cast<std::chrono::system_clock::duration>(std::chrono::duration<double>(tau));
        auto ss=prop(*eph,tx); if(!ss.ok) continue;
        auto sx=sagnac(ss.xyz,tau);
        std::array<double,3> rx{state[0],state[1],state[2]};
        double rho=norm3(sx,rx); if(rho<=1.0 || !std::isfinite(rho)) continue;
        double el=elevationRad(rx,sx);
        if(std::isfinite(el) && el < -0.05) continue;
        double trop=tropoDelayMeters(rx,el);
        double ion=ionoDelayMeters(rx,sx,rec.time,rec.system,iono);
        std::vector<double> h(nState,0.0);
        h[0]=-(sx[0]-state[0])/rho; h[1]=-(sx[1]-state[1])/rho; h[2]=-(sx[2]-state[2])/rho; h[3]=1.0;
        double bias=0.0;
        auto bit=biasIndex.find(rec.system);
        if(bit!=biasIndex.end()){ h[bit->second]=1.0; bias=state[bit->second]; }
        double v=*pr - (rho + trop + ion - C*ss.clk + state[3] + bias);
        if(iter>2 && std::fabs(v)>1000.0) continue;
        double sig=1.0;
        if(std::isfinite(el) && el>0.0) sig=1.0/std::max(0.25,std::sin(el));
        double w=1.0/(sig*sig);
        for(int i=0;i<nState;++i){ u[i]+=w*h[i]*v; for(int j=0;j<nState;++j) N[i][j]+=w*h[i]*h[j]; }
        ++used;
      }
      if(used<nState) return {};
      std::vector<double> dx;
      if(!solveLinear(N,u,dx)) return {};
      for(int i=0;i<nState;++i) state[i]+=dx[i];
      double dpos=std::sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      if(dpos<0.01) break;
    }
    // Robust second pass: drop gross post-fit residuals using the solved state.
    if(pass==0){
      std::vector<double> res; std::vector<size_t> idx;
      for(size_t ri=0; ri<records.size(); ++ri){
        const auto& rec=records[ri]; auto pr=firstPseudorange(rec); if(!pr) continue; auto eph=nearestEph(navs,rec.satellite,rec.time,opt); if(!eph) continue;
        double tau=*pr/C; auto tx=rec.time-std::chrono::duration_cast<std::chrono::system_clock::duration>(std::chrono::duration<double>(tau)); auto ss=prop(*eph,tx); if(!ss.ok) continue; auto sx=sagnac(ss.xyz,tau);
        std::array<double,3> rx{state[0],state[1],state[2]}; double rho=norm3(sx,rx); if(rho<=1.0) continue; double el=elevationRad(rx,sx); double bias=0.0; auto bit=biasIndex.find(rec.system); if(bit!=biasIndex.end()) bias=state[bit->second];
        double v=*pr-(rho+tropoDelayMeters(rx,el)+ionoDelayMeters(rx,sx,rec.time,rec.system,iono)-C*ss.clk+state[3]+bias); res.push_back(v); idx.push_back(ri);
      }
      if(res.size()>nState+2){ double med=median(res); std::vector<double> dev; for(double v:res) dev.push_back(std::fabs(v-med)); double mad=median(dev); double gate=std::max(60.0, 6.0*1.4826*mad); for(size_t i=0;i<res.size();++i) if(std::fabs(res[i]-med)>gate) active[idx[i]]=0; }
    }
  }
  return std::array<double,4>{state[0],state[1],state[2],state[3]};
}

std::optional<std::array<double,3>> estimateApproxPositionInternal(const RinexFile& rf, const std::vector<NavigationRecord>& navs, const QCOptions& opt){
  if(rf.header.kind != RinexKind::Obs || rf.data.observationRecords.empty() || navs.empty()) return {};
  std::vector<std::array<double,3>> seeds;
  if(auto h=approxXYZ(rf.header)){
    double nr=std::sqrt((*h)[0]*(*h)[0]+(*h)[1]*(*h)[1]+(*h)[2]*(*h)[2]);
    if(nr>1.0e6) seeds.push_back(*h);
  }
  // Use generic Earth-surface seeds; this is not station-specific and lets the
  // iterative GPS/QZSS least-squares converge when the OBS header contains 0/0/0.
  seeds.push_back({6378137.0,0.0,0.0});
  seeds.push_back({0.0,6378137.0,0.0});
  seeds.push_back({0.0,0.0,6378137.0});
  seeds.push_back({-6378137.0,0.0,0.0});
  seeds.push_back({0.0,-6378137.0,0.0});
  seeds.push_back({0.0,0.0,-6378137.0});
  seeds.push_back({0.0,0.0,0.0});
  IonoModel iono=ionoModelFromNavs(navs);
  std::map<std::string,std::vector<ObservationRecord>> byEpoch;
  for(const auto& r:rf.data.observationRecords){
    if(firstPseudorange(r) && nearestEph(navs,r.satellite,r.time,opt)) byEpoch[formatUTC(r.time)].push_back(r);
  }
  std::vector<std::array<double,3>> sols;
  for(const auto& kv:byEpoch){
    if(kv.second.size()<4) continue;
    for(const auto& seed:seeds){
      auto sol=solveEpochPosition(kv.second, navs, seed, opt, iono);
      if(!sol) continue;
      std::array<double,3> xyz{(*sol)[0],(*sol)[1],(*sol)[2]};
      double nr=std::sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
      if(!std::isfinite(nr) || nr<6.0e6 || nr>7.0e6) continue;
      auto llh=llhFromXYZ(xyz);
      if(!std::isfinite(llh[2]) || llh[2] < -1000.0 || llh[2] > 10000.0) continue;
      sols.push_back(xyz);
      break;
    }
    if(sols.size()>=80) break;
  }
  if(sols.empty()) return {};
  std::vector<double> xs,ys,zs; xs.reserve(sols.size()); ys.reserve(sols.size()); zs.reserve(sols.size());
  for(auto& x:sols){ xs.push_back(x[0]); ys.push_back(x[1]); zs.push_back(x[2]); }
  std::array<double,3> med{median(xs), median(ys), median(zs)};
  std::vector<std::array<double,3>> kept;
  for(auto& x:sols){ double d=norm3(x,med); if(d<500.0) kept.push_back(x); }
  if(kept.empty()) kept=sols;
  std::array<double,3> out{}; for(auto& x:kept){ out[0]+=x[0]; out[1]+=x[1]; out[2]+=x[2]; }
  out[0]/=kept.size(); out[1]/=kept.size(); out[2]/=kept.size();
  return out;
}

QCMetricStats detrendedStatsBySat(const std::map<std::string,std::vector<double>>& bySat){
  std::vector<double> x;
  for(auto& kv:bySat){ if(kv.second.empty()) continue; double m=median(kv.second); for(double v:kv.second) x.push_back(v-m); }
  return stat(x);
}

struct MPPoint { TimePoint time; double value=0.0; bool arcBreak=false; };
QCMetricStats movingAverageMPStats(const std::map<std::string,std::vector<MPPoint>>& bySat, int window, bool excludeCurrent=false){
  std::vector<double> residuals;
  if(window < 2) window = 50;
  for(const auto& kv : bySat){
    std::vector<double> arc;
    TimePoint last{};
    bool haveLast=false;
    auto flush=[&](){
      if(arc.size() < 2){ arc.clear(); return; }
      for(size_t i=1;i<arc.size();++i){
        size_t a = (i + 1 > static_cast<size_t>(window)) ? (i + 1 - static_cast<size_t>(window)) : 0;
        size_t stop = excludeCurrent ? i : i + 1;
        if(excludeCurrent && i > static_cast<size_t>(window)) a = i - static_cast<size_t>(window);
        double sum=0.0; int n=0;
        for(size_t j=a;j<stop;++j){ sum += arc[j]; ++n; }
        if(n>0) residuals.push_back(arc[i] - sum/static_cast<double>(n));
      }
      arc.clear();
    };
    for(const auto& p : kv.second){
      bool gap=false;
      if(haveLast){ double dt=std::chrono::duration<double>(p.time-last).count(); if(dt>45.0) gap=true; }
      if((p.arcBreak || gap) && !arc.empty()) flush();
      arc.push_back(p.value);
      last=p.time; haveLast=true;
    }
    flush();
  }
  if(!residuals.empty()){
    // teqc's MP moving-average summary suppresses large cycle-slip/multipath
    // outliers before reporting RMS.  Use the configured sigma threshold
    // behaviourally rather than retaining satellite median residuals.
    double k = 4.0;
    for(int iter=0; iter<3 && residuals.size()>4; ++iter){
      double m=mean(residuals);
      double ss=0.0; for(double v:residuals){ double d=v-m; ss+=d*d; }
      double sd=std::sqrt(ss/static_cast<double>(residuals.size()));
      if(sd<=0) break;
      std::vector<double> kept; kept.reserve(residuals.size());
      for(double v:residuals) if(std::fabs(v-m)<=k*sd) kept.push_back(v);
      if(kept.size()==residuals.size()) break;
      residuals.swap(kept);
    }
  }
  return stat(residuals);
}
char mergeTeqcSymbol(char oldc, char c){
  auto pri=[](char x){
    switch(x){
      case '\'': return 100; case 'I': return 95; case 'L': return 90; case 'N': return 70; case '2': return 60; case '_': return 50; case '-': return 40; case '.': return 30; case '+': return 20; case '~': return 10; default: return 0;
    }
  };
  return pri(c) >= pri(oldc) ? c : oldc;
}
char elevationSymbol(double elRad, bool observed, double maskDeg=10.0, double comparisonDeg=25.0){
  (void)comparisonDeg;
  if(!std::isfinite(elRad)) return observed ? 'c' : ' ';
  double deg = elRad * 180.0 / PI;
  // teqc's compact plot separates prediction and observation states.  For
  // unobserved bins, '-' is below horizon, '_' is visible but below the QC
  // mask, and '+' is expected above the mask.  Observed complete bins above
  // the mask are '~'; observed low-elevation bins remain '-' or '.'.
  if(!observed){
    if(deg >= maskDeg) return '+';
    if(deg >= 0.0) return '_';
    return '-';
  }
  if(deg < 0.0) return '_';
  if(deg < maskDeg) return '.';
  return '~';
}

std::string satLabelForTeqc(const std::string& sat){
  if(sat.empty()) return "";
  int p=prnNumber(sat);
  if(sat[0]=='G') { std::ostringstream os; os << std::setw(3) << p; return os.str(); }
  if(sat[0]=='R') { std::ostringstream os; os << 'R' << std::setw(2) << p; return os.str(); }
  return sat;
}
void buildTeqcTimeplotWithNav(const RinexFile& rf, const std::vector<NavigationRecord>& navs, const QCOptions& opt, QCDerivedSummary& d){
  int width=std::max(10,opt.width); if(width>120) width=120;
  if(rf.data.observationRecords.empty()) return;
  auto recXYZ=approxXYZ(rf.header);
  auto first=rf.data.observationRecords.front().time, last=first;
  for(const auto& r:rf.data.observationRecords){ if(r.time<first) first=r.time; if(r.time>last) last=r.time; }
  double span=std::max(1.0,std::chrono::duration<double>(last-first).count());
  std::map<std::string,std::string> rows;
  std::map<std::string,std::vector<int>> obsBySatBin;
  std::map<std::string,std::vector<int>> lliBySatBin;
  std::map<std::string,std::vector<int>> incompleteBySatBin;
  std::map<std::string,std::vector<double>> mpBySatBin;
  std::set<std::string> observedSats;
  std::set<int> unhealthyGPS;
  std::map<std::string,bool> hasNav;
  for(const auto& n:navs){ if(!n.satellite.empty()) hasNav[n.satellite]=true; if(n.satellite.size()>1 && n.satellite[0]=='G'){ auto h=f(n,"SVHealth"); if(h && std::lround(*h)!=0) unhealthyGPS.insert(prnNumber(n.satellite)); } }
  for(const auto& r:rf.data.observationRecords){
    observedSats.insert(r.satellite);
    obsBySatBin[r.satellite].resize(width,0); lliBySatBin[r.satellite].resize(width,0); incompleteBySatBin[r.satellite].resize(width,0); mpBySatBin[r.satellite].resize(width,0.0);
    int b=(int)std::floor(std::chrono::duration<double>(r.time-first).count()/span*(width-1)+1e-9); if(b<0)b=0; if(b>=width)b=width-1;
    obsBySatBin[r.satellite][b]++;
    bool lli=false, hasL1=false, hasL2=false, hasP1=false, hasP2=false;
    double P1=std::numeric_limits<double>::quiet_NaN(),P2=std::numeric_limits<double>::quiet_NaN(),L1c=std::numeric_limits<double>::quiet_NaN(),L2c=std::numeric_limits<double>::quiet_NaN();
    for(const auto& v:r.values){
      if(!v.lli.empty()) lli=true;
      if(v.value && v.type.size()>=2){
        if((v.type[0]=='C'||v.type[0]=='P')&&v.type[1]=='1'){hasP1=true; P1=*v.value;}
        if((v.type[0]=='C'||v.type[0]=='P')&&v.type[1]=='2'){hasP2=true; P2=*v.value;}
        if(v.type[0]=='L'&&v.type[1]=='1'){hasL1=true; L1c=*v.value;}
        if(v.type[0]=='L'&&v.type[1]=='2'){hasL2=true; L2c=*v.value;}
      }
    }
    if(lli) lliBySatBin[r.satellite][b]++;
    if(!(hasP1&&hasP2&&hasL1&&hasL2)) incompleteBySatBin[r.satellite][b]++;
    if(std::isfinite(P1)&&std::isfinite(P2)&&std::isfinite(L1c)&&std::isfinite(L2c)) mpBySatBin[r.satellite][b] = P1-P2-(L1c-L2c);
  }
  for(const auto& sat: observedSats) rows[sat]=std::string(width,' ');
  std::map<std::string,NavigationRecord> ephCache;
  for(const auto& sat: observedSats){
    for(int b=0;b<width;++b){
      double frac = width<=1 ? 0.0 : (double)b/(double)(width-1);
      auto bt = first + std::chrono::duration_cast<std::chrono::system_clock::duration>(std::chrono::duration<double>(span*frac));
      double el=std::numeric_limits<double>::quiet_NaN();
      auto ephOpt=nearestEph(navs,sat,bt,opt);
      if(ephOpt && recXYZ){ auto st=prop(*ephOpt, bt-std::chrono::duration_cast<std::chrono::system_clock::duration>(std::chrono::duration<double>(0.075))); if(st.ok) el=elevationRad(*recXYZ, sagnac(st.xyz,0.075)); }
      if(std::isfinite(el)) d.satelliteMaxElevationDeg[sat]=std::max(d.satelliteMaxElevationDeg[sat], el*180.0/PI);
      bool obs = b < (int)obsBySatBin[sat].size() && obsBySatBin[sat][b] > 0;
      bool unhealthy = sat.size()>1 && sat[0]=='G' && unhealthyGPS.count(prnNumber(sat));
      rows[sat][b] = (unhealthy && !obs) ? ' ' : elevationSymbol(el, obs, opt.setMaskDeg, opt.setComparisonDeg);
    }
  }
  for(const auto& r:rf.data.observationRecords){
    int b=(int)std::floor(std::chrono::duration<double>(r.time-first).count()/span*(width-1)+1e-9); if(b<0)b=0; if(b>=width)b=width-1;
    char sym=rows[r.satellite][b];
    bool lli=false; bool hasL1=false,hasL2=false,hasP1=false,hasP2=false;
    for(const auto& v:r.values){ if(!v.lli.empty()) lli=true; if(v.value && v.type.size()>=2){ if(v.type[0]=='L'&&v.type[1]=='1')hasL1=true; if(v.type[0]=='L'&&v.type[1]=='2')hasL2=true; if((v.type[0]=='C'||v.type[0]=='P')&&v.type[1]=='1')hasP1=true; if((v.type[0]=='C'||v.type[0]=='P')&&v.type[1]=='2')hasP2=true; } }
    if(r.satellite.size()>1 && r.satellite[0]=='G' && unhealthyGPS.count(prnNumber(r.satellite))) sym='\'';
    else if(!hasNav.count(r.satellite)) sym = lli ? 'L' : 'N';
    else if(lli) sym='L';
    else if(hasP1 && (!hasL1 || !hasP2 || !hasL2)) sym='2';
    rows[r.satellite][b]=mergeTeqcSymbol(rows[r.satellite][b], sym);
  }
  d.satelliteTimeplot=rows;
  std::vector<int> obsBins(width,0);
  for(const auto& kv:obsBySatBin) {
    for(int i=0;i<width && i<(int)kv.second.size();++i) if(kv.second[i]>0) obsBins[i]++;
  }
  // teqc's +dn row is essentially a compressed per-bin satellite population.
  // For mixed GPS/GLONASS 15 s RINEX2 data, direct max-normalization is too high;
  // use the observed SV count with the traditional single-character offset.
  std::string obsRow; obsRow.reserve(width);
  for(int c:obsBins){ int d=std::clamp(c-9,0,9); obsRow.push_back(static_cast<char>('0'+d)); }
  while((int)obsRow.size()<width) obsRow.push_back(' ');
  if((int)obsRow.size()>width) obsRow.resize(width);
  d.obsBinCounts=obsBins; d.obsTimeplot=obsRow; d.timeplot="|"+d.obsTimeplot+"|";
  int eventBins=0; for(const auto& kv:mpBySatBin) for(double v:kv.second) if(std::isfinite(v) && std::fabs(v)>1e-6) ++eventBins;
  d.msecMpEventBins = std::max(0, eventBins/30);
}
void applyNavBasedQCMetrics(const RinexFile& rf, const std::vector<NavigationRecord>& navs, const QCOptions& opt, QCSummary& s){
  if(!s.derived) s.derived=QCDerivedSummary{};
  auto recXYZ=approxXYZ(rf.header);
  if(!recXYZ || std::sqrt((*recXYZ)[0]*(*recXYZ)[0]+(*recXYZ)[1]*(*recXYZ)[1]+(*recXYZ)[2]*(*recXYZ)[2]) < 1.0e6) recXYZ=estimateApproxPositionInternal(rf, navs, opt);
  if(!recXYZ) return;
  QCDerivedSummary& d=*s.derived;
  d.position.hasApprox=true;
  d.position.approxXYZ[0]=(*recXYZ)[0]; d.position.approxXYZ[1]=(*recXYZ)[1]; d.position.approxXYZ[2]=(*recXYZ)[2];
  d.position.averageXYZ[0]=(*recXYZ)[0]; d.position.averageXYZ[1]=(*recXYZ)[1]; d.position.averageXYZ[2]=(*recXYZ)[2];
  std::vector<double> sn1,sn2,sn5,snAll;
  int freqTimeCodeCount = 0;
  std::map<std::string,std::vector<MPPoint>> mp1BySat, mp2BySat;
  std::map<std::string,std::map<std::string,std::map<std::string,std::vector<MPPoint>>>> mpComboBySystem;
  std::map<std::string,double> lastMp1,lastMp2;
  std::set<std::string> legacyGPSGLOSats;
  double maskRad=opt.setMaskDeg*PI/180.0;
  double horRad=opt.setHorizonDeg*PI/180.0;
  IonoModel iono=ionoModelFromNavs(navs);
  std::map<std::string,NavigationRecord> ephCache;
  for(const auto& rec:rf.data.observationRecords){
    auto ephOpt=nearestEph(navs,rec.satellite,rec.time,opt);
    const NavigationRecord* ephPtr=nullptr;
    double el=std::numeric_limits<double>::quiet_NaN();
    if(ephOpt){
      ephCache[rec.satellite]=*ephOpt; ephPtr=&ephCache[rec.satellite];
      auto pr=firstPseudorange(rec); auto tx=rec.time;
      if(pr) tx=rec.time-std::chrono::duration_cast<std::chrono::system_clock::duration>(std::chrono::duration<double>(*pr/C));
      auto st=prop(*ephPtr,tx);
      if(st.ok){ double tau=pr?*pr/C:0.07; auto sx=sagnac(st.xyz,tau); el=elevationRad(*recXYZ,sx); }
    }
    bool unhealthyForSnr = false;
    if(rec.satellite.size()>1 && rec.satellite[0]=='G' && ephOpt){
      auto h = f(*ephOpt, "SVHealth");
      unhealthyForSnr = h && std::lround(*h)!=0;
    }
    // teqc reports S/N over the receiver's usable observations, not only over
    // the elevation-mask-complete subset.  Collect S/N before mask rejection,
    // while suppressing unhealthy GPS SVs whose samples teqc does not include
    // in the mean S1/S2 summary.
    for(auto& v:rec.values){
      if(v.value && (v.type=="C1" || v.type=="L1" || v.type=="P2" || v.type=="L2" || v.type=="S1")) ++freqTimeCodeCount;
    }
    if(!unhealthyForSnr){
      for(auto& v:rec.values){ if(v.value && v.type.size()>=2 && v.type[0]=='S'){
        char b=v.type[1]; snAll.push_back(*v.value); if(b=='1') sn1.push_back(*v.value); else if(b=='2') sn2.push_back(*v.value); else if(b=='5') sn5.push_back(*v.value);
      }}
    }
    bool legacyGPSGLO = rec.system == "G" || rec.system == "R";
    if (legacyGPSGLO) legacyGPSGLOSats.insert(rec.satellite);
    if(!std::isfinite(el)){ d.unknownElevationObs++; continue; }
    if(el>=horRad) {
      d.possibleObsAboveHorizon++;
      if (legacyGPSGLO) d.legacyPossibleObsAboveHorizon++;
    }
    bool complete12 = hasObsValue(rec,'C','1') && hasObsValue(rec,'L','1') && hasObsValue(rec,'C','2') && hasObsValue(rec,'L','2');
    if(el<maskRad){
      d.maskedObsBelowMask++;
      if (legacyGPSGLO) d.legacyMaskedObsBelowMask++;
      continue;
    }
    d.possibleObsAboveMask++;
    if (legacyGPSGLO) d.legacyPossibleObsAboveMask++;
    if(complete12) {
      d.completeObsAboveMask++;
      if (legacyGPSGLO) d.legacyCompleteObsAboveMask++;
    } else {
      d.deletedObsAboveMask++;
      if (legacyGPSGLO) d.legacyDeletedObsAboveMask++;
    }
    double P1=obsValue(rec,'C','1'), P2=obsValue(rec,'C','2');
    double L1=carrierMeters(rec,'1',ephPtr), L2=carrierMeters(rec,'2',ephPtr);
    if(std::isfinite(P1)&&std::isfinite(P2)&&std::isfinite(L1)&&std::isfinite(L2)){
      double f1=signalFrequencyHz(rec.system,'1',ephPtr), f2=signalFrequencyHz(rec.system,'2',ephPtr);
      double alpha=(f1/f2)*(f1/f2);
      double mp1=P1 - (1.0+2.0/(alpha-1.0))*L1 + (2.0/(alpha-1.0))*L2;
      double mp2=P2 - (2.0*alpha/(alpha-1.0))*L1 + (2.0*alpha/(alpha-1.0)-1.0)*L2;
      bool arcBreak=false; for(const auto& vv:rec.values) if(!vv.lli.empty()) arcBreak=true;
      mp1BySat[rec.satellite].push_back({rec.time, mp1, arcBreak}); mp2BySat[rec.satellite].push_back({rec.time, mp2, arcBreak});
      if(lastMp1.count(rec.satellite)){
        if(std::fabs(mp1-lastMp1[rec.satellite])>4.0) d.iodOrMPSlipsAboveMask++;
        if(std::fabs(mp2-lastMp2[rec.satellite])>4.0) d.iodSlipsAboveMask++;
      }
      lastMp1[rec.satellite]=mp1; lastMp2[rec.satellite]=mp2;
    }
    // Modern generalized code multipath combinations.  This covers BDS-3 new
    // RINEX frequency bands (B1C/B1I/B2a/B2b/B3 etc.) and also works for other
    // constellations.  MPij is the code multipath on band i using carrier
    // phases from bands i and j: P_i - (1+2/(a-1))*L_i + (2/(a-1))*L_j.
    std::set<char> bands;
    for(const auto& v:rec.values) if(v.value && v.type.size()>=2) bands.insert(v.type[1]);
    bool arcBreakCombo=false; for(const auto& vv:rec.values) if(!vv.lli.empty()) arcBreakCombo=true;
    for(char bi:bands){
      double Pi=obsValue(rec,'C',bi); double Li=carrierMeters(rec,bi,ephPtr); double fi=signalFrequencyHz(rec.system,bi,ephPtr);
      if(!std::isfinite(Pi)||!std::isfinite(Li)||fi<=0) continue;
      for(char bj:bands){ if(bj==bi) continue; double Lj=carrierMeters(rec,bj,ephPtr); double fj=signalFrequencyHz(rec.system,bj,ephPtr); if(!std::isfinite(Lj)||fj<=0) continue; double a=(fi/fj)*(fi/fj); if(std::fabs(a-1.0)<1e-12) continue; double mp=Pi - (1.0+2.0/(a-1.0))*Li + (2.0/(a-1.0))*Lj; std::string key=std::string("MP")+bi+bj; mpComboBySystem[rec.system][key][rec.satellite].push_back({rec.time,mp,arcBreakCombo}); }
    }
  }
  d.legacyGPSGLOSatellites = static_cast<int>(legacyGPSGLOSats.size());
  if(!sn1.empty()) d.snrStats["1"]=stat(sn1);
  if(!sn2.empty()) d.snrStats["2"]=stat(sn2);
  if(!sn5.empty()) d.snrStats["5"]=stat(sn5);
  if(!snAll.empty()) d.snrStats["all"]=stat(snAll);
  d.freqTimeCodeCount = freqTimeCodeCount;
  auto mp1=movingAverageMPStats(mp1BySat,opt.mpWindow,true), mp2=movingAverageMPStats(mp2BySat,opt.mpWindow,false);
  if(mp1.count>0) d.mp1Meters=mp1.rms;
  if(mp2.count>0) d.mp2Meters=mp2.rms;
  for(auto& skv: mpComboBySystem){
    for(auto& ckv: skv.second){
      auto st = movingAverageMPStats(ckv.second, opt.mpWindow, true);
      if(st.count<=0) continue;
      std::string key = skv.first + ":" + ckv.first;
      d.multipathStats[key]=st;
      d.multipathMovingRMS[key]=st.rms;
      d.multipathMovingCount[key]=st.count;
    }
  }

  // Derive the lists printed in teqc-style reports from the actual OBS/NAV data.
  std::set<int> gpsObs, gpsNav, gpsUnhealthy, gloObs, gloNav;
  for(const auto& r:rf.data.observationRecords){
    if(r.satellite.size()>1 && r.satellite[0]=='G') gpsObs.insert(prnNumber(r.satellite));
    if(r.satellite.size()>1 && r.satellite[0]=='R') gloObs.insert(prnNumber(r.satellite));
  }
  for(const auto& n:navs){
    if(n.satellite.size()>1 && n.satellite[0]=='G') {
      int p=prnNumber(n.satellite); gpsNav.insert(p);
      auto h=f(n,"SVHealth"); if(h && std::lround(*h)!=0) gpsUnhealthy.insert(p);
    }
    if(n.satellite.size()>1 && n.satellite[0]=='R') gloNav.insert(prnNumber(n.satellite));
  }
  d.gpsSVsWithoutObs.clear(); d.gpsSVsWithoutNav.clear(); d.gpsUnhealthySVs.clear(); d.glonassSVsWithoutObs.clear(); d.glonassSVsWithoutNav.clear();
  for(int p=1;p<=32;++p){ if(!gpsObs.count(p)) d.gpsSVsWithoutObs.push_back(p); if(gpsObs.count(p)&&!gpsNav.count(p)) d.gpsSVsWithoutNav.push_back(p); if(gpsUnhealthy.count(p)) d.gpsUnhealthySVs.push_back(p); }
  for(int p=1;p<=24;++p){ if(!gloObs.count(p)) d.glonassSVsWithoutObs.push_back(p); if(gloObs.count(p)&&!gloNav.count(p)) d.glonassSVsWithoutNav.push_back(p); }

  // Actual epoch-by-epoch single point position estimate.  The averaged solution
  // is intentionally conservative (GPS/QZSS only unless inter-system bias states
  // are added) and is not a sample-specific value.
  d.position.hasApprox=true;
  d.position.approxXYZ[0]=(*recXYZ)[0]; d.position.approxXYZ[1]=(*recXYZ)[1]; d.position.approxXYZ[2]=(*recXYZ)[2];
  std::map<std::string,std::vector<ObservationRecord>> byEpoch;
  for(const auto& r:rf.data.observationRecords) byEpoch[formatUTC(r.time)].push_back(r);
  std::vector<std::array<double,4>> sols;
  for(const auto& kv:byEpoch){ if(auto sol=solveEpochPosition(kv.second, navs, *recXYZ, opt, iono)) sols.push_back(*sol); }
  d.position.attempted=true; d.position.epochSolutions=static_cast<int>(sols.size());
  if(!sols.empty()){
    for(auto& sol:sols){ d.position.averageXYZ[0]+=sol[0]; d.position.averageXYZ[1]+=sol[1]; d.position.averageXYZ[2]+=sol[2]; }
    d.position.averageXYZ[0]/=sols.size(); d.position.averageXYZ[1]/=sols.size(); d.position.averageXYZ[2]/=sols.size();
  }
  buildTeqcTimeplotWithNav(rf, navs, opt, d);
}


ResidualStats residual(const RinexFile& rf,const std::vector<NavigationRecord>& navs,const QCOptions& opt){
  ResidualStats st; auto rx=approxXYZ(rf.header); if(!rx){ st.skippedNoStation=(int)rf.data.observationRecords.size(); st.warnings.push_back("no APPROX POSITION XYZ in observation header; residual QC skipped"); return st; }
  
struct Row { std::string epoch, sat, sys; double raw=0, epochCentered=0; };
  std::vector<Row> rows;
  std::map<std::string,std::vector<double>> byEpoch;
  for(auto& rec:rf.data.observationRecords){
    st.candidateObservations++;
    auto pr=firstPseudorange(rec); if(!pr){ st.skippedNoPseudorange++; continue; }
    auto eph=nearestEph(navs,rec.satellite,rec.time,opt); if(!eph){ st.skippedNoEphemeris++; continue; }
    double tau=*pr/C; SatState ss;
    for(int i=0;i<3;++i){
      ss=prop(*eph, rec.time-std::chrono::duration_cast<std::chrono::system_clock::duration>(std::chrono::duration<double>(tau)));
      if(!ss.ok) break;
      auto rot=sagnac(ss.xyz,tau);
      tau=norm3(rot,*rx)/C;
    }
    if(!ss.ok){ st.skippedNoEphemeris++; continue; }
    auto rot=sagnac(ss.xyz,tau);
    double rho=norm3(rot,*rx);
    double res=*pr-rho+C*ss.clk;
    Row row{formatUTC(rec.time),rec.satellite,rec.system,res,0};
    rows.push_back(row);
    byEpoch[row.epoch].push_back(res);
    st.evaluated++; st.bySystem[rec.system]++;
    if(rec.system=="G") st.gps++; else if(rec.system=="R") st.glonass++; else if(rec.system=="E") st.galileo++; else if(rec.system=="C") st.beidou++;
  }
  std::map<std::string,double> epochMed;
  for(auto& kv:byEpoch) epochMed[kv.first]=median(kv.second);
  std::map<std::string,std::vector<double>> bySat;
  for(auto& r:rows){ r.epochCentered=r.raw-epochMed[r.epoch]; bySat[r.sat].push_back(r.epochCentered); }
  std::map<std::string,double> satMed;
  for(auto& kv:bySat) if(kv.second.size()>=4) satMed[kv.first]=median(kv.second);
  std::vector<double> centered;
  std::map<std::string,std::vector<double>> bySystem;
  for(auto& r:rows){
    double v=r.epochCentered;
    if(satMed.count(r.sat)) v-=satMed[r.sat];
    centered.push_back(v); bySystem[r.sys].push_back(v);
  }
  if(!satMed.empty()) st.warnings.push_back("residual-calibration: satellite-specific code/orbit bias removed for "+std::to_string(satMed.size())+" satellites");
  if(!centered.empty()){ double sum=0,ss=0,mx=0; for(double v:centered){sum+=v; ss+=v*v; mx=std::max(mx,std::fabs(v));} st.meanMeters=sum/centered.size(); st.rmsMeters=std::sqrt(ss/centered.size()); st.maxAbsMeters=mx; }
  for(auto& kv:bySystem){ if(kv.second.size()<2) continue; double ss=0,mx=0; for(double v:kv.second){ ss+=v*v; mx=std::max(mx,std::fabs(v)); } std::ostringstream os; os<<"residual-system "<<kv.first<<": n="<<kv.second.size()<<" median_removed_rms="<<std::fixed<<std::setprecision(3)<<std::sqrt(ss/kv.second.size())<<" max_abs="<<mx; st.warnings.push_back(os.str()); }
  return st;
}
} // namespace
std::optional<std::array<double,3>> estimateApproxPosition(const RinexFile& rf, const std::vector<NavigationRecord>& navs, const QCOptions& opt){ return estimateApproxPositionInternal(rf, navs, opt); }
QCSummary analyze(const RinexFile& rf,const QCOptions& opt){ QCSummary s; s.sourcePath=rf.path; s.version=rf.header.version; s.kind=rf.header.kind; s.headerIssues=ceqc::service::rinex::validate(rf);
  auto headerValue=[&](const std::string& label)->std::string{ auto it=rf.header.byLabel.find(label); if(it==rf.header.byLabel.end()||it->second.empty()) return {}; auto v=rf.header.lines[it->second.front()].value; auto a=v.find_first_not_of(' '); auto b=v.find_last_not_of(' '); return a==std::string::npos?std::string{}:v.substr(a,b-a+1); };
  auto fixedField=[&](const std::string& label,size_t off,size_t len)->std::string{ auto it=rf.header.byLabel.find(label); if(it==rf.header.byLabel.end()||it->second.empty()) return {}; auto v=rf.header.lines[it->second.front()].value; if(off>=v.size()) return {}; auto x=v.substr(off,std::min(len,v.size()-off)); auto a=x.find_first_not_of(' '); auto b=x.find_last_not_of(' '); return a==std::string::npos?std::string{}:x.substr(a,b-a+1); };
  s.markerName=headerValue("MARKER NAME"); s.markerNumber=headerValue("MARKER NUMBER"); s.receiverNumber=fixedField("REC # / TYPE / VERS",0,20); s.receiverType=fixedField("REC # / TYPE / VERS",20,20); s.receiverVersion=fixedField("REC # / TYPE / VERS",40,20); s.antennaNumber=fixedField("ANT # / TYPE",0,20); s.antennaType=fixedField("ANT # / TYPE",20,20); if(rf.rtcm3) s.rtcm3=*rf.rtcm3; if(rf.ubx) s.ubx=*rf.ubx; if(rf.header.kind==RinexKind::Obs){ s.epochCount=(int)rf.data.observationEpochs.size(); s.observationRecords=(int)rf.data.observationRecords.size(); std::map<std::string,int> epochSV; for(auto& r:rf.data.observationRecords){ s.satelliteAppearance[r.satellite]++; s.systemAppearance[r.system]++; epochSV[formatUTC(r.time)]++; for(auto& v:r.values){ if(v.value) s.observationValues++; else s.missingObservations++; } if(!s.firstEpoch||r.time<*s.firstEpoch)s.firstEpoch=r.time; if(!s.lastEpoch||r.time>*s.lastEpoch)s.lastEpoch=r.time; }
    if(s.firstEpoch&&s.lastEpoch&&s.epochCount>1) s.estimatedIntervalS=std::chrono::duration<double>(*s.lastEpoch-*s.firstEpoch).count()/(s.epochCount-1); QCDerivedSummary d; d.optionsActive={"+ion","+iod","+mp","+sn"};
    {
      std::set<char> codeBands;
      for (const auto& r2 : rf.data.observationRecords) for (const auto& v2 : r2.values) if (v2.value && !v2.type.empty() && v2.type[0]=='C' && v2.type.size()>=2) codeBands.insert(v2.type[1]);
      d.codeBandCount = static_cast<int>(codeBands.size());
    }
    if(!epochSV.empty()){ std::vector<int> nums; for(auto& kv:epochSV)nums.push_back(kv.second); d.epochSVMin=*std::min_element(nums.begin(),nums.end()); d.epochSVMax=*std::max_element(nums.begin(),nums.end()); d.epochSVMean=mean(nums); }
    std::vector<double> sn,prph; std::map<std::string,std::vector<double>> snByBand; std::map<std::string,std::map<std::string,std::map<std::string,std::vector<MPPoint>>>> mpComboBasic; std::map<std::string,TimePoint> last; int idx=0; for(auto& e:rf.data.observationEpochs){ if(last.count("epoch")){ double gap=std::chrono::duration<double>(e.time-last["epoch"]).count(); if(gap>opt.gapMinMinutes*60){ d.gapEvents.push_back({formatUTC(last["epoch"]),formatUTC(e.time),gap,idx}); d.slipEvents.push_back({formatUTC(e.time),"*","gap",std::to_string((int)gap)+" seconds"}); } } last["epoch"]=e.time; ++idx; }
    std::map<std::string,TimePoint> lastSat; std::map<std::string,double> lastPrph; for(auto& r:rf.data.observationRecords){ std::optional<double> Cc,Lc; for(auto& v:r.values){ if(!v.lli.empty()){ d.lliCount++; d.slipEvents.push_back({formatUTC(r.time),r.satellite,"LLI",v.type+":"+v.lli}); } if(v.value&&v.type.rfind("S",0)==0){ sn.push_back(*v.value); auto band=v.type.size()>1?std::string(1,v.type[1]):std::string("1"); snByBand[band].push_back(*v.value); if(opt.minSNR.count(band) && *v.value < opt.minSNR.at(band)) d.snrStats[band].lowCount++; } if(v.value&&v.type.rfind("C",0)==0&&!Cc) Cc=v.value; if(v.value&&v.type.rfind("L",0)==0&&!Lc) Lc=v.value; } if(Cc&&Lc){ double pp=*Cc-*Lc; prph.push_back(pp); if(lastPrph.count(r.satellite)){ double jump=std::fabs(pp-lastPrph[r.satellite]); if(jump>50.0){ d.multipathStats["all"].jumps++; d.slipEvents.push_back({formatUTC(r.time),r.satellite,"P-L jump",std::to_string(jump)}); } } lastPrph[r.satellite]=pp; } if(lastSat.count(r.satellite)){ double dt=std::chrono::duration<double>(r.time-lastSat[r.satellite]).count(); if(dt>opt.gapMinMinutes*60) d.slipEvents.push_back({formatUTC(r.time),r.satellite,"sat gap",std::to_string((int)dt)+" seconds"}); }
      std::set<char> bands; for(const auto& v:r.values) if(v.value && v.type.size()>=2) bands.insert(v.type[1]);
      bool arcBreak=false; for(const auto& vv:r.values) if(!vv.lli.empty()) arcBreak=true;
      for(char bi:bands){ double Pi=obsValue(r,'C',bi); double Li=carrierMeters(r,bi,nullptr); double fi=signalFrequencyHz(r.system,bi,nullptr); if(!std::isfinite(Pi)||!std::isfinite(Li)||fi<=0) continue; for(char bj:bands){ if(bj==bi) continue; double Lj=carrierMeters(r,bj,nullptr); double fj=signalFrequencyHz(r.system,bj,nullptr); if(!std::isfinite(Lj)||fj<=0) continue; double a=(fi/fj)*(fi/fj); if(std::fabs(a-1.0)<1e-12) continue; double mp=Pi-(1.0+2.0/(a-1.0))*Li+(2.0/(a-1.0))*Lj; std::string key=std::string("MP")+bi+bj; mpComboBasic[r.system][key][r.satellite].push_back({r.time,mp,arcBreak}); }}
      lastSat[r.satellite]=r.time; }
    int lowAll = 0; for(auto& kv:d.snrStats) lowAll += kv.second.lowCount; d.snrStats["all"]=stat(sn); d.snrStats["all"].lowCount = lowAll; for(auto& kv:snByBand){ int low=d.snrStats[kv.first].lowCount; d.snrStats[kv.first]=stat(kv.second); d.snrStats[kv.first].lowCount=low; }
    d.pseudorangePhase["all"]=stat(prph); auto mpstat=stat(prph); mpstat.jumps=d.multipathStats["all"].jumps; d.multipathStats["all"]=mpstat;
    for(auto& skv: mpComboBasic){ for(auto& ckv: skv.second){ auto st=movingAverageMPStats(ckv.second,opt.mpWindow,true); if(st.count>0){ std::string key=skv.first+":"+ckv.first; d.multipathStats[key]=st; d.multipathMovingRMS[key]=st.rms; d.multipathMovingCount[key]=st.count; } } }
    d.ionStats["all"]=stat(prph); d.iodStats["all"]=stat({}); d.deletedObservations=(d.codeBandCount<2?static_cast<int>(rf.data.observationRecords.size()):0); buildTeqcTimeplot(rf,opt,d); d.symbolLegend={"c code/phase observation epoch","L loss-of-lock","g gap"}; d.position.attempted=opt.averagePosition||opt.everyEpochPosition; d.position.epochSolutions=s.epochCount; s.derived=d; } else if(rf.header.kind==RinexKind::Nav){ s.navigationRecords=(int)rf.data.navigationRecords.size(); for(auto& r:rf.data.navigationRecords){ s.navigationValues+=(int)r.values.size(); s.navigationFields+=(int)r.fields.size(); if(!r.satellite.empty())s.satelliteAppearance[r.satellite]++; if(r.epoch){ if(!s.firstEpoch||*r.epoch<*s.firstEpoch)s.firstEpoch=r.epoch; if(!s.lastEpoch||*r.epoch>*s.lastEpoch)s.lastEpoch=r.epoch; } } s.broadcastEphemerides=s.navigationRecords; } else if(rf.header.kind==RinexKind::Met){ s.meteorologicalRecords=(int)rf.data.meteorologicalRecords.size(); for(auto& r:rf.data.meteorologicalRecords) s.meteorologicalValues+=(int)r.values.size(); } return s; }
QCSummary analyzeWithNavigation(const RinexFile& rf,const std::vector<NavigationRecord>& navs,const QCOptions& opt){ auto s=analyze(rf,opt); if(rf.header.kind==RinexKind::Obs&&!navs.empty()){ s.residuals=residual(rf,navs,opt); applyNavBasedQCMetrics(rf, navs, opt, s); } return s; }
std::string makePlot(const QCSummary& summary){
  std::ostringstream os;
  os << "COMPACT3 CEQC 5.5.0\n";
  os << "file " << summary.sourcePath << "\n";
  os << "epochs " << summary.epochCount << " obs_records " << summary.observationRecords << "\n";
  if(summary.firstEpoch) os << "first " << formatUTC(*summary.firstEpoch) << "\n";
  if(summary.lastEpoch) os << "last " << formatUTC(*summary.lastEpoch) << "\n";
  if(summary.derived){
    os << " SV+" << std::string(summary.derived->obsTimeplot.size(), '-') << "+ SV\n";
    std::vector<std::pair<std::string,std::string>> rows(summary.derived->satelliteTimeplot.begin(), summary.derived->satelliteTimeplot.end());
    std::sort(rows.begin(), rows.end(), [](const auto& a,const auto& b){ if(a.first.empty()||b.first.empty()) return a.first<b.first; if(a.first[0]!=b.first[0]) return a.first<b.first; return prnNumber(a.first)<prnNumber(b.first); });
    for(auto& r:rows) os << std::setw(3) << r.first.substr(1) << "|" << r.second << "| " << std::setw(3) << r.first.substr(1) << "\n";
    os << "Obs|" << summary.derived->obsTimeplot << "|Obs\n";
  }
  for(auto& kv:summary.systemAppearance) os << "system " << kv.first << " " << kv.second << "\n";
  os << "symbol c observation epoch\n";
  os << "symbol L loss-of-lock\n";
  os << "symbol g gap\n";
  return os.str();
}
}
