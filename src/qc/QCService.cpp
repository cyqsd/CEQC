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
#include <iomanip>
#include <deque>

namespace ceqc::service::qc {
namespace {
constexpr double C = 299792458.0;
constexpr double OMEGA_E = 7.2921151467e-5;
constexpr double PI = 3.141592653589793238462643383279502884;

// Internal observation epochs are read on the observation file time scale.  The
// uploaded RTCM/RINEX test set declares GPS time for OBS, while BDS D1/D2
// ephemerides use BDT and GLONASS FDMA ephemerides use UTC(SU)/UTC-like epochs.
// For broadcast propagation and Toe/Toc matching we must compare on the
// constellation navigation time scale; otherwise BDS is propagated 14 s too late
// and GLONASS about 18 s too late, which produces kilometre-level residuals.
constexpr double BDS_MINUS_GPS_SEC = -14.0;  // BDT = GPST - 14 s.
TimePoint shiftSeconds(TimePoint t, double sec){
  return t + std::chrono::duration_cast<TimePoint::duration>(std::chrono::duration<double>(sec));
}
int gpsMinusUtcSeconds(TimePoint gpsLikeTime){
  struct LeapEntry { TimePoint effective; int offset; };
  static const LeapEntry leaps[] = {
    {makeUTC(1981,7,1,0,0,0.0),1}, {makeUTC(1982,7,1,0,0,0.0),2},
    {makeUTC(1983,7,1,0,0,0.0),3}, {makeUTC(1985,7,1,0,0,0.0),4},
    {makeUTC(1988,1,1,0,0,0.0),5}, {makeUTC(1990,1,1,0,0,0.0),6},
    {makeUTC(1991,1,1,0,0,0.0),7}, {makeUTC(1992,7,1,0,0,0.0),8},
    {makeUTC(1993,7,1,0,0,0.0),9}, {makeUTC(1994,7,1,0,0,0.0),10},
    {makeUTC(1996,1,1,0,0,0.0),11}, {makeUTC(1997,7,1,0,0,0.0),12},
    {makeUTC(1999,1,1,0,0,0.0),13}, {makeUTC(2006,1,1,0,0,0.0),14},
    {makeUTC(2009,1,1,0,0,0.0),15}, {makeUTC(2012,7,1,0,0,0.0),16},
    {makeUTC(2015,7,1,0,0,0.0),17}, {makeUTC(2017,1,1,0,0,0.0),18},
  };
  int out=0;
  for(const auto& e:leaps) if(gpsLikeTime>=e.effective) out=e.offset;
  return out;
}
TimePoint navScaleTimeFromObsTime(TimePoint obsTime, const std::string& system){
  if(system=="C") return shiftSeconds(obsTime, BDS_MINUS_GPS_SEC);
  if(system=="R") return shiftSeconds(obsTime, -static_cast<double>(gpsMinusUtcSeconds(obsTime)));
  return obsTime;
}

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

std::vector<int> histogramLinear(const std::vector<double>& values, int bins, double lo=std::numeric_limits<double>::quiet_NaN(), double hi=std::numeric_limits<double>::quiet_NaN()){
  if(bins<=0) bins=18; if(bins>200) bins=200;
  std::vector<int> h(static_cast<size_t>(bins),0);
  std::vector<double> v; v.reserve(values.size());
  for(double x:values) if(std::isfinite(x)) v.push_back(x);
  if(v.empty()) return h;
  if(!std::isfinite(lo)) lo=*std::min_element(v.begin(),v.end());
  if(!std::isfinite(hi)) hi=*std::max_element(v.begin(),v.end());
  if(!(hi>lo)){ h[0]=static_cast<int>(v.size()); return h; }
  for(double x:v){
    int b=static_cast<int>(std::floor((x-lo)/(hi-lo)*bins));
    if(b<0) b=0; if(b>=bins) b=bins-1;
    h[static_cast<size_t>(b)]++;
  }
  return h;
}
std::string histogramText(const std::vector<int>& h){
  std::ostringstream os; os << "[";
  for(size_t i=0;i<h.size();++i){ if(i) os << ","; os << h[i]; }
  os << "]"; return os.str();
}
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
  if(rf.data.observationRecords.empty()){
    d.timeplot="|"+std::string(width,' ')+"|";
    d.obsTimeplot=std::string(width,'0');
    d.navTimeplot=std::string(width,' ');
    d.positionTimeplot=std::string(width,' ');
    return;
  }
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
  d.navTimeplot=std::string(width,'n');      // OBS-only: no NAV QC was available.
  d.positionTimeplot=std::string(width,'s'); // OBS-only: position QC skipped, not solved.
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
bool orbitUsable(const NavigationRecord& n){
  auto explicitFlag=n.fields.find("OrbitUsable");
  if(explicitFlag!=n.fields.end() && (!std::isfinite(explicitFlag->second.value) || explicitFlag->second.value <= 0.5)) return false;
  if(n.system=="G" || n.system=="J" || n.system=="E" || n.system=="C") {
    auto a=f(n,"SqrtA"), toe=f(n,"Toe"), m0=f(n,"M0"), ecc=f(n,"Eccentricity");
    if(!a || !toe || !m0 || !ecc) return false;
    if(!std::isfinite(*a) || *a < 5000.0 || *a > 7000.0) return false;
    if(!std::isfinite(*toe) || *toe < 0.0 || *toe >= 604800.0) return false;
    if(!std::isfinite(*ecc) || *ecc < 0.0 || *ecc >= 1.0) return false;
  }
  return true;
}
std::optional<NavigationRecord> nearestEph(const std::vector<NavigationRecord>& navs,const std::string& sat,TimePoint t,const QCOptions& opt){
  double best=std::numeric_limits<double>::infinity();
  std::optional<NavigationRecord> out;
  std::string sys=sat.empty()?"":sat.substr(0,1);
  if(opt.noOrbitSystems.count(sys) && opt.noOrbitSystems.at(sys)) return {};
  TimePoint navTime = navScaleTimeFromObsTime(t, sys);
  for(auto& n:navs){
    if(n.satellite!=sat || !n.epoch || n.fields.empty() || !orbitUsable(n)) continue;
    TimePoint cmpTime = navScaleTimeFromObsTime(t, n.system);
    double refdt=std::fabs(std::chrono::duration<double>(cmpTime-*n.epoch).count());
    double dt=refdt;
    if(auto toe=f(n,"Toe")){
      double sow=0.0;
      if(n.system=="E") sow=std::fmod(std::max(0.0,std::chrono::duration<double>(cmpTime-makeUTC(1999,8,22,0,0,0.0)).count()),604800.0);
      else if(n.system=="C") sow=std::fmod(std::max(0.0,std::chrono::duration<double>(cmpTime-makeUTC(2006,1,1,0,0,0.0)).count()),604800.0);
      else sow=std::fmod(std::max(0.0,std::chrono::duration<double>(cmpTime-makeUTC(1980,1,6,0,0,0.0)).count()),604800.0);
      double toeDt=sow-*toe;
      while(toeDt>302400.0) toeDt-=604800.0;
      while(toeDt<-302400.0) toeDt+=604800.0;
      dt=std::min(refdt,std::fabs(toeDt));
    }
    if(dt<best){best=dt; out=n;}
  }
  (void)navTime;
  if(best>6*3600) return {};
  return out;
}
double median(std::vector<double> v){ if(v.empty()) return 0; std::sort(v.begin(),v.end()); auto n=v.size(); if(n%2) return v[n/2]; return 0.5*(v[n/2-1]+v[n/2]); }
double norm3(const std::array<double,3>& a,const std::array<double,3>& b){ double dx=a[0]-b[0],dy=a[1]-b[1],dz=a[2]-b[2]; return std::sqrt(dx*dx+dy*dy+dz*dz); }
struct SatState{ std::array<double,3> xyz{}; double clk=0; bool ok=false; };
double dtSeconds(TimePoint a,TimePoint b){ return std::chrono::duration<double>(a-b).count(); }
std::array<double,3> sagnac(const std::array<double,3>& x,double tau);
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
SatState propagateForPseudorange(const NavigationRecord& e,
                                 TimePoint rxTime,
                                 double pseudorange,
                                 const std::array<double,3>& rxXYZ,
                                 double& tauOut){
  double tau = pseudorange / C;
  SatState ss;
  TimePoint navRxTime = navScaleTimeFromObsTime(rxTime, e.system);
  for(int i=0;i<3;++i){
    auto tx = navRxTime - std::chrono::duration_cast<std::chrono::system_clock::duration>(
      std::chrono::duration<double>(tau));
    ss = prop(e, tx);
    if(!ss.ok) break;
    auto rot = sagnac(ss.xyz, tau);
    tau = norm3(rot, rxXYZ) / C;
  }
  tauOut = tau;
  return ss;
}
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
    // GLONASS FDMA wavelengths depend on the satellite frequency channel.  If
    // no navigation record/header slot is available, do not guess k=0 because
    // that creates false tens-of-metres MP values in OBS-only QC.
    if(!eph) return 0.0;
    int k=0; auto fk=f(*eph,"FrequencyNumber"); if(fk) k=(int)std::lround(*fk);
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
  for(const auto& rec:records){ if(opt.noPositionSystems.count(rec.system) && opt.noPositionSystems.at(rec.system)) continue; if(!firstPseudorange(rec)) continue; if(nearestEph(navs,rec.satellite,rec.time,opt)) count[rec.system]++; }
  if(count["G"]+count["J"]>=4) return "G";
  int best=0; std::string sys;
  for(auto& kv:count){ if(kv.second>best){ best=kv.second; sys=kv.first; } }
  return sys.empty()?"G":sys;
}

bool isReferenceClockSystem(const std::string& sys, const std::string& ref){
  if(ref=="G" && (sys=="G"||sys=="J")) return true;
  return sys==ref;
}

bool systemDisabledForPosition(const QCOptions& opt, const std::string& sys){
  auto it=opt.noPositionSystems.find(sys);
  return it!=opt.noPositionSystems.end() && it->second;
}

struct EpochPositionSolution {
  std::array<double,4> state{};
  double postfitRmsM=std::numeric_limits<double>::infinity();
  double maxResidualM=std::numeric_limits<double>::infinity();
  int usedSVs=0;
  std::map<std::string,std::set<std::string>> usedBySystem;
};

struct SystemScreenResult {
  std::map<std::string,double> centeredRmsBySystem;
  std::map<std::string,int> countBySystem;
  std::map<std::string,bool> excludeFromPosition;
};

SystemScreenResult screenPositionSystems(const RinexFile& rf, const std::vector<NavigationRecord>& navs, const std::array<double,3>& rx, const QCOptions& opt){
  SystemScreenResult out;
  IonoModel iono=ionoModelFromNavs(navs);
  // Raw UBX pseudoranges can carry a large, time-varying receiver clock common
  // mode.  A single median over the whole session is therefore not a valid
  // system-quality screen: it falsely rejects otherwise healthy GPS/QZSS data
  // when the receiver clock drifts by kilometres.  Screen each constellation on
  // epoch-centred residuals instead; this keeps satellite/orbit/code scatter and
  // removes only the common receiver-clock term that SPP will estimate anyway.
  std::map<std::string,std::map<std::string,std::vector<double>>> bySysEpoch;
  for(const auto& rec:rf.data.observationRecords){
    if(systemDisabledForPosition(opt, rec.system)) continue;
    auto pr=firstPseudorange(rec); if(!pr) continue;
    auto eph=nearestEph(navs,rec.satellite,rec.time,opt); if(!eph) continue;
    double tau=0.0;
    auto ss=propagateForPseudorange(*eph, rec.time, *pr, rx, tau); if(!ss.ok) continue;
    auto sx=sagnac(ss.xyz,tau);
    double rho=norm3(sx,rx); if(!std::isfinite(rho) || rho<=1.0) continue;
    double el=elevationRad(rx,sx);
    double res=*pr - (rho + tropoDelayMeters(rx,el) + ionoDelayMeters(rx,sx,rec.time,rec.system,iono) - C*ss.clk);
    if(std::isfinite(res)) bySysEpoch[rec.system][formatUTC(rec.time)].push_back(res);
  }
  for(auto& skv:bySysEpoch){
    std::vector<double> centered;
    int rawCount=0;
    for(auto& ekv:skv.second){
      rawCount += static_cast<int>(ekv.second.size());
      if(ekv.second.empty()) continue;
      // Need at least two satellites to remove a common-mode clock without
      // destroying all information.  With one satellite, keep it out of the
      // screen instead of making it look perfect.
      if(ekv.second.size()<2) continue;
      double med=median(ekv.second);
      for(double v:ekv.second) centered.push_back(v-med);
    }
    out.countBySystem[skv.first]=rawCount;
    if(centered.size()<20) continue;
    auto st=stat(centered);
    double maxAbs=std::max(std::fabs(st.min), std::fabs(st.max));
    out.centeredRmsBySystem[skv.first]=st.rms;
    // GPS/QZSS are retained unless they are wildly broken.  Other systems are
    // admitted to the SPP only after their broadcast propagation/code model has
    // metre-to-tens-of-metres scatter.  This prevents kilometre-level BDS/GLO
    // model errors from being hidden inside receiver clock or ISB estimates.
    double rmsGate = (skv.first=="G" || skv.first=="J") ? 150.0 : 80.0;
    double maxGate = (skv.first=="G" || skv.first=="J") ? 800.0 : 500.0;
    if(st.rms > rmsGate || maxAbs > maxGate) out.excludeFromPosition[skv.first]=true;
  }
  return out;
}

std::optional<EpochPositionSolution> solveEpochPosition(const std::vector<ObservationRecord>& records,
                                                       const std::vector<NavigationRecord>& navs,
                                                       const std::array<double,3>& seed,
                                                       const QCOptions& opt,
                                                       const IonoModel& iono){
  int eligibleSVs=0;
  for(const auto& rec:records){
    if(opt.noPositionSystems.count(rec.system) && opt.noPositionSystems.at(rec.system)) continue;
    if(!firstPseudorange(rec)) continue;
    if(!nearestEph(navs,rec.satellite,rec.time,opt)) continue;
    ++eligibleSVs;
  }
  if(eligibleSVs < std::max(4, opt.minSVs)) return {};
  std::string refSys=sppReferenceSystem(records,navs,opt);
  std::map<std::string,int> biasIndex;
  for(const auto& rec:records){
    if(opt.noPositionSystems.count(rec.system) && opt.noPositionSystems.at(rec.system)) continue;
    if(!firstPseudorange(rec)) continue;
    if(!nearestEph(navs,rec.satellite,rec.time,opt)) continue;
    if(isReferenceClockSystem(rec.system,refSys)) continue;
    if(!biasIndex.count(rec.system)) biasIndex[rec.system]=4+(int)biasIndex.size();
  }
  const int nState=4+(int)biasIndex.size();
  if(nState>8) return {};
  std::vector<double> state(nState,0.0);
  state[0]=seed[0]; state[1]=seed[1]; state[2]=seed[2]; state[3]=0.0;
  {
    std::vector<double> clkGuess;
    std::array<double,3> rx{state[0],state[1],state[2]};
    for(const auto& rec:records){
      if(opt.noPositionSystems.count(rec.system) && opt.noPositionSystems.at(rec.system)) continue;
      auto pr=firstPseudorange(rec); if(!pr) continue;
      auto eph=nearestEph(navs,rec.satellite,rec.time,opt); if(!eph) continue;
      double tau=0.0;
      auto ss=propagateForPseudorange(*eph, rec.time, *pr, rx, tau); if(!ss.ok) continue;
      auto sx=sagnac(ss.xyz,tau);
      double rho=norm3(sx,rx); if(rho<=1.0 || !std::isfinite(rho)) continue;
      double el=elevationRad(rx,sx);
      double trop=tropoDelayMeters(rx,el);
      double ion=ionoDelayMeters(rx,sx,rec.time,rec.system,iono);
      clkGuess.push_back(*pr - (rho + trop + ion - C*ss.clk));
    }
    if(clkGuess.size() >= 4) state[3]=median(clkGuess);
  }
  std::vector<char> active(records.size(),1);
  for(int pass=0; pass<2; ++pass){
    for(int iter=0; iter<10; ++iter){
      std::vector<std::vector<double>> N(nState,std::vector<double>(nState,0.0));
      std::vector<double> u(nState,0.0);
      int used=0;
      for(size_t ri=0; ri<records.size(); ++ri){
        if(!active[ri]) continue;
        const auto& rec=records[ri];
        if(opt.noPositionSystems.count(rec.system) && opt.noPositionSystems.at(rec.system)) continue;
        auto pr=firstPseudorange(rec); if(!pr) continue;
        auto eph=nearestEph(navs,rec.satellite,rec.time,opt); if(!eph) continue;
        std::array<double,3> rx{state[0],state[1],state[2]};
        double tau=0.0;
        auto ss=propagateForPseudorange(*eph, rec.time, *pr, rx, tau); if(!ss.ok) continue;
        auto sx=sagnac(ss.xyz,tau);
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
        const auto& rec=records[ri]; if(systemDisabledForPosition(opt, rec.system)) continue; auto pr=firstPseudorange(rec); if(!pr) continue; auto eph=nearestEph(navs,rec.satellite,rec.time,opt); if(!eph) continue;
        std::array<double,3> rx{state[0],state[1],state[2]};
        double tau=0.0; auto ss=propagateForPseudorange(*eph, rec.time, *pr, rx, tau); if(!ss.ok) continue; auto sx=sagnac(ss.xyz,tau);
        double rho=norm3(sx,rx); if(rho<=1.0) continue; double el=elevationRad(rx,sx); double bias=0.0; auto bit=biasIndex.find(rec.system); if(bit!=biasIndex.end()) bias=state[bit->second];
        double v=*pr-(rho+tropoDelayMeters(rx,el)+ionoDelayMeters(rx,sx,rec.time,rec.system,iono)-C*ss.clk+state[3]+bias); res.push_back(v); idx.push_back(ri);
      }
      if(res.size()>nState+2){ double med=median(res); std::vector<double> dev; for(double v:res) dev.push_back(std::fabs(v-med)); double mad=median(dev); double gate=std::max(60.0, 6.0*1.4826*mad); for(size_t i=0;i<res.size();++i) if(std::fabs(res[i]-med)>gate) active[idx[i]]=0; }
    }
  }

  EpochPositionSolution sol;
  sol.state={state[0],state[1],state[2],state[3]};
  std::vector<double> finalResiduals;
  for(size_t ri=0; ri<records.size(); ++ri){
    if(!active[ri]) continue;
    const auto& rec=records[ri];
    if(systemDisabledForPosition(opt, rec.system)) continue;
    auto pr=firstPseudorange(rec); if(!pr) continue;
    auto eph=nearestEph(navs,rec.satellite,rec.time,opt); if(!eph) continue;
    std::array<double,3> rx{state[0],state[1],state[2]};
    double tau=0.0;
    auto ss=propagateForPseudorange(*eph, rec.time, *pr, rx, tau); if(!ss.ok) continue;
    auto sx=sagnac(ss.xyz,tau);
    double rho=norm3(sx,rx); if(rho<=1.0 || !std::isfinite(rho)) continue;
    double el=elevationRad(rx,sx);
    if(std::isfinite(el) && el < -0.05) continue;
    double bias=0.0; auto bit=biasIndex.find(rec.system); if(bit!=biasIndex.end()) bias=state[bit->second];
    double v=*pr-(rho+tropoDelayMeters(rx,el)+ionoDelayMeters(rx,sx,rec.time,rec.system,iono)-C*ss.clk+state[3]+bias);
    if(!std::isfinite(v)) continue;
    finalResiduals.push_back(v);
    sol.usedBySystem[rec.system].insert(rec.satellite);
  }
  sol.usedSVs=static_cast<int>(finalResiduals.size());
  if(sol.usedSVs<nState) return {};
  auto rst=stat(finalResiduals);
  sol.postfitRmsM=rst.rms;
  sol.maxResidualM=std::max(std::fabs(rst.min), std::fabs(rst.max));
  return sol;
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
      std::array<double,3> xyz{sol->state[0],sol->state[1],sol->state[2]};
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
QCMetricStats movingAverageMPStats(const std::map<std::string,std::vector<MPPoint>>& bySat, int window, bool excludeCurrent=false, double sigmaGate=4.0){
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
    double k = std::isfinite(sigmaGate) && sigmaGate > 0.0 ? sigmaGate : 4.0;
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
  std::vector<int> navCandidateBins(width,0);
  std::vector<int> navMissingBins(width,0);
  std::vector<int> positionEligibleBins(width,0);
  std::set<std::string> observedSats;
  std::set<int> unhealthyGPS;
  std::map<std::string,bool> hasUsableNav;
  for(const auto& n:navs){
    if(!n.satellite.empty() && orbitUsable(n)) hasUsableNav[n.satellite]=true;
    if(n.satellite.size()>1 && n.satellite[0]=='G'){ auto h=f(n,"SVHealth"); if(h && std::lround(*h)!=0) unhealthyGPS.insert(prnNumber(n.satellite)); }
  }
  for(const auto& r:rf.data.observationRecords){
    observedSats.insert(r.satellite);
    obsBySatBin[r.satellite].resize(width,0); lliBySatBin[r.satellite].resize(width,0); incompleteBySatBin[r.satellite].resize(width,0); mpBySatBin[r.satellite].resize(width,0.0);
    int b=(int)std::floor(std::chrono::duration<double>(r.time-first).count()/span*(width-1)+1e-9); if(b<0)b=0; if(b>=width)b=width-1;
    obsBySatBin[r.satellite][b]++;
    if(firstPseudorange(r)){
      if(nearestEph(navs,r.satellite,r.time,opt)) {
        navCandidateBins[b]++;
        if(!(opt.noPositionSystems.count(r.system) && opt.noPositionSystems.at(r.system))) positionEligibleBins[b]++;
      } else {
        navMissingBins[b]++;
      }
    }
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
    else if(!hasUsableNav.count(r.satellite)) sym = lli ? 'L' : 'N';
    else if(lli) sym='L';
    else if(hasP1 && (!hasL1 || !hasP2 || !hasL2)) sym='2';
    rows[r.satellite][b]=mergeTeqcSymbol(rows[r.satellite][b], sym);
  }
  d.satelliteTimeplot=rows;
  std::vector<int> obsBins(width,0);
  for(const auto& kv:obsBySatBin) {
    for(int i=0;i<width && i<(int)kv.second.size();++i) if(kv.second[i]>0) obsBins[i]++;
  }
  // Build a composite QC timeplot from the per-satellite states.  It must not
  // collapse to all "9" merely because the position solution succeeded; gaps,
  // missing navigation, incomplete dual-frequency observations and low elevation
  // bins are preserved with teqc-like priority.
  std::string obsRow; obsRow.reserve(width);
  auto pri=[&](char x){
    switch(x){ case 'L': return 90; case '2': return 70; case 'N': return 0; case '.': return 60; case '_': return 55; case '-': return 50; case '~': return 40; case '+': return 30; case 'c': return 20; default: return 0; }
  };
  for(int b=0;b<width;++b){
    if(b >= (int)obsBins.size() || obsBins[b] == 0){ obsRow.push_back('g'); continue; }
    char best='c'; int bestp=0;
    for(const auto& kv:rows){ if(b < (int)kv.second.size()){ char c=kv.second[b]; int p=pri(c); if(p>bestp){ bestp=p; best=c; } } }
    obsRow.push_back(best);
  }
  std::string navRow; navRow.reserve(width);
  std::string posRow; posRow.reserve(width);
  int minForPos=std::max(4,opt.minSVs);
  for(int b=0;b<width;++b){
    if(b >= (int)obsBins.size() || obsBins[b] == 0){ navRow.push_back('g'); posRow.push_back('g'); continue; }
    if(navCandidateBins[b] > 0) navRow.push_back('+');
    else if(navMissingBins[b] > 0) navRow.push_back('N');
    else navRow.push_back('n');
    if(positionEligibleBins[b] >= minForPos) posRow.push_back('P'); else posRow.push_back('s');
  }
  d.obsBinCounts=obsBins; d.obsTimeplot=obsRow; d.navTimeplot=navRow; d.positionTimeplot=posRow; d.timeplot="|"+d.obsTimeplot+"|";
  int eventBins=0; for(const auto& kv:mpBySatBin) for(double v:kv.second) if(std::isfinite(v) && std::fabs(v)>1e-6) ++eventBins;
  d.msecMpEventBins = std::max(0, eventBins/30);
}
void applyNavBasedQCMetrics(const RinexFile& rf, const std::vector<NavigationRecord>& navs, const QCOptions& opt, QCSummary& s){
  if(!s.derived) s.derived=QCDerivedSummary{};
  QCDerivedSummary& d=*s.derived;
  const double obsOnlyMp1Meters = d.mp1Meters;
  const double obsOnlyMp2Meters = d.mp2Meters;
  const auto obsOnlyMultipathStats = d.multipathStats;
  const auto obsOnlyMultipathMovingRMS = d.multipathMovingRMS;
  const auto obsOnlyMultipathMovingCount = d.multipathMovingCount;
  QCOptions posOpt=opt;
  auto recXYZ=approxXYZ(rf.header);
  if(!recXYZ || std::sqrt((*recXYZ)[0]*(*recXYZ)[0]+(*recXYZ)[1]*(*recXYZ)[1]+(*recXYZ)[2]*(*recXYZ)[2]) < 1.0e6){
    // Bootstrap missing/zero station coordinates conservatively.  Do not let
    // unverified BDS/GLO models pull the seed hundreds or thousands of metres.
    QCOptions seedOpt=opt;
    seedOpt.noPositionSystems["R"]=true;
    seedOpt.noPositionSystems["C"]=true;
    recXYZ=estimateApproxPositionInternal(rf, navs, seedOpt);
  }
  if(!recXYZ){
    // Still build observation/navigation timeplots when NAV is present; only the
    // elevation/position rows are degraded.  This matches teqc's behaviour of
    // not hiding OBS/NAV availability just because station coordinates are absent.
    d.position.candidateEpochs=s.epochCount;
    d.position.warnings.push_back("no station coordinates; position row skipped");
    buildTeqcTimeplotWithNav(rf, navs, posOpt, d);
    return;
  }

  auto sysScreen=screenPositionSystems(rf,navs,*recXYZ,opt);
  for(const auto& kv:sysScreen.excludeFromPosition){
    if(kv.second){
      posOpt.noPositionSystems[kv.first]=true;
      std::ostringstream w;
      w << "position-system-screen: " << kv.first << " excluded from SPP";
      auto rit=sysScreen.centeredRmsBySystem.find(kv.first);
      auto cit=sysScreen.countBySystem.find(kv.first);
      if(rit!=sysScreen.centeredRmsBySystem.end()) w << " epoch_centered_rms=" << std::fixed << std::setprecision(2) << rit->second << " m";
      if(cit!=sysScreen.countBySystem.end()) w << " samples=" << cit->second;
      d.position.warnings.push_back(w.str());
    }
  }

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
      auto pr=firstPseudorange(rec);
      if(pr){
        double tau=0.0;
        auto st=propagateForPseudorange(*ephPtr,rec.time,*pr,*recXYZ,tau);
        if(st.ok){ auto sx=sagnac(st.xyz,tau); el=elevationRad(*recXYZ,sx); }
      }
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
  auto mp1=movingAverageMPStats(mp1BySat,opt.mpWindow,true,opt.mpSigmas), mp2=movingAverageMPStats(mp2BySat,opt.mpWindow,false,opt.mpSigmas);
  if(mp1.count>0){ d.mp1Meters=mp1.rms; d.multipathMovingCount["MP1"]=mp1.count; }
  if(mp2.count>0){ d.mp2Meters=mp2.rms; d.multipathMovingCount["MP2"]=mp2.count; }
  for(auto& skv: mpComboBySystem){
    for(auto& ckv: skv.second){
      auto st = movingAverageMPStats(ckv.second, opt.mpWindow, true, opt.mpSigmas);
      if(st.count<=0) continue;
      std::string key = skv.first + ":" + ckv.first;
      d.multipathStats[key]=st;
      d.multipathMovingRMS[key]=st.rms;
      d.multipathMovingCount[key]=st.count;
    }
  }
  // MP is an observation-domain code/carrier combination.  NAV is used above
  // for elevation/rise-set/position diagnostics only; do not let availability of
  // broadcast ephemerides change the default MP RMS reported by +qcq.
  d.mp1Meters = obsOnlyMp1Meters;
  d.mp2Meters = obsOnlyMp2Meters;
  d.multipathStats = obsOnlyMultipathStats;
  d.multipathMovingRMS = obsOnlyMultipathMovingRMS;
  d.multipathMovingCount = obsOnlyMultipathMovingCount;
  if(opt.multipath){
    auto checkMp=[&](const std::string& name, double rmsMeters){
      std::string key=name;
      if(key.rfind("MP",0)==0) key=key.substr(2);
      auto it=opt.mpRMSCM.find(key);
      if(it!=opt.mpRMSCM.end() && std::isfinite(rmsMeters) && rmsMeters*100.0 > it->second){
        std::ostringstream w; w << "multipath-rms: " << name << " rms=" << std::fixed << std::setprecision(2) << rmsMeters*100.0 << " cm exceeds " << it->second << " cm";
        d.thresholdWarnings.push_back(w.str());
      }
    };
    if(d.mp1Meters>0) checkMp("MP12", d.mp1Meters);
    if(d.mp2Meters>0) checkMp("MP21", d.mp2Meters);
    for(const auto& kv:d.multipathStats){
      auto pos=kv.first.find(":MP");
      if(pos!=std::string::npos) checkMp(kv.first.substr(pos+1), kv.second.rms);
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

  // Actual epoch-by-epoch single point position estimate.  Numeric convergence is
  // not enough: CEQC accepts an epoch only after post-fit residual, height and
  // continuity gates pass.  Systems rejected by the residual screen above are
  // not allowed to contaminate the mixed-constellation solution.
  d.position.hasApprox=true;
  d.position.approxXYZ[0]=(*recXYZ)[0]; d.position.approxXYZ[1]=(*recXYZ)[1]; d.position.approxXYZ[2]=(*recXYZ)[2];
  d.position.residualRmsGateM=50.0;
  d.position.residualMaxGateM=300.0;
  d.position.jumpGateM=std::min(opt.positionJumpM, 250.0);
  std::map<std::string,std::vector<ObservationRecord>> byEpoch;
  for(const auto& r:rf.data.observationRecords) byEpoch[formatUTC(r.time)].push_back(r);
  std::vector<std::array<double,4>> acceptedSols;
  d.position.candidateEpochs=static_cast<int>(byEpoch.size());
  std::optional<std::array<double,3>> previousAccepted;
  std::map<std::string,std::string> epochStatus;
  for(const auto& kv:byEpoch){
    int eligible=0;
    for(const auto& rec:kv.second){
      if(systemDisabledForPosition(posOpt, rec.system)) continue;
      if(firstPseudorange(rec) && nearestEph(navs,rec.satellite,rec.time,posOpt)) ++eligible;
    }
    if(eligible < std::max(4,posOpt.minSVs)){
      d.position.skippedInsufficientSVs++;
      epochStatus[kv.first]="INSUFFICIENT_SVS";
      if(opt.everyEpochPosition){ QCEpochPosition ep; ep.time=kv.first; ep.usedSVs=eligible; ep.status="INSUFFICIENT_SVS"; d.epochPositions.push_back(ep); }
      continue;
    }
    auto sol=solveEpochPosition(kv.second, navs, *recXYZ, posOpt, iono);
    if(!sol){
      d.position.rejectedBadGeometry++;
      d.position.rejectedByStatus["NO_SOLUTION"]++;
      epochStatus[kv.first]="NO_SOLUTION";
      if(opt.everyEpochPosition){ QCEpochPosition ep; ep.time=kv.first; ep.usedSVs=eligible; ep.status="NO_SOLUTION"; d.epochPositions.push_back(ep); }
      continue;
    }
    d.position.epochNumericSolutions++;
    auto llh=llhFromXYZ({sol->state[0],sol->state[1],sol->state[2]});
    std::string status="OK";
    if(!std::isfinite(sol->postfitRmsM) || sol->postfitRmsM>d.position.residualRmsGateM || sol->maxResidualM>d.position.residualMaxGateM){
      status="REJECT_BAD_RESIDUAL";
      d.position.rejectedBadResidual++;
    } else if(!std::isfinite(llh[2]) || llh[2]<opt.positionHMinM || llh[2]>opt.positionHMaxM){
      status="REJECT_BAD_HEIGHT";
      d.position.rejectedBadHeight++;
    } else if(previousAccepted){
      std::array<double,3> xyz{sol->state[0],sol->state[1],sol->state[2]};
      double jump=norm3(xyz,*previousAccepted);
      if(jump>d.position.jumpGateM){
        status="REJECT_BAD_JUMP";
        d.position.rejectedBadJump++;
      }
    }

    if(status=="OK"){
      acceptedSols.push_back(sol->state);
      previousAccepted=std::array<double,3>{sol->state[0],sol->state[1],sol->state[2]};
      for(const auto& ukv:sol->usedBySystem) d.position.usedSVsBySystem[ukv.first]=std::max(d.position.usedSVsBySystem[ukv.first], (int)ukv.second.size());
    } else {
      d.position.rejectedByStatus[status]++;
    }
    epochStatus[kv.first]=status;
    if(opt.everyEpochPosition){
      QCEpochPosition ep;
      ep.time=kv.first; ep.x=sol->state[0]; ep.y=sol->state[1]; ep.z=sol->state[2];
      ep.latDeg=llh[0]*180.0/PI; ep.lonDeg=llh[1]*180.0/PI; ep.heightM=llh[2]; ep.clockBiasM=sol->state[3];
      ep.postfitRmsM=sol->postfitRmsM; ep.maxResidualM=sol->maxResidualM; ep.usedSVs=sol->usedSVs; ep.status=status;
      d.epochPositions.push_back(ep);
    }
  }
  d.position.attempted=true; d.position.skippedNoNavigation=false; d.position.epochSolutions=static_cast<int>(acceptedSols.size());
  if(!acceptedSols.empty()){
    d.position.averageXYZ[0]=d.position.averageXYZ[1]=d.position.averageXYZ[2]=0.0;
    for(auto& sol:acceptedSols){ d.position.averageXYZ[0]+=sol[0]; d.position.averageXYZ[1]+=sol[1]; d.position.averageXYZ[2]+=sol[2]; }
    d.position.averageXYZ[0]/=acceptedSols.size(); d.position.averageXYZ[1]/=acceptedSols.size(); d.position.averageXYZ[2]/=acceptedSols.size();
  }
  buildTeqcTimeplotWithNav(rf, navs, posOpt, d);
  // Replace the eligibility-only position row with the quality-gated result.
  if(!epochStatus.empty() && !rf.data.observationRecords.empty()){
    int width=std::max(10,opt.width); if(width>120) width=120;
    auto first=rf.data.observationRecords.front().time, last=first;
    for(const auto& r:rf.data.observationRecords){ if(r.time<first) first=r.time; if(r.time>last) last=r.time; }
    double span=std::max(1.0,std::chrono::duration<double>(last-first).count());
    std::string gated(width,'s');
    for(const auto& kv:byEpoch){
      if(kv.second.empty()) continue;
      int b=(int)std::floor(std::chrono::duration<double>(kv.second.front().time-first).count()/span*(width-1)+1e-9);
      if(b<0)b=0; if(b>=width)b=width-1;
      auto sit=epochStatus.find(kv.first);
      if(sit==epochStatus.end()) continue;
      if(sit->second=="OK") gated[b]='P';
      else if(sit->second=="INSUFFICIENT_SVS") gated[b]='s';
      else gated[b]='x';
    }
    for(int i=0;i<width && i<(int)d.positionTimeplot.size();++i){ if(d.positionTimeplot[i]=='g') gated[i]='g'; }
    d.positionTimeplot=gated;
  }
  for(auto& ev : d.riseSetEvents){ auto it=d.satelliteMaxElevationDeg.find(ev.satellite); if(it!=d.satelliteMaxElevationDeg.end()){ ev.maxElevationDeg=it->second; ev.hasEphemeris=true; } }
}


ResidualStats residual(const RinexFile& rf,const std::vector<NavigationRecord>& navs,const QCOptions& opt){
  ResidualStats st;
  auto rx=approxXYZ(rf.header);
  if(!rx){
    st.skippedNoStation=(int)rf.data.observationRecords.size();
    st.warnings.push_back("no APPROX POSITION XYZ in observation header; residual QC skipped");
    return st;
  }

  struct Row { std::string epoch, sat, sys; double raw=0, epochCentered=0, biasRemoved=0; };
  std::vector<Row> rows;
  std::map<std::string,std::vector<double>> byEpoch;
  std::map<std::string,std::vector<double>> rawBySystem;

  for(auto& rec:rf.data.observationRecords){
    st.candidateObservations++;
    auto pr=firstPseudorange(rec);
    if(!pr){ st.skippedNoPseudorange++; st.skippedNoPseudorangeBySystem[rec.system]++; continue; }
    auto eph=nearestEph(navs,rec.satellite,rec.time,opt);
    if(!eph){ st.skippedNoEphemeris++; st.skippedNoEphemerisBySystem[rec.system]++; continue; }
    double tau=0.0;
    SatState ss = propagateForPseudorange(*eph, rec.time, *pr, *rx, tau);
    if(!ss.ok){ st.skippedNoEphemeris++; st.skippedNoEphemerisBySystem[rec.system]++; continue; }
    auto rot=sagnac(ss.xyz,tau);
    double rho=norm3(rot,*rx);
    // Raw code-minus-range residual.  Corrections intentionally remain minimal here;
    // the bias-removed series below is reported separately so callers can see both
    // the absolute broadcast/code behaviour and the de-biased QC scatter.
    double res=*pr-rho+C*ss.clk;
    Row row{formatUTC(rec.time),rec.satellite,rec.system,res,0,0};
    rows.push_back(row);
    byEpoch[row.epoch].push_back(res);
    rawBySystem[row.sys].push_back(res);
    st.evaluated++; st.bySystem[rec.system]++;
    if(rec.system=="G") st.gps++; else if(rec.system=="R") st.glonass++; else if(rec.system=="E") st.galileo++; else if(rec.system=="C") st.beidou++;
  }

  std::vector<double> rawValues;
  rawValues.reserve(rows.size());
  for(const auto& r:rows) rawValues.push_back(r.raw);
  auto rawStat=stat(rawValues);
  st.rawMeanMeters=rawStat.mean; st.rawRmsMeters=rawStat.rms; st.rawMaxAbsMeters=std::max(std::fabs(rawStat.min),std::fabs(rawStat.max));
  for(auto& kv:rawBySystem) st.rawBySystem[kv.first]=stat(kv.second);

  std::map<std::string,double> epochMed;
  for(auto& kv:byEpoch) epochMed[kv.first]=median(kv.second);
  std::map<std::string,std::vector<double>> bySat;
  std::vector<double> epochCenteredValues;
  std::map<std::string,std::vector<double>> epochCenteredBySystem;
  for(auto& r:rows){
    r.epochCentered=r.raw-epochMed[r.epoch];
    bySat[r.sat].push_back(r.epochCentered);
    epochCenteredValues.push_back(r.epochCentered);
    epochCenteredBySystem[r.sys].push_back(r.epochCentered);
  }
  auto estat=stat(epochCenteredValues);
  st.meanMeters=estat.mean;
  st.rmsMeters=estat.rms;
  st.maxAbsMeters=std::max(std::fabs(estat.min),std::fabs(estat.max));
  for(auto& kv:epochCenteredBySystem) st.biasRemovedBySystem[kv.first]=stat(kv.second);

  // Diagnostic-only satellite-bias removal.  Do not use these values as the
  // primary QC residual because they can hide wrong constellation models.
  std::map<std::string,double> satMed;
  for(auto& kv:bySat) if(kv.second.size()>=4) satMed[kv.first]=median(kv.second);
  std::vector<double> satBiasRemoved;
  std::map<std::string,std::vector<double>> satBiasRemovedBySystem;
  for(auto& r:rows){
    double v=r.epochCentered;
    if(satMed.count(r.sat)) v-=satMed[r.sat];
    r.biasRemoved=v;
    satBiasRemoved.push_back(v);
    satBiasRemovedBySystem[r.sys].push_back(v);
  }
  if(!satMed.empty()) st.warnings.push_back("diagnostic-only: satellite-specific residual bias estimated for "+std::to_string(satMed.size())+" satellites; not used for primary QC");
  auto sb=stat(satBiasRemoved);
  st.satBiasRemovedMeanMeters=sb.mean;
  st.satBiasRemovedRmsMeters=sb.rms;
  st.satBiasRemovedMaxAbsMeters=std::max(std::fabs(sb.min),std::fabs(sb.max));
  for(auto& kv:satBiasRemovedBySystem) st.satBiasRemovedBySystem[kv.first]=stat(kv.second);
  return st;
}
} // namespace
std::optional<std::array<double,3>> estimateApproxPosition(const RinexFile& rf, const std::vector<NavigationRecord>& navs, const QCOptions& opt){ return estimateApproxPositionInternal(rf, navs, opt); }
QCSummary analyze(const RinexFile& rf,const QCOptions& opt){ QCSummary s; s.sourcePath=rf.path; s.version=rf.header.version; s.kind=rf.header.kind; s.headerIssues=ceqc::service::rinex::validate(rf);
  auto headerValue=[&](const std::string& label)->std::string{ auto it=rf.header.byLabel.find(label); if(it==rf.header.byLabel.end()||it->second.empty()) return {}; auto v=rf.header.lines[it->second.front()].value; auto a=v.find_first_not_of(' '); auto b=v.find_last_not_of(' '); return a==std::string::npos?std::string{}:v.substr(a,b-a+1); };
  auto fixedField=[&](const std::string& label,size_t off,size_t len)->std::string{ auto it=rf.header.byLabel.find(label); if(it==rf.header.byLabel.end()||it->second.empty()) return {}; auto v=rf.header.lines[it->second.front()].value; if(off>=v.size()) return {}; auto x=v.substr(off,std::min(len,v.size()-off)); auto a=x.find_first_not_of(' '); auto b=x.find_last_not_of(' '); return a==std::string::npos?std::string{}:x.substr(a,b-a+1); };
  s.markerName=headerValue("MARKER NAME"); s.markerNumber=headerValue("MARKER NUMBER"); s.receiverNumber=fixedField("REC # / TYPE / VERS",0,20); s.receiverType=fixedField("REC # / TYPE / VERS",20,20); s.receiverVersion=fixedField("REC # / TYPE / VERS",40,20); s.antennaNumber=fixedField("ANT # / TYPE",0,20); s.antennaType=fixedField("ANT # / TYPE",20,20); if(rf.rtcm3) s.rtcm3=*rf.rtcm3; if(rf.ubx) s.ubx=*rf.ubx; if(rf.header.kind==RinexKind::Obs){ s.epochCount=(int)rf.data.observationEpochs.size(); s.observationRecords=(int)rf.data.observationRecords.size(); std::map<std::string,int> epochSV; std::map<std::string,TimePoint> satFirst, satLast; std::map<std::string,int> satObsCount; std::map<std::string,std::map<std::string,std::vector<double>>> svCodeValues; int completeRecords=0, partialRecords=0, yCodeObservations=0; for(auto& r:rf.data.observationRecords){ s.satelliteAppearance[r.satellite]++; s.systemAppearance[r.system]++; epochSV[formatUTC(r.time)]++; if(!satFirst.count(r.satellite)||r.time<satFirst[r.satellite]) satFirst[r.satellite]=r.time; if(!satLast.count(r.satellite)||r.time>satLast[r.satellite]) satLast[r.satellite]=r.time; satObsCount[r.satellite]++; int missingInRecord=0; for(auto& v:r.values){ if(v.value){ s.observationValues++; if((!v.type.empty() && (v.type[0]=='C'||v.type[0]=='P'))) svCodeValues[r.satellite][v.type].push_back(*v.value); if(r.system=="G" && (v.type=="P1"||v.type=="P2"||v.type=="C1W"||v.type=="C2W"||v.type=="L1W"||v.type=="L2W"||v.type=="S1W"||v.type=="S2W")) yCodeObservations++; } else { s.missingObservations++; missingInRecord++; } } if(missingInRecord==0) completeRecords++; else partialRecords++; if(!s.firstEpoch||r.time<*s.firstEpoch)s.firstEpoch=r.time; if(!s.lastEpoch||r.time>*s.lastEpoch)s.lastEpoch=r.time; }
    if(s.firstEpoch&&s.lastEpoch&&s.epochCount>1) s.estimatedIntervalS=std::chrono::duration<double>(*s.lastEpoch-*s.firstEpoch).count()/(s.epochCount-1); QCDerivedSummary d; d.ionEnabled=opt.ion; d.iodEnabled=opt.iod; d.multipathEnabled=opt.multipath; d.snrEnabled=opt.snr; d.lliEnabled=opt.lli; d.pseudorangePhaseEnabled=opt.pseudorangePhase; d.clockSlipsEnabled=opt.clockSlips; d.plotEnabled=opt.plot; d.symbolCodesEnabled=opt.symbolCodes; d.allSymbolsEnabled=opt.allSymbols; d.riseSetEnabled=opt.riseSet; d.ssvEnabled=opt.ssv; d.svprEnabled=opt.svpr; d.dataIndicatorsEnabled=opt.dataIndicators; d.ceqcExtensionEnabled=opt.ceqcExtension; d.yCodeEnabled=opt.yCode; d.everyEpochXYZ=opt.everyEpochXYZ; d.everyEpochGeodetic=opt.everyEpochGeodetic; d.everyEpochDecimal=opt.everyEpochDecimal; d.minSVsUsed=opt.minSVs; if(opt.ion) d.optionsActive.push_back("+ion"); if(opt.iod) d.optionsActive.push_back("+iod"); if(opt.multipath) d.optionsActive.push_back("+mp"); if(opt.snr) d.optionsActive.push_back("+sn"); if(opt.lli) d.optionsActive.push_back("+lli"); if(opt.pseudorangePhase) d.optionsActive.push_back("+pl"); if(opt.clockSlips) d.optionsActive.push_back("+cl"); if(opt.everyEpochPosition) d.optionsActive.push_back("+eep"); if(opt.positionOnly) d.optionsActive.push_back("+pos"); if(opt.symbolCodes) d.optionsActive.push_back("+sym"); if(opt.allSymbols) d.optionsActive.push_back("++sym"); if(!opt.riseSet) d.optionsActive.push_back("-rs"); if(opt.ssv) d.optionsActive.push_back("+ssv"); if(opt.svpr) d.optionsActive.push_back("+svpr"); if(opt.mpRaw) d.optionsActive.push_back("+mp_raw"); if(opt.yCode) d.optionsActive.push_back("+Y-code"); if(opt.dataIndicators) d.optionsActive.push_back("+data"); if(opt.ceqcExtension) d.optionsActive.push_back("+ceqc_ext");
    d.dataCompleteness.completeRecords=completeRecords; d.dataCompleteness.partialRecords=partialRecords; d.dataCompleteness.missingValues=s.missingObservations; d.dataCompleteness.yCodeObservations=yCodeObservations;
    for(const auto& kv:satFirst){ QCRiseSetEvent ev; ev.satellite=kv.first; ev.first=formatUTC(kv.second); ev.last=formatUTC(satLast[kv.first]); ev.durationHours=std::max(0.0,std::chrono::duration<double>(satLast[kv.first]-kv.second).count()/3600.0); ev.obsCount=satObsCount[kv.first]; d.riseSetEvents.push_back(ev); }
    std::sort(d.riseSetEvents.begin(), d.riseSetEvents.end(), [](const auto& a,const auto& b){ if(a.satellite.empty()||b.satellite.empty()) return a.satellite<b.satellite; if(a.satellite[0]!=b.satellite[0]) return a.satellite<b.satellite; return prnNumber(a.satellite)<prnNumber(b.satellite); });
    for(const auto& satkv:svCodeValues){ for(const auto& codekv:satkv.second){ d.svPseudorangeStats[satkv.first][codekv.first]=stat(codekv.second); } }
    
    {
      std::set<char> codeBands;
      for (const auto& r2 : rf.data.observationRecords) for (const auto& v2 : r2.values) if (v2.value && !v2.type.empty() && v2.type[0]=='C' && v2.type.size()>=2) codeBands.insert(v2.type[1]);
      d.codeBandCount = static_cast<int>(codeBands.size());
    }
    if(!epochSV.empty()){ std::vector<int> nums; for(auto& kv:epochSV)nums.push_back(kv.second); d.epochSVMin=*std::min_element(nums.begin(),nums.end()); d.epochSVMax=*std::max_element(nums.begin(),nums.end()); d.epochSVMean=mean(nums); }
    std::vector<double> sn,prph,ionVals,mpVals; std::map<std::string,std::vector<double>> snByBand; std::map<std::string,std::vector<double>> snBySystemCode; std::map<std::string,std::map<std::string,std::map<std::string,std::vector<MPPoint>>>> mpComboBasic; std::map<std::string,TimePoint> last; int idx=0; for(auto& e:rf.data.observationEpochs){ if(last.count("epoch")){ double gap=std::chrono::duration<double>(e.time-last["epoch"]).count(); if(gap>opt.gapMinMinutes*60){ d.gapEvents.push_back({formatUTC(last["epoch"]),formatUTC(e.time),gap,idx}); if(opt.clockSlips) d.slipEvents.push_back({formatUTC(e.time),"*","gap",std::to_string((int)gap)+" seconds"}); } } last["epoch"]=e.time; ++idx; }
    std::map<std::string,TimePoint> lastSat; std::map<std::string,double> lastPrph; for(auto& r:rf.data.observationRecords){ std::optional<double> Cc,Lc; for(auto& v:r.values){ if(opt.lli && !v.lli.empty()){ d.lliCount++; if(opt.clockSlips) d.slipEvents.push_back({formatUTC(r.time),r.satellite,"LLI",v.type+":"+v.lli}); } if(opt.snr && v.value&&v.type.rfind("S",0)==0){ sn.push_back(*v.value); auto band=v.type.size()>1?std::string(1,v.type[1]):std::string("1"); snByBand[band].push_back(*v.value); std::string codeKey=r.system+":"+v.type; snBySystemCode[codeKey].push_back(*v.value); if(opt.minSNR.count(band) && *v.value < opt.minSNR.at(band)){ d.snrStats[band].lowCount++; d.snrStats[codeKey].lowCount++; } } if(v.value&&v.type.rfind("C",0)==0&&!Cc) Cc=v.value; if(v.value&&v.type.rfind("L",0)==0&&!Lc) Lc=v.value; } if(Cc&&Lc){ double pp=*Cc-*Lc; prph.push_back(pp); if(lastPrph.count(r.satellite)){ double jump=std::fabs(pp-lastPrph[r.satellite]); double jumpCM=jump*100.0; if(opt.ion && jumpCM>opt.ionJumpCM) d.ionStats["all"].jumps++; if(opt.iod && lastSat.count(r.satellite)){ double dtmin=std::max(1e-9,std::chrono::duration<double>(r.time-lastSat[r.satellite]).count()/60.0); if(jumpCM/dtmin>opt.iodJumpCMPerMin) d.iodStats["all"].jumps++; } double ms=C*1.0e-3; double rem=std::fabs(jump-std::round(jump/ms)*ms); if(opt.clockSlips && rem<=opt.msecTol) d.clockSlipCount++; if(jump>50.0){ if(opt.multipath) d.multipathStats["all"].jumps++; if(opt.clockSlips) d.slipEvents.push_back({formatUTC(r.time),r.satellite,"P-L jump",std::to_string(jump)}); } } lastPrph[r.satellite]=pp; } if(lastSat.count(r.satellite)){ double dt=std::chrono::duration<double>(r.time-lastSat[r.satellite]).count(); if(opt.clockSlips && dt>opt.gapMinMinutes*60) d.slipEvents.push_back({formatUTC(r.time),r.satellite,"sat gap",std::to_string((int)dt)+" seconds"}); }
      std::set<char> bands; for(const auto& v:r.values) if(v.value && v.type.size()>=2) bands.insert(v.type[1]);
      if(opt.ion){
        for(char bi:bands){ double Pi=obsValue(r,'C',bi); if(!std::isfinite(Pi)) continue; for(char bj:bands){ if(bj<=bi) continue; double Pj=obsValue(r,'C',bj); if(std::isfinite(Pj)) ionVals.push_back(Pi-Pj); }}
      }
      bool arcBreak=false; for(const auto& vv:r.values) if(!vv.lli.empty()) arcBreak=true;
      if(opt.multipath) for(char bi:bands){ double Pi=obsValue(r,'C',bi); double Li=carrierMeters(r,bi,nullptr); double fi=signalFrequencyHz(r.system,bi,nullptr); if(!std::isfinite(Pi)||!std::isfinite(Li)||fi<=0) continue; for(char bj:bands){ if(bj==bi) continue; double Lj=carrierMeters(r,bj,nullptr); double fj=signalFrequencyHz(r.system,bj,nullptr); if(!std::isfinite(Lj)||fj<=0) continue; double a=(fi/fj)*(fi/fj); if(std::fabs(a-1.0)<1e-12) continue; double mp=Pi-(1.0+2.0/(a-1.0))*Li+(2.0/(a-1.0))*Lj; mpVals.push_back(mp); std::string key=std::string("MP")+bi+bj; mpComboBasic[r.system][key][r.satellite].push_back({r.time,mp,arcBreak}); }}
      lastSat[r.satellite]=r.time; }
    if(opt.snr){ int lowAll = 0; for(auto& kv:d.snrStats) lowAll += kv.second.lowCount; d.snrStats["all"]=stat(sn); d.snrStats["all"].lowCount = lowAll; for(auto& kv:snByBand){ int low=d.snrStats[kv.first].lowCount; d.snrStats[kv.first]=stat(kv.second); d.snrStats[kv.first].lowCount=low; } for(auto& kv:snBySystemCode){ int low=d.snrStats[kv.first].lowCount; d.snrStats[kv.first]=stat(kv.second); d.snrStats[kv.first].lowCount=low; } }
    if(opt.pseudorangePhase) d.pseudorangePhase["all"]=stat(prph);
    if(opt.multipath){ auto mpstat=stat(mpVals); mpstat.jumps=d.multipathStats["all"].jumps; d.multipathStats["all"]=mpstat; for(auto& skv: mpComboBasic){ for(auto& ckv: skv.second){ auto st=movingAverageMPStats(ckv.second,opt.mpWindow,true,opt.mpSigmas); if(st.count>0){ std::string key=skv.first+":"+ckv.first; d.multipathStats[key]=st; d.multipathMovingRMS[key]=st.rms; d.multipathMovingCount[key]=st.count; } } } int mp1c=0, mp2c=0; double mp1ss=0.0, mp2ss=0.0; for(const auto& kv:d.multipathStats){ if(kv.first.size()>=4 && kv.first.find("MP12")!=std::string::npos){ mp1c+=kv.second.count; mp1ss+=kv.second.rms*kv.second.rms*kv.second.count; } if(kv.first.size()>=4 && kv.first.find("MP21")!=std::string::npos){ mp2c+=kv.second.count; mp2ss+=kv.second.rms*kv.second.rms*kv.second.count; } } if(mp1c>0){ d.multipathMovingCount["MP1"]=mp1c; if(d.mp1Meters<=0) d.mp1Meters=std::sqrt(mp1ss/mp1c); } if(mp2c>0){ d.multipathMovingCount["MP2"]=mp2c; if(d.mp2Meters<=0) d.mp2Meters=std::sqrt(mp2ss/mp2c); } }
    if(opt.ion){ auto st=stat(ionVals); st.jumps=d.ionStats["all"].jumps; d.ionStats["all"]=st; d.histogramSamples["ion"]=(int)ionVals.size(); if(!ionVals.empty()) d.histograms["ion"]=histogramLinear(ionVals,opt.ionBins); } if(opt.iod){ auto st=stat({}); st.jumps=d.iodStats["all"].jumps; d.iodStats["all"]=st; } if(opt.snr){ d.histogramSamples["snr"]=(int)sn.size(); if(!sn.empty()) d.histograms["snr"]=histogramLinear(sn,opt.snBins,0.0,60.0); } if(opt.multipath){ d.histogramSamples["mp"]=(int)mpVals.size(); if(!mpVals.empty()) d.histograms["mp"]=histogramLinear(mpVals,opt.mpBins); } if(opt.pseudorangePhase){ d.histogramSamples["pseudorange_phase"]=(int)prph.size(); if(!prph.empty()) d.histograms["pseudorange_phase"]=histogramLinear(prph,opt.bins); } if(opt.pseudorangePhase && !prph.empty()){ auto pst=stat(prph); double sd=sampleStdDev(pst); if(sd>0 && std::fabs(pst.mean)>opt.codeSigmas*sd){ std::ostringstream w; w<<"code-sigma: mean exceeds "<<opt.codeSigmas<<" sigma gate"; d.thresholdWarnings.push_back(w.str()); } } d.deletedObservations=(d.codeBandCount<2?static_cast<int>(rf.data.observationRecords.size()):0); buildTeqcTimeplot(rf,opt,d); if(opt.symbolCodes||opt.allSymbols) d.symbolLegend={"c code/phase observation epoch","L loss-of-lock","g gap","N observed without navigation","2 incomplete dual-frequency","~ observed above mask","+ expected above mask","_ below mask","- below horizon"}; d.position.candidateEpochs=s.epochCount; d.position.skippedNoNavigation=(opt.averagePosition||opt.everyEpochPosition); d.position.attempted=false; d.position.epochSolutions=0; s.derived=d; } else if(rf.header.kind==RinexKind::Nav){ s.navigationRecords=(int)rf.data.navigationRecords.size(); for(auto& r:rf.data.navigationRecords){ s.navigationValues+=(int)r.values.size(); s.navigationFields+=(int)r.fields.size(); if(!r.satellite.empty())s.satelliteAppearance[r.satellite]++; if(r.epoch){ if(!s.firstEpoch||*r.epoch<*s.firstEpoch)s.firstEpoch=r.epoch; if(!s.lastEpoch||*r.epoch>*s.lastEpoch)s.lastEpoch=r.epoch; } } s.broadcastEphemerides=s.navigationRecords; } else if(rf.header.kind==RinexKind::Met){ s.meteorologicalRecords=(int)rf.data.meteorologicalRecords.size(); for(auto& r:rf.data.meteorologicalRecords) s.meteorologicalValues+=(int)r.values.size(); } return s; }
QCSummary analyzeWithNavigation(const RinexFile& rf,const std::vector<NavigationRecord>& navs,const QCOptions& opt){ auto s=analyze(rf,opt); if(rf.header.kind==RinexKind::Obs&&!navs.empty()){ s.residuals=residual(rf,navs,opt); applyNavBasedQCMetrics(rf, navs, opt, s); } return s; }
std::string makePlot(const QCSummary& summary){
  auto axis = [](size_t width){
    if (width == 72) return std::string("-----------|-----------|-----------|-----------|-----------|-----------|");
    std::string r(width, '-');
    for (size_t i = 11; i < width; i += 12) r[i] = '|';
    return r;
  };
  auto label = [](const std::string& sat){
    if (sat.empty()) return std::string();
    int p = prnNumber(sat);
    if (sat[0] == 'G') { std::ostringstream os; os << std::setw(3) << p; return os.str(); }
    if (sat[0] == 'R') { std::ostringstream os; os << 'R' << std::setw(2) << p; return os.str(); }
    return sat;
  };
  auto firstNonSpaceLocal = [](const std::string& row){ for(size_t i=0;i<row.size();++i) if(row[i]!=' ') return (int)i; return 9999; };
  auto lastNonSpaceLocal = [](const std::string& row){ for(int i=(int)row.size()-1;i>=0;--i) if(row[(size_t)i]!=' ') return i; return -1; };
  auto complexity = [](const std::string& row){ int c=0; for(char ch: row) if(ch!=' ' && ch!='~') ++c; return c; };

  std::ostringstream os;
  size_t width = 72;
  std::string obsLine(width, ' ');
  if(summary.derived && !summary.derived->obsTimeplot.empty()) { width = summary.derived->obsTimeplot.size(); obsLine = summary.derived->obsTimeplot; }
  std::string ax = axis(width);
  os << " SV+" << ax << "+ SV\n";
  if(summary.derived){
    std::vector<std::pair<std::string,std::string>> rows(summary.derived->satelliteTimeplot.begin(), summary.derived->satelliteTimeplot.end());
    std::sort(rows.begin(), rows.end(), [&](const auto& a,const auto& b){
      if(a.first.empty()||b.first.empty()) return a.first<b.first;
      char sa=a.first[0], sb=b.first[0];
      if(sa!=sb) return sa<sb;
      auto key=[&](const auto& r){
        const std::string& row=r.second;
        bool unhealthy=row.find("'")!=std::string::npos;
        int first=firstNonSpaceLocal(row), last=lastNonSpaceLocal(row);
        int complex=complexity(row);
        int prn=prnNumber(r.first);
        int sysWeight = r.first.empty()?9:(r.first[0]=='G'?0:(r.first[0]=='R'?1:2));
        int group=5;
        if(unhealthy) group=0; else if(first>0) group=3; else if(!row.empty() && row[0]=='2') group=1; else if(complex<=2) group=2; else group=4;
        double mxEl = summary.derived->satelliteMaxElevationDeg.count(r.first) ? summary.derived->satelliteMaxElevationDeg.at(r.first) : -999.0;
        int elevKey = static_cast<int>(std::lround(-mxEl*10.0));
        return std::tuple<int,int,int,int,int,int,int>(sysWeight, group, unhealthy?last:first, elevKey, complex, last, prn);
      };
      return key(a)<key(b);
    });
    for(auto& r: rows){
      std::string lab=label(r.first);
      os << std::setw(3) << lab << '|' << r.second << '|' << std::setw(3) << lab << "\n";
    }
  }
  bool navAssisted = (summary.derived && summary.derived->position.epochSolutions > 0) || !summary.navInputFiles.empty();
  if(navAssisted){
    os << "-dn|" << std::string(width, ' ') << "|-dn\n";
    os << "+dn|" << obsLine << "|+dn\n";
    std::string plus10 = obsLine;
    for(char& c: plus10){ if(c>='0' && c<='8') c='9'; else if(c=='9') c='a'; }
    os << "+10|" << plus10 << "|+10\n";
    std::string pos = (summary.derived && !summary.derived->positionTimeplot.empty()) ? summary.derived->positionTimeplot : std::string(width,' ');
    for(char& c: pos) if(c=='P') c='o';
    os << "Pos|" << pos << "|Pos\n";
  } else {
    os << "Obs|" << obsLine << "|Obs\n";
  }
  os << "Clk|" << std::string(width, ' ') << "|Clk\n";
  os << "   +" << ax << "+   \n";
  if(summary.firstEpoch && summary.lastEpoch){
    os << "first: " << formatUTC(*summary.firstEpoch) << "\n";
    os << "last : " << formatUTC(*summary.lastEpoch) << "\n";
  }
  os << "symbols: c code/phase, L loss-of-lock, N no-nav, 2 incomplete dual-frequency, ~ observed above mask, + expected, _ below mask, - below horizon\n";
  return os.str();
}}
