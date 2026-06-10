#include "ceqc/translate/Translator.hpp"
#include "ceqc/rinex/RinexService.hpp"
#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace ceqc::service::translator {
std::shared_ptr<Translator> makeRTCM3();
std::shared_ptr<Translator> makeUBX();
namespace { std::string lower(std::string s){for(auto&c:s)c=std::tolower((unsigned char)c);return s;} std::vector<unsigned char> prefix(const std::string& p){ std::ifstream f(p,std::ios::binary); std::vector<unsigned char> b(512); f.read((char*)b.data(),b.size()); b.resize((size_t)std::max<std::streamsize>(0,f.gcount())); return b; }
class RinexTranslator final: public Translator { public: FormatInfo format() const override{return {"rinex",{"rnx","obs","nav","met"},"RINEX OBS/NAV/MET",true};} bool probe(const std::string& path,const std::vector<unsigned char>& p) const override{ std::string s((const char*)p.data(),p.size()); auto lp=lower(path); return s.find("RINEX VERSION / TYPE")!=std::string::npos || lp.ends_with(".rnx") || lp.ends_with("o") || lp.ends_with("n") || lp.ends_with("p");} std::vector<RinexFile> decode(const std::string& path) const override{return {ceqc::service::rinex::readFile(path)};} };
class Placeholder final: public Translator { FormatInfo i_; public: Placeholder(std::string n,std::vector<std::string>a,std::string d):i_{std::move(n),std::move(a),std::move(d),false}{} FormatInfo format() const override{return i_;} bool probe(const std::string& path,const std::vector<unsigned char>&) const override{auto lp=lower(path); if(lp.find(i_.name)!=std::string::npos)return true; for(auto&a:i_.aliases) if(lp.find(lower(a))!=std::string::npos)return true; return false;} std::vector<RinexFile> decode(const std::string&) const override{throw std::runtime_error("translator registered but decoder not implemented: "+i_.name);} };
}
Registry::Registry(){
  items_.push_back(std::make_shared<RinexTranslator>());
  items_.push_back(makeRTCM3());
  items_.push_back(makeUBX());
  items_.push_back(std::make_shared<Placeholder>("trimble-dat",std::vector<std::string>{"trimble","dat","t00","t01","t02"},"Trimble DAT/T0x download files"));
  items_.push_back(std::make_shared<Placeholder>("trimble-rt17",std::vector<std::string>{"rt17"},"Trimble RT17 stream"));
  items_.push_back(std::make_shared<Placeholder>("trimble-tsip",std::vector<std::string>{"tsip"},"Trimble TSIP stream"));
  items_.push_back(std::make_shared<Placeholder>("ashtech",std::vector<std::string>{"ashtech-b","ashtech-e","ashtech-s","ashtech-d","ashtech-r","ashtech-u","ash"},"Ashtech download and real-time streams"));
  items_.push_back(std::make_shared<Placeholder>("leica-lb2",std::vector<std::string>{"lb2","leica"},"Leica LB2"));
  items_.push_back(std::make_shared<Placeholder>("leica-mdb",std::vector<std::string>{"mdb"},"Leica MDB"));
  items_.push_back(std::make_shared<Placeholder>("leica-ds",std::vector<std::string>{"ds"},"Leica DS file set"));
  items_.push_back(std::make_shared<Placeholder>("topcon-tps",std::vector<std::string>{"tps","jps","javad","topcon"},"Topcon TPS / Javad JPS"));
  items_.push_back(std::make_shared<Placeholder>("septentrio-sbf",std::vector<std::string>{"sbf","septentrio"},"Septentrio SBF"));
  items_.push_back(std::make_shared<Placeholder>("navcom",std::vector<std::string>{"navcom"},"NavCom binary"));
  items_.push_back(std::make_shared<Placeholder>("binex",std::vector<std::string>{"bnx"},"BINEX"));
  items_.push_back(std::make_shared<Placeholder>("nmea",std::vector<std::string>{"nmea0183"},"NMEA 0183 text stream"));
  items_.push_back(std::make_shared<Placeholder>("motorola-oncore",std::vector<std::string>{"oncore"},"Motorola Oncore"));
  items_.push_back(std::make_shared<Placeholder>("rockwell-zodiac",std::vector<std::string>{"zodiac"},"Rockwell Zodiac"));
  items_.push_back(std::make_shared<Placeholder>("argo",std::vector<std::string>{"argo-dat","argo-orb"},"ARGO DAT/ORB"));
  items_.push_back(std::make_shared<Placeholder>("ti4100",std::vector<std::string>{"gesar","bepp","core","ti-rom"},"TI-4100 formats"));
  items_.push_back(std::make_shared<Placeholder>("rogue-turbo",std::vector<std::string>{"conanbinary","turbobinary","turborogue","turbostar"},"ConanBinary / TurboBinary"));
}
const Translator* Registry::find(const std::string& name) const{auto n=lower(name); for(auto& t:items_){auto f=t->format(); if(lower(f.name)==n)return t.get(); for(auto&a:f.aliases) if(lower(a)==n)return t.get();} return nullptr;}
const Translator* Registry::detect(const std::string& path) const{auto p=prefix(path); for(auto& t:items_) if(t->probe(path,p)) return t.get(); return nullptr;}
std::vector<FormatInfo> Registry::formats() const{std::vector<FormatInfo> f; for(auto&t:items_)f.push_back(t->format()); return f;}
}
