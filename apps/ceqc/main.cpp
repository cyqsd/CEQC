#include "ceqc/app/Application.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
int main(int argc,char**argv){
  std::vector<std::string> args;
  std::unique_ptr<std::ofstream> outFile, errFile;
  for(int i=1;i<argc;++i){
    std::string a=argv[i];
    if((a=="+out"||a=="++out") && i+1<argc){ outFile=std::make_unique<std::ofstream>(argv[++i], a=="++out"?std::ios::app:std::ios::out); continue; }
    if((a=="+err"||a=="++err") && i+1<argc){ errFile=std::make_unique<std::ofstream>(argv[++i], a=="++err"?std::ios::app:std::ios::out); continue; }
    args.emplace_back(std::move(a));
  }
  std::ostream& out = (outFile && *outFile) ? *outFile : std::cout;
  std::ostream& err = (errFile && *errFile) ? *errFile : std::cerr;
  ceqc::app::Application app(out,err); return app.run(args);
}
