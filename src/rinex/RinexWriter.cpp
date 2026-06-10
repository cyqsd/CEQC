#include "ceqc/rinex/RinexService.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
namespace ceqc::service::rinex {
void write(std::ostream& os, const RinexFile& rf){
  // Anubis 3.11 has a known fragile RINEX2 reader path that aborts on some
  // LF-only RINEX2 continuation records even when all records are 80 columns.
  // CRLF is valid text transport for RINEX and keeps the numerical content and
  // fixed columns unchanged, so emit RINEX2 files with CRLF.  RTKLIB/RTKPLOT is
  // Windows-centric and RTKCONV-EX golden files use CRLF, so +rtklib/+rtkplot
  // emits CRLF for all requested versions. Strict RINEX3/4 remain LF.
  const char* eol = (writerOptions().rtklibCompat || (rf.header.version > 0.0 && rf.header.version < 3.0)) ? "\r\n" : "\n";
  for(const auto& h: rf.header.lines) os<<h.raw<<eol;
  for(const auto& b: rf.body) os<<b<<eol;
}
void writeFile(const std::string& path,const RinexFile& rf){
  // Open in binary mode even though the payload is text.  On Windows, text-mode
  // streams translate every '\n' to '\r\n'.  Because write() deliberately
  // emits RTKCONV/RTKPLOT-compatible CRLF for +rtklib/+rtkplot, text mode would
  // turn those records into "\r\r\n", which RTKPlot handles poorly and which
  // no longer matches RTKCONV-EX golden files.  Binary mode preserves the exact
  // requested newline sequence on Linux and Windows.
  std::ofstream f(path, std::ios::binary);
  if(!f) throw std::runtime_error("cannot create "+path);
  write(f,rf);
}
}
