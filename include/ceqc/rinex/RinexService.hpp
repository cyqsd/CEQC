#pragma once
#include "ceqc/core/Command.hpp"
#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace ceqc::service::rinex {
using namespace ceqc::model;

struct ParserOptions { bool relax=false; bool extendRinex2=false; bool reformat=false; };
struct WriterOptions { bool rtklibCompat=false; };
void setParserOptions(ParserOptions options);
ParserOptions parserOptions();
void setWriterOptions(WriterOptions options);
WriterOptions writerOptions();
RinexFile readFile(const std::string& path);
RinexFile readStream(const std::string& path, std::istream& is);
void parseContent(RinexFile& rf);
std::vector<ValidationIssue> validate(const RinexFile& rf);
void write(std::ostream& os, const RinexFile& rf);
void writeFile(const std::string& path, const RinexFile& rf);
void applyHeaderEdits(RinexFile& rf, const std::map<std::string,std::string>& edits);
void applyObsTypeFilter(RinexFile& rf, const std::string& obsList, bool exclude, bool renameOnly=false);
void applyMetTypeFilter(RinexFile& rf, const std::string& obsList, bool exclude, bool renameOnly=false);
RinexFile windowObservation(const RinexFile& rf, const std::optional<TimePoint>& start, const std::optional<TimePoint>& end);
RinexFile decimate(const RinexFile& rf, const DecimationSpec& spec);
RinexFile merge(const std::vector<RinexFile>& files, RinexKind kind, double targetVersion);
std::string normalizeHeaderLabel(const std::string& opt);
}
