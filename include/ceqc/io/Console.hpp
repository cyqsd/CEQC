#pragma once
#include "ceqc/qc/QC.hpp"
#include <iosfwd>
namespace ceqc::view {
void printHelp(std::ostream& os);
void printVersion(std::ostream& os);
void printIssues(std::ostream& os, const std::string& path, const std::vector<ceqc::model::ValidationIssue>& issues);
void printQC(std::ostream& os, const ceqc::model::QCSummary& s, bool quiet=false, bool teqcCompat=false);
}
