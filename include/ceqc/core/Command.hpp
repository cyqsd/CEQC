#pragma once
#include "ceqc/qc/QC.hpp"
#include <chrono>
#include <map>
#include <optional>
#include <string>
#include <vector>
namespace ceqc::model {
struct DecimationSpec { bool enabled=false; std::chrono::seconds interval{0}; std::chrono::seconds offset{0}; std::string raw; };
struct SVSelection { bool specified=false; bool includeMode=true; bool all=true; std::vector<int> prns; };
struct Operation {
  bool showHelp=false, showVersion=false, showID=false, showFormats=false, showConfig=false, showAllConfig=false, showBCF=false, verifyOnly=false, qc=false, quietQC=false, teqcCompat=false;
  std::string outputObs, outputNav, outputMet, outputBinex, translatorName; double targetVersion=0;
  std::vector<std::string> configFiles; std::string teqcGolden, teqcDiff, teqcEOL="lf";
  std::map<std::string,std::string> obsHeaderEdits, navHeaderEdits, metHeaderEdits, headerEdits, binexMetadata;
  std::vector<std::string> extractLabels; std::optional<TimePoint> windowStart, windowEnd;
  DecimationSpec obsDecimation, navDecimation, metDecimation; std::string obsSummaryTarget; bool obsSummaryAppend=false;
  std::map<std::string,SVSelection> svSelections; int maxRxChannels=12; int maxRxSVs=12; bool maxRxSVsSpecified=false; std::map<std::string,int> maxExpectedSVs; bool useAllChannels=true; bool allowNaNObs=true; bool svDuplicates=false; bool orderSVByPRN=false; bool rtklibCompat=false; bool relaxRinex=false; bool extendRinex2=false; bool reformatRinex=false;
  std::vector<std::string> ephIncludeTypes; std::vector<std::string> ephExcludeTypes;
  std::vector<std::string> inputFiles; QCOptions qcOptions;
};
}
