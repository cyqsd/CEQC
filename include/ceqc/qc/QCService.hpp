#pragma once
#include "ceqc/qc/QC.hpp"
#include <string>
#include <vector>
#include <optional>
#include <array>
namespace ceqc::service::qc {
using namespace ceqc::model;
QCSummary analyze(const RinexFile& rf, const QCOptions& options = {});
QCSummary analyzeWithNavigation(const RinexFile& rf, const std::vector<NavigationRecord>& navs, const QCOptions& options = {});
std::optional<std::array<double,3>> estimateApproxPosition(const RinexFile& rf, const std::vector<NavigationRecord>& navs, const QCOptions& options = {});
std::string makePlot(const QCSummary& summary);
}
