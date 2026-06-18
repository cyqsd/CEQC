#pragma once
#include "ceqc/qc/QC.hpp"
#include <iosfwd>

namespace ceqc::view {
void writeQCJson(std::ostream& os, const ceqc::model::QCSummary& s);
}
