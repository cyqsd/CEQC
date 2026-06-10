#pragma once
#include "ceqc/core/Command.hpp"
#include <string>
#include <vector>
namespace ceqc::cli {
ceqc::model::Operation parseArgs(const std::vector<std::string>& args);
}
