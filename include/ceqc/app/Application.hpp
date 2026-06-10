#pragma once
#include <iosfwd>
#include <string>
#include <vector>
namespace ceqc::app {
class Application {
  std::ostream& out_; std::ostream& err_;
public:
  Application(std::ostream& out, std::ostream& err): out_(out), err_(err) {}
  int run(const std::vector<std::string>& args);
};
}
