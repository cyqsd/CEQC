#pragma once
#include "ceqc/rinex/Rinex.hpp"
#include <memory>
#include <string>
#include <vector>

namespace ceqc::service::translator {
using namespace ceqc::model;
struct FormatInfo { std::string name; std::vector<std::string> aliases; std::string description; bool implemented=false; };
class Translator {
public:
  virtual ~Translator()=default;
  virtual FormatInfo format() const=0;
  virtual bool probe(const std::string& path, const std::vector<unsigned char>& prefix) const=0;
  virtual std::vector<RinexFile> decode(const std::string& path) const=0;
};
class Registry {
  std::vector<std::shared_ptr<Translator>> items_;
public:
  Registry();
  const Translator* find(const std::string& name) const;
  const Translator* detect(const std::string& path) const;
  std::vector<FormatInfo> formats() const;
};
}
