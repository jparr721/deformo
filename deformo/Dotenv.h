#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace config_keys {
const static std::string kVertexShaderLocation = "VERTEX_SHADER_LOCATION";
const static std::string kFragmentShaderLocation = "FRAGMENT_SHADER_LOCATION";
} // namespace config_keys

class Configuration {
  public:
    Configuration();

    auto Load(const std::string& env_path) -> bool;
    auto Value(const std::string& key) -> std::string;

  private:
    std::unordered_map<std::string, std::string> values_;
};

static const auto env_config = std::make_unique<Configuration>();
