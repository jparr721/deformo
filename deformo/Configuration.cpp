#include "Configuration.h"
#include "Utils.h"

#include <cassert>
#include <filesystem>

Configuration::Configuration() {
    assert(Load(".env") && "FAILED TO LOAD ENV CONFIG");
}

auto Configuration::Value(const std::string& key) -> std::string {
    return values_[key];
}

auto Configuration::Load(const std::string& env_path) -> bool {
    assert(std::filesystem::exists(env_path) && "FILE DOES NOT EXIST");

    std::ifstream file = utils::OpenFile(env_path);
    for (std::string line; std::getline(file, line);) {
        const auto kv = utils::Split(line, "=");
        assert(line.length() > 3 && "IMPROPER ENV CONFIG");
        values_.insert(kv);
    }

    return true;
}
