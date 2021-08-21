#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

template <typename T> class CsvFile {
  public:
    CsvFile(const std::string& filename, const std::vector<std::string>& keys)
        : fs_(), n_keys_(keys.size()) {
        fs_.exceptions(std::ios::failbit | std::ios::badbit | std::ios::out |
                       std::ios::app);
        fs_.open(filename);

        for (int i = 0; i < keys.size(); ++i) {
            const auto key = keys.at(i);
            if (i < keys.size() - 1) {
                fs_ << key << ",";
            } else {
                fs_ << key;
            }
        }

        fs_ << std::endl;
    }

    ~CsvFile() {
        fs_.flush();
        fs_.close();
    }

    CsvFile& operator<<(const std::vector<T>& values) {
        assert(values.size() == n_keys_);
        for (int i = 0; i < values.size(); ++i) {
            const auto value = values.at(i);
            if (i < values.size() - 1) {
                fs_ << value << ",";
            } else {
                fs_ << value;
            }
        }

        fs_ << std::endl;

        return *this;
    }

  private:
    unsigned int n_keys_ = 0;

    std::ofstream fs_;
};
