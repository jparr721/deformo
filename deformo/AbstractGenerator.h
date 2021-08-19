#pragma once

#include "DeformoAssert.h"
#include "Numerics.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

template <typename T = std::vector<Real>> class AbstractGenerator {
  public:
    std::string dataset_path;
    std::unordered_map<std::string_view, std::vector<T>> dataset;

    virtual ~AbstractGenerator() = default;

    virtual auto Capture() -> void = 0;
    virtual auto WriteDatasetToCSVPath() -> bool final {
        _CsvFile c(dataset_path);
        return c.Absorb(dataset);
    }

  private:
    class _CsvFile {
      public:
        _CsvFile(const std::string& filename) : fs_() {
            fs_.exceptions(std::ios::failbit | std::ios::badbit);
            fs_.open(filename);
        }

        ~_CsvFile() {
            fs_.flush();
            fs_.close();
        }

        template <typename _T> _CsvFile& operator<<(const _T& value) {
            fs_ << value << ",";
            return *this;
        }

        template <typename _T>
        auto Absorb(
            const std::unordered_map<std::string_view, std::vector<_T>>& data)
            -> bool {
            // First, drop the keys in.
            std::vector<std::string_view> keys;
            int n_keys = data.size();
            int i = 0;
            for (const auto& [k, _] : data) {
                if (i < n_keys - 1) {
                    fs_ << k << ",";
                } else {
                    fs_ << k;
                }
                keys.emplace_back(k);
                ++i;
            }

            End();

            // Now, go line-by-line.
            const int row_len = data.begin()->second.size();
            for (i = 0; i < row_len; ++i) {
                for (int k = 0; k < keys.size(); ++k) {
                    const auto key = keys.at(k);
                    if (k < keys.size() - 1) {
                        fs_ << data.at(key)[i] << ",";
                    } else {
                        fs_ << data.at(key)[i];
                    }
                }
                End();
            }

            return true;
        }

      private:
        std::ofstream fs_;

        auto End() -> void { fs_ << std::endl; }
    };

    DeformoAssertion assertion_;
};
