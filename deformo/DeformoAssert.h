#pragma once

#include <cassert>
#include <sstream>
#include <string>
#include <iostream>

class DeformoAssertion {
  public:
    enum class ErrorLevel {
        kDebug = 0x00,
        kInfo,
        kWarn,
        kCritical,
    };

    ErrorLevel error_level;

    explicit DeformoAssertion(ErrorLevel level = ErrorLevel::kCritical)
        : error_level(level) {}

    template <typename... Context>
    void Assert(bool condition, const std::string& function_name,
                const std::string& file, int line, Context... context) {
        if (!condition) {
            std::ostringstream oss;
            LogAll(oss, context...);

            const std::string assert_caller = __FUNCTION__;

            std::cout << "[" << ErrorLevelToString() << "] " << assert_caller
                      << "() had an error in: " << function_name << "() in "
                      << file << ":" << line << ":" << std::endl
                      << oss.str() << std::endl;

            assert(condition && error_level == ErrorLevel::kCritical);
        }
    }

  private:
    template <typename T> void LogAll(std::ostream& o, T val) {
        o << val << " ";
    }

    template <typename T, typename... Context>
    void LogAll(std::ostream& o, T val, Context... context) {
        LogAll(o, val);
        LogAll(o, context...);
    }

    auto ErrorLevelToString() -> std::string {
        switch (error_level) {
        case ErrorLevel::kDebug:
            return "Debug";
        case ErrorLevel::kInfo:
            return "Info";
        case ErrorLevel::kWarn:
            return "Warning";
        case ErrorLevel::kCritical:
            return "Critical";
        }
    }
};
