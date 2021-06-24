#include "QTUtils.h"

namespace utils::qt {
std::string QStringToString(const QString& string) {
    return string.toUtf8().constData();
}
} // namespace utils::qt
