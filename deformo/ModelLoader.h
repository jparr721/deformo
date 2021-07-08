#pragma once

#include "Tensorflow.h"
#include <string>

class ModelLoader {
  public:
    explicit ModelLoader(const std::string& path);
};
