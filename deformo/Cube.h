#pragma once

#include "Numerics.h"

class Cube {
  public:
    Cube();
    Cube(const Vector3r& position);

    auto Positions() -> VectorXr;

  private:
    VectorXr positions_;
};
