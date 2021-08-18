#pragma once

#include <string>
#include <vector>
#include "Numerics.h"
#include "ThirdParty/cppflow/model.h"
#include "ThirdParty/cppflow/ops.h"

class NeuralNetwork {
  public:
    explicit NeuralNetwork(const std::string& model_path);

    auto Predict(const MatrixXr& input) -> std::vector<Real>;
    auto Predict(const Tensor3r& input) -> std::vector<Real>;
    auto Predict(const cppflow::tensor& input) -> std::vector<Real>;

    private:
    cppflow::model m_;
};
