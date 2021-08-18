#include "NeuralNetwork.h"
#include <iostream>

NeuralNetwork::NeuralNetwork(const std::string& model_path)
    : m_(cppflow::model(model_path)) {}

auto NeuralNetwork::Predict(const cppflow::tensor& input) -> std::vector<Real> {
    auto output = m_(input);

    std::cout << output << std::endl;

    const auto values = output.get_data<Real>();

    return values;
}
