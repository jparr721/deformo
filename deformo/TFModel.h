#pragma once

#include "Tensorflow.h"
#include <string>
#include <memory>

class TfModel {
  public:
    explicit TfModel(const std::string& model_path);
    auto Predict() -> int;

  private:
    const std::string model_path_;

    std::shared_ptr<TF_Graph> graph_;
    std::shared_ptr<TF_Session> session_;
    //std::unique_ptr<TF_Status> status_;

    auto ReadModelFile(const std::string& model_path) const -> TF_Buffer*;
};
