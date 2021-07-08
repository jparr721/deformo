#include "TFModel.h"
#include <cassert>
#include <filesystem>
#include <fstream>

TfModel::TfModel(const std::string& model_path) : model_path_(model_path) {}

auto TfModel::ReadModelFile(const std::string& model_path) const -> TF_Buffer* {
    assert(std::filesystem::exists(model_path) &&
           "MODEL FILE NOT FOUND AT PATH");

    const int model_size = std::filesystem::file_size(model_path);

    assert(model_size > 0 && "MODEL FILE EMPTY");

    std::string model_data;
    std::ifstream in(model_path);

    assert(in.good() && "FAILED TO READ MODEL PATH");

    std::getline(in, model_data,
                 std::string::traits_type::to_char_type(
                     std::string::traits_type::eof()));

    return TF_NewBufferFromString(model_data.c_str(), model_size);
}
