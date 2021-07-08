#pragma once

#include <iostream>
#include <memory>
#include <tensorflow/c/c_api.h>

//inline auto TfOkay(const std::shared_ptr<TF_Status>& status) {
//    if (TF_GetCode(status.get()) != TF_OK) {
//        std::cerr << "Error: " << TF_Message(status.get()) << std::endl;
//        return 0;
//    }
//
//    return 1;
//}
