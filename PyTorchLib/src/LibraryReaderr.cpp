//
// Created by anichols on 9/17/23.
//

#include "LibraryReaderr.h"

#include <torch/torch.h>
#include <iostream>
#include <torch/nn/modules/loss.h>
#include <torch/optim/adam.h>

#include <QDebug>

void LibraryReaderr::test() {
    torch::Tensor tensor = torch::rand({2, 3});
    std::cout << tensor << std::endl;
}
