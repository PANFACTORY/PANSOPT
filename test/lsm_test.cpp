#include "../src/lsm.h"

#include <gtest/gtest.h>

using namespace PANSOPT;

TEST(LSMTest, LSMTest1) {
    auto getg = [](const std::vector<double> &s) {
        return std::accumulate(s.begin(), s.end(), 0.0) - 10.0;
    };

    std::vector<double> s(12, 1), dfds(12, -1), dgds(12, 1);
    LSM<double> optimizer(3, 4, 0.03, 0.05);

    for (int i = 0; i < 100; ++i) {
        optimizer.UpdateVariables(s, dfds, getg(s), dgds);
        std::cout << getg(s) << std::endl;
    }
}