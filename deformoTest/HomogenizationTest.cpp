#include "pch.h"

#include <Eigen/SparseLU>
#include <iostream>
#include <memory>

#include "../deformo/Numerics.h"
#include "../deformo/Homogenization.cpp"
#include "../deformo/Homogenization.h"
#include "../deformo/Rve.cpp"
#include "../deformo/Rve.h"

const auto rve = std::make_shared<Rve>(10, 10, 10, 0.1, 0.1, 0.1);

TEST(TestHomogenization, TestHexahedron) {
	ASSERT_TRUE(rve.get() != nullptr);
	
	const auto homogenization = std::make_shared<Homogenization>(rve);
	ASSERT_TRUE(homogenization.get() != nullptr);

	const auto hexahedron =
            homogenization->ComputeHexahedron(0.5, 0.5, 0.5);
	ASSERT_EQ(hexahedron.size(), 4);
}

TEST(TestHomogenization, TestComputeDegreesOfFreedom) {
	ASSERT_TRUE(rve.get() != nullptr);
	
	const auto homogenization = std::make_shared<Homogenization>(rve);
	ASSERT_TRUE(homogenization.get() != nullptr);

	homogenization->ComputeDegreesOfFreedom(1000);
}
