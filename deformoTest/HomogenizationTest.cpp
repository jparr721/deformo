#include "pch.h"

#include <Eigen/SparseLU>
#include <iostream>
#include <memory>

#include "../deformo/Numerics.h"
#include "../deformo/Homogenization.cpp"
#include "../deformo/Homogenization.h"
#include "../deformo/Rve.cpp"
#include "../deformo/Rve.h"

TEST(TestHomogenization, TestHexahedron) {
	const auto rve = std::make_shared<Rve>(10, 10, 10, 0.1, 0.1, 0.1);
	ASSERT_TRUE(rve.get() != nullptr);
	
	const auto homogenization = std::make_shared<Homogenization>(rve);
	ASSERT_TRUE(homogenization.get() != nullptr);

	const auto hexahedron =
            homogenization->ComputeHexahedron(0.5, 0.5, 0.5);

	for (const auto v : hexahedron) {
	  std::cout << "Rows: " << v.rows() << " Cols: " << v.cols() << std::endl;
	  std::cout << v << std::endl;
	}
}
