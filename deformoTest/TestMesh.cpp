
#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pch.h"
#include "CppUnitTest.h"

#include "../deformo/EigenTypes.h"
#include "../deformo/Mesh.h"
#include "../deformo/Mesh.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace SimulationTest
{
	TEST_CLASS(TestMesh)
	{
	public: 
		TEST_METHOD(MeshConstructsProperly) {
			Eigen::VectorXd vertices(12);
			vertices << 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1;

		    const auto mesh = std::make_unique<Mesh>(vertices);

			const int nrows = mesh->vertices.rows();
			const int nmatches = mesh->indices.size();

			Assert::AreEqual(nrows, 12);
			Assert::AreEqual(nmatches, 6);
		}
	};
}
