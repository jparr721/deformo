#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../deformo/EigenTypes.h"
#include "pch.h"
#include "CppUnitTest.h"
#include "../deformo/Simulation.h"
#include "../deformo/Simulation.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace SimulationTest
{
	TEST_CLASS(TestSimulation)
	{
	public:
		TEST_METHOD(SimulationConstructsProperly)
		{
			Eigen::VectorXd particles(10);
			particles << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

			const auto fixed_indices =
			std::vector<Eigen::Vector3d>{Eigen::Vector3d(1, 2, 3)};
			const auto sim = new Simulation(particles, fixed_indices);

			Assert::IsNotNull(sim);

			delete sim;
		}
	};
}
