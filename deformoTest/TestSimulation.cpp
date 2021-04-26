#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../deformo/EigenTypes.h"
#include "pch.h"
#include "CppUnitTest.h"
#include "../deformo/Simulation.h"
#include "../deformo/Simulation.cpp"
#include "../deformo/Integrators.h"
#include "../deformo/Integrators.cpp"

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

			const auto sim = std::make_unique<Simulation>(1., 210e6, 0.3, particles);

			Assert::AreEqual(40000., 1 / std::pow(0.005, 2));
			Assert::IsNotNull(sim.get());
		}
	};
}
