#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pch.h"
#include "CppUnitTest.h"

#include "../deformo/EigenTypes.h"
#include "../deformo/Integrators.h"
#include "../deformo/Integrators.cpp"
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

			const auto sim = std::make_unique<Simulation>(1., 210e6, 0.3, particles);

			Assert::IsNotNull(sim.get());
		}

		TEST_METHOD(MatlabBookParameterComparison) 
		{
		  Eigen::VectorXd displacements(6);
		  displacements << 0, 0, 0.5, 0, 0.5, 0.25;

		  const auto sim = std::make_unique<Simulation>(1., 210e6, 0.3, displacements);
		  //Eigen::Matrix66d kk1;
		  //kk1 << 2.0192, 0, 0, -1.0096, -2.0192, 1.0096; 
		  //kk1 << 0, 5.7692, -0.8654, 0, 0.8654, -5.7692; 
		  //kk1 << 0, -0.8654, 1.4423, 0, -1.4423, 0.8654; 
		  //kk1 << -1.0096, 0, 0, 0.5048, 1.0096, -0.5048;
		  //kk1 << -2.0192, 0.8654, -1.4423, 1.0096, 3.4615, -1.8750;
		  //kk1 << 1.0096, -5.7692, 0.8654, -0.5048, -1.8750, 6.2740;

		  //Assert::AreEqual(kk1, std::get<0>(sim->k[0]));
		  Assert::IsNotNull(sim.get());
		}
	};
}
