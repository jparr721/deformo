#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pch.h"
#include "CppUnitTest.h"

#include "../deformo/EigenTypes.h"
#include "../deformo/Input.h"
#include "../deformo/Input.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace SimulationTest
{
	TEST_CLASS(TestInput)
	{
	public:
		TEST_METHOD(InputConstructsProperly)
		{
			const auto input = std::make_unique<Input>();
			Assert::IsNotNull(input.get());
		}

		TEST_METHOD(MatlabBookParameterComparison) 
		{
		}
	};
}
