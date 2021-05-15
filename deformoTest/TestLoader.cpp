#include "pch.h"
#include "CppUnitTest.h"
#include <filesystem>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>

#include "../deformo/Loader.h"
#include "../deformo/Loader.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace deformoTest
{
	TEST_CLASS(TestLoader)
	{
	public:
		TEST_METHOD(LoaderLoadsTetgen){
			const std::string node_data =
R"(2  3  0  0
   1    0  0  0
   2    1  0  0
#)";

			const std::string ele_data =
R"(2  4  0
    1       8     3     1     4
    2       1     6     7     5
#)";

			const std::string face_data =
R"(2  4  0
    1       8     3     1    -1 
    2       1     6     7    -1 
#)";


            const std::string cdir = std::filesystem::current_path().string();

			const std::filesystem::path node_path = std::filesystem::path(cdir + "/v.node");
			const std::filesystem::path ele_path = std::filesystem::path(cdir + "/v.ele");
                        const std::filesystem::path face_path =
                            std::filesystem::path(cdir + "/v.face");

			std::ofstream node_file_ptr;
			node_file_ptr.open(node_path);
			node_file_ptr << node_data;
			node_file_ptr.close();

			std::ofstream face_file_ptr;
			face_file_ptr.open(face_path);
			face_file_ptr << face_data;
			face_file_ptr.close();

			std::ofstream ele_file_ptr;
			ele_file_ptr.open(ele_path);
			ele_file_ptr << ele_data;
			ele_file_ptr.close();

			Eigen::MatrixXf V;
			Eigen::MatrixXf F;
			Eigen::MatrixXf T;

			loader::ReadTetgenVertexFile(V, node_path.string());
			loader::ReadTetgenFaceFile(F, face_path.string());
			loader::ReadTetgenEleFile(T, ele_path.string());

			// Clean up after the tests complete
			std::filesystem::remove(node_path);
			std::filesystem::remove(ele_path);
			std::filesystem::remove(face_path);

			Assert::AreEqual(static_cast<int>(V.rows()), 2);
			Assert::AreEqual(static_cast<int>(V.cols()), 3);
			Assert::AreEqual(V(1, 0), 1.f);

			Assert::AreEqual(static_cast<int>(T.rows()), 2);
			Assert::AreEqual(static_cast<int>(T.cols()), 4);
			Assert::AreEqual(T(1, 0), 1.f);
			Assert::AreEqual(T(1, 1), 6.f);
			Assert::AreEqual(T(1, 2), 7.f);
			Assert::AreEqual(T(1, 3), 5.f);

			Assert::AreEqual(static_cast<int>(F.rows()), 2);
			Assert::AreEqual(static_cast<int>(F.cols()), 3);
			Assert::AreEqual(T(1, 0), 1.f);
			Assert::AreEqual(T(1, 1), 6.f);
			Assert::AreEqual(T(1, 2), 7.f);
		}
	};
}
