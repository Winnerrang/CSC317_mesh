#include "write_obj.h"
#include <fstream>
#include <cassert>
#include <iostream>

// asume input matrix starts from 0
// therefore we need to add 1 to it for each indices
bool write_obj(
  const std::string & filename,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & UV,
  const Eigen::MatrixXi & UF,
  const Eigen::MatrixXd & NV,
  const Eigen::MatrixXi & NF)
{
  assert((F.size() == 0 || F.cols() == 3 || F.cols() == 4) && "F must have 3 or 4 columns");
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here:
  
  std::ofstream ofs(filename, std::ios_base::out);

  // Vertex
  
  for (int row = 0; row < V.rows(); row++) {
	  ofs << "v";
	  for (int col = 0; col < V.cols(); col++) {
		  ofs << " " << V(row, col);
	  }

	  ofs << std::endl;
  }
  
  // texture coordinates
  for (int row = 0; row < UV.rows(); row++) {
	  ofs << "vt";
	  for (int col = 0; col < UV.cols(); col++) {
		  ofs << " " << UV(row, col);
	  }

	  ofs << std::endl;
  }

  // normal vectors
  for (int row = 0; row < NV.rows(); row++) {
	  ofs << "vn";
	  for (int col = 0; col < NV.cols(); col++) {
		  ofs << " " << NV(row, col);
	  }

	  ofs << std::endl;
  }

  // faces
  assert(F.rows() == UF.rows() && F.rows() == NF.rows() && "Number of faces are not equal in three matrices\n");
  assert(F.cols() == UF.cols() && F.cols() == NF.cols() && "Number of cols are not equal in three matrices\n");

  for (int row = 0; row < F.rows(); row++) {
	  ofs << "f";
	  for (int col = 0; col < F.cols(); col++) {

		  //matrix index starts from 0, but .obj file's index start from 1
		  ofs << " " << F(row, col) + 1 << "/" << UF(row, col) + 1 << "/" << NF(row, col) + 1;

	  }

	  ofs << std::endl;
  }
  ////////////////////////////////////////////////////////////////////////////
  return true;
}
