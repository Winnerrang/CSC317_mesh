#include "vertex_triangle_adjacency.h"
#include <assert.h>
bool notIn(std::vector<int>& VF, int face);

void vertex_triangle_adjacency(
  const Eigen::MatrixXi & F,
  const int num_vertices,
  std::vector<std::vector<int> > & VF)
{
  VF.resize(num_vertices);
  ////////////////////////////////////////////////////////////////////////////
  // loop through each faces and add the faces to the corresponding vertex
  
  for (int f = 0; f < F.rows(); f++) {
	  assert(F(f, 0) >= 0 && F(f, 1) < VF.size() && "Vertex index out of boundary\n");
	  assert(F(f, 1) >= 0 && F(f, 1) < VF.size() && "Vertex index out of boundary\n");
	  assert(F(f, 2) >= 0 && F(f, 2) < VF.size() && "Vertex index out of boundary\n");

	  //check for duplicate
	  assert(notIn(VF[F(f, 0)], f));
	  assert(notIn(VF[F(f, 1)], f));
	  assert(notIn(VF[F(f, 2)], f));

	  VF[F(f, 0)].push_back(f);
	  VF[F(f, 1)].push_back(f);
	  VF[F(f, 2)].push_back(f);
  }


  ////////////////////////////////////////////////////////////////////////////
}

bool notIn(std::vector<int>& VF, int face) {
	for (auto i : VF) {
		if (i == face) return false;
	}

	return true;
}