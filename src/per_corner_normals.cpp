#include "per_corner_normals.h"
#include "triangle_area_normal.h"
// Hint:
#include "vertex_triangle_adjacency.h"
#include <iostream>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>

void per_corner_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double corner_threshold,
  Eigen::MatrixXd & N)
{
  N = Eigen::MatrixXd::Zero(F.rows()*3,3);
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here:
  
  // create data structure for vertex -> adjacent faces
  std::vector<std::vector<int>> VF;
  vertex_triangle_adjacency(F, V.rows(), VF);

  for (int f = 0; f < F.rows(); f++) {
	  auto refNormal = triangle_area_normal(V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)));
	  refNormal.normalize();
	  
	  // calculate the corner normal for each vertex of the triangle
	  for (int v = 0; v < F.cols(); v++) {
		  int vertexIndex = F(f, v);

		  assert(vertexIndex >= 0 && vertexIndex < VF.size() && "Your vertex_triangle_adjacency is wrong!");
		  
		  for (int i = 0; i < VF[vertexIndex].size(); i++) {
			  int adjFace = VF[vertexIndex][i];
			  
			  
			  auto adjAreaNormal = triangle_area_normal(V.row(F(adjFace, 0)), V.row(F(adjFace, 1)), V.row(F(adjFace, 2)));
			  auto unitAdjNormal = adjAreaNormal.normalized();

			  if (refNormal.dot(unitAdjNormal) > std::cos(corner_threshold / 180 * M_PI)) {
				  N.row(f * 3 + v) += adjAreaNormal;
			  }
		  }

	  }
  }

  for (int row = 0; row < N.rows(); row++) {

	  assert(N.row(row) != Eigen::RowVector3d::Zero());

	  N.row(row).normalize();
  }
  ////////////////////////////////////////////////////////////////////////////
}
