#include "per_vertex_normals.h"
#include "triangle_area_normal.h"

void per_vertex_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & N)
{
  N = Eigen::MatrixXd::Zero(V.rows(),3);
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here:
  
  // find the total area normal for each vertex
  for (int f = 0; f < F.rows(); f++) {
	  auto normal = triangle_area_normal(V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)));

	  // add area normal to these three vertices of the triangle
	  N.row(F(f, 0)) += normal;
	  N.row(F(f, 1)) += normal;
	  N.row(F(f, 2)) += normal;
  }

  // normalize every normal for each vertex
  for (int v = 0; v < N.rows(); v++) {
	  assert(N.row(v) != Eigen::RowVector3d::Zero() && "One of the vertex is 0 vector");

	  N.row(v).normalize();
  }
  ////////////////////////////////////////////////////////////////////////////
}
