#include "per_face_normals.h"
#include "triangle_area_normal.h"

void per_face_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & N)
{
  ////////////////////////////////////////////////////////////////////////////
  // Replace with your code:

	N.resize(F.rows(), 3);
	for (int f = 0; f < F.rows(); f++) {
		auto normal = triangle_area_normal(V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)));
		N.row(f) = normal.normalized();
	}

  ////////////////////////////////////////////////////////////////////////////
}
