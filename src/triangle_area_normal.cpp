#include "triangle_area_normal.h"
#include <Eigen/Geometry>
#include <iostream>
Eigen::RowVector3d triangle_area_normal(
  const Eigen::RowVector3d & a, 
  const Eigen::RowVector3d & b, 
  const Eigen::RowVector3d & c)
{
  ////////////////////////////////////////////////////////////////////////////
  // Replace with your code:
  ////////////////////////////////////////////////////////////////////////////


	return 0.5 * (b - a).cross(c - a);
}
