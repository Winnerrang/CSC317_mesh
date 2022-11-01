#include "sphere.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>

void sphere(
  const int num_faces_u,
  const int num_faces_v,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & UV,
  Eigen::MatrixXi & UF,
  Eigen::MatrixXd & NV,
  Eigen::MatrixXi & NF)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here:
	
	const double r = 1;
	
	V.resize(num_faces_u * (num_faces_v - 1) + 2, 3);
	F.resize(num_faces_u * num_faces_v, 4);

	UV.resize((num_faces_u + 1) * (num_faces_v + 1), 2);
	UF.resize(num_faces_u * num_faces_v, 4);
	NV.resize(num_faces_u * (num_faces_v - 1) + 2, 3);
	NF.resize(num_faces_u * num_faces_v, 4);

	// stored the vertice' index location in F, UF and NF
	// will be used for generate the connectivity matrix

	auto vertexIndices = Eigen::MatrixXi(num_faces_v + 1, num_faces_u);
	int row = 0;
	int textureRow = 0;

	// the reason v is <= is because the end point is not the start point
	// but for u, 0 degree is equal to 360 degree, thus < is used
	for (int vStep = 0; vStep <= num_faces_v; vStep++) {
		
		bool cornerCase = false;
		for (int uStep = 0; uStep < num_faces_u; uStep++) {
			if (cornerCase) break;
			double x, y, z, phi, theta;


			// find the location of the points
			if (vStep == 0) {
				z = 0;
				x = 0;
				y = -1;
				cornerCase = true;
			}
			else if (vStep == num_faces_v) {
				z = 0;
				x = 0;
				y = 1;
				cornerCase = true;
			}
			else {
				// start from the bottom
				phi = M_PI - (M_PI / num_faces_v) * vStep;
				theta = (2 * M_PI / num_faces_u) * uStep;

				
				z = r * std::sin(phi) * std::cos(theta);
				x = r * std::sin(phi) * std::sin(theta);
				y = r * std::cos(phi);
			}

			auto normal = Eigen::Vector3d(x, y, z).normalized();
			
			double u = ((double)uStep) / num_faces_u;
			double v = ((double)vStep) / num_faces_v;

			assert(row < V.rows() && "Row out of boundaries\n");
			V(row, 0) = x;
			V(row, 1) = y;
			V(row, 2) = z;

			UV(row, 0) = u;
			UV(row, 1) = v;

			NV(row, 0) = normal.x();
			NV(row, 1) = normal.y();
			NV(row, 2) = normal.z();

			if (vStep == 0 || vStep == num_faces_v) {
				vertexIndices.row(vStep) = Eigen::RowVectorXi::Ones(num_faces_u);
				vertexIndices.row(vStep) *= row;
			}
			else {
				vertexIndices(vStep, uStep) = row;
			}
			row++;
		}
	}

	row = 0;
	auto textureIndices = Eigen::MatrixXi(num_faces_v + 1, num_faces_u + 1);
	for (int vStep = 0; vStep <= num_faces_v; vStep++) {

		for (int uStep = 0; uStep <= num_faces_u; uStep++) {

			double u = ((double)uStep) / num_faces_u;
			double v = ((double)vStep) / num_faces_v;

			assert(row < UV.rows() && "Row out of boundaries\n");


			UV(row, 0) = u;
			UV(row, 1) = v;


			textureIndices(vStep, uStep) = row;

			row++;
		}
	}


	row = 0;
	for (int r = 0; r < vertexIndices.rows() - 1; r++) {
		for (int c = 0; c < vertexIndices.cols(); c++) {
			auto indices = Eigen::RowVector4i();
			indices(0) = vertexIndices(r, c);
			indices(1) = vertexIndices(r, (c + 1) % vertexIndices.cols());
			indices(2) = vertexIndices(r + 1, (c + 1) % vertexIndices.cols());
			indices(3) = vertexIndices(r + 1, c);



			assert(row < F.rows() && "Row out of boundary for F\n");
			F.row(row) = indices;
			NF.row(row) = indices;
			row++;
		}
	}

	row = 0;
	for (int r = 0; r < textureIndices.rows() - 1; r++) {
		for (int c = 0; c < textureIndices.cols() - 1; c++) {
			auto indices = Eigen::RowVector4i();
			indices(0) = textureIndices(r, c);
			indices(1) = textureIndices(r, c + 1);
			indices(2) = textureIndices(r + 1, c + 1);
			indices(3) = textureIndices(r + 1, c);



			assert(row < UF.rows() && "Row out of boundary for F\n");
			UF.row(row) = indices;
			
			row++;
		}
	}


  ////////////////////////////////////////////////////////////////////////////
}
