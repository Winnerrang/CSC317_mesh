#include "catmull_clark.h"
#include <unordered_map>
#include <utility>
#include <functional>
#include <assert.h>
#include <vector>
#include <iostream>

typedef std::pair<int, int> pair;

//ref: https://stackoverflow.com/questions/32685540/why-cant-i-compile-an-unordered-map-with-a-pair-as-key
struct pair_hash {
	template <class T1, class T2>
	std::size_t operator () (const std::pair<T1, T2>& p) const {
		auto h1 = std::hash<T1>{}(p.first);
		auto h2 = std::hash<T2>{}(p.second);
		
		if (h1 != h2) {
			return h1 ^ h2;
		}
		else {
			return h1;
		}
		
	}
};

void generateFPandEdge(Eigen::MatrixXd& SV, 
	Eigen::MatrixXi& SF, 
	Eigen::VectorXi& FFP, 
	Eigen::MatrixXi& EV, 
	std::unordered_map<pair, int, pair_hash>& VE,
	std::vector<std::vector<int>>& VNF, 
	std::vector<std::vector<int>>& ENF);

void generateEPandVNE(Eigen::MatrixXd& SV,
	Eigen::MatrixXi& SF,
	Eigen::MatrixXi& EV, 
	std::vector<std::vector<int>>& VNE, 
	Eigen::VectorXi& FFP, 
	Eigen::MatrixXi& FEP,
	std::vector<std::vector<int>>& ENF);

void updateOriginalPoint(Eigen::MatrixXd& SV,
	int origNumVertex, 
	Eigen::MatrixXi& EV,
	std::vector<std::vector<int>>& VNE, 
	std::vector<std::vector<int>>& VNF, 
	Eigen::VectorXi& FFP);

void generateFaces(Eigen::MatrixXi& SF,
	Eigen::VectorXi& FFP, 
	Eigen::MatrixXi& FEP, 
	Eigen::MatrixXi& newFace);

void catmull_clark(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int num_iters,
  Eigen::MatrixXd & SV,
  Eigen::MatrixXi & SF)
{
  ////////////////////////////////////////////////////////////////////////////
  // Replace with your code here:
	
	SV = V;
	SF = F;

	for (int i = 0; i < num_iters; i++) {
		
		assert(SF.cols() == 4 && "It is not a quad surface");
		int origNumVertex = SV.rows();
		auto FFP = Eigen::VectorXi(SF.rows());						// face to face point
		auto FEP = Eigen::MatrixXi(SF.rows(), SF.cols());			// face to edge point, same order as SF
		auto EV = Eigen::MatrixXi();								// edges to vertex
		auto VE = std::unordered_map<pair, int, pair_hash>();		// pair of vertex to edge
		auto VNF = std::vector<std::vector<int>>(SV.rows());		// near faces for each vertex
		auto VNE = std::vector<std::vector<int>>(SV.rows());		// near edges for each vertex
		auto ENF = std::vector<std::vector<int>>();					// near faces for each edge
		// create new face point
		generateFPandEdge(SV, SF, FFP, EV, VE, VNF, ENF);

		// create new edge point
		generateEPandVNE(SV, SF, EV, VNE, FFP, FEP, ENF);
		
		
		// update original vertex point
		updateOriginalPoint(SV, origNumVertex, EV, VNE, VNF, FFP);
		
		// generate new faces
		Eigen::MatrixXi newFace(SF.rows() * 4, 4);
		generateFaces(SF, FFP, FEP, newFace);
		SF = newFace;
	}
	
  
  ////////////////////////////////////////////////////////////////////////////
}


void generateFPandEdge(Eigen::MatrixXd& SV,
	Eigen::MatrixXi& SF,
	Eigen::VectorXi& FFP,
	Eigen::MatrixXi& EV,
	std::unordered_map<pair, int, pair_hash>& VE,
	std::vector<std::vector<int>>& VNF,
	std::vector<std::vector<int>>& ENF) {

	
	// using the euler formula to find the number of edges
	// ignore nesis for now, resize if needed
	EV.resize(SV.rows() + SF.rows() - 2, 2);

	ENF.resize(EV.rows());

	int newEdgeIndex = 0;


	int newVertexIndex = SV.rows();
	SV.conservativeResize(SF.rows() + SV.rows(), SV.cols());

	
	for (int f = 0; f < SF.rows(); f++) {
		Eigen::RowVector3d facePoint = Eigen::RowVector3d::Zero();

		for (int v = 0; v < SF.cols(); v++) {

			facePoint += SV.row(SF(f, v));

			// two point of edge
			int start = SF(f, v);
			int end = SF(f, (v + 1) % 4);

			auto vertexPair = pair(start, end);
			auto iter = VE.find(vertexPair);

			// find the reverse order
			if (iter == VE.end()) {
				vertexPair = pair(end, start);
				iter = VE.find(vertexPair);
			}
			

			int edgeIndex;
			if (iter == VE.end()) {

				assert(newEdgeIndex <= EV.rows());

				// there is hole in the surface
				if (newEdgeIndex == EV.rows()) {
					EV.conservativeResize(EV.rows() + 2, EV.cols());
					ENF.resize(EV.rows());
				}

				EV.row(newEdgeIndex) = Eigen::RowVector2i(start, end);
				VE[vertexPair] = newEdgeIndex;
				edgeIndex = newEdgeIndex;
				newEdgeIndex++;
			}
			else {
				edgeIndex = iter->second;
			}
			
			VNF[start].push_back(f);
			ENF[edgeIndex].push_back(f);
			
		}

		facePoint /= 4;
		assert(newVertexIndex < SV.rows() && "Out of boundary in vertex matrix");
		SV.row(newVertexIndex) = facePoint;
		FFP(f) = newVertexIndex;
		newVertexIndex++;
	}

	assert(newEdgeIndex == EV.rows() && "number of edges is too small");
	assert(newVertexIndex == SV.rows() && "number of new face point is too small");
}


void generateEPandVNE(Eigen::MatrixXd& SV,
	Eigen::MatrixXi& SF,
	Eigen::MatrixXi& EV,
	std::vector<std::vector<int>>& VNE,
	Eigen::VectorXi& FFP,
	Eigen::MatrixXi& FEP,
	std::vector<std::vector<int>>& ENF) {


	int newVertexIndex = SV.rows();
	SV.conservativeResize(SV.rows() + EV.rows(), SV.cols());

	for (int eIndex = 0; eIndex < EV.rows(); eIndex++) {
		
		Eigen::RowVector3d edgePoint = Eigen::RowVector3d::Zero();
		
		for (int v = 0; v < EV.cols(); v++) {
			VNE[EV(eIndex, v)].push_back(eIndex);
			edgePoint += SV.row(EV(eIndex, v));
		}

		// the order of the edge in two near faces
		Eigen::Vector2i order;

		assert(ENF[eIndex].size() == 2);
		for (int f = 0; f < ENF[eIndex].size(); f++) {
			int nFace = ENF[eIndex][f];
			edgePoint += SV.row(FFP(nFace));

			for (int v = 0; v < SF.cols(); v++) {
				int vIndex = SF(nFace, v);
				int nvIndex = SF(nFace, (v + 1) % 4);
				if ((vIndex == EV(eIndex, 0) && nvIndex == EV(eIndex, 1)) || (vIndex == EV(eIndex, 1) && nvIndex == EV(eIndex, 0))) {
					order[f] = v;
				}
			}

		}

		edgePoint /= 4;
		assert(newVertexIndex < SV.rows() && "Out of boundary in vertex matrix");
		SV.row(newVertexIndex) = edgePoint;
		FEP(ENF[eIndex][0], order[0]) = newVertexIndex;
		FEP(ENF[eIndex][1], order[1]) = newVertexIndex;
		newVertexIndex++;

		
	}
}

void updateOriginalPoint(Eigen::MatrixXd& SV,
	int origNumVertex,
	Eigen::MatrixXi& EV,
	std::vector<std::vector<int>>& VNE,
	std::vector<std::vector<int>>& VNF,
	Eigen::VectorXi& FFP) {

	for (int v = 0; v < origNumVertex; v++) {
		
		Eigen::RowVector3d edgeTotal(Eigen::RowVector3d::Zero()), faceTotal(Eigen::RowVector3d::Zero()), vertexTotal(Eigen::RowVector3d::Zero());
		
		for (auto e : VNE[v]) {
			auto start = EV(e ,0);
			auto end = EV(e, 1);
			edgeTotal += (SV.row(start) + SV.row(end)) / 2;
		}

		for (auto f : VNF[v]) {
			faceTotal += SV.row(FFP(f));
		}

		SV.row(v) = (faceTotal / VNF[v].size() + 2 * edgeTotal / VNE[v].size() + (VNF[v].size() - 3) * SV.row(v)) / VNF[v].size();

	}

}

void generateFaces(Eigen::MatrixXi& SF,
	Eigen::VectorXi& FFP,
	Eigen::MatrixXi& FEP,
	Eigen::MatrixXi& newFace) {
	newFace.resize(SF.rows() * 4, SF.cols());

	int newFaceIndex = 0;
	for (int f = 0; f < SF.rows(); f++) {
		for (int v = 0; v < SF.cols(); v++) {
 			newFace.row(newFaceIndex) = Eigen::RowVector4i(SF(f, v), FEP(f, v), FFP(f), FEP(f, (v + 3) % 4));
			newFaceIndex++;
		}
	}

	assert(newFaceIndex == newFace.rows());
}