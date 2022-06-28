#ifndef SIMULATION_H
#define SIMULATION_H

#include <lclab2.h>

/*
	Problems:
	For some reason the accelerations/forces are incredibly massive...
	Masses have been changed to be really massive to try and counteract this
	effect but a real solution is necessary...
*/

struct Vertex {
	// Role is necessary to determine the mass of each node
	// M(bulk) = 10^-9 g
	// M(surface) = 5*10^-10 g
	// M(edge) = 2.5*10^-10 g
	// M(corner) = 1.25*10^-10 g

	enum Role { Bulk = BIT(0), Surface = BIT(1), Glass = BIT(2), Top = BIT(3), Side = BIT(4), Edge = BIT(4), Corner = BIT(5) };
	Eigen::Vector3d position, director, velocity, acceleration, force, velocity_next;
	LC::scalar mass;
	Role role;
};

typedef std::array<unsigned int, 8> VolumeElement;
typedef Eigen::Vector3d Force;

class Simulation {
	
	struct ParameterList {
		// Bulk modulus (Pa)
		LC::scalar kappa = 28e6;
		// Shear modulus (Pa)
		LC::scalar mu = 570e3;
		// Limiting extension of polymer chains
		int I_m = 30;
		// Estimated change in the order parameter from nematic to isotropic (Try 0.533 as well)
		LC::scalar UdS = 80e3; // Pascals
		// Step time (seconds)
		LC::scalar dT = 1e-7;
		bool PBCs = false;

		// Initial cell dimensions
		std::array<LC::scalar, 3> cell_dims;
		std::array<int, 3> voxels;

		// Corner mass (kg)
		LC::scalar cornerMass = 1.25e-13;
		// Edge mass (kg)
		LC::scalar edgeMass = 2.5e-13;
		// Surface mass (kg)
		LC::scalar surfaceMass = 5.0e-13;
		// Bulk mass (kg)
		LC::scalar bulkMass = 1.0e-12;

		// Light direction (coming from top of cell)
		Eigen::Vector3d lightDirection{ 0., 0., -1. };
		// Material axis
		Eigen::Vector3d materialAxis{ 0., 0., 1. };
		// Director axis (eventually replace with dynamic director field)
		Eigen::Vector3d directorAxis{ 0., 1., 1. };
	};

public:

	/*
		Initialization of vertices and elements from structured director data
		
		Calling with no parameters initializes with 10x10x10 voxels, 10um x 10um x 10um cell dimensions
		with the director field aligned with (0,0,1)
	*/
	void InitVerticesAndElements(const LC::scalar* directors = 0,
		const std::array<int, 3>& voxels = { 16,16,4 },
		const std::array<LC::scalar, 3>& cell_dims = {1.6e-6,1.6e-6,0.4e-6});

	/*
		Routines
	*/

	// Update next velocity and acceleration for elemIdx
	void updateDynamicsAtElement(const unsigned int& elemIdx);

	// Update nodes
	void updateNodes();

	void setNextDynamicsToZero();

	/*
		Access data
	*/
	Vertex* getNodes();
	unsigned int getNumNodes();
	const VolumeElement* getElements();
	unsigned int getNumElements();
	ParameterList& Parameters();

private:
	ParameterList m_params;

	std::unique_ptr<Vertex[]> m_nodes;
	unsigned int m_nNodes;

	std::unique_ptr<VolumeElement[]> m_elements;
	unsigned int m_nElements;
};

void Simulation::InitVerticesAndElements(const LC::scalar* directors,
	const std::array<int, 3>& voxels,
	const std::array<LC::scalar, 3>& cell_dims) {

	std::array<int, 3> evoxels = voxels;
	m_params.cell_dims = cell_dims;
	m_params.voxels = voxels;

	bool PBC = m_params.PBCs;

	// periodic boundary conditions
	for (int i = 0 + 2 * PBC; i < 3; i++)
		--evoxels[i];


	m_nNodes = voxels[0] * voxels[1] * voxels[2];
	m_nElements = evoxels[0] * evoxels[1] * evoxels[2];

	if (!m_nNodes) return;

	m_nodes = std::unique_ptr<Vertex[]>(new Vertex[m_nNodes]);
	m_elements = std::unique_ptr<VolumeElement[]>(new VolumeElement[m_nElements]);

	auto index = [voxels, PBC](int i, int j, int k) {
		if (PBC) {
			i = i >= voxels[0] ? i - voxels[0] : i < 0 ? voxels[0] - 1 - i : i;
			j = j >= voxels[1] ? j - voxels[1] : j < 0 ? voxels[1] - 1 - j : j;
			k = k >= voxels[2] ? k - voxels[2] : k < 0 ? voxels[2] - 1 - k : k;
		}
		return (unsigned int)(i + j * voxels[0] + k * voxels[0] * voxels[1]);
	};

	auto eindex = [evoxels, voxels](int i, int j, int k) {
		return (unsigned int)(i + j * evoxels[0] + k * evoxels[0] * evoxels[1]);
	};

	std::array<LC::scalar, 3> dr = { cell_dims[0] / (voxels[0] - 1 + PBC),
		cell_dims[1] / (voxels[1] - 1 + PBC),
		cell_dims[2] / (voxels[2] - 1) };


	// Initialize nodes and elements
	for (int i = 0; i < voxels[0]; i++) {
		for (int j = 0; j < voxels[1]; j++) {
			for (int k = 0; k < voxels[2]; k++) {

				unsigned int idx = index(i, j, k);

				m_nodes[idx].acceleration = Eigen::Vector3d{ 0.,0.,0. };
				m_nodes[idx].velocity = Eigen::Vector3d{ 0.,0.,0. };
				m_nodes[idx].velocity_next = Eigen::Vector3d{ 0.,0.,0. };

				unsigned int elem_idx = eindex(i, j, k);
				std::array<int, 3> ri = { i, j, k };

				m_nodes[idx].role = Vertex::Role::Bulk;
				m_nodes[idx].mass = m_params.bulkMass;

				for (int d = 0; d < 3; d++) {

					// Define the position of the node
					if (d < 2)
						m_nodes[idx].position[d] = -0.5 * cell_dims[d] + dr[d] * (ri[d]+ 0.5 * PBC);
					else
						m_nodes[idx].position[d] = -0.5 * cell_dims[d] + dr[d] * ri[d];
					// Classify the node as surface node
					if (ri[d] % evoxels[d] == 0) {
						if (PBC && d < 2) {
							// lateral xy dims go on forever hence bulk nodes
						}
						else {
							// Note that if side node + PBC is also on the top or bottom, it is still a surface node
							m_nodes[idx].role = Vertex::Role::Surface;
							m_nodes[idx].mass = m_params.surfaceMass;
						}
					}

				}

				if (m_nodes[idx].role == Vertex::Role::Surface) {
					// Check if top or glass node
					if (k == 0)
						m_nodes[idx].role = Vertex::Role(m_nodes[idx].role | Vertex::Role::Glass);
					else if (k == evoxels[2])
						m_nodes[idx].role = Vertex::Role(m_nodes[idx].role | Vertex::Role::Top);

					if ((i % evoxels[0] == 0 || j % evoxels[1] == 0) && !PBC)
						m_nodes[idx].role = Vertex::Role(m_nodes[idx].role | Vertex::Role::Side);

					// Determine if edge
					// (0, 1) -> xy
					// (0, 2) -> xz
					// (1, 2) -> yz

					// There can only be edges and corners without PBCs
					if (!PBC) {

						for (int a = 0; a < 3; a++) {
							for (int b = a + 1; b < 3; b++) {
								if (ri[a] % evoxels[a] == 0 && ri[b] & evoxels[b] == 0) {
									m_nodes[idx].role = Vertex::Role(m_nodes[idx].role | Vertex::Role::Edge);
									m_nodes[idx].mass = m_params.edgeMass;
								}

							}
						}

						// Determine if corner
						for (int a = 0; a < 2; a++) {
							for (int b = 0; b < 2; b++) {
								for (int c = 0; c < 2; c++) {
									if (i == a * evoxels[0] && j == b * evoxels[1] && k == c * evoxels[2]) {
										m_nodes[idx].role = Vertex::Role(m_nodes[idx].role | Vertex::Role::Corner);
										m_nodes[idx].mass = m_params.cornerMass;
									}
								}
							}
						}
					}

				}

				if (directors) {
					m_nodes[idx].director[0] = directors[idx];
					m_nodes[idx].director[1] = directors[idx + m_nNodes];
					m_nodes[idx].director[2] = directors[idx + 2 * m_nNodes];
				}
				else {
					m_params.directorAxis.normalize();
					m_nodes[idx].director = m_params.directorAxis;
				}

				// Initialize the vertex indices defining each element

				// 0 -> (0,0,0)
				// 1 -> (1,0,0)
				// 2 -> (0,1,0)
				// 3 -> (1,1,0)
				// 4 -> (0,0,1)
				// 5 -> (1,0,1)
				// 6 -> (0,1,1)
				// 7 -> (1,1,1)
				if (i < evoxels[0] && j < evoxels[1] && k < evoxels[2]) {
					m_elements[elem_idx][0] = index(i, j, k);
					m_elements[elem_idx][1] = index(i + 1, j, k);
					m_elements[elem_idx][2] = index(i, j + 1, k);
					m_elements[elem_idx][3] = index(i + 1, j + 1, k);
					m_elements[elem_idx][4] = index(i, j, k + 1);
					m_elements[elem_idx][5] = index(i + 1, j, k + 1);
					m_elements[elem_idx][6] = index(i, j + 1, k + 1);
					m_elements[elem_idx][7] = index(i + 1, j + 1, k + 1);
				}
			}
		}
	}

}

void Simulation::setNextDynamicsToZero() {
	for (unsigned int idx = 0; idx < m_nNodes; idx++) {
		m_nodes[idx].force = Eigen::Vector3d{ 0.,0.,0. };
		m_nodes[idx].velocity_next = Eigen::Vector3d{ 0.,0.,0. };
	}
}

/*
	Compute the force at the volume element index listed
*/
void Simulation::updateDynamicsAtElement(const unsigned int& elemIdx) {
	/*
		Steps to compute the force
		1. Compute deformation gradient tensor (lambda)
		2. Compute gradient terms
		3. Compute total force
	*/

	// Deformation gradient tensor
	std::array<Eigen::Vector3d, 8> r_i;
	LC::scalar I_m = m_params.I_m;
	bool PBC = m_params.PBCs;

	// Fill r_i from surrounding vertices
	for (int i = 0; i < 8; i++) {
		unsigned int sub_i = m_elements[elemIdx][i];
		unsigned int sub_0 = m_elements[elemIdx][0];
		// m_elements[elemIdx][i] ranges over local indices of elemIdx
		// defined during initialization of mesh
		// Position defined in voxel coordinate space (dimensionless)
		if (PBC) {
			// Translate to the opposite side of the cell if on one edge
			for (int d = 0; d < 2; d++) {
				LC::scalar dr = m_nodes[sub_i].position[d] - m_nodes[sub_0].position[d];
				if (dr > 0.5 * m_params.cell_dims[d]) {
					dr = dr - m_params.cell_dims[d];
				}
				if (dr < -0.5 * m_params.cell_dims[d]) {
					dr = dr + m_params.cell_dims[d];
				}
				r_i[i][d] = 2. * (m_nodes[sub_0].position[d] + dr);
			}
			// No PBCs in z-direction
			r_i[i][2] = 2. * m_nodes[sub_i].position[2];
		}
		else
		// Note each element has a volume of 8 integrated over (a,b,g) space -> Jacobian transformation
		// requires multiplication by 2 to renormalize volume
			r_i[i] = 2. * m_nodes[sub_i].position;
	}

	std::array<LC::scalar, 3> dr = { m_params.cell_dims[0] / (m_params.voxels[0] - 1 + PBC),
			m_params.cell_dims[1] / (m_params.voxels[1] - 1 + PBC),
			m_params.cell_dims[2] / (m_params.voxels[2] - 1) };

	std::array<LC::scalar, 3> ra, rb, rg, rab, rbg, rag, rabg;

	for (int d = 0; d < 3; d++) {
		// This is fixed by multiplying each component to be twice the size that is ordinarily expected

		// Each component is rescaled to be dimensionless
		ra[d] = 1. / 8. * (-r_i[0][d] + r_i[1][d] - r_i[2][d] + r_i[3][d] - r_i[4][d] + r_i[5][d] - r_i[6][d] + r_i[7][d]) / dr[d];
		rb[d] = 1. / 8. * (-r_i[0][d] - r_i[1][d] + r_i[2][d] + r_i[3][d] - r_i[4][d] - r_i[5][d] + r_i[6][d] + r_i[7][d]) / dr[d];
		rg[d] = 1. / 8. * (-r_i[0][d] - r_i[1][d] - r_i[2][d] - r_i[3][d] + r_i[4][d] + r_i[5][d] + r_i[6][d] + r_i[7][d]) / dr[d];
		rab[d] = 1. / 8. * (r_i[0][d] - r_i[1][d] - r_i[2][d] + r_i[3][d] + r_i[4][d] - r_i[5][d] - r_i[6][d] + r_i[7][d]) / dr[d];
		rag[d] = 1. / 8. * (r_i[0][d] - r_i[1][d] + r_i[2][d] - r_i[3][d] - r_i[4][d] + r_i[5][d] - r_i[6][d] + r_i[7][d]) / dr[d];
		rbg[d] = 1. / 8. * (r_i[0][d] + r_i[1][d] - r_i[2][d] - r_i[3][d] - r_i[4][d] - r_i[5][d] + r_i[6][d] + r_i[7][d]) / dr[d];
		rabg[d] = 1. / 8. * (-r_i[0][d] + r_i[1][d] + r_i[2][d] - r_i[3][d] + r_i[4][d] - r_i[5][d] - r_i[6][d] + r_i[7][d]) / dr[d];
	}

	auto drda = [ra, rab, rag, rabg](LC::scalar b, LC::scalar g) {
		Eigen::Vector3d result;
		for (int d = 0; d < 3; d++) {
			result[d] = ra[d] + b * rab[d] + g * rag[d] + g * b * rabg[d];
		}
		return result;
	};
	auto drdb = [rb, rab, rbg, rabg](LC::scalar a, LC::scalar g) {
		Eigen::Vector3d result;
		for (int d = 0; d < 3; d++) {
			result[d] = rb[d] + a * rab[d] + g * rbg[d] + g * a * rabg[d];
		}
		return result;
	};
	auto drdg = [rg, rag, rbg, rabg](LC::scalar a, LC::scalar b) {
		Eigen::Vector3d result;
		for (int d = 0; d < 3; d++) {
			result[d] = rg[d] + a * rag[d] + b * rbg[d] + b * a * rabg[d];
		}
		return result;
	};

	// Compute the invariants of the deformation gradient tensor (Gaussian integration)
	LC::scalar gInt = 1./ sqrt(3.);
	LC::scalar a, b, g, bg, ag, ab;

	Eigen::Vector3d gradI, gradJ, gradQeps;
	// Iterate over the eight grid points (+-gInt, +-gInt, +-gInt) computing the force

	for (int idx = 0; idx< 8; idx++) {

		// Current node in element elemIdx
		unsigned int sub_i = m_elements[elemIdx][idx];
				
		Force force = { 0., 0., 0. };

		// 0 -> (-1,-1,-1)
		// 1 -> (1,-1,-1)
		// 2 -> (-1,1,-1)
		// 3 -> (1,1,-1)
		// 4 -> (-1,-1,1)
		// 5 -> (1,-1,1)
		// 6 -> (-1,1,1)
		// 7 -> (1,1,1)

		if (idx == 0) {
			a = -gInt;
			b = -gInt;
			g = -gInt;
			bg = -(1 - b) * (1 - g);
			ag = -(1 - a) * (1 - g);
			ab = -(1 - a) * (1 - b);
		}
		else if (idx == 1) {
			a = gInt;
			b = -gInt;
			g = -gInt;
			bg = (1 - b) * (1 - g);
			ag = -(1 + a) * (1 - g);
			ab = -(1 + a) * (1 - b);
		}
		else if (idx == 2) {
			a = -gInt;
			b = gInt;
			g = -gInt;
			bg = -(1 + b) * (1 - g);
			ag = (1 - a) * (1 - g);
			ab = -(1 - a) * (1 + b);
		}
		else if (idx == 3) {
			a = gInt;
			b = gInt;
			g = -gInt;
			bg = (1 + b) * (1 - g);
			ag = (1 + a) * (1 - g);
			ab = -(1 + a) * (1 + b);
		}
		else if (idx == 4) {
			a = -gInt;
			b = -gInt;
			g = gInt;
			bg = -(1 - b) * (1 + g);
			ag = -(1 - a) * (1 + g);
			ab = (1 - a) * (1 - b);
		}
		else if (idx == 5) {
			a = gInt;
			b = -gInt;
			g = gInt;
			bg = (1 - b) * (1 + g);
			ag = -(1 + a) * (1 + g);
			ab = (1 + a) * (1 - b);
		}
		else if (idx == 6) {
			a = -gInt;
			b = gInt;
			g = gInt;
			bg = -(1 + b) * (1 + g);
			ag = (1 - a) * (1 + g);
			ab = (1 - a) * (1 + b);
		}
		else if (idx == 7) {
			a = gInt;
			b = gInt;
			g = gInt;
			bg = (1 + b) * (1 + g);
			ag = (1 + a) * (1 + g);
			ab = (1 + a) * (1 + b);
		}

		// Deformation gradient tensor terms

		auto drda_eval = drda(b, g);
		auto drdb_eval = drdb(a, g);
		auto drdg_eval = drdg(a, b);

		LC::scalar J = drda_eval.dot(drdb_eval.cross(drdg_eval));
		LC::scalar Jmax = pow(1.2, 3);

		//LC_INFO("J = {0}", J);
				
		// Do not update if the volume elem is negative or I1 greater than I_m
		// Evalate terms needed for force at each node
		LC::scalar I1 = drda_eval.squaredNorm() + drdb_eval.squaredNorm() + drdg_eval.squaredNorm();
		if (J < 0. || I1 >= I_m) {
			continue;
		}

		LC::scalar J1_3 = pow(J, 1. / 3.);
		LC::scalar J2_3 = J1_3 * J1_3;
		LC::scalar J5_3 = J2_3 * J2_3 * J1_3;

		// Derive strain tensor
		Eigen::Matrix3d lambda;
		for (int d = 0; d < 3; d++) {
			lambda(d, 0) = drda_eval[d];
			lambda(d, 1) = drdb_eval[d];
			lambda(d, 2) = drdg_eval[d];
		}

		// Right Cauchy-Green strain tensor
		Eigen::Matrix3d eps = lambda.transpose() * lambda;

		//LC_INFO("J = {0}", J);

		// Defined by pnas.1811823115
		LC::scalar UdS = m_params.UdS; // Pascals

		Eigen::Vector3d n = m_nodes[sub_i].director;

		gradI = 0.25 * (bg * drda_eval + ag * drdb_eval + ab * drdg_eval);
		gradJ = 1./8. * (bg * drdb_eval.cross(drdg_eval) + ag * drdg_eval.cross(drda_eval) + ab * drda_eval.cross(drdb_eval));
		gradQeps = 0.25 * (n[0] * bg + n[1] * ag + n[2] * ab) * (n[0] * drda_eval + n[1] * drdb_eval + n[2] * drdg_eval);

		// Compute the force

		LC::scalar n_an_beps_ab = 0.0;
		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < 3; jj++)
				n_an_beps_ab += n[ii] * n[jj] * eps(ii, jj);

		// Material resistance to shear deformations
		Force term1 = -0.5 * m_params.mu * (I_m - 3) / (I_m - I1) * gradI;
		// Material resistance to changes in volume
		Force term2 = -m_params.kappa * log(J) / J * gradJ;
		// Coupling of material strain to nematic ordering
		Force term3 = UdS * (1. / J2_3 * gradQeps - gradI / 3. / J2_3 - 2. * (n_an_beps_ab - I1 / 3.) * gradJ / 3. / J5_3);

		force = 
			term1+ 
			term2+
			term3
			;

		// Eliminate invalid results
		if (force != force)
			force = { 0., 0., 0. };

		// To convert to actual force, multiply by dr^2
		for (int d = 0; d < 3; d++) {
			force[d] *= dr[d]* dr[d];
		}

		m_nodes[sub_i].force += force;

	}

}

void Simulation::updateNodes() {
	// Set next velocity, acceleration at node to zero
	setNextDynamicsToZero();

	// Compute acceleration/velocity contributions over all nodes using elements
	for (unsigned int i = 0; i < m_nElements; i++) {
		updateDynamicsAtElement(i);
	}

	for (unsigned int i = 0; i < m_nNodes; i++) {

		// Compute the next position (only if not glass)
		// ***Note this can be made more complex by updating the x,y components
		// of the glass position and fixing the z position.


		// Set next acceleration, velocity
		// Compute the velocity with damping param eta = 0.1
		LC::scalar eta = 0.1;
		Eigen::Vector3d accel = m_nodes[i].force / m_nodes[i].mass;

		//LC_INFO("Acceleration = {0}", accel.norm());

		m_nodes[i].velocity_next = (1 - eta) * m_nodes[i].velocity + 0.5 * (m_nodes[i].acceleration + accel) * m_params.dT;

		m_nodes[i].velocity = m_nodes[i].velocity_next;

		auto newPosition =  m_nodes[i].position + (m_nodes[i].velocity * m_params.dT
			+ accel * pow(m_params.dT, 2));

		// Don't update lateral components if it's a side surface
		//if (!(m_nodes[i].role & Vertex::Role::Side)) {
		m_nodes[i].position[0] = newPosition[0];
		m_nodes[i].position[1] = newPosition[1];
		//}

		// Don't update the z position of glass nodes
		if (!(m_nodes[i].role & Vertex::Role::Glass))
			m_nodes[i].position[2] = newPosition[2];
		
			// Don't let nodes pass through bottom of glass
		if (m_nodes[i].position[2] < -0.5 * m_params.cell_dims[2])
			m_nodes[i].position[2] = -0.5 * m_params.cell_dims[2];

		// Don't allow node to pass through sides of cell
		//for (int d = 0; d < 2; d++) {
		//	if (m_nodes[i].position[d] < -0.5 * m_params.voxels[d])
		//		m_nodes[i].position[d] = -0.5 * m_params.voxels[d];
		//	if (m_nodes[i].position[d] > 0.5 * m_params.voxels[d])
		//		m_nodes[i].position[d] = 0.5 * m_params.voxels[d];
		//}

		
	}

}

// Simple functions to access data
Vertex* Simulation::getNodes() {
	return m_nodes.get();
}
unsigned int Simulation::getNumNodes() {
	return m_nNodes;
}
const VolumeElement* Simulation::getElements() {
	return m_elements.get();
}
unsigned int Simulation::getNumElements() {
	return m_nElements;
}

Simulation::ParameterList& Simulation::Parameters() {
	return m_params;
}

#endif