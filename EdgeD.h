#ifndef EDGED_H
#define EDGED_H


struct EdgeD
{
	int edge_i;
	int edge_a_i_global;
	int edge_b_i_global;
	int edge_a_i_local;
	int edge_b_i_local;

	std::vector<Eigen::MatrixXd*> stiff_K_ijs;
	std::vector<int> K_belongs;
	Eigen::Vector3d Stress;

};


#endif