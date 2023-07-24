//==========================================================================//
//@author KatsuyaKikuchi
//@brief Implementation of the ParticleD.h
//==========================================================================//

#include "ParticleD.h"

//==========================================================================//
ParticleD::ParticleD(const Eigen::Vector3d& grid)
//各頂点の位置や速度,加速度,力のほか隣接点などの情報をもったParticleDクラス
	:grid(grid), exp_pos(grid), expago_pos(grid), prime_pos(grid), initial_pos(grid),
	deltax_in_model(Eigen::Vector3d::Zero()),
	deltaxago_in_model(Eigen::Vector3d::Zero()),
	force(Eigen::Vector3d::Zero()),
	displace(Eigen::Vector3d::Zero()),
	velocity(Eigen::Vector3d::Zero()),
	acceleration(Eigen::Vector3d::Zero()),
	fixed(false),
	drawgrid(grid.x() + 100.0, grid.y(), grid.z() + 200.0),
	imgrid(grid*cipher)
{
}
ParticleD::~ParticleD() {}
//==========================================================================//
//==========================================================================//
//	@start		   				ループ設定									//
//==========================================================================//
void ParticleD::Draw()const {
	DrawCircle(int(grid.x()), int(grid.y()), 3, RED, true);
}
void ParticleD::Update(const Eigen::Vector3d &next_grid) {
	acceleration = (((next_grid - this->grid) / TIME_STEP)- velocity) / TIME_STEP;
	velocity = (next_grid - this->grid) / TIME_STEP;
	this->grid = next_grid;
}
void ParticleD::Update_Grid(const Eigen::Vector3d &next_grid) {
	this->grid = next_grid;
}
void ParticleD::Update_Draw_Grid(const Eigen::Vector3d &next_grid) {
	this->drawgrid = next_grid;
}
//初期速度を入れる様
void ParticleD::Update_Velocity(const Eigen::Vector3d &init_velocity) {
	this->velocity = init_velocity;
}//==========================================================================//
//	@end		   				ループ設定									//
//==========================================================================//
//==========================================================================//
void ParticleD::Set_Force(const Eigen::Vector3d &f) {
	this->force = f;
}
void ParticleD::Set_Fixed(const bool& fixed) {
	this->fixed = fixed;
}
// 
void ParticleD::Set_Exp_Pos(const Eigen::Vector3d &pos) {
	this->exp_pos = pos;
}
//
void ParticleD::Set_ExpAgo_Pos(const Eigen::Vector3d &pos) {
	this->expago_pos = pos;
}
// 
void ParticleD::Set_Prime_Pos(const Eigen::Vector3d &pos) {
	this->prime_pos = pos;
}
//初期位置にした(OldFEMで使う)
void ParticleD::Set_Initial_Pos(const Eigen::Vector3d &pos) {
	this->initial_pos = pos;
}
void ParticleD::Set_Deltax_In_Model(const Eigen::Vector3d &pos) {
	this->deltax_in_model = pos;
}
void ParticleD::Set_Velocity_In_Model(const Eigen::Vector3d& vel) {
	this->velocity = vel;
}
void ParticleD::Set_DeltaxAgo_In_Model(const Eigen::Vector3d &pos) {
	this->deltaxago_in_model = pos;
}
void ParticleD::Set_IM_Grid(int a) {
	this->imgrid = this->grid;
}
void ParticleD::Set_Mass(double a) {
	this->p_mass = a;
}
//==========================================================================//
//==========================================================================//
void ParticleD::Add_Force(const Eigen::Vector3d &f) {
	this->force += f;
}

//重複してるかどうかを確認し、重複していたらEdgeに追加しない
void ParticleD::Add_Edge_Info(const int global_i, const int global_j, const int local_i, const int local_j, const Eigen::MatrixXd kij, const int tetra_group_id) {
	EdgeD temppEdge;
	temppEdge.edge_a_i_global = global_i;
	temppEdge.edge_b_i_global = global_j;
	temppEdge.edge_a_i_local = local_i;
	temppEdge.edge_b_i_local = local_j;
	temppEdge.K_belongs.push_back(tetra_group_id);
	temppEdge.stiff_K_ijs.push_back(new Eigen::MatrixXd(kij));

	bool findsame = false;
	for (unsigned int pei = 0; pei < p_Edges.size(); pei++) {
		if ((p_Edges[pei].edge_a_i_local == temppEdge.edge_a_i_local) && (p_Edges[pei].edge_b_i_local == temppEdge.edge_b_i_local)) {
			findsame = true;
			bool beloned = false;
			for (unsigned int bli = 0; bli < p_Edges[pei].K_belongs.size(); bli++) {
				if (tetra_group_id == p_Edges[pei].K_belongs[bli]) {
					(*p_Edges[pei].stiff_K_ijs[bli]) = (*p_Edges[pei].stiff_K_ijs[bli]) + kij;
					beloned = true;
					break;
				}
			}
			if (!beloned) {
				p_Edges[pei].stiff_K_ijs.push_back(new Eigen::MatrixXd(kij));
				p_Edges[pei].K_belongs.push_back(tetra_group_id);
			}
			break;
		}
	}
	//同じものがなかった辺情報を追加する
	if (!findsame) {
		temppEdge.edge_i = p_Edges.size();
		p_Edges.push_back(temppEdge);
	}
}
//一つ前(k)の予測位置と現在の予測位置(k+1)の差を計算する
double ParticleD::Add_convergence_iteration(double convite) {
	return convite + (this->exp_pos - this->expago_pos).squaredNorm();
}
//一つ前(k)の予測位置と現在の予測位置(k+1)の差を計算する(差を利用)
double ParticleD::Add_convergence_iteration2(double convite) {
	return convite + (this->deltax_in_model - this->deltaxago_in_model).squaredNorm();
}
//==========================================================================//
//==========================================================================//
const Eigen::Vector3d& ParticleD::Get_Displace()const {
	return this->displace;
}
const Eigen::Vector3d& ParticleD::Get_Force()const {
	return this->force;
}
//ParticleDのグローバル座標を取得
const Eigen::Vector3d& ParticleD::Get_Grid()const {
	return this->grid;
}
const Eigen::Vector3d& ParticleD::Get_Vel()const {
	return this->velocity;
}
const Eigen::Vector3d& ParticleD::Get_Acc()const {
	return this->acceleration;
}
const Eigen::Vector3d& ParticleD::Get_Exp_Pos()const {
	return this->exp_pos;
}
const Eigen::Vector3d& ParticleD::Get_Prime_Pos()const {
	return this->prime_pos;
}
const Eigen::Vector3d& ParticleD::Get_Initial_Pos()const {
	return this->initial_pos;
}
const Eigen::Vector3d& ParticleD::Get_ExpAgo_Pos()const {
	return this->expago_pos;
}
const Eigen::Vector3d& ParticleD::Get_Deltax_In_Model()const {
	return this->deltax_in_model;
}
const Eigen::Vector3d& ParticleD::Get_Mean_Vel()const {
	return this->mean_vel;
}
const Eigen::Vector3d& ParticleD::Get_DeltaxAgo_In_Model()const {
	return this->deltaxago_in_model;
}
const Eigen::Vector3d& ParticleD::Get_Draw_Grid()const {
	return this->drawgrid;
}
const Eigen::Vector3d& ParticleD::Get_IM_Grid()const {
	return this->imgrid;
}
const double& ParticleD::Get_Mass()const {
	return this->p_mass;
}
//==========================================================================//
//==========================================================================//
bool ParticleD::Is_Fixed()const {
	return this->fixed;
}
//==========================================================================//