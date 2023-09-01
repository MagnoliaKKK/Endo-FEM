//==========================================================================//
//@author KatsuyaKikuchi
//@brief 各パーティクル
//==========================================================================//

#ifndef _PARTICLED
#define _PARTICLED

#include "Debug.h"
#include "EdgeD.h"

class ParticleD {
	//各頂点の位置や速度,加速度,力のほか隣接点などの情報をもったParticleDクラス

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW	//Eigenをメンバ変数に使うためのマクロ定義

	ParticleD(const Eigen::Vector3d& grid);
	ParticleD(const Eigen::Vector3d& grid, int a);
	~ParticleD();

	void Draw()const;
	void Update(const Eigen::Vector3d &next_grid);
	void Update_Grid(const Eigen::Vector3d &next_grid);
	void Update_Draw_Grid(const Eigen::Vector3d &next_grid);
	void Update_Velocity(const Eigen::Vector3d &next_grid);

	void Set_Force(const Eigen::Vector3d &f);
	void Set_Fixed(const bool& fixed);
	void Set_Exp_Pos(const  Eigen::Vector3d &pos);
	void Set_Prime_Pos(const  Eigen::Vector3d &pos);
	void Set_Initial_Pos(const  Eigen::Vector3d &pos);
	void Set_ExpAgo_Pos(const  Eigen::Vector3d &pos);

	//この節点のFEMにおける解ベクトルをセット
	void Set_Deltax_In_Model(const  Eigen::Vector3d &pos);
	void Set_Velocity_In_Model(const  Eigen::Vector3d &vel);
	
	void Set_DeltaxAgo_In_Model(const  Eigen::Vector3d &pos);

	void Set_IM_Grid(int a);
	void Set_Mass(double a);

	void Add_Force(const Eigen::Vector3d &f);
	//重複してるかどうかを確認し、重複していたらEdgeに追加しない
	void Add_Edge_Info(const int global_i, const int global_j, const int local_i, const int local_j, const Eigen::MatrixXd kij, const int tetrea_group_id);

	const Eigen::Vector3d& Get_Displace()const;
	const Eigen::Vector3d& Get_Force()const;
	const Eigen::Vector3d& Get_Grid()const;
	const Eigen::Vector3d& Get_Vel()const;
	const Eigen::Vector3d& Get_Acc()const;
	const Eigen::Vector3d& Get_Exp_Pos()const;
	const Eigen::Vector3d& Get_Prime_Pos()const;
	const Eigen::Vector3d& Get_Initial_Pos()const;
	const Eigen::Vector3d& Get_ExpAgo_Pos()const;
	const Eigen::Vector3d& Get_Deltax_In_Model()const;
	const Eigen::Vector3d& Get_Mean_Vel()const;
	const Eigen::Vector3d& Get_DeltaxAgo_In_Model()const;
	const Eigen::Vector3d& ParticleD::Get_Draw_Grid()const;
	const Eigen::Vector3d& ParticleD::Get_IM_Grid()const;
	const double& ParticleD::Get_Mass()const;
	

	bool Is_Fixed()const;

	//addition by tei
	int p_id;
	std::vector<EdgeD> p_Edges;//辺情報
	double p_mass;	//ParticleDの質量

	std::vector<int> p_belong_TetraGroup_ids; //p_belong_TetraGroup_ids; 点pが所属するグループのidの配列

	double Add_convergence_iteration(double convite);//一つ前(k)の予測位置と現在の予測位置(k+1)の差を計算する
	double Add_convergence_iteration2(double convite);//一つ前(k)の予測位置と現在の予測位置(k+1)の差を計算する(差を利用)
private:
	Eigen::Vector3d grid;			//グローバル座標
	Eigen::Vector3d force;			//力
	Eigen::Vector3d displace;
	Eigen::Vector3d velocity;		//速度
	Eigen::Vector3d acceleration;	//加速度
	Eigen::Vector3d exp_pos;		//予測位置（都度変更される）
	Eigen::Vector3d expago_pos;		//予測位置（都度変更される,その一つ前を記録）
	Eigen::Vector3d prime_pos;		//弾性力以外による変位した位置（1ステップで固定）
	Eigen::Vector3d initial_pos;	//初期位置
	Eigen::Vector3d deltax_in_model;//線形方程式の解ベクトル(モデルで共有)
	Eigen::Vector3d mean_vel;
	Eigen::Vector3d deltaxago_in_model;//線形方程式の解ベクトル(モデルで共有)(1つ前を記録)
	bool fixed;						//固定されているかどうか
	Eigen::Vector3d drawgrid;	//描画用のグローバル座標
	Eigen::Vector3d imgrid;	//分割用のグローバル座標
};

#endif
#pragma once

