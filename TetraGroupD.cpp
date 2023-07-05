//------------------------------------------
//@author KatsuyaKikuchi
//@brief Implementation of the TetraGroupD.h
//------------------------------------------

#include "TetraGroupD.h"
#include <typeinfo>
#include "UTClapack.h"



//==========================================================================//
TetraGroupD::TetraGroupD(std::vector< TetraElementD* > elements, ObjectData data, double mass, int tetra_Group_id)
	:elements(elements),
	rotate_matrix(Eigen::Matrix3d::Identity(3, 3)),
	quaternion(Eigen::Quaterniond(Eigen::Matrix3d::Identity(3, 3))),
	center_grid(Eigen::Vector3d::Zero()),
	origin_center_grid(Eigen::Vector3d::Zero()),
	mass(mass), f_damping(0.35), v_damping(0.8), data(data),
	rotate_matrix_trans(Eigen::Matrix3d::Identity(3, 3)),
	tetra_group_id(tetra_Group_id),
	particles(Create_Particles(elements)),
	particle_num(particles.size())
{
	Init();
	BuildMatrix();
}
TetraGroupD::~TetraGroupD() {
	//particlesに格納されているポインタは外部で定義。
	//TetraGroupD内では解放しないように注意。
}
//==========================================================================//
//==========================================================================//
//	@start		   				初期設定									//
//==========================================================================//
void TetraGroupD::Init() {
	Set_Size_para(particle_num);				 // 各変数の要素数を決定する
	Set_Size_para2(particles);					 // 各変数の行列のサイズを決定する
}
void TetraGroupD::BuildMatrix() {
	//==========================================================================//
	//	@start		   			剛性行列の計算									//
	//==========================================================================//

	//各四面体要素の質量行列を作成する
	for (auto _e : elements) {
		_e->Create_M_Matrix(data.density);
	}

	//グループの質量行列を作成する
	Create_M_Matrix();
	std::cout << "Create Object M Matrix Of Group" << tetra_group_id << std::endl;

	//particleの質量を質量行列から計算する
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		m_In_Group[pi] = M_Matrix(3 * pi, 3 * pi);
	}

	//グループの対角化質量行列の作成
	//Create_Diag_M_Matrix();
	//グループの対角化質量SuM行列の作成
	//重心を求めるときに使う
	//Create_SUM_M_Matrix();
	//質量行列の対角成分の和が質量と同じであることの確かめ
	//double cluence;
	//for (int pi = 0; pi < particle_num; pi++) {
	//	cluence += M_Matrix(3 * pi, 3 * pi);
	//}
	//std::cout << "cluence is " << cluence << "true mass is " << mass << std::endl;

	//Create_Center_Grid();

	//各四面体要素の剛性行列の計算
	//1と2がある
	/*
	for (auto _e : elements) {
		_e->Create_Stiffness_Matrix2(origin_center_grid, data.young, data.poisson);
	}
	*/
	//グループの剛性行列を作成する
	/*Create_Stiffness_Matrix();
	std::cout << "Create Object K Matrix Of Group" << tetra_group_id << std::endl;
	*/

	//重さベクトルの生成(単位はNニュートン)
	Eigen::VectorXd g(3 * particles.size());
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		size_t index = std::distance(particles.begin(), it);
		g.block(3 * index, 0, 3, 1) = Eigen::Vector3d(0.0, Gravity, 0.0);//重力加速度は9.8m/s2としている
	}
	this->gravity = g;

	//各particleに隣接するparticleの情報を加える
	//Find_Edges();
	//std::cout << "Find Edge Of Group" << tetra_group_id << std::endl;

	//CRS形式で剛性行列を格納する
	//Create_Local_Stiffmatrix_Array();
	//std::cout << "Create_Local_Stiffmatrix_Array" << tetra_group_id << std::endl;
	//std::cout << "Particles of "<<tetra_group_id << " is " << std::endl;

	//グループに入っている節点の順番を確認
	for (auto _p:particles) {
		std::cout << _p->p_id << std::endl;
	}
	std::cout << std::endl;
}
//グループの剛性行列を作成する
void TetraGroupD::Create_Stiffness_Matrix() {
	stiffness_matrix = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	for (auto _e : elements) {
		for (auto p1_it = particles.begin(); p1_it != particles.end(); ++p1_it) {
			size_t p1_index = std::distance(particles.begin(), p1_it);

			for (auto p2_it = particles.begin(); p2_it != particles.end(); ++p2_it) {
				size_t p2_index = std::distance(particles.begin(), p2_it);

				stiffness_matrix.block(3 * p1_index, 3 * p2_index, 3, 3) =
					stiffness_matrix.block(3 * p1_index, 3 * p2_index, 3, 3) + _e->Get_K_Submatrix(*p1_it, *p2_it);
			}
		}
	}
}

//グループの質量行列を作成する
void TetraGroupD::Create_M_Matrix() {
	M_Matrix = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	for (auto _e : elements) {
		for (auto p1_it = particles.begin(); p1_it != particles.end(); ++p1_it) {
			size_t p1_index = std::distance(particles.begin(), p1_it);
			for (auto p2_it = particles.begin(); p2_it != particles.end(); ++p2_it) {
				size_t p2_index = std::distance(particles.begin(), p2_it);

				M_Matrix.block(3 * p1_index, 3 * p2_index, 3, 3) =
					M_Matrix.block(3 * p1_index, 3 * p2_index, 3, 3) + _e->Get_M_Submatrix(*p1_it, *p2_it);
			}
		}
	}
	//質量行列の逆行列を作成
	//std::cout << M_Matrix << std::endl;
	inv_M_Matrix = M_Matrix.inverse();
}
//グループの対角化質量行列を作成する
void TetraGroupD::Create_Diag_M_Matrix() {
	ParticleD* pit;
	Diag_M_Matrix = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	for (unsigned int pi = 0; pi < particle_num;pi++) {
		pit = particles[pi];
		Diag_M_Matrix.block(3 * pi, 3 * pi, 3, 3) = pit->p_mass * Eigen::Matrix3d::Identity();
	}
	//std::cout << "Diag_M_Matrix" <<std::endl;
	//std::cout << Diag_M_Matrix << std::endl;
}

//グループの対角化SUM質量行列を作成する
void TetraGroupD::Create_SUM_M_Matrix() {
	//ParticleD* pit;
	SUM_M_Matrix = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	Eigen::MatrixXd SUMsub_M_Matrix = Eigen::MatrixXd::Zero(3 , 3 * particle_num);
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		//SUMsub_M_Matrix.block(0 , 3 * pi, 3, 3) = (M_Matrix_C(3 * pi,3 * pi)/Group_Mass) * Eigen::Matrix3d::Identity();
		SUMsub_M_Matrix.block(0, 3 * pi, 3, 3) = (m_In_Group[pi] / mass) * Eigen::Matrix3d::Identity();
	}
	//std::cout << "SUMsub_M_Matrix" << std::endl;
	//std::cout << SUMsub_M_Matrix << std::endl;
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		SUM_M_Matrix.block(3 * pi, 0, 3, 3 * particle_num) = SUMsub_M_Matrix;
	}
	/*std::cout << "SUM_M_Matrix" << std::endl;
	std::cout << SUM_M_Matrix << std::endl;*/
	std::cout << "Create SUM_M_Matrix Of Group " << tetra_group_id<< std::endl;
}
//グループの減衰行列を作成する
void TetraGroupD::Create_Damping_Matrix(){
	//Damping_Matrix = B_damping * Eigen::MatrixXd::Identity(3 * particle_num, 3 * particle_num);
	//質量減衰のみ
	Damping_Matrix = M_damping * TIME_STEP * M_Matrix_C;

	//剛性減衰のみ
	//Damping_Matrix = M_damping * TIME_STEP * stiffness_matrix;

	//レイリー近似
	//Damping_Matrix = TIME_STEP * (M_damping * M_Matrix_C + K_damping * stiffness_matrix);

	MassDamInv_Matrix = (M_Matrix_C + Damping_Matrix).inverse();

	Eigen::MatrixXd MassDamInvMassT = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	MassDamInvMassT = MassDamInv_Matrix * M_Matrix_C * TIME_STEP;
	DammingT_Matrix_Sparse.setZero();
	DammingT_Matrix_Sparse = MassDamInvMassT.sparseView();

	//質量中心行列のSparse作成(JacobiやConstの計算で使用)
	Eigen::MatrixXd MassCondi = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	MassCondi = Eigen::MatrixXd::Identity(3 * particles.size(), 3 * particles.size()) - SUM_M_Matrix;//(I-Mj,cm)
	Eigen::MatrixXd Damm_Matrix = M_Matrix_C + Damping_Matrix;//(Mj+Cj')
	MassCondi_Sparse.setZero();
	Damm_Matrix_Sparse.setZero();
	MassCondi_Sparse = MassCondi.sparseView();
	Damm_Matrix_Sparse = Damm_Matrix.sparseView();
	std::cout << "Create Damping Matrix of Group " << tetra_group_id << std::endl;
}

//TetraElementに使用されているパーティクルを重複なくparticlesに格納する
std::vector< ParticleD* > TetraGroupD::Create_Particles(std::vector< TetraElementD* > elements) {
	std::vector< ParticleD* > p;
	for (auto _e : elements) {
		std::vector< ParticleD* > temp_p = _e->Get_Particle();
		for (auto _p : temp_p) {
			p.push_back(_p);
		}
	}
	//重複要素を削除
	std::sort(p.begin(), p.end());
	p.erase(std::unique(p.begin(), p.end()), p.end());
	return p;
}
//静止時の重心やローカルベクトルの生成
//節点質量
void TetraGroupD::Create_Center_Grid() {
	origin_center_grid = Eigen::Vector3d::Zero();
	center_grid = Eigen::Vector3d::Zero();
	for (unsigned int oi = 0; oi < particle_num; oi++) {
		origin_center_distance[oi] = Eigen::Vector3d::Zero();
		center_distance[oi] = Eigen::Vector3d::Zero();
		origin_local_grid[oi] = Eigen::Vector3d::Zero();
	}
	//静止状態のグループの位置ベクトル生成
	Eigen::VectorXd PosVector = Eigen::VectorXd::Zero(3 * particle_num);
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		PosVector.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Grid();
	}

	//力がかかってないとき(time=0)における重心ベクトルの作成
	Eigen::VectorXd x(3 * particles.size());
	Eigen::VectorXd centerog(3 * particles.size());
	//質量行列が対角成分のみでないのとき
	if (mdiag == FALSE) {
		for (auto it = particles.begin(), end = particles.end(); it != end; it++) {
			size_t index = std::distance(particles.begin(), it);
			x.block(3 * index, 0, 3, 1) = (*it)->Get_Grid();
		}
		centerog = (M_Matrix / mass)*x;
		for (auto it = particles.begin(), end = particles.end(); it != end; it++) {
			size_t index = std::distance(particles.begin(), it);
			origin_center_grid += centerog.block(3 * index, 0, 3, 1);
		}
	}
	//質量行列が対角成分のみのとき
	else {
		/*
		for (unsigned int oi = 0; oi < particle_num; oi++) {
			origin_center_grid += M_Matrix_C(3 * oi, 3 * oi) * particles[oi]->Get_Grid();
		}
		origin_center_grid = origin_center_grid / Group_Mass;
		*/
		origin_center_grid = SUM_M_Matrix.block(0, 0, 3, 3 * particle_num) * PosVector;
	}
	//シミュレーションによって重心の位置は変化するので、変化する重心ベクトルを設定
	center_grid = origin_center_grid;
	std::cout << "center_grid is " << std::endl;
	std::cout << origin_center_grid << std::endl;

	//各particleにおける初めの重心からの距離(変化しない),
	//現在の重心からの距離(変化する),各particleのローカル座標初期位置(変化しない)の設定
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		//ParticleD* pit = particles[pi];
		origin_center_distance[pi] = particles[pi]->Get_Grid() - origin_center_grid;
		center_distance[pi] = particles[pi]->Get_Grid() - origin_center_grid;
		//各particleのローカル座標初期位置は最初は回転していない(rotate_matrix=単位行列)
		//origin_local_grid[pi] = rotate_matrix *(pit->Get_Grid() - origin_center_grid);
		origin_local_grid[pi] = particles[pi]->Get_Grid() - origin_center_grid;
		//std::cout << "origin_local_grid "<< particles[pi]->p_id <<" is"<< std::endl;
		//std::cout << origin_local_grid[pi] << std::endl;
	}

	//Constの計算で使うベクトルの生成
	OrigineVector = Eigen::VectorXd::Zero(3 * particles.size());
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		OrigineVector.block(3 * pi, 0, 3, 1) = origin_local_grid[pi];
	}
	//std::cout << "OrigineVector" << OrigineVector << std::endl;
	/*
	Eigen::VectorXd OrigineVector2 = Eigen::VectorXd::Zero(3 * particles.size());
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		OrigineVector2.block(3 * pi, 0, 3, 1) = PosVector.block(3 * pi, 0, 3, 1) - origin_center_grid;
	}
	*/
	//std::cout << "OrigineVector" << OrigineVector-OrigineVector2 << std::endl;
}
//局所剛性行列の作成
void TetraGroupD::Create_Local_Stiffness_Matrix() {
	//各四面体要素の剛性行列の計算
	//1と2がある
	
	for (auto _e : elements) {
		/*形状関数はdetだからたぶんどこ原点でも大丈夫
		Eigen::Vector3d inner_ele_p = Eigen::Vector3d::Zero();
		for (auto _p:_e->Get_Particle()) {
			inner_ele_p += _p->Get_Grid();
		}
		inner_ele_p = inner_ele_p / 4.0;
	//要素の中心(重心ではない)を基準にするならこっち
	_e->Create_Stiffness_Matrix2(inner_ele_p, data.young, data.poisson);
	*/
	//グループの重心を基準にするならこっち
	_e->Create_Stiffness_Matrix2(origin_center_grid, data.young, data.poisson);
	}

	//1のほう(なにかしらの計算ミスがあるかも？)
	/*
	for (auto _e : elements) {
		_e->Create_Stiffness_Matrix(origin_center_grid, data.young, data.poisson);
	}
	*/

	//グループの剛性行列を作成する
	Create_Stiffness_Matrix();

	StiffnessTT_Matrix_Sparse.setZero();
	StiffnessTT_Matrix_Sparse = (stiffness_matrix * TIME_STEP * TIME_STEP).sparseView();

	std::cout << "Create Object K Matrix Of Group" << tetra_group_id << std::endl;
	std::cout << StiffnessTT_Matrix_Sparse << std::endl;
}
//節点情報などを付加する
void TetraGroupD::Create_Information() {
	//各particleに隣接するparticleの情報を加える(菊池さんか丁さんが書いたからようわからん)
	Find_Edges();
	std::cout << "Find Edge Of Group" << tetra_group_id << std::endl;
}

//==========================================================================//
//	@end		   				初期設定									//
//==========================================================================//

//==========================================================================//
//	@start		   				ループ設定									//
//==========================================================================//
void TetraGroupD::Draw()const {
	for (auto _e : elements) {
		_e->Draw();
	}
}
void TetraGroupD::Set_Gravity() {
}

//各particleに隣接する点を記録する
void TetraGroupD::Find_Edges() {
	for (auto _e : elements) {
		std::vector<ParticleD*> temp_ps = _e->Get_Particle();
		for (auto it = particles.begin(); it != particles.end(); ++it) {
			size_t index1 = std::distance(particles.begin(), it);

			for (auto it2 = particles.begin(); it2 != particles.end(); ++it2) {
				size_t index2 = std::distance(particles.begin(), it2);
				size_t p1_index = -1, p2_index = -1;

				for (auto tempit = temp_ps.begin(); tempit != temp_ps.end(); ++tempit) {
					size_t index = std::distance(temp_ps.begin(), tempit);
					if (*it == *tempit) { p1_index = index; }
					if (*it2 == *tempit) { p2_index = index; }
				}
				if (-1 == p1_index || -1 == p2_index) {
					continue;
				}
				//重複してるかどうかを確認し、重複していたらEdgeに追加しない
				//(モデルにおける1のid,モデルにおける1のid,グループにおける1のid,グループにおける2のid,剛性行列の値,グループのid)
				(*it)->Add_Edge_Info((*it)->p_id, (*it2)->p_id, index1, index2, _e->Get_K_Submatrix(*it, *it2), this->tetra_group_id);
			}
		}
	}
}

//行列を値が存在するところのみ格納する
void TetraGroupD::Create_Local_Stiffmatrix_Array() {

	//隣接点を記録する(隣接する点とのみ剛性行列は値をもつ)
	for (auto _e : elements) {
		std::vector<ParticleD*> temp_ps = _e->Get_Particle();
		for (auto it = particles.begin(); it != particles.end(); ++it) {
			size_t index1 = std::distance(particles.begin(), it);
			for (auto it2 = it; it2 != particles.end(); ++it2) {
				size_t index2 = std::distance(particles.begin(), it2);
				size_t p1_index = -1, p2_index = -1;

				for (auto tempit = temp_ps.begin(); tempit != temp_ps.end(); ++tempit) {
					size_t index = std::distance(temp_ps.begin(), tempit);
					if (*it == *tempit) { p1_index = index; }
					if (*it2 == *tempit) { p2_index = index; }
				}
				if (-1 == p1_index || -1 == p2_index) {
					continue;
				}
				stiffmatrix_valued_list[index1].push_back(index2);
				if (index1 != index2) {
					stiffmatrix_valued_list_sym[index2].push_back(index1);
				}
			}
		}
	}

	//重複して記録してある点があるので消去する
	for (unsigned int vi = 0; vi < stiffmatrix_valued_list.size(); vi++) {
		sort(stiffmatrix_valued_list[vi].begin(), stiffmatrix_valued_list[vi].end());
		stiffmatrix_valued_list[vi].erase(unique(stiffmatrix_valued_list[vi].begin(), stiffmatrix_valued_list[vi].end()), stiffmatrix_valued_list[vi].end());
	}
	//重複して記録してある点があるので消去する
	for (unsigned int vi = 0; vi < stiffmatrix_valued_list_sym.size(); vi++) {
		sort(stiffmatrix_valued_list_sym[vi].begin(), stiffmatrix_valued_list_sym[vi].end());
		stiffmatrix_valued_list_sym[vi].erase(unique(stiffmatrix_valued_list_sym[vi].begin(), stiffmatrix_valued_list_sym[vi].end()), stiffmatrix_valued_list_sym[vi].end());
	}

	//隣接リストから剛性行列を抜き出しリストにする
	Eigen::MatrixXi rotated_indexesmatrix = Eigen::MatrixXi::Zero(3 * particles.size(), 3 * particles.size());
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		size_t index1 = std::distance(particles.begin(), it);
		ParticleD* pit = particles[index1];
		for (unsigned int vi = 0; vi < stiffmatrix_valued_list[index1].size(); vi++) {
			//隣接リストから剛性行列を抜き出す
			Eigen::Matrix3d* local = new Eigen::Matrix3d((stiffness_matrix.block(3 * index1, 3 * stiffmatrix_valued_list[index1][vi], 3, 3)));
			rotated_stiffmatrixs.push_back(local);

			//同じ頂点ではないとき
			if (index1 != stiffmatrix_valued_list[index1][vi]) {
				rotated_indexesmatrix(index1, stiffmatrix_valued_list[index1][vi]) = rotated_stiffmatrixs.size() - 1;
			}
		}
		//グループのparticleの数の分、配列をとる
		local_xs.push_back(new Eigen::Vector3d(Eigen::Vector3d::Zero()));
	}
	//値がある場所を記録する
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		size_t index1 = std::distance(particles.begin(), it);
		for (unsigned int vi = 0; vi < stiffmatrix_valued_list_sym[index1].size(); vi++) {
			valued_sym_rotated_indexes.push_back(rotated_indexesmatrix(stiffmatrix_valued_list_sym[index1][vi], index1));
		}
	}
}

//予測位置(弾性力を考慮しない)を事前計算する
void TetraGroupD::Calc_Exp_Pos() {
	f_Local = Eigen::VectorXd::Zero(3 * particles.size());
	x_Local = Eigen::VectorXd::Zero(3 * particles.size());
	v_Local = Eigen::VectorXd::Zero(3 * particles.size());
	for (auto it = particles.begin(), end = particles.end(); it != end; it++) {
		size_t index = std::distance(particles.begin(), it);
		f_Local.block(3 * index, 0, 3, 1) = (*it)->Get_Force();
		x_Local.block(3 * index, 0, 3, 1) = (*it)->Get_Grid();
		v_Local.block(3 * index, 0, 3, 1) = (*it)->Get_Vel();
	}

	//重力を加え速度を更新する
	for (auto it = particles.begin(), end = particles.end(); it != end; it++) {
		size_t index = std::distance(particles.begin(), it);
		v_Local.block(3 * index, 0, 3, 1) = v_Local.block(3 * index, 0, 3, 1) + (f_Local.block(3 * index, 0, 3, 1) / (*it)->Get_Mass()) * TIME_STEP + gravity.block(3 * index, 0, 3, 1) * TIME_STEP;
	}

	//弾性力を除いた時の位置の更新
	Eigen::VectorXd p = Eigen::VectorXd::Zero(3 * particles.size());
	p = x_Local + v_Local * TIME_STEP;
	for (auto it = particles.begin(), end = particles.end(); it != end; it++) {
		size_t index = std::distance(particles.begin(), it);
		if (!((*it)->Is_Fixed())) {
			(*it)->Set_Exp_Pos(p.block(3 * index, 0, 3, 1));
			(*it)->Set_Prime_Pos(p.block(3 * index, 0, 3, 1));
			//(*it)->Set_Initial_Pos(p.block(3 * index, 0, 3, 1));
		}
		//固定点の場合はずっと同じ値が入っている気がする
		//(だから下みたいなことをしなくてもいいはず)
		else {
			(*it)->Set_Exp_Pos((*it)->Get_Grid());
			(*it)->Set_Prime_Pos((*it)->Get_Grid());
			//(*it)->Set_Initial_Pos((*it)->Get_Grid());
		}
		//debug
		/*
		std::cout << (*it)->p_id << "particle" << std::endl;
		std::cout<<(*it)->Get_Exp_Pos()<< std::endl;
		std::cout << (*it)->Get_Grid() << std::endl;
		std::cout << "Finished EXP " << std::endl;
		*/
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Calc_Exp" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);
}
// 予測位置(弾性力を考慮しない)を事前計算する(damperをいれる)
void TetraGroupD::Calc_Exp_Pos2() {
	f_Local = Eigen::VectorXd::Zero(3 * particles.size());
	x_Local = Eigen::VectorXd::Zero(3 * particles.size());
	v_Local = Eigen::VectorXd::Zero(3 * particles.size());
	
	/*
	for (auto it = particles.begin(), end = particles.end(); it != end; it++) {
		size_t index = std::distance(particles.begin(), it);
		f_Local.block(3 * index, 0, 3, 1) = (*it)->Get_Force();
		x_Local.block(3 * index, 0, 3, 1) = (*it)->Get_Grid();
		v_Local.block(3 * index, 0, 3, 1) = (*it)->Get_Vel();
	}
	*/
	for (unsigned int pi = 0; pi < particle_num;pi++) {
		f_Local.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Force();
		x_Local.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Grid();
		v_Local.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Vel();
	}

	//重力を加え速度を更新する
	/*
	for (auto it = particles.begin(), end = particles.end(); it != end; it++) {
		size_t index = std::distance(particles.begin(), it);
		v_Local.block(3 * index, 0, 3, 1) = v_Local.block(3 * index, 0, 3, 1) + (f_Local.block(3 * index, 0, 3, 1) / (*it)->Get_Mass()) * TIME_STEP + gravity.block(3 * index, 0, 3, 1) * TIME_STEP;
	}
	*/
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		v_Local.block(3 * pi, 0, 3, 1) = v_Local.block(3 * pi, 0, 3, 1) + (f_Local.block(3 * pi, 0, 3, 1) / M_Matrix_C(3*pi,3*pi)) * TIME_STEP + gravity.block(3 * pi, 0, 3, 1) * TIME_STEP;
	}

	//弾性力を除いた時の位置の更新(予測位置)
	Eigen::VectorXd mp = Eigen::VectorXd::Zero(3 * particles.size());
	if (useSparse) {
		mp = x_Local + DammingT_Matrix_Sparse * v_Local;
	}
	else {
		mp = x_Local + MassDamInv_Matrix * M_Matrix_C * v_Local * TIME_STEP;
		//std::cout << "MassDamInv_Matrix * M_Matrix_C" << std::endl;
		//std::cout << MassDamInv_Matrix * M_Matrix_C << std::endl;
	}
	/*
	for (auto it = particles.begin(), end = particles.end(); it != end; it++) {
		size_t index = std::distance(particles.begin(), it);
		if (!((*it)->Is_Fixed())) {
			(*it)->Set_Exp_Pos(mp.block(3 * index, 0, 3, 1));
			(*it)->Set_Prime_Pos(mp.block(3 * index, 0, 3, 1));
		}
		else {
			(*it)->Set_Exp_Pos((*it)->Get_Grid());
			(*it)->Set_Prime_Pos((*it)->Get_Grid());
		}
	}*/
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		//固定されていなければ予測位置更新
		if (!(particles[pi]->Is_Fixed())) {
			particles[pi]->Set_Exp_Pos(mp.block(3 * pi, 0, 3, 1));
			particles[pi]->Set_Prime_Pos(mp.block(3 * pi, 0, 3, 1));
		}
		//固定されていたら一応GRIDをいれておく
		else {
			particles[pi]->Set_Exp_Pos(particles[pi]->Get_Grid());
			particles[pi]->Set_Prime_Pos(particles[pi]->Get_Grid());
		}
		//std::cout<<"particle"<< particles[pi]->p_id<< std::endl;
		//std::cout << particles[pi]->Get_Exp_Pos() << std::endl;
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Calc_Exp" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);

	//Constの計算で使うベクトルの生成
	PrimeVector = Eigen::VectorXd::Zero(3 * particles.size());
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		PrimeVector.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Prime_Pos();
	}
}
// 予測位置(弾性力を考慮しない)を事前計算する(damperをいれる+グループごとの予測位置をもつ)
void TetraGroupD::Calc_Exp_Pos3() {
	//初期化
	f_Local = Eigen::VectorXd::Zero(3 * particles.size());
	x_Local = Eigen::VectorXd::Zero(3 * particles.size());
	v_Local = Eigen::VectorXd::Zero(3 * particles.size());

	//現時点での位置や力、速度のベクトル化
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		f_Local.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Force();
		x_Local.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Grid();
		v_Local.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Vel();
	}
	//重力を加え速度を更新する
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		v_Local.block(3 * pi, 0, 3, 1) = v_Local.block(3 * pi, 0, 3, 1) + (f_Local.block(3 * pi, 0, 3, 1) / M_Matrix_C(3 * pi, 3 * pi)) * TIME_STEP + gravity.block(3 * pi, 0, 3, 1) * TIME_STEP;
		//v_Local.block(3 * pi, 0, 3, 1) = v_Local.block(3 * pi, 0, 3, 1) + (f_Local.block(3 * pi, 0, 3, 1) / m_In_Group[pi]) * TIME_STEP + gravity.block(3 * pi, 0, 3, 1) * TIME_STEP;
	}

	//弾性力を除いた時の位置の更新(予測位置)
	//Constの計算や回転行列の推定で使うベクトルの生成
	PrimeVector = Eigen::VectorXd::Zero(3 * particles.size());
	if (useSparse) {
		PrimeVector = x_Local + DammingT_Matrix_Sparse * v_Local;
	}
	else {
		PrimeVector = x_Local + MassDamInv_Matrix * M_Matrix_C * v_Local * TIME_STEP;
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Calc_Exp" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);
	//固定点は静止時の値に修正しておく(3x3単位なら上でやったほうがいい)
	/*
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->Is_Fixed()) {
			//PrimeVector.block(3 * pi, 0, 3, 1) = InitialVector.block(3 * pi, 0, 3, 1);
			PrimeVector.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Grid();
		}
	}
	*/
	//std::cout << "calc EXP" << std::endl;
}
// グループごとの予測位置の計算
// Calculate the predicted position for each group
void TetraGroupD::Calc_Exp_Pos_Group() {
	//初期化
	f_Local = Eigen::VectorXd::Zero(3 * particles.size());
	x_Local = Eigen::VectorXd::Zero(3 * particles.size());
	v_Local = Eigen::VectorXd::Zero(3 * particles.size());
	
	feclearexcept(FE_ALL_EXCEPT);

	//現時刻での位置や力、速度のベクトル化
	//vectorization of position, force, and velocity at current time
	for (unsigned int pi = 0; pi < particle_num; pi++) {

		f_Local.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Force();
		x_Local.block(3 * pi, 0, 3, 1) = GroupGridVector.block(3 * pi, 0, 3, 1);
		v_Local.block(3 * pi, 0, 3, 1) = GroupVelVector.block(3 * pi, 0, 3, 1);
		//初期に力をかける場合は下のようにする（ただし力をかけてどうなるのかは知らない）
		//If you want to apply a force early on, do the following (but don't know what happens when you apply the force)
		if (UseIniForce) {
			//particles[pi]->Set_Force(Eigen::Vector3d::Zero());
			if (tetra_group_id==1 && particles[pi]->p_id == 61&& countup<1000) {
				f_Local.block(3 * pi, 0, 3, 1) = Eigen::Vector3d(0, 1.0e+4, 0);
			}
		}

		/*feclearexcept(FE_ALL_EXCEPT);
		if (fetestexcept(FE_ALL_EXCEPT)) {
			std::cout << "P0, exception" << std::endl;
		}*/
	}
	//重力を加え速度を更新する
	//Apply gravity and update velocity
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		v_Local.block(3 * pi, 0, 3, 1) = v_Local.block(3 * pi, 0, 3, 1) + (f_Local.block(3 * pi, 0, 3, 1) / M_Matrix_C(3 * pi, 3 * pi)) * TIME_STEP + gravity.block(3 * pi, 0, 3, 1) * TIME_STEP;
		Eigen::VectorXd tempA = Eigen::VectorXd::Zero(3 * particles.size());
		Eigen::VectorXd tempB = Eigen::VectorXd::Zero(3 * particles.size());
		Eigen::VectorXd tempC = Eigen::VectorXd::Zero(3 * particles.size());
		

		//tempA.block(3 * pi, 0, 3, 1) = v_Local.block(3 * pi, 0, 3, 1);
		//if (fetestexcept(FE_INVALID)) {
		//	std::cout << "FE_INVALID Calc_Exp2" << std::endl;
		//}
		//feclearexcept(FE_ALL_EXCEPT);
		//tempB.block(3 * pi, 0, 3, 1) = (f_Local.block(3 * pi, 0, 3, 1) / M_Matrix_C(3 * pi, 3 * pi)) * TIME_STEP;
		//if (fetestexcept(FE_INVALID)) {
		//	std::cout << "FE_INVALID Calc_Exp3" << std::endl;
		//}
		//feclearexcept(FE_ALL_EXCEPT);
		//tempC.block(3 * pi, 0, 3, 1) = gravity.block(3 * pi, 0, 3, 1) * TIME_STEP;
		//if (fetestexcept(FE_INVALID)) {
		//	std::cout << "FE_INVALID Calc_Exp5" << std::endl;
		//}
		//feclearexcept(FE_ALL_EXCEPT);

		//v_Local.block(3 * pi, 0, 3, 1) = tempA + tempB + tempC;
		//if (fetestexcept(FE_INVALID)) {
		//	std::cout << "FE_INVALID Calc_Exp4" << std::endl;
		//}
		//feclearexcept(FE_ALL_EXCEPT);

		/*feclearexcept(FE_ALL_EXCEPT);
		if (fetestexcept(FE_ALL_EXCEPT)) {
			std::cout << "P1, exception" << std::endl;
		}*/

	}
	/*if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Calc_Exp000" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);*/
	//弾性力を除いた時の位置の更新(予測位置)
	//Update position when elastic forces are excluded (predicted position)
	//Constの計算や回転行列の推定で使うベクトルの生成
	//Generate vectors for use in calculating Vector b(: Ax=b ) and estimating rotation matrices
	PrimeVector = Eigen::VectorXd::Zero(3 * particles.size());
	if (useSparse) {
		PrimeVector = x_Local + DammingT_Matrix_Sparse * v_Local;
	}
	else {
		PrimeVector = x_Local + MassDamInv_Matrix * M_Matrix_C * v_Local * TIME_STEP;
		if (fetestexcept(FE_INVALID)) {
			std::cout << "FE_INVALID Calc_Exp111" << std::endl;
		}
		feclearexcept(FE_ALL_EXCEPT);
	}
	//std::cout << "Prime Vector" << std::endl << PrimeVector;
	/*if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Calc_Exp" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);*/

	//固定点の扱い方はいまだによくわかっていない
	//I still don't really understand how to handle fixed points.

	//固定点は静止時の値に修正しておく(3x3単位なら上でやったほうがいい)
	/*
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->Is_Fixed()) {
			PrimeVector.block(3 * pi, 0, 3, 1) = InitialVector.block(3 * pi, 0, 3, 1);
			//PrimeVector.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Grid();
		}
	}
	*/
	//std::cout << "calc EXP" << std::endl;
}

//CRSによって行列計算をうまく行いながらFEM計算をする
void TetraGroupD::Calc_CRSFEM() {
	ParticleD* pit;
	Eigen::Vector3d xcm = Eigen::Vector3d::Zero();

	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		xcm += pit->p_mass * particles[pi]->Get_Exp_Pos();
	}
	xcm = xcm / mass;

	int _rsmi = 0;
	int _rsm_symi = 0;

	for (unsigned int pi = 0; pi < particle_num; pi++) {
		(*local_xs[pi]) = rotate_matrix.transpose()*(particles[pi]->Get_Exp_Pos() - xcm) - origin_local_grid[pi];
	}

	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];

		Eigen::Vector3d u = Eigen::Vector3d::Zero();

		//上三角行列の部分の計算
		unsigned int svl_size = stiffmatrix_valued_list[pi].size();
		for (unsigned int smvi = 0; smvi <svl_size; smvi++) {
			unsigned int p_local_i = stiffmatrix_valued_list[pi][smvi];

			u += (*rotated_stiffmatrixs[_rsmi]) * (*local_xs[p_local_i]);
			_rsmi++;
		}
		//下三角行列の部分の計算
		unsigned int svlsy_size = stiffmatrix_valued_list_sym[pi].size();

		for (unsigned int smyi = 0; smyi < svlsy_size; smyi++) {
			unsigned int p_local_i = stiffmatrix_valued_list_sym[pi][smyi];

			u += (*rotated_stiffmatrixs[valued_sym_rotated_indexes[_rsm_symi]]).transpose() * (*local_xs[p_local_i]);
			_rsm_symi++;
		}

		//外力のみで計算した位置に弾性力を考慮した項を加える
		if (!((pit)->Is_Fixed())) {
			u = rotate_matrix * u / pit->p_mass * TIME_STEP * TIME_STEP * f_damping;
			//バネダンパ下記途中
			pit->Set_Exp_Pos(pit->Get_Prime_Pos() - u);
			if(M_damping!=0){
				pit->Set_Exp_Pos(pit->Get_Exp_Pos() - (pit->Get_Vel() / pit->p_mass) * TIME_STEP * TIME_STEP * M_damping);
			}
			pit->Set_Initial_Pos(pit->Get_Prime_Pos() - u);
		}
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID FEM" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);
}
//CRS形式を使わずに疎行列FEM計算をする
void TetraGroupD::Calc_FEM() {
	ParticleD* pit;

	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		pLocal.block(3 * pi, 0, 3, 1) = pit->Get_Exp_Pos();
	}

	Eigen::Vector3d xcm = Eigen::Vector3d::Zero();

	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		xcm += pit->p_mass * particles[pi]->Get_Exp_Pos();
	}
	xcm = xcm / mass;

	for (unsigned int pi = 0; pi < particle_num; pi++) {
		xgLocal.block(3 * pi, 0, 3, 1) = rotate_matrix.transpose()*(particles[pi]->Get_Exp_Pos() - xcm) - origin_local_grid[pi];
	}
	//上と下とで計算は同じだが、下のほうが速い
	//f_inLocal = R_Matrix * stiffness_matrix * trans_R_Matrix * (pLocal - xgLocal);
	f_inLocal = stiffness_matrix * xgLocal;

	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		if (!(pit)->Is_Fixed()) {
			Eigen::Vector3d u = rotate_matrix*f_inLocal.block(3 * pi, 0, 3, 1) * TIME_STEP * TIME_STEP*f_damping / pit->p_mass;
			pit->Set_Exp_Pos(pit->Get_Prime_Pos() - u);
			pit->Set_Initial_Pos(pit->Get_Prime_Pos() - u);
		}
	}
}

//回転行列の更新
void TetraGroupD::Update_Rotate() {
	//Eigen::Matrix3d Agorotation = Eigen::Matrix3d::Zero();
	//Agorotation = rotate_matrix;
	if (useAPD) {
		Create_Rotate_Matrix_APD();
	}
	else {
		Create_Rotate_Matrix();
	}
	//APDをはかった
	//Create_Rotate_Matrix_APD_Debug(Agorotation);
}

//svdによる極分解によって回転行列を計算する
void TetraGroupD::Create_Rotate_Matrix() {
	MicroSecondTimer mtRotationApq;
	MicroSecondTimer mtRotationSVD;
	MicroSecondTimer mtRotationSparse;
	mtRotationApq.setid(20);
	mtRotationSVD.setid(22);
	mtRotationSparse.setid(23);

	mtRotationApq.startMyTimer();
	//線形変換AはApqとAqqの積で表現できる
	Eigen::Matrix3d Apq = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d tempA = Eigen::Matrix3d::Zero();

	//Apqの計算
	if (useSparse) {
		//現在の重心の計算(x_cm^j)
		center_grid = Eigen::Vector3d::Zero();
		for (unsigned int pi = 0; pi < particle_num; pi++) {
			center_grid[0] += SUM_M_Matrix(0, 3 * pi) * PrimeVector[3 * pi];
			center_grid[1] += SUM_M_Matrix(1, 3 * pi) * PrimeVector[3 * pi + 1];
			center_grid[2] += SUM_M_Matrix(2, 3 * pi) * PrimeVector[3 * pi + 2];
		}

		//Apqの計算
		for (unsigned int pi = 0; pi < particle_num; pi++) {
			tempA = (PrimeVector.block(3 * pi, 0, 3, 1) - center_grid) * origin_center_distance[pi].transpose();
			//Apq += tempA;
			Apq += M_Matrix_C(3 * pi, 3 * pi) * tempA;
			//std::cout << "origin_center_distance" << std::endl << origin_center_distance[pi] << std::endl;

			//Apq += m_In_Group[pi] * tempA;
		}
		//std::cout << "Apq" << std::endl << Apq << std::endl;
		//std::cout << "PrimeVector" << std::endl << PrimeVector << std::endl;
		
	}
	else {
		//現在の重心の計算(x_cm^j)
		Eigen::VectorXd center_grid_Vector = Eigen::VectorXd::Zero(3 * particle_num, 3 * particle_num);
		center_grid_Vector = SUM_M_Matrix * PrimeVector;

		//Apqの計算
		for (unsigned int pi = 0; pi < particle_num; pi++) {
			tempA = (PrimeVector.block(3 * pi, 0, 3, 1) - center_grid_Vector.block(3 * pi, 0, 3, 1)) * origin_center_distance[pi].transpose();
			Apq += M_Matrix_C(3 * pi, 3 * pi) * tempA;
			//Apq += m_In_Group[pi] * tempA;
		}
	}
	mtRotationApq.endMyTimer();

	mtRotationSVD.startMyTimer();
	//回転以外の成分行列SはApqを極分解して得られる
	Eigen::JacobiSVD <Eigen::MatrixXd> svd(Apq, Eigen::ComputeFullU | Eigen::ComputeFullV);

	//この時刻におけるグループの回転行列Rと逆行列
	rotate_matrix = svd.matrixU() * svd.matrixV().transpose();
	mtRotationSVD.endMyTimer();

	//std::cout << rotate_matrix << std::endl;
	mtRotationSparse.startMyTimer();
	rotate_matrix_trans = rotate_matrix.transpose();
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID R" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);

	Eigen::MatrixXd rotate_matrix3N = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	//Rn_Matrix_Sparse.setZero();
	//Rn_MatrixTR_Sparse.setZero();
	// Calc rotate_matrix3N
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		rotate_matrix3N.block(3 * pi, 3 * pi, 3, 3) = rotate_matrix;
	}
	Rn_Matrix_Sparse = rotate_matrix3N.sparseView();
	Rn_MatrixTR_Sparse = (rotate_matrix3N.transpose()).sparseView();
	mtRotationSparse.endMyTimer();

	//Debug
	//かかる時間を計測する
	Apqtime += mtRotationApq.getDt();
	APDtime += mtRotationSVD.getDt();
	APDSparsetime += mtRotationSparse.getDt();
	if (APDcount == 200) {
		std::cout << "ApqTime " << tetra_group_id << " is " << std::setprecision(10) << Apqtime / 200.0 << std::endl;
		std::cout << "SVDTime " << tetra_group_id << " is " << std::setprecision(10) << APDtime / 200.0 << std::endl;
		std::cout << "SparseTime " << tetra_group_id << " is " << std::setprecision(10) << APDSparsetime / 200.0 << std::endl;
	}
	APDcount++;
}
//極分解によって回転行列を計算する
void TetraGroupD::Create_Rotate_Matrix_APD() {
	MicroSecondTimer mtRotationApq;
	MicroSecondTimer mtRotationAPD;
	MicroSecondTimer mtRotationSparse;
	mtRotationApq.setid(20);
	mtRotationAPD.setid(21);
	mtRotationSparse.setid(23);

	mtRotationApq.startMyTimer();
	//線形変換AはApqとAqqの積で表現できる
	Eigen::Matrix3d Apq = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d tempA = Eigen::Matrix3d::Zero();

	//Apqの計算
	
	if (useSparse) {
		//現在の重心の計算(x_cm^j)
		center_grid = Eigen::Vector3d::Zero();
		for (unsigned int pi = 0; pi < particle_num; pi++) {
			center_grid[0] += SUM_M_Matrix(0, 3 * pi) * PrimeVector[3 * pi];
			center_grid[1] += SUM_M_Matrix(0, 3 * pi) * PrimeVector[3 * pi + 1];
			center_grid[2] += SUM_M_Matrix(0, 3 * pi) * PrimeVector[3 * pi + 2];
		}

		//Apqの計算
		for (unsigned int pi = 0; pi < particle_num; pi++) {
			tempA = (PrimeVector.block(3*pi,0,3,1) - center_grid) * origin_center_distance[pi].transpose();
			Apq += M_Matrix_C(3*pi,3*pi) * tempA;
			//Apq += m_In_Group[pi] * tempA;
		}
	}
	else {
		//現在の重心の計算(x_cm^j)
		Eigen::VectorXd center_grid_Vector = Eigen::VectorXd::Zero(3*particle_num, 3 * particle_num);
		center_grid_Vector = SUM_M_Matrix * PrimeVector;

		//Apqの計算
		for (unsigned int pi = 0; pi < particle_num; pi++) {
			tempA = (PrimeVector.block(3 * pi, 0, 3, 1) - center_grid_Vector.block(3 * pi, 0, 3, 1)) * origin_center_distance[pi].transpose();
			Apq += M_Matrix_C(3 * pi, 3 * pi) * tempA;
			//Apq += m_In_Group[pi] * tempA;
		}
	}
	/*
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		ParticleD* pit = particles[pi];
		center_grid += (pit->Get_Mass() / Group_Mass) * pit->Get_Prime_Pos();
	}
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		ParticleD* pit = particles[pi];
		tempA = (pit->Get_Prime_Pos() - center_grid) * origin_center_distance[pi].transpose();
		Apq += pit->Get_Mass() * tempA;
	}
	*/
	mtRotationApq.endMyTimer();


	mtRotationAPD.startMyTimer();
	//回転以外の成分行列SはApqを極分解して得られる
	Eigen::Vector3d omega = Eigen::Vector3d::Identity();
	double tau = 0.1e-5;
	Eigen::Vector3d gradR = Eigen::Vector3d::Zero();
	Eigen::Matrix3d HesseR = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d S = Eigen::Matrix3d::Zero();
	int countupAPD = 0;
	//回転以外の成分行列SはApqを極分解して得られる
	//std::cout << quaternion.w() << "," << quaternion.x() << "," << quaternion.y() << "," << quaternion.z() << "," << std::endl;
	
	for (unsigned int ci = 0; ci < 20; ci++) {
		Eigen::Matrix3d R = quaternion.matrix();
		if (fetestexcept(FE_INVALID)) {
			std::cout << "FE_INVALID R" << std::endl;
		}
		feclearexcept(FE_ALL_EXCEPT);
		S = R.transpose() * Apq;
		gradR = axlAPD(S);
		//std::cout << gradR << std::endl;
		HesseR = S.trace()* Eigen::Matrix3d::Identity() - (S + S.transpose()) * 0.5;
		//std::cout << HesseR << std::endl;
		omega = -1 * HesseR.inverse() * gradR;
		//std::cout << omega<< std::endl;
		double w = omega.norm();
		if (w < 1.0e-9)
			break;
		//std::cout << "w is " << w << std::endl;
		omega = clamp2(omega, -1 * PI, PI);
		//quaternion = Eigen::Quaterniond(Eigen::AngleAxisd(w, (1.0 / w) * omega)) * quaternion;
		quaternion = quaternion * Exp2(omega);
		countupAPD++;
	}
	rotate_matrix = quaternion.matrix();
	mtRotationAPD.endMyTimer();

	mtRotationSparse.startMyTimer();
	//std::cout << "countupAPD " << countupAPD << std::endl;
	//std::cout << rotate_matrix << std::endl;
	Eigen::MatrixXd rotate_matrix3N = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	Rn_Matrix_Sparse.setZero();
	Rn_MatrixTR_Sparse.setZero();
	// Calc rotate_matrix3N
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		rotate_matrix3N.block(3 * pi, 3 * pi, 3, 3) = rotate_matrix;
	}
	Rn_Matrix_Sparse = rotate_matrix3N.sparseView();
	Rn_MatrixTR_Sparse = (rotate_matrix3N.transpose()).sparseView();
	mtRotationSparse.endMyTimer();

	//Debug
	//かかる時間を計測する
	Apqtime += mtRotationApq.getDt();
	APDtime += mtRotationAPD.getDt();
	APDSparsetime += mtRotationSparse.getDt();
	if (APDcount == 300) {
		std::cout << "ApqTime " << tetra_group_id << " is " << std::setprecision(10) << Apqtime / 300.0 << std::endl;
		std::cout << "APDTime " << tetra_group_id << " is " << std::setprecision(10) << APDtime / 300.0 << std::endl;
		std::cout << "SparseTime " << tetra_group_id << " is " << std::setprecision(10) << APDSparsetime / 300.0 << std::endl;
	}
	APDcount++;
}
void TetraGroupD::Create_Rotate_Matrix_APD_Debug(Eigen::Matrix3d temp) {
	//線形変換AはApqとAqqの積で表現できる
	Eigen::Matrix3d Apq = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d tempA = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d APDrotate_matrix = Eigen::Matrix3d::Zero();
	Eigen::Quaterniond quaterniond = Eigen::Quaterniond(Eigen::Matrix3d::Identity(3, 3));

	//Debug用
	Eigen::Matrix3d InitialR = Eigen::Matrix3d::Zero();
	InitialR(0,0) = 1.0; InitialR(1, 1) = cos(PI/4); InitialR(2, 2) = cos(PI / 4);
	InitialR(1, 2) = -1 * sin(PI / 4); InitialR(2, 1) = sin(PI / 4);
	//std::cout << InitialR << std::endl;
	quaterniond = Eigen::Quaterniond(InitialR);
	//quaterniond = Eigen::Quaterniond(temp);

	//Apqの計算
	//codeでは正規化(Mullerでやるために)した.Apq/Group\Mass さらに最大値で全体を割った
	//現在の重心の計算(x_cm^j)
	if (useSparse) {
		//現在の重心の計算(x_cm^j)
		center_grid = Eigen::Vector3d::Zero();
		for (unsigned int pi = 0; pi < particle_num; pi++) {
			center_grid[0] += SUM_M_Matrix(0, 3 * pi) * PrimeVector[3 * pi];
			center_grid[1] += SUM_M_Matrix(0, 3 * pi) * PrimeVector[3 * pi + 1];
			center_grid[2] += SUM_M_Matrix(0, 3 * pi) * PrimeVector[3 * pi + 2];
		}

		//Apqの計算
		for (unsigned int pi = 0; pi < particle_num; pi++) {
			tempA = (PrimeVector.block(3 * pi, 0, 3, 1) - center_grid) * origin_center_distance[pi].transpose();
			//Apq += tempA;
			Apq += M_Matrix_C(3 * pi, 3 * pi) * tempA;
			//Apq += (pit->Get_Mass() / Group_Mass) * tempA;
			//Apq += m_In_Group[pi] * tempA;
		}
	}
	else {
		//現在の重心の計算(x_cm^j)
		Eigen::VectorXd center_grid_Vector = Eigen::VectorXd::Zero(3 * particle_num, 3 * particle_num);
		center_grid_Vector = SUM_M_Matrix * PrimeVector;

		//Apqの計算
		for (unsigned int pi = 0; pi < particle_num; pi++) {
			tempA = (PrimeVector.block(3 * pi, 0, 3, 1) - center_grid_Vector.block(3 * pi, 0, 3, 1)) * origin_center_distance[pi].transpose();
			Apq += M_Matrix_C(3 * pi, 3 * pi) * tempA;
			//Apq += m_In_Group[pi] * tempA;
		}
	}
	//double maxApq = Apq.maxCoeff();
	//Apq = Apq / maxApq;
	//回転以外の成分行列SはApqを極分解して得られる

	Eigen::Vector3d omega = Eigen::Vector3d::Identity();
	double tau = 0.1e-5;
	Eigen::Vector3d gradR = Eigen::Vector3d::Zero();
	Eigen::Matrix3d HesseR = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d S = Eigen::Matrix3d::Zero();
	int countupAPD = 0;
	//回転以外の成分行列SはApqを極分解して得られる
	//std::cout << quaternion.w() << "," << quaternion.x() << "," << quaternion.y() << "," << quaternion.z() << "," << std::endl;

	//KKB18の方法
	
	for (unsigned int ci = 0; ci < 20; ci++) {
		Eigen::Matrix3d R = quaterniond.matrix();
		InsertAPD(rotate_matrix, R, ci);
		S = R.transpose() * Apq;
		gradR = axlAPD(S);
		HesseR = S.trace()* Eigen::Matrix3d::Identity() - (S + S.transpose()) * 0.5;
		omega = -1 * HesseR.inverse() * gradR;
		double w = omega.norm();
		if (w < 1.0e-9)
			break;
		//std::cout << "w is " << w << std::endl;
		omega = clamp2(omega, -1 * PI, PI);
		quaterniond = Eigen::Quaterniond(Eigen::AngleAxisd(w, (1.0 / w) * omega)) * quaterniond;
		//quaterniond = quaternion * Exp2(omega);
	}
	

	//Muller[MBCM16]の方法
	/*
	for (unsigned int countupAPD = 0; countupAPD < 20; countupAPD++) {
		Eigen::Matrix3d R = quaterniond.matrix();
		InsertAPD(rotate_matrix,R, countupAPD);
		Eigen::Vector3d omega = (R.col(0).cross(Apq.col(0)) + R.col
		(1).cross(Apq.col(1)) + R.col(2).cross(Apq.col(2))
			) * (1.0 / fabs(R.col(0).dot(Apq.col(0)) + R.col
			(1).dot(Apq.col(1)) + R.col(2).dot(Apq.col(2))) +
				1.0e-9);
		double w = omega.norm();
		if (w < 1.0e-9)
			break;
		quaterniond = Eigen::Quaterniond(Eigen::AngleAxisd(w, (1.0 / w) * omega)) * quaterniond;
		//quaterniond = quaterniond * Exp2(-1 * omega);
	}
	*/
	
	//std::cout << "MaxRotationVector"<< tetra_group_id <<"Group" << std::setprecision(10) << MaxRotationVector << std::endl;
	//std::cout << "MinRotationVector"<< tetra_group_id <<"Group" << std::setprecision(10) << MinRotationVector << std::endl;
	
	if(APDcount == 200){
		std::cout << "MaxRotationVector" << tetra_group_id << "Group" << std::setprecision(10) << MaxRotationVector << std::endl;
		std::cout << "MinRotationVector"<< tetra_group_id <<"Group" << std::setprecision(10) << MinRotationVector << std::endl;
	}
	int tert=1;
	if (APDcount == 201) {
		while (tert) {
		}
	}
	
	//std::cout << "countupAPD " << countupAPD << std::endl;
}
void TetraGroupD::InsertAPD(Eigen::Matrix3d p,Eigen::Matrix3d q, int c) {
	double tempK = 0.0;

	//角距離の計算
	double Theta = ((p.transpose() * q).trace() - 1) / 2;
	//std::cout << "Theta" << std::setprecision(10) << Theta<< std::endl;
	tempK = acos(Theta);
	//std::cout << "Insert" << std::setprecision(10) << tempK << std::endl;
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Angler" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);

	//代入
	if (MaxRotationVector[c] < tempK) {
		MaxRotationVector[c] = tempK;
		//std::cout << "Insert" << std::setprecision(4) << tempK<< std::endl;
	}
	if (MinRotationVector[c] > tempK) {
		MinRotationVector[c] = tempK;
		//std::cout << "Insert" << std::setprecision(4) << tempK << std::endl;
	}
}
Eigen::Vector3d TetraGroupD::axlAPD(Eigen::Matrix3d a) {
	Eigen::Vector3d g = Eigen::Vector3d::Zero();
	g[0] = a(1, 2) - a(2, 1);
	g[1] = a(2, 0) - a(0, 2);
	g[2] = a(0, 1) - a(1, 0);
	return g;
}
void TetraGroupD::Update_Rotate2() {
	//線形変換AはApqとAqqの積で表現できる
	Eigen::Matrix3d Apq = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d tempA = Eigen::Matrix3d::Zero();

	//Apqの計算
	if (useSparse) {
		//Apqの計算
		for (unsigned int pi = 0; pi < particle_num; pi++) {
			tempA = PrimeVector.block(3 * pi, 0, 3, 1) * InitialVector.block(3 * pi, 0, 3, 1).transpose();
			//Apq += tempA;
			Apq += M_Matrix_C(3 * pi, 3 * pi) * tempA;
			//Apq += m_In_Group[pi] * tempA;
		}
	}
	//回転以外の成分行列SはApqを極分解して得られる
	Eigen::Vector3d omega = Eigen::Vector3d::Identity();
	double tau = 0.1e-5;
	Eigen::Vector3d gradR = Eigen::Vector3d::Zero();
	Eigen::Matrix3d HesseR = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d S = Eigen::Matrix3d::Zero();
	int countupAPD = 0;
	//回転以外の成分行列SはApqを極分解して得られる
	//std::cout << quaternion.w() << "," << quaternion.x() << "," << quaternion.y() << "," << quaternion.z() << "," << std::endl;

	for (unsigned int ci = 0; ci < 20; ci++) {
		Eigen::Matrix3d R = quaternion.matrix();
		if (fetestexcept(FE_INVALID)) {
			std::cout << "FE_INVALID R" << std::endl;
		}
		feclearexcept(FE_ALL_EXCEPT);
		S = R.transpose() * Apq;
		gradR = axlAPD(S);
		//std::cout << gradR << std::endl;
		HesseR = S.trace()* Eigen::Matrix3d::Identity() - (S + S.transpose()) * 0.5;
		//std::cout << HesseR << std::endl;
		omega = -1 * HesseR.inverse() * gradR;
		//std::cout << omega<< std::endl;
		double w = omega.norm();
		if (w < 1.0e-9)
			break;
		//std::cout << "w is " << w << std::endl;
		omega = clamp2(omega, -1 * PI, PI);
		//quaternion = Eigen::Quaterniond(Eigen::AngleAxisd(w, (1.0 / w) * omega)) * quaternion;
		quaternion = quaternion * Exp2(omega);
		countupAPD++;
	}
	rotate_matrix = quaternion.matrix();

	//std::cout << "countupAPD " << countupAPD << std::endl;
	//std::cout << rotate_matrix << std::endl;
	Eigen::MatrixXd rotate_matrix3N = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	Rn_Matrix_Sparse.setZero();
	Rn_MatrixTR_Sparse.setZero();
	// Calc rotate_matrix3N
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		rotate_matrix3N.block(3 * pi, 3 * pi, 3, 3) = rotate_matrix;
	}
	Rn_Matrix_Sparse = rotate_matrix3N.sparseView();
	Rn_MatrixTR_Sparse = (rotate_matrix3N.transpose()).sparseView();
}
//==========================================================================//
//	@end		   				ループ設定									//
//==========================================================================//
//==========================================================================//
double TetraGroupD::Get_Volume() {
	double volume = 0.0;
	for (auto _e : elements) {
		volume += _e->Get_Volume();
	}
	return volume;
}

Eigen::Vector3d TetraGroupD::Get_Prime_Pos(int pid) {
	int element_p_id = -1;
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid)
			temp = particles[pi]->Get_Prime_Pos();
	}
	return temp;
}
Eigen::Vector3d TetraGroupD::Get_Initial_Pos(int pid) {
	int element_p_id = -1;
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid)
			temp = particles[pi]->Get_Initial_Pos();
	}
	return temp;
}
//idからグループにおけるその節点の予測位置を取得する
Eigen::Vector3d TetraGroupD::Get_Exp_Pos(int pid) {
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			temp = particles[pi]->Get_Exp_Pos();
			//std::cout << particles[pi]->Get_Exp_Pos() << std::endl;
			//std::cout << temp << std::endl;
			return temp;
		}
	} 
	return temp;
}
//idからその節点のグループでの位置を取得する
Eigen::Vector3d TetraGroupD::Get_X_In_Group(int pid) {
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			temp = x_In_Group.block(3 * pi, 0, 3, 1);
			return temp;
		}
	}
	return temp;
}
//idからその節点のグループでの解を取得する
Eigen::Vector3d TetraGroupD::Get_Deltax_In_Group(int pid) {
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			temp = Deltax_In_Group.block(3 * pi, 0, 3, 1);
			return temp;
		}
	}
	std::cout << "ERROR1153 NOT FOUND" << std::endl;
	return temp;
}
Eigen::Vector3d TetraGroupD::Get_DeltaxNew(int pid) {
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			temp = DeltaxNew.block(3 * pi, 0, 3, 1);
			return temp;
		}
	}
	std::cout << "ERROR1153 NOT FOUND" << std::endl;
	return temp;
}
//idからその節点のグループでの予測位置を取得する
Eigen::Vector3d TetraGroupD::Get_Exp_In_Group(int pid) {
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			temp = PrimeVector.block(3 * pi, 0, 3, 1);
			return temp;
		}
	}
	std::cout << "ERROR1165 NOT FOUND" << std::endl;
	return temp;
}
//idからその節点のグループでの位置を取得する
Eigen::Vector3d TetraGroupD::Get_Grid_In_Group(int pid) {
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			temp = GroupGridVector.block(3 * pi, 0, 3, 1);
			return temp;
		}
	}
	std::cout << "ERROR1177 NOT FOUND" << std::endl;
	return temp;
}
//idからその節点のグループでの速度を取得する
Eigen::Vector3d TetraGroupD::Get_Vel_In_Group(int pid) {
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			temp = GroupVelVector.block(3 * pi, 0, 3, 1);
			return temp;
		}
	}
	std::cout << "ERROR1189 NOT FOUND" << std::endl;
	return temp;
}

//idからグループにおけるその節点の質量を取得する
double TetraGroupD::Get_GMass_In_Group(int pid) {
	double tgmass = 0.0;
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			tgmass = m_In_Group[pi];
			return tgmass;
		}
	}
	return tgmass;
}
//idから全体を考慮したその節点の質量を取得する
double TetraGroupD::Get_Mass_In_Group(int pid) {
	double tmass = 0.0;
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			tmass = particles[pi]->Get_Mass();
			return tmass;
		}
	}
	return tmass;
}

//節点idからグループにおけるその節点のindexを取得する
unsigned int TetraGroupD::Get_Group_Index(int pid){
	int element_p_id = -1;
	unsigned int index = 0;
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			index = pi;
			return index;
		}
	}
	return index;
}

void TetraGroupD::Set_Group_Mass(double a) {
	this->Group_Mass = a;
}
double TetraGroupD::Get_Group_Mass() {
	return this->Group_Mass;
}
/*
//idからexp_posを変数のベクトル値に変更する
void TetraGroupD::Update_CoExp_Pos(int pid,Eigen::Vector3d a) {
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			particles[pi]->Set_Exp_Pos(a);
		}
	}
}
*/

void TetraGroupD::ReSet_Fbind_Pos() {
	//std::cout << bind_force_iterative << std::endl;
	bind_force_iterative = Eigen::VectorXd::Zero(3* particles.size());
}

//グループのparticleを取得
std::vector<ParticleD*> TetraGroupD::Get_Particle()const {
	return particles;
}
//各particleにおける初めの重心からの距離を取得
std::vector<Eigen::Vector3d> TetraGroupD::Get_origin_center_distance() {
	return origin_center_distance;
}
//各particleのローカル座標初期位置を取得
std::vector<Eigen::Vector3d> TetraGroupD::Get_origin_local_grid() {
	return origin_local_grid;
}
//各particleのローカル座標初期位置を取得
Eigen::VectorXd TetraGroupD::Get_bind_force() {
	return bind_force_iterative;
}
void TetraGroupD::Write_bind_force() {
	std::cout << bind_force_iterative << std::endl;
}
//==========================================================================//

//各変数の要素数を決定する(exp_pos,origin_center_distance,center_distance,origin_local_grid)
void TetraGroupD::Set_Size_para(int particle_num) {
	/*exp_pos_g.resize(particle_num);
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		exp_pos_g[pi] = particles[pi]->Get_Exp_Pos();
	}
	initial_pos_g.resize(particle_num);
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		initial_pos_g[pi] = particles[pi]->Get_Initial_Pos();
	}
	prime_pos_g.resize(particle_num);
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		prime_pos_g[pi] = particles[pi]->Get_Prime_Pos();
	}*/
	origin_center_distance.resize(particle_num);
	center_distance.resize(particle_num);
	origin_local_grid.resize(particle_num);
}
//各変数の行列のサイズを決定する(R_Matrix,trans_R_Matrix,M_Matrix)
void TetraGroupD::Set_Size_para2(std::vector< ParticleD* > particles) {
	R_Matrix = Eigen::MatrixXd::Identity(3 * particles.size(), 3 * particles.size());
	trans_R_Matrix = Eigen::MatrixXd::Identity(3 * particles.size(), 3 * particles.size());
	M_Matrix = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	M_Matrix_C = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	stiffness_matrix = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	SUM_M_Matrix = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	Damping_Matrix = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	MassDamInv_Matrix = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());

	stiffmatrix_valued_list.resize(particles.size());
	stiffmatrix_valued_list_sym.resize(particles.size());

	pLocal = Eigen::VectorXd(3 * particles.size());
	f_inLocal = Eigen::VectorXd(3 * particles.size());
	xgLocal = Eigen::VectorXd(3 * particles.size());
	centerLocal = Eigen::VectorXd(3 * particles.size());
	f_Local = Eigen::VectorXd(3 * particles.size());
	v_Local = Eigen::VectorXd(3 * particles.size());
	x_Local = Eigen::VectorXd(3 * particles.size());

	bind_force_iterative = Eigen::VectorXd::Zero(3 * particles.size());
	x_In_Group = Eigen::VectorXd::Zero(3 * particles.size());
	m_In_Group = Eigen::VectorXd::Zero(particles.size());
	Deltax_In_Group = Eigen::VectorXd::Zero(3 * particles.size());
	Deltax_CoFEM = Eigen::VectorXd::Zero(3 * particle_num);
	Deltax_Bind = Eigen::VectorXd::Zero(3 * particle_num);

	Jacobi_Matrix = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	DiagFEM_Matrix_iteration = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	F_FEM_Matrix_iteration = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	Constant_term_iteration = Eigen::VectorXd::Zero(3 * particles.size());
	PrimeVector = Eigen::VectorXd::Zero(3 * particles.size());
	OrigineVector = Eigen::VectorXd::Zero(3 * particles.size());
	InitialVector = Eigen::VectorXd::Zero(3 * particles.size());
	iterativeVector = Eigen::VectorXd::Zero(3 * particles.size());

	MaxRotationVector = Eigen::VectorXd::Zero(20);
	MinRotationVector = Eigen::VectorXd::Zero(20);
	for (int i = 0; i < 20;i++) {
		MinRotationVector[i] = 100.0;
	}
	//std::cout << MinRotationVector << std::endl;
}

//極分解でのニュートン法関連の関数
Eigen::Vector3d TetraGroupD::clamp2(Eigen::Vector3d x, double y, double z) {
	if (x.norm() < y) {
		return (y / x.norm()) * x;
	}
	else if (x.norm() > z) {
		return (z / x.norm()) * x;
	}
	else {
		return x;
	}
}
Eigen::Quaterniond TetraGroupD::Cay2(Eigen::Vector3d a) {
	double s = 0.25*a.transpose()*a;
	double x = (1.0 / (1.0 + s))*a.x();
	double y = (1.0 / (1.0 + s))*a.y();
	double z = (1.0 / (1.0 + s))*a.z();
	Eigen::Quaterniond qq = Eigen::Quaterniond((1.0 - s)/(1.0 + s), x, y, z);
	return  qq;
}
Eigen::Quaterniond TetraGroupD::Cay(Eigen::Vector3d a) {
	double x = (2 / (1 + a.transpose()*a))*a.x();
	double y = (2 / (1 + a.transpose()*a))*a.y();
	double z = (2 / (1 + a.transpose()*a))*a.z();
	Eigen::Quaterniond qq = Eigen::Quaterniond((1 - a.transpose()*a), x, y, z);
	return  qq;
}
Eigen::Quaterniond TetraGroupD::Exp(Eigen::Vector3d a) {
	double s = sin(a.norm());
	double x = s * a.x();
	double y = s * a.y();
	double z = s * a.z();
	Eigen::Quaterniond qq = Eigen::Quaterniond(cos(a.norm()), x, y, z);
	return  qq;
}
Eigen::Quaterniond TetraGroupD::Exp2(Eigen::Vector3d a) {
	double s = sin((a * 0.5).norm());
	double x = s * a.x() / a.norm();
	double y = s * a.y() / a.norm();
	double z = s * a.z() / a.norm();
	Eigen::Quaterniond qq = Eigen::Quaterniond(cos((a * 0.5).norm()), x, y, z);
	return  qq;
}
//implicitの計算
//CRS形式を使わずに計算する
//省略法の反復におけるLocal制約の計算
//fulpivotで連立方程式を解き反復する
//弾性力以外の力による位置を引いてから計算する
//解いた解は位置
void TetraGroupD::Calc_iterative_FEM_Fbind_pivot() {
	//k回目の変位ベクトルの生成
	Eigen::VectorXd  vector_u = Eigen::VectorXd::Zero(3 * particle_num);

	//FEM計算
	Eigen::FullPivLU<Eigen::MatrixXd> lu(Jacobi_Matrix);
	vector_u = lu.solve(Constant_term_iteration + bind_force_iterative);

	for (unsigned int pi = 0; pi < particle_num; pi++) {
		x_In_Group.block(3 * pi, 0, 3, 1) = vector_u.block(3 * pi, 0, 3, 1);
	}

	//中身の確認debug
	//OP_debug_iterative1(vector_u);
}
//implicitの計算
//CRS形式を使わずに計算する
//省略法の反復におけるLocal制約の計算
//fulpivotで連立方程式を解き反復する
//弾性力以外の力による位置を引いてから計算する
//解いた解は変位
void TetraGroupD::Calc_iterative_FEM_Fbind_pivot2() {
	//k回目の変位ベクトルの生成
	Eigen::VectorXd  vector_u = Eigen::VectorXd::Zero(3 * particle_num);

	//FEM計算
	//GMRESなどの反復法を用いて反復法(Local)を解く
	if (useGMRES) {
		//GMRES

		//前処理なし
		
		Eigen::GMRES< Eigen::SparseMatrix<double>> gmresFEM;
		gmresFEM.setMaxIterations(1);//外部反復の設定
		gmresFEM.set_restart(1);//内部反復の設定
		gmresFEM.compute(Jacobi_Matrix_Sparse);// Compute 

		//初期値0
		//vector_u = gmresFEM.solve(Constant_term_iteration + bind_force_iterative);
		//初期値は一つ前の値から計算したものをいれる
		vector_u = gmresFEM.solveWithGuess(Constant_term_iteration + bind_force_iterative, iterativeVector);
		//std::cout << "GMRES前なし：#iterations：" << gmresFEM.iterations() << "、推定エラー：" << gmresFEM.error() << std::endl;
		

		//前処理あり
		/*
		Eigen::GMRES< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > gmresFEM2;
		gmresFEM2.preconditioner().setFillfactor(7);  //Get the reference of the preconditioner and set properties
		//gmresFEM2.setTolerance(10e-3);//許容値の設定
		gmresFEM2.setMaxIterations(1);//外部反復の設定
		gmresFEM2.set_restart(1);//内部反復の設定
		gmresFEM2.compute(Jacobi_Matrix_Sparse);// Compute the ILUT factorization

		//初期値0
		//vector_u = gmresFEM2.solve(Constant_term_iteration + bind_force_iterative);
		//初期値は一つ前の値から計算したものをいれる
		vector_u = gmresFEM2.solveWithGuess(Constant_term_iteration + bind_force_iterative, iterativeVector);
		//std::cout << "GMRES前あり：#iterations：" << gmresFEM2.iterations() << "、推定エラー：" << gmresFEM2.error() << std::endl;
		*/
	}
	//LU分解を用いて反復法(Local)を解く
	else {
		Eigen::FullPivLU<Eigen::MatrixXd> lu(Jacobi_Matrix);
		vector_u = lu.solve(Constant_term_iteration + bind_force_iterative);
		//計算上はタイムステップをかける必要はない
		//vector_u = lu.solve(Constant_term_iteration + bind_force_iterative * TIME_STEP * TIME_STEP);
	}
	


	//std::cout << "オブジェクトの変位ベクトル is " << std::endl;
	//std::cout << x_In_Group << std::endl;

	//変位ベクトルを足す
	//本当はvectoruとDeltax_IN_Groupを同じにしてもよい

	for (unsigned int pi = 0; pi < particle_num; pi++) {
		x_In_Group.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Prime_Pos() + vector_u.block(3 * pi, 0, 3, 1);
		Deltax_In_Group.block(3 * pi, 0, 3, 1) = vector_u.block(3 * pi, 0, 3, 1);
	}
	
	//中身の確認debug
	//OP_debug_iterative1(vector_u);
}
//implicitの計算
//CRS形式を使わずに計算する
//反復におけるLocal制約の計算
//fulpivotとGMRESで連立方程式を解き反復する
//弾性力以外の力による位置を引いてから計算する
//解いた解は変位
void TetraGroupD::Calc_iterative_LocalFEM() {
	//k回目の変位ベクトルの生成
	//Eigen::VectorXd vector_u = Eigen::VectorXd::Zero(3 * particle_num);
	Deltax_In_Group = Eigen::VectorXd::Zero(3 * particle_num);
	//FEM計算
	//GMRESなどの反復法を用いて反復法(Local)を解く
	if (useGMRES) {
		//GMRES

		//前処理なし
		if (!usePreIte) {
			Eigen::GMRES< Eigen::SparseMatrix<double>> gmresFEM;
			gmresFEM.setMaxIterations(outerGMRES);//外部反復の設定
			gmresFEM.set_restart(innerGMRES);//内部反復の設定
			gmresFEM.compute(Jacobi_Matrix_Sparse);// Compute 

			//初期値0
			//vector_u = gmresFEM.solve(Constant_term_iteration + bind_force_iterative);

			//初期値は一つ前の値から計算したものをいれる
			MicroSecondTimer mtGMRESReal;
			mtGMRESReal.setid(32);
			mtGMRESReal.startMyTimer();
			Deltax_In_Group = gmresFEM.solveWithGuess(Constant_term_iteration + bind_force_iterative, iterativeVector);
			mtGMRESReal.endMyTimer();
			//std::cout << "Real GMRES is " << mtGMRESReal.getDt() << std::endl;
			//std::cout << "GMRES前なし：#iterations：" << gmresFEM.iterations() << "、推定エラー：" << gmresFEM.error() << std::endl;
		}
		else {
			//前処理あり
			
			Eigen::GMRES< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > gmresFEM;
			gmresFEM.preconditioner().setFillfactor(7);  //Get the reference of the preconditioner and set properties
			//gmresFEM.setTolerance(10e-3);//許容値の設定
			gmresFEM.setMaxIterations(outerGMRES);//外部反復の設定
			gmresFEM.set_restart(innerGMRES);//内部反復の設定
			gmresFEM.compute(Jacobi_Matrix_Sparse);// Compute the ILUT factorization

			//初期値0
			//vector_u = gmresFEM2.solve(Constant_term_iteration + bind_force_iterative);

			//初期値は一つ前の値から計算したものをいれる
			Deltax_In_Group = gmresFEM.solveWithGuess(Constant_term_iteration + bind_force_iterative, iterativeVector);
			//std::cout << "GMRES前あり：#iterations：" << gmresFEM2.iterations() << "、推定エラー：" << gmresFEM2.error() << std::endl;
		}
	}
	//LU分解を用いて反復法(Local)を解く
	else {
		Eigen::FullPivLU<Eigen::MatrixXd> lu(Jacobi_Matrix);
		Deltax_In_Group = lu.solve(Constant_term_iteration + bind_force_iterative);
	}

	//std::cout << "vectoru" << Deltax_In_Group << std::endl;

	//変位ベクトルを足す
	//本当はvectoruとDeltax_IN_Groupを同じにしてもよい

	/*
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		//x_In_Group.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Prime_Pos() + vector_u.block(3 * pi, 0, 3, 1);
		Deltax_In_Group.block(3 * pi, 0, 3, 1) = vector_u.block(3 * pi, 0, 3, 1);
	}
	std::cout << "vectoru" << vector_u << std::endl;
	*/

	//中身の確認debug
	//OP_debug_iterative1(vector_u);
}
//GMRESの前処理を行う
void TetraGroupD::Calc_GMRES_Pre() {
	Deltax_CoFEM = Eigen::VectorXd::Zero(3 * particle_num);
	//FEM計算
	//GMRESなどの反復法を用いて反復法(Local)を解く
	if (useGMRES) {
		//GMRES
		//前処理なし
		if (!usePreIte) {
			gmresFEM_Pre.setMaxIterations(outerGMRES);//外部反復の設定
			gmresFEM_Pre.set_restart(innerGMRES);//内部反復の設定
			gmresFEM_Pre.compute(Jacobi_Matrix_Sparse);// Compute
			Deltax_CoFEM = gmresFEM_Pre.solve(Constant_term_iteration);
		}
		else {
			//前処理あり
			gmresFEM_Pre2.preconditioner().setFillfactor(7);  //Get the reference of the preconditioner and set properties
			//gmresFEM.setTolerance(10e-3);//許容値の設定
			gmresFEM_Pre2.setMaxIterations(outerGMRES);//外部反復の設定
			gmresFEM_Pre2.set_restart(innerGMRES);//内部反復の設定
			gmresFEM_Pre2.compute(Jacobi_Matrix_Sparse);// Compute the ILUT factorization
			Deltax_CoFEM = gmresFEM_Pre2.solve(Constant_term_iteration);
		}
	}
}
//implicitの計算
//CRS形式を使わずに計算する
//反復におけるLocal制約の計算
//弾性力以外の力による位置を引いてから計算する
//解いた解は変位
//前処理済
void TetraGroupD::Calc_GMRES_FEM() {
	//解ベクトル
	Deltax_In_Group = Eigen::VectorXd::Zero(3 * particle_num);
	Deltax_Bind = Eigen::VectorXd::Zero(3 * particle_num);
	if (useGMRES) {
		//GMRES
		//前処理なし
		if (!usePreIte) {
			MicroSecondTimer mtGMRESReal;
			mtGMRESReal.setid(32);
			mtGMRESReal.startMyTimer();
			Deltax_Bind = gmresFEM_Pre.solve(bind_force_iterative);
			//Deltax_In_Group = gmresFEM_Pre.solveWithGuess(Constant_term_iteration + bind_force_iterative, iterativeVector);
			Deltax_In_Group = Deltax_CoFEM + Deltax_Bind;
			
			mtGMRESReal.endMyTimer();
		}
		else {
			//前処理あり
			MicroSecondTimer mtGMRESReal;
			mtGMRESReal.setid(32);
			mtGMRESReal.startMyTimer();
			Deltax_Bind = gmresFEM_Pre2.solve(bind_force_iterative);
			//Deltax_Bind = gmresFEM_Pre2.solveWithGuess(bind_force_iterative, iterativeVector);
			//Deltax_In_Group = gmresFEM_Pre2.solveWithGuess(Constant_term_iteration + bind_force_iterative, iterativeVector);
			Deltax_In_Group = Deltax_CoFEM + Deltax_Bind;
			//std::cout << "Deltax in group" << std::endl << Deltax_In_Group << std::endl;
			mtGMRESReal.endMyTimer();
		}
	}
}
//implicitの計算
//CRS形式を使わずに計算する
//省略法の反復におけるLocal制約の計算
//ヤコビ法で連立方程式を解き反復する
//弾性力以外の力による位置を引いてから計算する
//解いた解は変位
void TetraGroupD::Calc_iterative_FEM_Fbind_Jacobi2() {
	ParticleD* pit;
	Eigen::VectorXd ExpJacobiVector = Eigen::VectorXd::Zero(3 * particle_num);

	//k回目のDelta_EXPベクトル生成
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		ExpJacobiVector.block(3 * pi, 0, 3, 1) = pit->Get_Exp_Pos()- pit->Get_Prime_Pos();
	}

	//k回目の変位ベクトルの生成
	Eigen::VectorXd  vector_u = Eigen::VectorXd::Zero(3 * particle_num);
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		vector_u.block(3 * pi, 0, 3, 1) = pit->Get_Exp_Pos() - pit->Get_Prime_Pos();
	}

	//固有値の計算
	Eigen::EigenSolver< Eigen::MatrixXd > s(F_FEM_Matrix_iteration + DiagFEM_Matrix_iteration);
	//固有値
	//std::cout << "固有値\n" << s.eigenvalues() << std::endl;
	//std::cout << "固有ベクトル\n" << s.eigenvectors() << std::endl;
	//Detのチェック
	//std::cout << "Det A = "<< (F_FEM_Matrix_iteration + DiagFEM_Matrix_iteration).determinant() << std::endl;
	//対角優位かチェック
	/*for (unsigned int pii = 0; pii < 3 * particle_num; pii++) {
		double cccout = 0;
		for (unsigned int pij = 0; pij < 3 * particle_num; pij++) {
			cccout = cccout + abs(F_FEM_Matrix_iteration(pii, pij));
		}
		if (abs(DiagFEM_Matrix_iteration(pii, pii))> cccout ) {
			std::cout << "diagonal advantage" << std::endl;
		}
		else {
			std::cout << "no" << std::endl;
		}
		cccout = 0;
	}*/

	//FEM計算
	//Eigen::FullPivLU<Eigen::MatrixXd> lu(Jacobi_Matrix);
	//vector_u = lu.solve(Constant_term_iteration + bind_force_iterative * TIME_STEP * TIME_STEP);

	//vector_u = DiagFEM_Matrix_iteration.inverse()* (Constant_term_iteration + bind_force_iterative * TIME_STEP * TIME_STEP - F_FEM_Matrix_iteration * ExpJacobiVector);


	//ここでも反復したほうがいいのかな？
	for (int pi = 0; pi < 10; pi++) {
		vector_u = DiagFEM_Matrix_iteration.inverse()* (Constant_term_iteration + bind_force_iterative * TIME_STEP * TIME_STEP - F_FEM_Matrix_iteration * vector_u);
		//std::cout << "Jacobi itera is " << std::endl;
		//std::cout << DiagFEM_Matrix_iteration.inverse()* F_FEM_Matrix_iteration * vector_u << std::endl;
	}
	

	//変位ベクトルを足す
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		x_In_Group.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Prime_Pos() + vector_u.block(3 * pi, 0, 3, 1);
	}


	//中身の確認debug
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		//std::cout << "Constant D of  u" << "particle " << pit->p_id << "group is " << tetra_group_id << std::endl;
		//std::cout << x_In_Group.block(3 * pi, 0, 3, 1) << std::endl;
		//std::cout << vector_u.block(3 * pi, 0, 3, 1) << std::endl;
	}
}
//Fbindを更新する
void TetraGroupD::Update_Fbind_Pos() {
	ParticleD* pit;
	Eigen::VectorXd ExpVector = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::VectorXd ExpPreVector = Eigen::VectorXd::Zero(3 * particle_num);

	//expベクトルの生成
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		ExpVector.block(3 * pi, 0, 3, 1) = pit->Get_Exp_Pos();
		ExpPreVector.block(3 * pi, 0, 3, 1) = pit->Get_Prime_Pos();
	}
	/*
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		if (!(pit)->Is_Fixed()) {
			//bind_force_iterative.block(3 * pi, 0, 3, 1) = bind_force_iterative.block(3 * pi, 0, 3, 1) + F_bind_coeff * (pit->Get_Exp_Pos() - x_In_Group.block(3 * pi, 0, 3, 1));
			//bind_force_iterative.block(3 * pi, 0, 3, 1) = bind_force_iterative.block(3 * pi, 0, 3, 1) + Jacobi_Matrix * (pit->Get_Exp_Pos() - x_In_Group.block(3 * pi, 0, 3, 1));
		}
	}
	*/
	//加算していく方法
	bind_force_iterative = bind_force_iterative + F_bind_coeff * Jacobi_Matrix * (ExpVector - x_In_Group);
	//bind_force_iterative = bind_force_iterative + F_bind_coeff * Jacobi_Matrix * (DeltaxVector - Deltax_In_Group);
	
	//Delta xを次の反復法の初期値にする
	iterativeVector = ExpVector - ExpPreVector;
	//iterativeVector = DeltaxVector ;

	//一気に計算していく方法
	//bind_force_iterative = (Jacobi_Matrix * ExpVector) - Constant_term_iteration;
	//std::cout << "bindforce is " << std::endl;
	//std::cout << bind_force_iterative << std::endl;
}
//Fbindを更新する(差を利用)
void TetraGroupD::Update_Fbind_Pos2() {
	Eigen::VectorXd Deltax_In_Model_Vector = Eigen::VectorXd::Zero(3 * particle_num);
	//Deltax_modelベクトルの生成
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		Deltax_In_Model_Vector.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Deltax_In_Model();
	}
	//加算していく方法
	if (useSparse) {
		bind_force_iterative += F_bind_coeff * Jacobi_Matrix_Sparse * (Deltax_In_Model_Vector - Deltax_In_Group);
	}
	else {
		bind_force_iterative += F_bind_coeff * Jacobi_Matrix * (Deltax_In_Model_Vector - Deltax_In_Group);
	}
	if (!useUpBind) {
		bind_force_iterative = Eigen::VectorXd::Zero(3 * particle_num);
	}
	//bind_force_iterative = bind_force_iterative + F_bind_coeff * Jacobi_Matrix * (DeltaxVector - Deltax_In_Group);
	//Delta xを次の反復法の初期値にする
	iterativeVector = Deltax_In_Model_Vector;
	//std::cout << "Deltax_In_Model_Vector" << std::endl;
	//std::cout << Deltax_In_Model_Vector << std::endl;
}
//Fbindを更新する(差を利用,予測位置はグループで異なる)
void TetraGroupD::Update_Fbind_Pos3() {
	Eigen::VectorXd Deltax_In_Model_Vector = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::VectorXd Exp_In_Model_Vector = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::VectorXd Ago_Vector = Eigen::VectorXd::Zero(3 * particle_num);
	//Deltax_modelベクトルの生成
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		Deltax_In_Model_Vector.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Deltax_In_Model();
	}
	//予測位置ベクトルの生成(平均)
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		Exp_In_Model_Vector.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Exp_Pos();
	}
	//現在の位置ベクトルの生成
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		Ago_Vector.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Grid();
	}
	//加算
	if (useSparse) {
		bind_force_iterative += F_bind_coeff * Jacobi_Matrix_Sparse * (Deltax_In_Model_Vector - Deltax_In_Group + Exp_In_Model_Vector - PrimeVector);
		//bind_force_iterative += (F_bind_damping / TIME_STEP)*(Deltax_In_Model_Vector+ Exp_In_Model_Vector - Ago_Vector);
		//加算じゃないほう
		//bind_force_iterative = F_bind_coeff * Jacobi_Matrix_Sparse * (Deltax_In_Model_Vector - Deltax_In_Group + Exp_In_Model_Vector - PrimeVecto)+(F_bind_damping / TIME_STEP)*(Deltax_In_Model_Vector + Exp_In_Model_Vector - Ago_Vector);
	}
	else {
		bind_force_iterative += F_bind_coeff * Jacobi_Matrix * (Deltax_In_Model_Vector - Deltax_In_Group + Exp_In_Model_Vector - PrimeVector);
	}
	if (!useUpBind) {
		bind_force_iterative = Eigen::VectorXd::Zero(3 * particle_num);
	}
	if (tetra_group_id==2) {
		//std::cout << (Deltax_In_Model_Vector + Exp_In_Model_Vector - Ago_Vector) << std::endl;
		//std::cout << bind_force_iterative << std::endl;
	}
	//Delta xを次の反復法の初期値にする
	iterativeVector = Deltax_In_Model_Vector + Exp_In_Model_Vector - PrimeVector;
}
void TetraGroupD::Update_Fbind_Pos4() {
	for (unsigned int pi = 0; pi < particle_num;pi++) {
		//共有節点かどうか
		if ((particles[pi]->p_belong_TetraGroup_ids.size())>1) {
			//固定されていない点
			if (!(particles[pi]->Is_Fixed())) {
				Eigen::Vector3d Conv = Eigen::Vector3d::Zero();
				//(n)(Exp + Deltax)
				Conv = particles[pi]->p_belong_TetraGroup_ids.size() * (PrimeVector.block(3 * pi, 0, 3, 1) + Deltax_In_Group.block(3 * pi, 0, 3, 1));
				//std::cout << particles[pi]->p_belong_TetraGroup_ids.size() << std::endl;
				//(n+1)(Exp + Deltax) - n(MeanExp + MeanDeltax)
				Conv += -1 * particles[pi]->p_belong_TetraGroup_ids.size() * (particles[pi]->Get_Exp_Pos() + particles[pi]->Get_Deltax_In_Model());
				bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;

				//速度の計算
				//bind_force_iterative.block(3 * pi, 0, 3, 1) += (F_bind_damping / TIME_STEP) * Conv;
			}
			//固定点の場合
			else {
				Eigen::Vector3d Conv = Eigen::Vector3d::Zero();
				//(Exp + Deltax - Grid)
				Conv = PrimeVector.block(3 * pi, 0, 3, 1)+ Deltax_In_Group.block(3 * pi, 0, 3, 1)- particles[pi]->Get_Exp_Pos();
				bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;

				//速度の計算
				//bind_force_iterative.block(3 * pi, 0, 3, 1) += (F_bind_damping / TIME_STEP) * Conv;
			}
		}
		//固定点の場合
		else {
			if ((particles[pi]->Is_Fixed())) {
				Eigen::Vector3d Conv = Eigen::Vector3d::Zero();
				//(Exp + Deltax)
				Conv = PrimeVector.block(3 * pi, 0, 3, 1) + Deltax_In_Group.block(3 * pi, 0, 3, 1) - particles[pi]->Get_Exp_Pos();
				bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;

				//速度の計算
				//bind_force_iterative.block(3 * pi, 0, 3, 1) += (F_bind_damping / TIME_STEP) * Conv;
			}
			else {
				bind_force_iterative.block(3 * pi, 0, 3, 1) = Eigen::Vector3d::Zero();
			}
		}
		
	}
	if (!useUpBind) {
		bind_force_iterative = Eigen::VectorXd::Zero(3 * particle_num);
	}
	std::cout <<"Bind" <<tetra_group_id<<" is "<<bind_force_iterative << std::endl;
}
//予測位置が同じとき
void TetraGroupD::Update_Fbind_Pos5() {
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		//共有節点かどうか
		//本当は分離しないといけない
		if ((particles[pi]->p_belong_TetraGroup_ids.size()) > 1) {
				Eigen::Vector3d Conv = Eigen::Vector3d::Zero();
				//(n)(Exp + Deltax)
				Conv = particles[pi]->p_belong_TetraGroup_ids.size() * Deltax_In_Group.block(3 * pi, 0, 3, 1);
				//std::cout << particles[pi]->p_belong_TetraGroup_ids.size() << std::endl;
				Conv = Conv - particles[pi]->p_belong_TetraGroup_ids.size() * particles[pi]->Get_Deltax_In_Model();
				bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;
				//std::cout << "Bind" << particles[pi]->p_id << "of " << tetra_group_id << " is " << std::endl;
				//std::cout<< bind_force_iterative.block(3 * pi, 0, 3, 1) << std::endl;
		}
	}
	if (!useUpBind) {
		bind_force_iterative = Eigen::VectorXd::Zero(3 * particle_num);
	}
	//std::cout << "Bind" << tetra_group_id << " is " << bind_force_iterative << std::endl;
}
//Debug用
//節点も同じように処理
void TetraGroupD::Update_Fbind_Pos6() {
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		Eigen::Vector3d Conv = Eigen::Vector3d::Zero();
		//共有節点かどうか
		//本当は分離しないといけない
		if ((particles[pi]->p_belong_TetraGroup_ids.size()) > 1) {
			//固定されていない点
			if (!(particles[pi]->Is_Fixed())) {
				//(n)(Exp + Deltax)
				Conv = particles[pi]->p_belong_TetraGroup_ids.size() * (PrimeVector.block(3 * pi, 0, 3, 1) + Deltax_In_Group.block(3 * pi, 0, 3, 1));
				//std::cout << particles[pi]->p_belong_TetraGroup_ids.size() << std::endl;
				Conv = Conv - particles[pi]->p_belong_TetraGroup_ids.size() *  (particles[pi]->Get_Exp_Pos() + particles[pi]->Get_Deltax_In_Model());
				bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;
				//std::cout << "Bind" << particles[pi]->p_id << "of " << tetra_group_id << " is " << std::endl;
				//std::cout<< bind_force_iterative.block(3 * pi, 0, 3, 1) << std::endl;
			}
		}
		//固定点の場合
		if ((particles[pi]->Is_Fixed())) {
			//(Exp + Deltax)
			Conv = (PrimeVector.block(3 * pi, 0, 3, 1) + Deltax_In_Group.block(3 * pi, 0, 3, 1));
			//std::cout << particles[pi]->p_belong_TetraGroup_ids.size() << std::endl;
			Conv = Conv - (particles[pi]->Get_Exp_Pos());
			bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;
			//std::cout << "Bind" << particles[pi]->p_id << "of " << tetra_group_id << " is " << std::endl;
			//std::cout << bind_force_iterative.block(3 * pi, 0, 3, 1) << std::endl;
		}
		if ( ((particles[pi]->p_belong_TetraGroup_ids.size()) > 1) ) {
			if (Conv.squaredNorm() < 10e-3) {
				//std::cout << "Converce " << particles[pi]->p_id << std::endl;
			}
		}
	}
	//std::cout << "Bind" << tetra_group_id << " is " << bind_force_iterative << std::endl;
}

void TetraGroupD::Update_Fbind_Pos8() {
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		Eigen::Vector3d Conv = Eigen::Vector3d::Zero();
		//共有節点かどうか
		//本当は分離しないといけない
		if ((particles[pi]->p_belong_TetraGroup_ids.size()) > 1) {
			//固定されていない点
			if (!(particles[pi]->Is_Fixed())) {
				//(n)(Exp + Deltax)
				Conv = particles[pi]->p_belong_TetraGroup_ids.size() * (PrimeVector.block(3 * pi, 0, 3, 1) + DeltaxNew.block(3 * pi, 0, 3, 1));
				//std::cout << particles[pi]->p_belong_TetraGroup_ids.size() << std::endl;
				Conv = Conv - particles[pi]->p_belong_TetraGroup_ids.size() * (particles[pi]->Get_Exp_Pos() + particles[pi]->Get_Deltax_In_Model());
				bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;
				if (fetestexcept(FE_INVALID)) {
					std::cout << "FE_INVALID bindF" << std::endl;
				}
				feclearexcept(FE_ALL_EXCEPT);
				//std::cout << "Bind" << particles[pi]->p_id << "of " << tetra_group_id << " is " << std::endl;
				//std::cout<< bind_force_iterative.block(3 * pi, 0, 3, 1) << std::endl;
			}
		}
		//固定点の場合
		if ((particles[pi]->Is_Fixed())) {
			//(Exp + Deltax)
			Conv = (PrimeVector.block(3 * pi, 0, 3, 1) + Deltax_In_Group.block(3 * pi, 0, 3, 1));
			//std::cout << particles[pi]->p_belong_TetraGroup_ids.size() << std::endl;
			Conv = Conv - (particles[pi]->Get_Exp_Pos());
			bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID BindF1" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
			//std::cout << "Bind" << particles[pi]->p_id << "of " << tetra_group_id << " is " << std::endl;
			//std::cout << bind_force_iterative.block(3 * pi, 0, 3, 1) << std::endl;
		}
		if (((particles[pi]->p_belong_TetraGroup_ids.size()) > 1)) {
			if (Conv.squaredNorm() < 10e-3) {
				//std::cout << "Converce " << particles[pi]->p_id << std::endl;
			}
		}
	}
}
Eigen::Vector3d TetraGroupD::Calc_Distance() {
	Eigen::Vector3d Conv = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == 19) {
			Conv = particles[pi]->p_belong_TetraGroup_ids.size() * (PrimeVector.block(3 * pi, 0, 3, 1) + Deltax_In_Group.block(3 * pi, 0, 3, 1));
			Conv = Conv - particles[pi]->p_belong_TetraGroup_ids.size() *  (particles[pi]->Get_Exp_Pos() + particles[pi]->Get_Deltax_In_Model());
		}
	}
	return Conv;
}
//Debug用
//節点は別処理
void TetraGroupD::Update_Fbind_Pos7() {
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		Eigen::Vector3d Conv = Eigen::Vector3d::Zero();
		//共有節点かどうか
		//本当は分離しないといけない
		if ((particles[pi]->p_belong_TetraGroup_ids.size()) > 1) {
			//固定されていない点
			if (!(particles[pi]->Is_Fixed())) {
				//(n)(Exp + Deltax)
				Conv = particles[pi]->p_belong_TetraGroup_ids.size() * (PrimeVector.block(3 * pi, 0, 3, 1) + Deltax_In_Group.block(3 * pi, 0, 3, 1));
				//std::cout << particles[pi]->p_belong_TetraGroup_ids.size() << std::endl;
				Conv = Conv - particles[pi]->p_belong_TetraGroup_ids.size() *  (particles[pi]->Get_Exp_Pos() + particles[pi]->Get_Deltax_In_Model());
				bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;
				//std::cout << "Bind" << particles[pi]->p_id << "of " << tetra_group_id << " is " << std::endl;
				//std::cout<< bind_force_iterative.block(3 * pi, 0, 3, 1) << std::endl;
			}
		}
		//固定点の場合
		if ((particles[pi]->Is_Fixed())) {
			//(Exp + Deltax)
			Conv = Deltax_In_Group.block(3 * pi, 0, 3, 1);
			bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;
			//std::cout << "Bind" << particles[pi]->p_id << "of " << tetra_group_id << " is " << std::endl;
			//std::cout << bind_force_iterative.block(3 * pi, 0, 3, 1) << std::endl;
		}
		if (((particles[pi]->p_belong_TetraGroup_ids.size()) > 1)) {
			if (Conv.squaredNorm() < 10e-3) {
				//std::cout << "Converce " << particles[pi]->p_id << std::endl;
			}
		}
	}
	//std::cout << "Bind" << tetra_group_id << " is " << bind_force_iterative << std::endl;
}
void TetraGroupD::Calc_Convergence() {
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		Eigen::Vector3d WER = Eigen::Vector3d::Zero();
		WER = PrimeVector.block(3 * pi, 0, 3, 1) + Deltax_In_Group.block(3 * pi, 0, 3, 1) - (particles[pi]->Get_Deltax_In_Model() + particles[pi]->Get_Exp_Pos());
		//std::cout << particles[pi]->p_id << " of Group "<< tetra_group_id <<" is "<< WER.squaredNorm() << std::endl;
		Eigen::Vector3d EXP1 = Eigen::Vector3d::Zero();
		EXP1 = PrimeVector.block(3 * pi, 0, 3, 1) - particles[pi]->Get_Exp_Pos();
		std::cout << particles[pi]->p_id << " of Group " << tetra_group_id << " is " << EXP1.squaredNorm() << std::endl;
	}
}
void TetraGroupD::Calc_Convergence2() {
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->Is_Fixed()) {
			Eigen::Vector3d WER = Eigen::Vector3d::Zero();
			WER = PrimeVector.block(3 * pi, 0, 3, 1) + Deltax_In_Group.block(3 * pi, 0, 3, 1) - particles[pi]->Get_Exp_Pos();
			std::cout << particles[pi]->p_id << " of Group " << tetra_group_id << " is " << WER.squaredNorm() << std::endl;
		}
		else {
			Eigen::Vector3d WER = Eigen::Vector3d::Zero();
			WER = Deltax_In_Group.block(3 * pi, 0, 3, 1) - particles[pi]->Get_Deltax_In_Model();
			std::cout << particles[pi]->p_id << " of Group " << tetra_group_id << " is " << WER.squaredNorm() << std::endl;
		}
	}
}
//一つ前の予測位置と現在の予測位置でどれだけ反復で更新するのか記憶
//節点で計算するようにしたからこの関数は使わない
double TetraGroupD::Add_convergence_iteration(double convite) {
	ParticleD* pit;
	Eigen::VectorXd ExpVector = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::VectorXd ExpAgoVector = Eigen::VectorXd::Zero(3 * particle_num);

	//expベクトルの生成
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		ExpVector.block(3 * pi, 0, 3, 1) = pit->Get_Exp_Pos();
		ExpAgoVector.block(3 * pi, 0, 3, 1) = pit->Get_ExpAgo_Pos();
	}
	return convite + (ExpVector - ExpAgoVector).squaredNorm();
}
void TetraGroupD::Calc_Jacobi_Matrix_iteration() {
	//初期化
	Jacobi_Matrix = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	Jacobi_Matrix_Sparse.setZero();

	Eigen::MatrixXd Ident = Eigen::MatrixXd::Identity(3 * particle_num, 3 * particle_num);

	//回転行列を3Nにする
	Eigen::MatrixXd rotate_matrix3N = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	// Calc rotate_matrix3N
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		rotate_matrix3N.block(3 * pi, 3 * pi, 3, 3) = rotate_matrix;
	}

	//python用の出力
	//OP_python_M();
	//OP_python_MC();
	//OP_python_R(rotate_matrix3N);
	//OP_python_SumM();
	//OP_python_I(Ident);
	//OP_python_Stiff();


	Jacobi_Matrix = M_Matrix_C + Damping_Matrix + rotate_matrix3N * stiffness_matrix * rotate_matrix3N.transpose()*(Ident - SUM_M_Matrix) * TIME_STEP * TIME_STEP;

	//Jacobi_Matrix = Ident + MassDamInv_Matrix * stiffness_matrix * TIME_STEP * TIME_STEP - MassDamInv_Matrix * stiffness_matrix * TIME_STEP * TIME_STEP * SUM_M_Matrix;

	
	//Sparse化
	Jacobi_Matrix_Sparse = Jacobi_Matrix.sparseView();
	//std::cout << "Jacobi_Matrix_Sparse" << Jacobi_Matrix_Sparse << std::endl;

	//python用の出力
	//OP_python_Jacobi1(rotate_matrix3N);
	//OP_python_Jacobi();

	//固有値の確認
	//OP_eigenvalue(Jacobi_Matrix);
	
	
	//対角優位かチェック
	//OP_diag_advantage();

	//対称行列かチェック!!
	//OP_Symetric(rotate_matrix3N);

	//可換かチェック
	//OP_CommutativeSR(rotate_matrix3N);
	//やはり可換ではない

	/*
	//対角行列の生成
	for (unsigned int pii = 0; pii < 3 * particle_num; pii++) {
		DiagFEM_Matrix_iteration.block(pii,pii , 1, 1) = Jacobi_Matrix.block(pii, pii, 1, 1);
	}
	//対角以外の行列の生成
	F_FEM_Matrix_iteration = Jacobi_Matrix - DiagFEM_Matrix_iteration ;
	//2020/12/14現在使ってはいない
	*/

	//対角優位かチェック
	//OP_diag_advantage2()

	//その他もろもろの出力がみたいとき
	//OP_OtherMatrix(rotate_matrix3N);
}
void TetraGroupD::LHS() {
	Jacobi_Matrix_Sparse.setZero();
	Eigen::SparseMatrix<double> sparse_identity(3 * particle_num, 3 * particle_num);
	sparse_identity.setIdentity();
	

	//Jacobi_Matrix = Damm_Matrix_Sparse + Rn_Matrix_Sparse * StiffnessTT_Matrix_Sparse * Rn_MatrixTR_Sparse * MassCondi_Sparse;

	//Sparse化
	//Jacobi_Matrix_Sparse = Jacobi_Matrix.sparseView();


	//現公
	Jacobi_Matrix_Sparse = Damm_Matrix_Sparse + Rn_Matrix_Sparse * StiffnessTT_Matrix_Sparse * Rn_MatrixTR_Sparse * MassCondi_Sparse;
	//Jacobi_Matrix_Sparse = sparse_identity + (MassDamInv_Matrix).sparseView() * StiffnessTT_Matrix_Sparse - (MassDamInv_Matrix).sparseView() * StiffnessTT_Matrix_Sparse * SUM_M_Matrix.sparseView();
	
	
}

void TetraGroupD::RHS(){
	Constant_term_iteration = Eigen::VectorXd::Zero(3 * particle_num);

	//計算
	//MassCondi = Eigen::MatrixXd::Identity(3 * particles.size(), 3 * particles.size()) - SUM_M_Matrix;//(I-Mj,cm)
	Constant_term_iteration = Rn_Matrix_Sparse * StiffnessTT_Matrix_Sparse * (OrigineVector - Rn_MatrixTR_Sparse * MassCondi_Sparse * PrimeVector) + bind_force_iterative;

	//Constant_term_iteration = StiffnessTT_Matrix_Sparse * (OrigineVector - Rn_MatrixTR_Sparse * MassCondi_Sparse * PrimeVector) + Rn_MatrixTR_Sparse * bind_force_iterative * TIME_STEP * TIME_STEP;

}
void TetraGroupD::LHS0() {
	//初期化
	Jacobi_Matrix = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	Jacobi_Matrix_Inv = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	Jacobi_Matrix_Sparse.setZero();

	Eigen::MatrixXd Ident = Eigen::MatrixXd::Identity(3 * particle_num, 3 * particle_num);

	//回転行列を3Nにする
	Eigen::MatrixXd rotate_matrix3N = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	// Calc rotate_matrix3N
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		rotate_matrix3N.block(3 * pi, 3 * pi, 3, 3) = rotate_matrix;
	}

	//Jacobi_Matrix = M_Matrix_C + Damping_Matrix + rotate_matrix3N * stiffness_matrix * rotate_matrix3N.transpose() * (Ident - SUM_M_Matrix) * TIME_STEP * TIME_STEP;
	Jacobi_Matrix = Ident + TIME_STEP * TIME_STEP * MassDamInv_Matrix * stiffness_matrix * (Ident - SUM_M_Matrix);
	Jacobi_Matrix_Inv = Jacobi_Matrix.inverse();
	
	//Jacobi_Matrix = Ident + MassDamInv_Matrix * stiffness_matrix * (Ident - SUM_M_Matrix) * TIME_STEP * TIME_STEP;

	//Sparse化
	//Jacobi_Matrix_Sparse = Jacobi_Matrix.sparseView();
}

void TetraGroupD::RHS0() {
	Constant_term_iteration = Eigen::VectorXd::Zero(3 * particle_num);
	
	Eigen::MatrixXd Ident2 = Eigen::MatrixXd::Identity(3 * particle_num, 3 * particle_num);

	

	//回転行列を3Nにする
	Eigen::MatrixXd rotate_matrix3N = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	// Calc rotate_matrix3N
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		rotate_matrix3N.block(3 * pi, 3 * pi, 3, 3) = rotate_matrix;
	}

	//計算
	Constant_term_iteration = TIME_STEP * TIME_STEP * MassDamInv_Matrix * (stiffness_matrix * (OrigineVector - rotate_matrix3N.transpose() * PrimeVector + rotate_matrix3N.transpose() * SUM_M_Matrix) + rotate_matrix3N.transpose() * bind_force_iterative);
}

void TetraGroupD::Calc_Jacobi_Matrix_iteration_Sparse() {
	//初期化
	Jacobi_Matrix_Sparse.setZero();
	Eigen::MatrixXd Ident = Eigen::MatrixXd::Identity(3 * particle_num, 3 * particle_num);

	//Jacobi_Matrix = Damm_Matrix_Sparse + Rn_Matrix_Sparse * StiffnessTT_Matrix_Sparse * Rn_MatrixTR_Sparse * MassCondi_Sparse;

	//Sparse化
	//Jacobi_Matrix_Sparse = Jacobi_Matrix.sparseView();


	//現公
	Jacobi_Matrix_Sparse = Damm_Matrix_Sparse + Rn_Matrix_Sparse * StiffnessTT_Matrix_Sparse * Rn_MatrixTR_Sparse * MassCondi_Sparse;

	//Eigen::MatrixXd MCInv = MassDamInv_Matrix.inverse();//(mass+damping)^-1
	//Eigen::SparseMatrix<double> MCInv_Sparse = MCInv.sparseView();//sparse 化
	//Jacobi_Matrix_Sparse = Rn_MatrixTR_Sparse + StiffnessTT_Matrix_Sparse * MCInv_Sparse * Rn_MatrixTR_Sparse * (Ident - MassCondi_Sparse);

	//Updated A 
	//Jacobi_Matrix_Sparse = Ident + Damm_Matrix_Sparse.inverse() * StiffnessTT_Matrix_Sparse - Damm_Matrix_Sparse.inverse() * StiffnessTT_Matrix_Sparse * SUM_M_Matrix
	
	//重心じゃないやつ用
	//Jacobi_Matrix_Sparse = Damm_Matrix_Sparse + Rn_Matrix_Sparse * StiffnessTT_Matrix_Sparse * Rn_MatrixTR_Sparse;
	
	
	//Old用
	//Jacobi_Matrix_Sparse = Damm_Matrix_Sparse + StiffnessTT_Matrix_Sparse * MassCondi_Sparse;
	


	//Jacobi_Matrix_Sparse = Damm_Matrix_Sparse + F_bind_damping*Rn_Matrix_Sparse * StiffnessTT_Matrix_Sparse * Rn_MatrixTR_Sparse * MassCondi_Sparse;
	/*
	Jacobi_Matrix_Sparse = Rn_MatrixTR_Sparse * MassCondi_Sparse;
	Jacobi_Matrix_Sparse = StiffnessTT_Matrix_Sparse * Jacobi_Matrix_Sparse;
	Jacobi_Matrix_Sparse = Rn_Matrix_Sparse * Jacobi_Matrix_Sparse;
	Jacobi_Matrix_Sparse += Damm_Matrix_Sparse;
	*/
	//std::cout << "Jacobi_Matrix_Sparse" << Jacobi_Matrix_Sparse << std::endl;
}
//OldFEMの係数行列作成
void TetraGroupD::Calc_Jacobi_Matrix_iteration_Old() {
	//初期化
	Jacobi_Matrix = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	Jacobi_Matrix_Sparse.setZero();
	if (useSparse) {
		Jacobi_Matrix_Sparse = Damm_Matrix_Sparse + StiffnessTT_Matrix_Sparse;
	}
	else {
		Jacobi_Matrix = M_Matrix_C + Damping_Matrix + TIME_STEP * TIME_STEP * stiffness_matrix;
		Jacobi_Matrix_Sparse = Jacobi_Matrix.sparseView();
	}
	//std::cout << "Jacobi_Matrix_Sparse" << Jacobi_Matrix_Sparse << std::endl;

	//python用の出力
	//OP_python_M();
	//OP_python_MC();
	//OP_python_I(Ident);
	//OP_python_Stiff();


	//Jacobi_Matrix = M_Matrix_C + Damping_Matrix +  stiffness_matrix * TIME_STEP * TIME_STEP;

	//Sparse化
	//Jacobi_Matrix_Sparse = Jacobi_Matrix.sparseView();

	

	//std::cout << "Jacobi_Matrix_Sparse" << Jacobi_Matrix_Sparse << std::endl;

	//python用の出力
	//OP_python_Jacobi1(rotate_matrix3N);
	//OP_python_Jacobi();

	//固有値の確認
	//OP_eigenvalue(Jacobi_Matrix);


	//対角優位かチェック
	//OP_diag_advantage();

	//対称行列かチェック!!
	//OP_Symetric(rotate_matrix3N);

	//可換かチェック
	//OP_CommutativeSR(rotate_matrix3N);
	//やはり可換ではない

	/*
	//対角行列の生成
	for (unsigned int pii = 0; pii < 3 * particle_num; pii++) {
		DiagFEM_Matrix_iteration.block(pii,pii , 1, 1) = Jacobi_Matrix.block(pii, pii, 1, 1);
	}
	//対角以外の行列の生成
	F_FEM_Matrix_iteration = Jacobi_Matrix - DiagFEM_Matrix_iteration ;
	//2020/12/14現在使ってはいない
	*/

	//対角優位かチェック
	//OP_diag_advantage2()

}
void TetraGroupD::Calc_Constant_term_iteration(){
	ParticleD* pit;
	//初期化
	Constant_term_iteration = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::VectorXd ExpVector = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::VectorXd LocalVector = Eigen::VectorXd::Zero(3 * particle_num);

	//expベクトルの生成
	//本来は1ステップごとにやる意味はない
	//事前計算でやるべき
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		pit = particles[pi];
		ExpVector.block(3 * pi, 0, 3, 1) = pit->p_mass * pit->Get_Prime_Pos();
		LocalVector.block(3 * pi, 0, 3, 1) = origin_local_grid[pi];
	}

	//回転行列を3Nにする
	Eigen::MatrixXd rotate_matrix3N = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);

	// Calc rotate_matrix3N
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		rotate_matrix3N.block(3 * pi, 3 * pi, 3, 3) = rotate_matrix;
	}

	//計算
	Constant_term_iteration = ExpVector + rotate_matrix3N * stiffness_matrix * LocalVector * TIME_STEP * TIME_STEP;
	
	//debug
	/*for (unsigned int pi = 0; pi < particle_num; pi++) {
		//節点のidを確認
		std::cout << "particle is" << std::endl;
		std::cout << particles[pi]->p_id << std::endl;
		//値
		std::cout << "Constant value C is" << std::endl;
		std::cout << Constant_term_iteration.block(3 * pi, 0, 3, 1) << std::endl;
	}*/
}
void TetraGroupD::Calc_Constant_term_iteration2() {
	//初期化
	Constant_term_iteration = Eigen::VectorXd::Zero(3 * particle_num);
	//Eigen::VectorXd ExpPreVector = Eigen::VectorXd::Zero(3 * particle_num);
	//Eigen::VectorXd LocalVector = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::MatrixXd Ident2 = Eigen::MatrixXd::Identity(3 * particle_num, 3 * particle_num);

	//expベクトルの生成
	//本来は1ステップごとにやる意味はない
	//事前計算でやるべき
	/*
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		ExpPreVector.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Prime_Pos();
		LocalVector.block(3 * pi, 0, 3, 1) = origin_local_grid[pi];
	}
	*/

	//回転行列を3Nにする
	Eigen::MatrixXd rotate_matrix3N = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	// Calc rotate_matrix3N
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		rotate_matrix3N.block(3 * pi, 3 * pi, 3, 3) = rotate_matrix;
	}

	//計算
	Constant_term_iteration = TIME_STEP * TIME_STEP * rotate_matrix3N * stiffness_matrix * ( (rotate_matrix3N.transpose() * ((SUM_M_Matrix -  Ident2) * PrimeVector)) + OrigineVector);
	//Constant_term_iteration = TIME_STEP * TIME_STEP * stiffness_matrix * ((rotate_matrix3N.transpose() * ((SUM_M_Matrix - Ident2) * PrimeVector)) + OrigineVector);
	//debug
	/*for (unsigned int pi = 0; pi < particle_num; pi++) {
		//節点のidを確認
		std::cout << "particle is" << std::endl; 
		std::cout << particles[pi]->p_id << std::endl;
		//値
		std::cout << "Constant value D is" << std::endl;
		std::cout << Constant_term_iteration.block(3 * pi, 0, 3, 1) << std::endl;
	}*/

	//python用の出力
	//OP_python_ConstantDelta();
}
void TetraGroupD::Calc_Constant_term_iteration_Sparse() {
	//初期化
	Constant_term_iteration = Eigen::VectorXd::Zero(3 * particle_num);

	//計算
	//MassCondi = Eigen::MatrixXd::Identity(3 * particles.size(), 3 * particles.size()) - SUM_M_Matrix;//(I-Mj,cm)
	Constant_term_iteration = Rn_Matrix_Sparse * StiffnessTT_Matrix_Sparse * (OrigineVector - Rn_MatrixTR_Sparse * MassCondi_Sparse * PrimeVector);
	
	//Eigen::MatrixXd MCInv = MassDamInv_Matrix.inverse();//(Mj+Cj')^-1
	//Eigen::SparseMatrix<double> MCInv_Sparse = MCInv.sparseView();//sparse 化
	//Constant_term_iteration = StiffnessTT_Matrix_Sparse * (OrigineVector - Rn_MatrixTR_Sparse * MassCondi_Sparse * PrimeVector);

	/*Constant_term_iteration = MassCondi_Sparse * -1 * PrimeVector;
	Constant_term_iteration = Rn_MatrixTR_Sparse * Constant_term_iteration;
	Constant_term_iteration += OrigineVector;
	Constant_term_iteration = Rn_Matrix_Sparse * StiffnessTT_Matrix_Sparse * Constant_term_iteration;*/
	

	//重心じゃないやつ用
	/*
	Constant_term_iteration = -1 * PrimeVector;
	Constant_term_iteration = Rn_MatrixTR_Sparse * Constant_term_iteration;
	Constant_term_iteration += InitialVector;
	Constant_term_iteration = Rn_Matrix_Sparse * StiffnessTT_Matrix_Sparse * Constant_term_iteration;
	*/

	//Old用
	/*
	Constant_term_iteration = MassCondi_Sparse * -1 * PrimeVector;
	Constant_term_iteration = Constant_term_iteration;
	Constant_term_iteration += OrigineVector;
	Constant_term_iteration = StiffnessTT_Matrix_Sparse * Constant_term_iteration;
	*/

	//Constant_term_iteration = Rn_Matrix_Sparse * Constant_term_iteration;
	//std::cout << Constant_term_iteration << std::endl;
}
//OldFEMの定数値作成
void TetraGroupD::Calc_Constant_term_iteration_Old() {
	//初期化
	Constant_term_iteration = Eigen::VectorXd::Zero(3 * particle_num);

	//計算
	//Constant_term_iteration = TIME_STEP * TIME_STEP *  stiffness_matrix * (InitialVector - PrimeVector);
	if (useSparse) {
		Constant_term_iteration = StiffnessTT_Matrix_Sparse * (InitialVector - PrimeVector);
	}
	else {
		Constant_term_iteration = TIME_STEP * TIME_STEP * stiffness_matrix * (InitialVector - PrimeVector);
	}
	

	//debug
	/*for (unsigned int pi = 0; pi < particle_num; pi++) {
		//節点のidを確認
		std::cout << "particle is" << std::endl;
		std::cout << particles[pi]->p_id << std::endl;
		//値
		std::cout << "Constant value D is" << std::endl;
		std::cout << Constant_term_iteration.block(3 * pi, 0, 3, 1) << std::endl;
	}*/

	//python用の出力
	//OP_python_ConstantDelta();
}
void TetraGroupD::OP_python_M() {
	
	std::cout << "M_Matrix is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
	//値
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
	//値
	std::cout << M_Matrix.block(pi, pj, 1, 1) << ",";
	}
	std::cout << M_Matrix.block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
	//値
	std::cout << M_Matrix.block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << M_Matrix.block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
	
	//python用の出力
	/*
	std::cout << "M_Matrix_emergency is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
	//値
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
	//値
	std::cout << M_Matrix_emergency.block(pi, pj, 1, 1) << ",";
	}
	std::cout << M_Matrix_emergency.block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
	//値
	std::cout << M_Matrix_emergency.block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << M_Matrix_emergency.block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
	*/
}
//python用の出力
//モデル全体を考えた質量行列を出力
void TetraGroupD::OP_python_MC() {
	std::cout << "M_Matrix_C is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
		//値
		std::cout << "[";
		for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
			//値
			std::cout << M_Matrix_C.block(pi, pj, 1, 1) << ",";
		}
		std::cout << M_Matrix_C.block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
		//値
		std::cout << M_Matrix_C.block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << M_Matrix_C.block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
}
void TetraGroupD::OP_python_R(Eigen::MatrixXd rotate_matrix3N) {
	
	std::cout << "rotate_matrix3N is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
	//値
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
	//値
	std::cout << rotate_matrix3N.block(pi, pj, 1, 1) << ",";
	}
	std::cout << rotate_matrix3N.block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
	//値
	std::cout << rotate_matrix3N.block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << rotate_matrix3N.transpose().block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
	//python用の出力
	std::cout << "rotate_matrix3N.transpose() is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
	//値
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
	//値
	std::cout << rotate_matrix3N.transpose().block(pi, pj, 1, 1) << ",";
	}
	std::cout << rotate_matrix3N.transpose().block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
	//値
	std::cout << rotate_matrix3N.transpose().block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << rotate_matrix3N.transpose().block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
}
void TetraGroupD::OP_python_SumM() {
	std::cout << "SUM_M_Matrix / Group_Mass is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
		//値
		std::cout << "[";
		for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
			//値
			std::cout << (SUM_M_Matrix / Group_Mass).block(pi, pj, 1, 1) << ",";
		}
		std::cout << (SUM_M_Matrix / Group_Mass).block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
		//値
		std::cout << (SUM_M_Matrix / Group_Mass).block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << (SUM_M_Matrix / Group_Mass).block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
}
void TetraGroupD::OP_python_I(Eigen::MatrixXd Ident) {
	std::cout << "Ident is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
		//値
		std::cout << "[";
		for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
			//値
			std::cout << Ident.block(pi, pj, 1, 1) << ",";
		}
		std::cout << Ident.block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
		//値
		std::cout << Ident.block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << Ident.block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
}
void TetraGroupD::OP_python_Stiff() {
	std::cout << "stiffness_matrix is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
		//値
		std::cout << "[";
		for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
			//値
			std::cout << stiffness_matrix.block(pi, pj, 1, 1) << ",";
		}
		std::cout << stiffness_matrix.block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
		//値
		std::cout << stiffness_matrix.block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << stiffness_matrix.block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
}
//モデル全体を考えた質量行列を出力
void TetraGroupD::OP_python_Damping() {
	std::cout << "Damping_Matrix is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
		//値
		std::cout << "[";
		for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
			//値
			std::cout << Damping_Matrix.block(pi, pj, 1, 1) << ",";
		}
		std::cout << Damping_Matrix.block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
		//値
		std::cout << Damping_Matrix.block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << Damping_Matrix.block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
}
void TetraGroupD::OP_python_Jacobi1(Eigen::MatrixXd rotate_matrix3N, Eigen::MatrixXd Ident) {
	std::cout << "Jacobi Debug is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
		//値
		std::cout << "[";
		for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
			//値
			std::cout << (stiffness_matrix * rotate_matrix3N.transpose()*(Ident - (SUM_M_Matrix / Group_Mass))).block(pi, pj, 1, 1) << ",";
		}
		std::cout << (stiffness_matrix * rotate_matrix3N.transpose()*(Ident - (SUM_M_Matrix / Group_Mass))).block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
		//値
		std::cout << (stiffness_matrix * rotate_matrix3N.transpose()*(Ident - (SUM_M_Matrix / Group_Mass))).block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << (stiffness_matrix * rotate_matrix3N.transpose()*(Ident - (SUM_M_Matrix / Group_Mass))).block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
}
void TetraGroupD::OP_python_Jacobi() {
	std::cout << "Jacobi_Matrix is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
		//値
		std::cout << "[";
		for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
			//値
			std::cout << Jacobi_Matrix.block(pi, pj, 1, 1) << ",";
		}
		std::cout << Jacobi_Matrix.block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
		//値
		std::cout << Jacobi_Matrix.block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << Jacobi_Matrix.block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
}
//固有値の確認
void TetraGroupD::OP_eigenvalue(Eigen::MatrixXd Matrix) {
	Eigen::EigenSolver< Eigen::MatrixXd > s(Matrix);
	//固有値出力
	std::cout << "固有値\n" << s.eigenvalues() << std::endl;
	std::cout << "固有ベクトル\n" << s.eigenvectors() << std::endl;
}
void TetraGroupD::OP_diag_advantage() {
	std::cout << "diag check" << std::endl;
	for (unsigned int pii = 0; pii < 3 * particle_num; pii++) {
	double cccout = 0;
	for (unsigned int pij = 0; pij < 3 * particle_num; pij++) {
		cccout = cccout + abs(Jacobi_Matrix(pii, pij));
	}
	cccout = cccout - abs(Jacobi_Matrix(pii, pii));
	if (abs(Jacobi_Matrix(pii, pii)) > cccout) {
		std::cout << "diagonal" << std::endl;
	}
	else {
		std::cout << "nooooooooooooooooooooooooooooooooooooooooooooooo" << std::endl;
	}
		cccout = 0;
	}
}
void TetraGroupD::OP_diag_advantage2() {
	//対角優位かチェック2
	for (unsigned int pii = 0; pii < 3 * particle_num; pii++) {
		double cccout = 0;
		for (unsigned int pij = 0; pij < 3 * particle_num; pij++) {
			cccout = cccout + abs(F_FEM_Matrix_iteration(pii, pij));
		}
		if (abs(DiagFEM_Matrix_iteration(pii, pii))> cccout) {
			std::cout << "diagonal advantage2" << std::endl;
		}
		else {
			std::cout << "no2" << std::endl;
		}
		cccout = 0;
	}
}
void TetraGroupD::OP_Symetric(Eigen::MatrixXd rotate_matrix3N) {
	std::cout << "M" << std::endl;
	std::cout << M_Matrix << std::endl;
	std::cout << "M.T" << std::endl;
	std::cout << M_Matrix.transpose() << std::endl;

	std::cout << "stiffness_matrix" << std::endl;
	std::cout << std::setprecision(16) << stiffness_matrix << std::endl;
	std::cout << "ONE diff" << std::endl;
	for (unsigned int qii = 0; qii< 3 * particle_num; qii++) {
	for (unsigned int qij = 0; qij < 3 * particle_num; qij++) {
	std::cout << (stiffness_matrix).block(qii, qij, 1, 1) - ((stiffness_matrix.transpose()).block(qii, qij, 1, 1)) << ",";
	}
	std::cout << "F" << std::endl;
	}

	//ヤコビ行列が対称行列かチェック
	
	std::cout << "Jacobi_Matrix" << std::endl;
	for (unsigned int qii = 0; qii< 3 * particle_num; qii++) {
	for (unsigned int qij = 0; qij < 3 * particle_num; qij++) {
	std::cout<<(Jacobi_Matrix).block(qii, qij, 1, 1) - ((Jacobi_Matrix.transpose()).block(qii, qij, 1, 1)) << ",";
	}
	std::cout << "F" << std::endl;
	}

	std::cout << "ONE Mat" << std::endl;
	std::cout << std::setprecision(6) <<rotate_matrix3N * stiffness_matrix * rotate_matrix3N.transpose() * TIME_STEP * TIME_STEP << std::endl;
	std::cout << "ONE.T" << std::endl;
	std::cout << std::setprecision(6) <<(rotate_matrix3N * stiffness_matrix * rotate_matrix3N.transpose() * TIME_STEP * TIME_STEP).transpose() << std::endl;
	std::cout << "ONE diff" << std::endl;
	for (unsigned int qii = 0;qii< 3 * particle_num; qii++) {
	for (unsigned int qij = 0; qij < 3 * particle_num; qij++) {
	std::cout << (rotate_matrix3N * stiffness_matrix * rotate_matrix3N.transpose() * TIME_STEP * TIME_STEP).block(qii,qij , 1, 1) - ((rotate_matrix3N * stiffness_matrix * rotate_matrix3N.transpose() * TIME_STEP * TIME_STEP).transpose()).block(qii, qij, 1, 1)<<",";
	}
	std::cout << "F" << std::endl;
	}
	std::cout << "stiffness_matrix" << std::endl;
	std::cout << std::setprecision(16) << stiffness_matrix << std::endl;
	std::cout << "stiffness_matrix_TT" << std::endl;
	std::cout << std::setprecision(16) << stiffness_matrix - stiffness_matrix.transpose() << std::endl;
}
void TetraGroupD::OP_CommutativeSR(Eigen::MatrixXd Matrix) {
	std::cout << "stiffness_matrix * rotatetion" << std::endl;
	std::cout << "stiffness_matrix_TT" << std::endl;
	std::cout << std::setprecision(16) << stiffness_matrix * Matrix.transpose() - Matrix.transpose() * stiffness_matrix << std::endl;
}
void TetraGroupD::OP_OtherMatrix(Eigen::MatrixXd rotate_matrix3N) {
	std::cout << "rotate_matrix3N" << std::endl;
	std::cout << std::setprecision(10) << rotate_matrix3N << std::endl;
	std::cout << "SUM_M_Matrix / mass" << std::endl;
	std::cout << std::setprecision(10) << SUM_M_Matrix / mass << std::endl;
	std::cout << "DiagFEM_Matrix_iteration" << std::endl;
	std::cout << DiagFEM_Matrix_iteration << std::endl;
	std::cout << "F_FEM_Matrix_iteration" << std::endl;
	std::cout << F_FEM_Matrix_iteration << std::endl;
}
void TetraGroupD::OP_python_ConstantDelta() {
	//x0(差分じゃないときとの比較用)
	/*
	std::cout << "primary x is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
	//値
	std::cout << ExpPreVector.block(pi, 0, 1, 1) << ",";
	}
	std::cout << ExpPreVector.block(3 * particle_num - 1, 0, 1, 1) << "])" << std::endl;
	*/
	std::cout << "Constant value D is" << tetra_group_id << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
		//値
		std::cout << Constant_term_iteration.block(pi, 0, 1, 1) << ",";
	}
	std::cout << Constant_term_iteration.block(3 * particle_num - 1, 0, 1, 1) << "])" << std::endl;
}
//線形方程式の解いた値などをみる
void TetraGroupD::OP_debug_iterative1(Eigen::VectorXd v) {
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		std::cout << "Constant D of  u" << "particle " << particles[pi]->p_id << "group is " << tetra_group_id << std::endl;
		std::cout << x_In_Group.block(3 * pi, 0, 3, 1) << std::endl;
		std::cout << v.block(3 * pi, 0, 3, 1) << std::endl;
		if (particles[pi]->p_id == 3) {
			std::cout << "Hello u" << "particle " << particles[pi]->p_id << "group is " << tetra_group_id << std::endl;
			std::cout << bind_force_iterative.block(3 * pi, 0, 3, 1) << std::endl;
			std::cout << v.block(3 * pi, 0, 3, 1) << std::endl;
			std::cout << x_In_Group.block(3 * pi, 0, 3, 1) << std::endl;
		}
	}
	//EXPはグループで同じ値
	//X_in_groupはグループで違う値のはず
	std::cout << "Before shared Groval " << "groupnumbe is " << tetra_group_id  <<std::endl;
	std::cout << Get_Exp_Pos(3) << std::setprecision(6) << std::endl;
	std::cout << "Before shared Local " << "groupnumbe is " << tetra_group_id << std::endl;
	std::cout << Get_X_In_Group(3) << std::setprecision(6) << std::endl;
}
void TetraGroupD::OP_python_File() {
	std::ofstream outputfile("MinvK.txt");
	outputfile << "test";
	outputfile.close();
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * particle_num - 1; pi++) {
		//値
		std::cout << "[";
		for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
			//値
			std::cout << stiffness_matrix.block(pi, pj, 1, 1) << ",";
		}
		std::cout << stiffness_matrix.block(pi, 3 * particle_num - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * particle_num - 1; pj++) {
		//値
		std::cout << stiffness_matrix.block(3 * particle_num - 1, pj, 1, 1) << ",";
	}
	std::cout << stiffness_matrix.block(3 * particle_num - 1, 3 * particle_num - 1, 1, 1) << "]])" << std::endl;
}
void TetraGroupD::Set_Deltax_In_Group(Eigen::VectorXd a) {
	Deltax_In_Group = Eigen::VectorXd::Zero(3 * particles.size());
	Deltax_In_Group = a;
}
/*
void TetraGroupD::Mode_Analytic() {
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	double omega;
	double omega2;
	double ModeT;
	A = M_Matrix_C.inverse() * stiffness_matrix;
	Eigen::EigenSolver< Eigen::MatrixXd > s(A);
	Eigen::VectorXcd ve = s.eigenvectors().col(0);
	//Eigen::VectorXd vee = ve.cast(double);
	std::cout << s.eigenvectors() << std::endl;
	std::cout << ve << std::endl;
	ModeT = (2 * PI) / omega;
	std::cout << "Mode Analytic " << tetra_group_id << ":" << omega << std::endl;
}
*/