//------------------------------------------
//@author KatsuyaKikuchi
//@brief Implementation of the TetraGroupD.h
//------------------------------------------

#include "TetraGroupD.h"
#include <typeinfo>
#include "UTClapack.h"
#include <iostream>
#include <fstream>



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
	int bbb = 1;
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
	SUM_M.setZero();
	SUM_M = SUM_M_Matrix.sparseView();
	//std::cout << "SUM_M_Matrix" << std::endl;
	//std::cout << SUM_M_Matrix << std::endl;
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
	MassDamInvSparse.setZero();
	MassCondi_Sparse = MassCondi.sparseView();
	Damm_Matrix_Sparse = Damm_Matrix.sparseView();
	MassDamInvSparse = MassDamInv_Matrix.sparseView();//Mass + damping inverse sparse
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
	StiffnessSparse.setZero();
	StiffnessSparse = stiffness_matrix.sparseView();

	/*std::cout << "Create Object K Matrix Of Group" << tetra_group_id << std::endl;
	std::cout << StiffnessTT_Matrix_Sparse << std::endl;*/
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

	rotate_matrix3N = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	Eigen::MatrixXd rotate_matrix3NTR = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	//Rn_Matrix_Sparse.setZero();
	//Rn_MatrixTR_Sparse.setZero();
	// Calc rotate_matrix3N
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		rotate_matrix3N.block(3 * pi, 3 * pi, 3, 3) = rotate_matrix;
	}
	rotate_matrix3NTR = rotate_matrix3N.transpose();
	Rn_Matrix_Sparse = rotate_matrix3N.sparseView();
	Rn_MatrixTR_Sparse = (rotate_matrix3N.transpose()).sparseView();
	mtRotationSparse.endMyTimer();
	res = rotate_matrix3N * rotate_matrix3NTR;


	//std::ofstream file("rotate.txt");
	//if (file.is_open()) {
	//	// 输出矩阵的每个元素到文件
	//	for (int i = 0; i < rotate_matrix3N.rows(); ++i) {
	//		for (int j = 0; j < rotate_matrix3N.cols(); ++j) {
	//			file << rotate_matrix3N(i, j) << " ";
	//		}
	//		file << std::endl;
	//	}
	//	// 关闭文件流
	//	file.close();
	//}
	//std::ofstream file("res.txt");

	//if (file.is_open()) {
	//	// 输出矩阵的每个元素到文件
	//	for (int i = 0; i < res.rows(); ++i) {
	//		for (int j = 0; j < res.cols(); ++j) {
	//			file << res(i, j) << " ";
	//		}
	//		file << std::endl;
	//	}
	//	// 关闭文件流
	//	file.close();
	//}




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

Eigen::Vector3d TetraGroupD::axlAPD(Eigen::Matrix3d a) {
	Eigen::Vector3d g = Eigen::Vector3d::Zero();
	g[0] = a(1, 2) - a(2, 1);
	g[1] = a(2, 0) - a(0, 2);
	g[2] = a(0, 1) - a(1, 0);
	return g;
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
Eigen::Vector3d TetraGroupD::Get_Velocity(int pid) {
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		if (particles[pi]->p_id == pid) {
			temp = GroupVelVector.block(3 * pi, 0, 3, 1);
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

void TetraGroupD::Set_Group_Mass(double a) {
	this->Group_Mass = a;
}

void TetraGroupD::ReSet_Fbind_Pos() {
	//std::cout << bind_force_iterative << std::endl;
	bind_force_iterative = Eigen::VectorXd::Zero(3* particles.size());
}

//グループのparticleを取得
std::vector<ParticleD*> TetraGroupD::Get_Particle()const {
	return particles;
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

Eigen::Quaterniond TetraGroupD::Exp2(Eigen::Vector3d a) {
	double s = sin((a * 0.5).norm());
	double x = s * a.x() / a.norm();
	double y = s * a.y() / a.norm();
	double z = s * a.z() / a.norm();
	Eigen::Quaterniond qq = Eigen::Quaterniond(cos((a * 0.5).norm()), x, y, z);
	return  qq;
}


//GMRESの前処理を行う

void TetraGroupD::CalcDeltax()
{

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

	solver.compute(Jacobi_Matrix_Sparse);
	DeltaxNew = solver.solve(Constant_term_iteration);
	DeltaxNew = Rn_MatrixTR_Sparse * DeltaxNew;

}


void TetraGroupD::Update_Fbind_Pos8() {
	for (unsigned int pi = 0; pi < particle_num; pi++) {
		Eigen::Vector3d Conv = Eigen::Vector3d::Zero();
		Eigen::Vector3d Conv1 = Eigen::Vector3d::Zero();
		//共有節点かどうか
		//本当は分離しないといけない
		if ((particles[pi]->p_belong_TetraGroup_ids.size()) > 1)  // 如果是共同点（点所属组数大于1
		{
			//固定されていない点
			if (!(particles[pi]->Is_Fixed())) {
				//(n)(Exp + Deltax)
				Conv = (PrimeVector.block(3 * pi, 0, 3, 1) + DeltaxNew.block(3 * pi, 0, 3, 1));
				Conv1 = GroupVelVector.block(3 * pi, 0, 3, 1);
				//std::cout << DeltaxNew.block(3 * pi, 0, 3, 1) << std::endl;
				//std::cout << particles[pi]->p_belong_TetraGroup_ids.size() << std::endl;
				Conv = Conv - (particles[pi]->Get_Exp_Pos() + particles[pi]->Get_Deltax_In_Model());
				Conv1 = Conv1 - particles[pi]->Get_Mean_Vel();
				//Conv = Conv -  (particles[pi]->Get_Exp_Pos() + particles[pi]->Get_Deltax_In_Model());
				bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_coeff * Conv;
				bind_force_iterative.block(3 * pi, 0, 3, 1) += F_bind_damping * Conv1;

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
			Conv = (PrimeVector.block(3 * pi, 0, 3, 1) + DeltaxNew.block(3 * pi, 0, 3, 1));
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
		//if (((particles[pi]->p_belong_TetraGroup_ids.size()) > 1)) {
		//	if (Conv.squaredNorm() < 10e-3) {
		//		std::cout << "Converce " << particles[pi]->p_id << std::endl;
		//	}
		//	else
		//		std::cout << "NotConvergence " << std::endl;
		//}
	}
}
//Debug用
//節点は別処理


void TetraGroupD::LHS0() {
	//初期化
	Jacobi_Matrix = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	Jacobi_Matrix_Inv = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	Jacobi_Matrix_Sparse.setZero();

	Eigen::MatrixXd Ident = Eigen::MatrixXd::Identity(3 * particle_num, 3 * particle_num);

	//回転行列を3Nにする

	//Jacobi_Matrix = M_Matrix_C + Damping_Matrix + rotate_matrix3N * stiffness_matrix * rotate_matrix3N.transpose() * (Ident - SUM_M_Matrix) * TIME_STEP * TIME_STEP;
	//Jacobi_Matrix = Ident + TIME_STEP * TIME_STEP * MassDamInv_Matrix * stiffness_matrix - TIME_STEP * TIME_STEP * MassDamInv_Matrix * stiffness_matrix * SUM_M_Matrix;
	Jacobi_Matrix = Ident + TIME_STEP * TIME_STEP * MassDamInvSparse * StiffnessSparse * MassCondi_Sparse;
	Jacobi_Matrix_Sparse = Jacobi_Matrix.sparseView();
	/*Jacobi_Matrix = Ident + TIME_STEP * TIME_STEP * MassDamInv_Matrix * stiffness_matrix * (Ident - SUM_M_Matrix);
	Jacobi_Matrix_Inv = Jacobi_Matrix.inverse();*/
	//std::ofstream file("Jacobi_Matrix.txt");

	//if (file.is_open()) {
	//	// 输出矩阵的每个元素到文件
	//	for (int i = 0; i < Jacobi_Matrix.rows(); ++i) {
	//		for (int j = 0; j < Jacobi_Matrix.cols(); ++j) {
	//			file << Jacobi_Matrix(i, j) << " ";
	//		}
	//		file << std::endl;
	//	}
	//	// 关闭文件流
	//	file.close();
	//}

	
	//Jacobi_Matrix = Ident + MassDamInv_Matrix * stiffness_matrix * (Ident - SUM_M_Matrix) * TIME_STEP * TIME_STEP;

	//Sparse化
	//Jacobi_Matrix_Sparse = Jacobi_Matrix.sparseView();
}

void TetraGroupD::RHS0() {
	Constant_term_iteration = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::MatrixXd Ident = Eigen::MatrixXd::Identity(3 * particle_num, 3 * particle_num);
	Eigen::MatrixXd tempA = Eigen::MatrixXd::Zero(3 * particle_num, 3 * particle_num);
	Eigen::VectorXd tempB = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::VectorXd tempC = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::VectorXd tempE = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::VectorXd tempF = Eigen::VectorXd::Zero(3 * particle_num);
	Eigen::SparseMatrix<double> tempD;
	tempD.setZero();

	//tempA = TIME_STEP * TIME_STEP * MassDamInv_Matrix * stiffness_matrix;
	//tempD = TIME_STEP * TIME_STEP * MassDamInvSparse * StiffnessSparse;
	//tempB = (OrigineVector - rotate_matrix3N.transpose() * (PrimeVector - SUM_M_Matrix * PrimeVector));
	//tempE = (OrigineVector - Rn_MatrixTR_Sparse * (PrimeVector - SUM_M * PrimeVector));
	//tempC = MassDamInv_Matrix * rotate_matrix3N.inverse() * bind_force_iterative;
	//tempF = MassDamInvSparse * Rn_MatrixTR_Sparse * bind_force_iterative;
	//Constant_term_iteration = tempD * tempE + tempF;
	Constant_term_iteration = (TIME_STEP * TIME_STEP * MassDamInvSparse * StiffnessSparse) * (OrigineVector - Rn_MatrixTR_Sparse * (PrimeVector - SUM_M * PrimeVector)) + MassDamInvSparse * Rn_MatrixTR_Sparse * 0.01 * 0.01* bind_force_iterative;

	//Constant_term_iteration = (TIME_STEP * TIME_STEP * MassDamInv_Matrix * stiffness_matrix) *
		//(OrigineVector - rotate_matrix3N.transpose() * (PrimeVector - SUM_M_Matrix * PrimeVector)) + MassDamInv_Matrix
		//* rotate_matrix3N.inverse() * bind_force_iterative;

	//計算
	/*Constant_term_iteration = (TIME_STEP * TIME_STEP * MassDamInvSparse * StiffnessSparse) * 
		(OrigineVector - Rn_MatrixTR_Sparse * (PrimeVector * (SUM_M - Ident))) + MassDamInvSparse
		* Rn_MatrixTR_Sparse * bind_force_iterative;*/


}

void TetraGroupD::Set_Deltax_In_Group(Eigen::VectorXd a) {
	Deltax_In_Group = Eigen::VectorXd::Zero(3 * particles.size());
	Deltax_In_Group = a;
}
