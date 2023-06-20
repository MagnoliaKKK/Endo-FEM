//===========================================================================//
//@author IKUTO ENDO
//@brief UseBlockObjectDouble.hの実装
//===========================================================================//
#include "UseBlockObjectDouble.h"
//オブジェクトをN(x),M(y),L(z)で分割する
//===========================================================================//
UseBlockObjectDouble::UseBlockObjectDouble(std::vector<ParticleD*> p, ObjectData data) //コンストラクタ
	: ObjectD(p, data)//親クラスはObjectクラス
{
	std::cout << "Create Object Using Local Stiffness Matrix with Group" << std::endl;
	Init();			//初期設定
}
UseBlockObjectDouble::~UseBlockObjectDouble() { //デコンストラクタ
	for (auto _g : groups) {
		delete _g;
	}
	for (auto _t : tetras) {
		delete _t;
	}
	groups.clear();
	tetras.clear();
}
//===========================================================================//
//==========================================================================//
//	@start		   				初期設定									//
//==========================================================================//
void UseBlockObjectDouble::Init() {
	Delaunay_Triangulation();  //ドロネー三角形分割をする
	Create_Groups();		   //グループを作る
	data_name = "Use_GroupBlock";//オブジェクトの名前をUse_GroupBlockとする
}
//==========================================================================//
//	@start		   			  グループ作成								    //
//==========================================================================//
void UseBlockObjectDouble::Create_Groups() {
	//----------------------------------------------------------------------//
	//-------------------------TetraSet 初期設定----------------------------//

	int tetra_group_id = 0;

	//四面体要素の共通面を探索する
	for (auto _t : tetras) {
		_t->Create_Faces();
		tetra_map.insert(std::map<TetraElementD*, bool>::value_type(_t, false));
	}
	for (unsigned int t1 = 0; t1 < tetras.size() - 1; ++t1) {
		for (unsigned int t2 = t1 + 1; t2 < tetras.size(); ++t2) {
			if (tetras[t1]->has_Common_Face(tetras[t2])) {
				tetras[t1]->Push_Side_Tetra(tetras[t2]);
				tetras[t2]->Push_Side_Tetra(tetras[t1]);
			}
		}
	}
	//----------------------------------------------------------------------//
	//-------------------------N個のGroupの作成----------------------------//
	std::vector< TetraElementD* > remain = this->Create_Group_Candidate();
	//----------------------------------------------------------------------//
	//N-1個の作成
	std::vector<TetraElementD*> temp_t;
	for (int ix = 0; ix < dividedbyn; ix++) {
		for (int iy = 0; iy < dividedbym; iy++) {
			for (int iz = 0; iz < dividedbyl; iz++) {
				for (auto _t : remain) {
					if (GroupdivideNML(_t, ix,iy,iz)) {
						temp_t.push_back(_t);
					}
				}
				Create_Group(temp_t, tetra_group_id);
				temp_t.clear();
				tetra_group_id++;
				remain = this->Create_Group_Candidate();
			}
		}
	}
	//N個目の作成
	Create_Group(remain, tetra_group_id);
	tetra_group_id++;
	remain = this->Create_Group_Candidate();
	if (remain.empty()) {
		std::cout << "All complete : " << std::endl;
	}
	//----------------------------------------------------------------------//
	std::cout << "success create All-Group. All-Group size : " << groups.size() << std::endl;
	//グループ分け完了

	//共有節点行列の作成
	if (groups.size() > 1) {
		//Node_Sharing_Matrixの作成
		Sum_particlenum = 0;
		for (auto _g : groups) {
			Sum_particlenum += _g->particle_num;
		}
		std::cout << "Sum_particlenum is " << Sum_particlenum << std::endl;
		Node_Sharing_Matrix = Eigen::MatrixXd::Zero(3 * Sum_particlenum, 3 * Sum_particlenum);
		int tempp1 = 0;
		int tempp2 = 0;
		int tempp3 = 0;
		int tempp4 = 0;
		int tempp5 = 0;
		for (unsigned int i = 0; i < groups.size(); i++) {
			for (unsigned int j = 0; j < groups[i]->particle_num; j++) {
				tempp1 = groups[i]->particles[j]->p_id;
				tempp5 = 0;
				tempp5 += groups[i]->particle_num;;
				for (unsigned int k = i + 1; k < groups.size(); k++) {
					for (unsigned int l = 0; l < groups[k]->particle_num; l++) {
						tempp2 = groups[k]->particles[l]->p_id;
						if (tempp1 == tempp2) {
							Node_Sharing_Matrix.block(3 * tempp3, 3 * tempp4 + 3 * j, 3, 3) = Eigen::Matrix3d::Identity();
							Node_Sharing_Matrix.block(3 * tempp3, 3 * tempp4 + 3 * tempp5 + 3 * l, 3, 3) = -1 * Eigen::Matrix3d::Identity();
							tempp3++;
							Share_particle_id.push_back(tempp2);
						}
					}
					tempp5 += groups[k]->particle_num;
				}
			}
			tempp4 += groups[i]->particle_num;
		}
		Share_particlenum = tempp3;
		//これでたぶん 3*tempp3 x 3* Sum_particlenum の行列ができるはず
		//std::cout << "オブジェクトの共有節点行列 is " << std::endl;
		//std::cout << Node_Sharing_Matrix.block(0, 0, 3 * Share_particlenum, 3 * Sum_particlenum) << std::endl;
		
		std::cout << "オブジェクトの共有数 is " << Share_particlenum << std::endl;
		std::cout << "オブジェクトの共有節点 is " << std::endl;
		for (int i = 0; i < Share_particlenum; i++) {
			std::cout << Share_particle_id[i] << std::endl;
		}
	}
	for (auto _g : groups) {
		AllParticlenum += _g->particle_num;
	}
	std::cout << "オブジェクトの全数 is " << AllParticlenum << std::endl;

	double tempM = 0.0;
	//モデルの体積計算
	double tempV = 0.0;
	for (auto _g : groups) {
		tempV += _g->Get_Volume();
	}
	//質量の決定
	if (mSys == 0) {
		//グループの質量の決定
		for (auto _g : groups) {
			_g->Set_Group_Mass(_g->Get_Volume() * this->data.density);
		}
		double Pmass = (tempV * this->data.density) / particles.size();
		for (auto _p : particles) {
			_p->Set_Mass(Pmass);
		}
		for (auto _p : particles) {
			std::cout << _p->Get_Mass() << std::endl;
		}
		//グループの対角質量行列の作成
		for (auto _g : groups) {
			for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
				_g->M_Matrix_C.block(3 * pi, 3 * pi, 3, 3) = Pmass * Eigen::Matrix3d::Identity(3, 3);
			}
			_g->Create_SUM_M_Matrix();
		}
	}
	else if(mSys == 1) {
		//質量で和をとる
		double Massp = 0.0;
		for (auto _p : particles) {
			Massp = 0.0;
			if (_p->p_belong_TetraGroup_ids.size() == 0) {
				std::cout << "ERROR_UseGrouptetraObject124" << std::endl;
			}
			for (auto _g : _p->p_belong_TetraGroup_ids) {
				TetraGroupD* tg = groups[_g];
				Massp += tg->Get_GMass_In_Group(_p->p_id);
				//std::cout << "tg" << tg->tetra_group_id << std::endl;
				//std::cout << "p_id" << p->p_id << std::endl;
				//std::cout << "delta_p" << delta_p << std::endl;
			}
			_p->Set_Mass(Massp);
		}
		for (auto _p : particles) {
			std::cout << _p->Get_Mass() << std::endl;
		}
		for (auto _g : groups) {
			tempM = 0.0;
			//モデル全体の質量から計算したM行列の作成
			for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
				_g->M_Matrix_C.block(3 * pi, 3 * pi, 3, 3) = _g->particles[pi]->Get_Mass() * Eigen::Matrix3d::Identity(3, 3);
				tempM += _g->particles[pi]->Get_Mass();
			}
			_g->Set_Group_Mass(tempM);
			_g->Create_SUM_M_Matrix();
		}
	}
	//質量行列はFEMのものを使う
	else {
		for (auto _g : groups) {
			//質量行列の代入
			_g->M_Matrix_C = Eigen::MatrixXd::Zero(3 * _g->particles.size(), 3 * _g->particles.size());
			_g->M_Matrix_C = _g->M_Matrix;
			std::ofstream outputfile("Mass_Matrix.txt", std::ios_base::app);
			outputfile << "Mass_Matrix "<< _g->tetra_group_id<<" is " << std::endl;
			outputfile << std::setprecision(3) << _g->M_Matrix_C << std::endl;
			outputfile.close();
			
			//グループの質量
			double tempM = 0.0;
			for (unsigned int pi = 0; pi < _g->particle_num;pi++) {
				tempM += _g->Get_GMass_In_Group(_g->particles[pi]->p_id);
			}
			//std::cout << "Group" << tetra_group_id << " Mass is " << tempM << std::endl;
			//std::cout << "Group" << tetra_group_id << " Mass is " << _g->mass << std::endl;
			_g->Set_Group_Mass(tempM);
			_g->Create_SUM_M_Matrix();
		}
	}
	//剛性行列の作成
	MicroSecondTimer mtStiffness;
	mtStiffness.setid(33);
	double StiffnessTime = 0.0;
	for (auto _g : groups) {
		mtStiffness.startMyTimer();
		_g->Create_Center_Grid();
		_g->Create_Local_Stiffness_Matrix();
		_g->Create_Damping_Matrix();
		mtStiffness.endMyTimer();
		StiffnessTime += mtStiffness.getDt();

		//節点情報などの計算
		_g->Create_Information();

		//各種counter
		_g->APDcount = 0;
		_g->Apqtime = 0.0;
		_g->APDtime = 0.0;
		//_g->Mode_Analytic();

		//Python用の出力
		//_g->OP_python_MC();
		//_g->OP_python_Stiff();
		//_g->OP_python_Damping();

		//初期座標のベクトル生成
		_g->InitialVector = Eigen::VectorXd::Zero(3 * _g->particle_num);
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			_g->InitialVector.block(3 * pi, 0, 3, 1) = _g->particles[pi]->Get_Initial_Pos();
		}
		//グループごとの座標を初期化

		//自由浮動させて、ゴースト力がでないことを実験
		//初期座標を回転させる場合
		/*
		Eigen::Matrix3d InitialR = Eigen::Matrix3d::Zero();
		InitialR(0, 0) = 1.0; InitialR(1, 1) = cos(PI / 4); InitialR(2, 2) = cos(PI / 4);
		InitialR(1, 2) = -1 * sin(PI / 4); InitialR(2, 1) = sin(PI / 4);
		*/

		//初期座標のベクトル生成
		_g->GroupGridVector = Eigen::VectorXd::Zero(3 * _g->particle_num);
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			//自由浮遊の実験では下を使用（回転行列InitialRを掛ける）
			//_g->GroupGridVector.block(3 * pi, 0, 3, 1) = InitialR * _g->particles[pi]->Get_Initial_Pos();
			_g->GroupGridVector.block(3 * pi, 0, 3, 1) = _g->particles[pi]->Get_Initial_Pos();
		}
		//グループごとの速度の初期化
		_g->GroupVelVector = Eigen::VectorXd::Zero(3 * _g->particle_num);
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			_g->GroupVelVector.block(3 * pi, 0, 3, 1) = Eigen::Vector3d::Zero();
			if (!(_g->particles[pi]->Is_Fixed())) {
				//速度に任意の値をいれることができる
				//_g->GroupVelVector.block(3 * pi, 0, 3, 1) += Eigen::Vector3d(0, 10, 0);
			}
		}
		//節点の質量の確認
		/*
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			std::cout << _g->M_Matrix_C(3 * pi, 3 * pi) - _g->Get_GMass_In_Group(_g->particles[pi]->p_id) << std::endl;
		}
		*/
		/*
		Eigen::VectorXd InitialVector2 = Eigen::VectorXd::Zero(3 * _g->particle_num);
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			InitialVector2.block(3 * pi, 0, 3, 1) = _g->particles[pi]->Get_Grid();
		}
		std::cout << _g->InitialVector - InitialVector2 << std::endl;
		*/
		//SUM_M_Matrixが重心を計算できているのか確認
		//std::cout << "centorgrid2" << std::endl;
		//std::cout << _g->SUM_M_Matrix *_g->InitialVector << std::endl;

		/*
		std::cout << "OrigineVector" << std::endl;
		std::cout << _g->OrigineVector << std::endl;
		std::cout << "OrigineVector2" << std::endl;

		//Eigen::VectorXd tempO = Eigen::MatrixXd::Identity(3 * _g->particle_num, 3 * _g->particle_num)*_g->InitialVector;
		//Eigen::MatrixXd tempO = Eigen::MatrixXd::Zero(3 * _g->particle_num, 3 * _g->particle_num);
		//tempO = -1 * _g->SUM_M_Matrix;
		// (unsigned int pi = 0; pi < 3 * _g->particle_num;pi++) {
		//	tempO(pi, pi) += 1.0 * _g->Get_Group_Mass();
		//}
		//tempO /= _g->Get_Group_Mass();
		//Eigen::VectorXd tempCM = -1 * _g->SUM_M_Matrix *_g->InitialVector;
		std::cout << Eigen::MatrixXd::Identity(3*_g->particle_num, 3*_g->particle_num)*_g->InitialVector- _g->SUM_M_Matrix *_g->InitialVector << std::endl;
		//std::cout << Eigen::MatrixXd::Identity(3 * _g->particle_num, 3 * _g->particle_num)*_g->InitialVector  << std::endl;
		//std::cout << -1* _g->SUM_M_Matrix *_g->InitialVector << std::endl;
		//std::cout << tempO + tempCM << std::endl;
		//std::cout << (tempO * _g->InitialVector) << std::endl;
		*/
	}
	//std::cout << "Stiffenss Time is " << StiffnessTime << std::endl;
	//Debug用
	GMREScount = 0;
	ConbiteGMRES = Eigen::VectorXd::Zero(20);
	//Debug出力
	std::ofstream outputfile("MinvK.txt");
	for (auto _g : groups) {
		Eigen::MatrixXd Define = Eigen::MatrixXd::Zero(_g->particle_num, _g->particle_num);
		Define = _g->M_Matrix_C.inverse() *_g->stiffness_matrix;
		outputfile << "np.array([";
		for (unsigned int pi = 0; pi < 3 * _g->particle_num - 1; pi++) {
			//値
			outputfile << "[";
			for (unsigned int pj = 0; pj < 3 * _g->particle_num - 1; pj++) {
				//値
				outputfile << Define.block(pi, pj, 1, 1) << ",";
			}
			outputfile << Define.block(pi, 3 * _g->particle_num - 1, 1, 1) << "]," << "\n";
		}
		outputfile << "[";
		for (unsigned int pj = 0; pj < 3 * _g->particle_num - 1; pj++) {
			//値
			outputfile << Define.block(3 * _g->particle_num - 1, pj, 1, 1) << ",";
		}
		outputfile << Define.block(3 * _g->particle_num - 1, 3 * _g->particle_num - 1, 1, 1) << "]])" << "\n";
		outputfile << "\n" << "\n" << "\n" << "\n" << "\n";
	}
	outputfile.close();
	std::cout << "OutPut Define" << std::endl;
}
//==========================================================================//
//	@end		   			  グループ作成								    //
//==========================================================================//


// 物体の中心座標
Eigen::Vector3d UseBlockObjectDouble::Get_Center_Grid(std::vector<ParticleD*> p) {
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	for (auto _p : p) {
		center += _p->Get_Grid();
	}
	return center / p.size();
}
// グループの体積を計算して、グループの質量、剛性行列を作成する
void UseBlockObjectDouble::Create_Group(std::vector<TetraElementD*> tetra_set, int tetra_group_id) {
	if (tetra_set.size() == 0) { return; }
	double volume = 0.0;
	for (auto _t : tetra_set) {
		volume += _t->Get_Volume();
		tetra_map[_t] = true;
	}
	TetraGroupD* g = new TetraGroupD(tetra_set, data, data.density * volume, tetra_group_id);
	groups.push_back(g);

	//頂点がどこのグループに所属したのか記録する
	std::vector<ParticleD*> temp_p = g->Get_Particle();
	for (unsigned int i = 0; i < temp_p.size(); i++) {
		temp_p[i]->p_belong_TetraGroup_ids.push_back(tetra_group_id);
	}
}
//==========================================================================//
//	@end		   				初期設定									//
//==========================================================================//
//==========================================================================//
//	@start		   				ループ設定									//
//==========================================================================//
void UseBlockObjectDouble::Update() {
	double Exptime = 0.0;
	double Rotatetime = 0.0;
	for (auto _g : groups) {
		//弾性力以外の力による位置更新を計算
		//std::cout << "Hello STEP" << countup << std::endl;
		//std::cout << std::fixed << std::setprecision(10) << particles[3]->Get_Grid() << std::endl;
		//初期速度代入
		//particles[3]->Update_Velocity(Eigen::Vector3d(0.0, 0.1, 0.0));

		mtCEPos.startMyTimer();
		//_g->Calc_Exp_Pos3();
		_g->Calc_Exp_Pos_Group();//Predicted position xj(i+1)
		mtCEPos.endMyTimer();
		Exptime += mtCEPos.getDt();
		//std::cout<<"Hello EXP"<<countup<<std::endl;
		//std::cout << std::fixed << std::setprecision(10) <<particles[3]->Get_Exp_Pos() <<std::endl;

		//回転を計算
		mtUpRotate.startMyTimer();
		if (useUpRotate) {
			_g->Update_Rotate();
			//_g->Update_Rotate2();
		}
		else {
			_g->rotate_matrix = Eigen::Matrix3d::Identity();
			Eigen::MatrixXd rotate_matrix3N = Eigen::MatrixXd::Identity(3 * particles.size(), 3 * particles.size());
			for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
				rotate_matrix3N.block(3 * pi, 3 * pi, 3, 3) = _g->rotate_matrix;
			}
			_g->Rn_Matrix_Sparse.setZero();
			_g->Rn_MatrixTR_Sparse.setZero();
			_g->Rn_Matrix_Sparse = rotate_matrix3N.sparseView();
			_g->Rn_MatrixTR_Sparse = (rotate_matrix3N.transpose()).sparseView();
		}
		mtUpRotate.endMyTimer();
		Rotatetime += mtUpRotate.getDt();
		//std::cout << "create Rotate" << std::endl;
	}
	//std::cout << "Exp Time is " << Exptime << std::endl;
	//std::cout << "Rotate Time is " << Rotatetime << std::endl;
	//予測位置の平均を節点が保持しておく(Objectがvectorとしてもっていてもいい)
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Set_Exp_Pos(Calc_New_Exp_Pos_Mean(_p));
		}
		else {
			//std::cout << "fixed:" << _p->p_id<<"tetra"<< _p->p_belong_TetraGroup_ids.size() << std::endl;
			_p->Set_Exp_Pos(_p->Get_Grid());
		}
		if (fetestexcept(FE_INVALID)) {
			std::cout << "FE_INVALID PBD 263" << std::endl;
		}
		feclearexcept(FE_ALL_EXCEPT);
	}
	//EXPがちゃんと計算できているか確認
	/*
	std::cout << "Initial" << std::endl;
	for (auto _g : groups) {
		_g->Calc_Convergence();
	}
	*/
	//std::cout << "OK" << std::endl;
	mtCconstr.startMyTimer();
	MicroSecondTimer mtUpadateJacobi;
	MicroSecondTimer mtPrecomp;
	MicroSecondTimer mtFbinditera;
	mtUpadateJacobi.setid(30);
	mtPrecomp.setid(32);
	mtFbinditera.setid(33);
	double JacobiTime = 0.0;
	double PrecompTime = 0.0;
	double FbinditeraTime = 0.0;
	Solve_Constraints13(PBD_LOOP);
	//差分法かLUかを選択する
	//if (whichmethodused) {
	//	for (auto _g : groups) {
	//		mtCP_3.startMyTimer();
	//		mtUpadateJacobi.startMyTimer();
	//		//初期化してもしなくてもいい
	//		_g->iterativeVector = Eigen::VectorXd::Zero(3 * _g->particles.size());
	//		if (useSparse) {
	//			_g->Calc_Jacobi_Matrix_iteration_Sparse();
	//			_g->Calc_Constant_term_iteration_Sparse();
	//		}
	//		else {
	//			_g->Calc_Jacobi_Matrix_iteration();
	//			_g->Calc_Constant_term_iteration2();
	//		}
	//		mtUpadateJacobi.endMyTimer();
	//		JacobiTime += mtUpadateJacobi.getDt();
	//		//Jacobi行列の前処理
	//		mtPrecomp.startMyTimer();
	//		_g->Calc_GMRES_Pre();
	//		mtPrecomp.endMyTimer();
	//		PrecompTime += mtPrecomp.getDt();
	//		mtCP_3.endMyTimer();
	//	}
	//	//std::cout << "Jacobi time is " << JacobiTime << std::endl;
	//	//std::cout << "Precomp time is " << PrecompTime << std::endl;
	//	//反復法
	//	//std::cout << "OK" << std::endl;
	//	mtFbinditera.startMyTimer();
	//	Solve_Constraints12(PBD_LOOP);
	//	mtFbinditera.endMyTimer();
	//	FbinditeraTime = mtFbinditera.getDt();
		//std::cout << "Fbinditera time is " << FbinditeraTime << std::endl;
	//}
	//else {
	//	for (auto _g : groups) {
	//		mtCP_3.startMyTimer();
	//		//初期化してもしなくてもいい
	//		//_g->iterativeVector = Eigen::VectorXd::Zero(3 * _g->particles.size());
	//		_g->Calc_Jacobi_Matrix_iteration();
	//		_g->Calc_Constant_term_iteration2();
	//		mtCP_3.endMyTimer();
	//	}
	//	//反復法
	//	Solve_Constraints10_LU(PBD_LOOP);
	//}
	mtCconstr.endMyTimer();
	std::ostringstream sstr;
	std::ostringstream sstr2;
	std::ostringstream sstr3;
	std::ostringstream sstr4;
	std::ostringstream sstr5;
	std::ostringstream sstr6;
	//各計算にかかった時間を出力
	sstr << std::fixed;
	unsigned int string_color = GetColor(255, 255, 255);
	sstr << "Rotate is " << std::setprecision(4) << mtUpRotate.getDt() << ",EXFis " << std::setprecision(4) << mtCEPos.getDt()
		<< ",Jacobi is " << std::setprecision(4) << mtCP_3.getDt()
		/*<< ",FEM is " << std::setprecision(4) << mtCP_2.getDt() */
		<< ",Const" << std::setprecision(4) << mtCconstr.getDt() << std::endl;
	//     回転の計算時間,弾性力以外の力による位置更新の計算時間(ms)
	//	   Jacobi行列の更新にかかる時間
	//	   一回目の有限要素法による位置更新の計算時間(グループの総和),制約条件による位置の修正にかかる時間
	DrawString(0, 15, sstr.str().data(), string_color);
	sstr.str("");

	//マウスの座標を出力(外力を働かせるようにもできる)
	/*
	int Mouse = GetMouseInput();
	if (Mouse & MOUSE_INPUT_LEFT) {
		int x, y;
		GetMousePoint(&x, &y);
		sstr2 << "MousePoint(" << x << ", " << y << ")" << std::endl;
		sstr3 << "Force is (" << Outofforce[0] << "," << Outofforce[1] << "," << Outofforce[2] << ")" << "N" << std::endl;
		sstr3 << "Force is (" << Outofforce[0] << "," << Outofforce[1] << "," << Outofforce[2] << ")" << "N" << std::endl;
	}
	*/
	unsigned int string_color2 = GetColor(255, 255, 1);
	DrawString(0, 49, sstr2.str().data(), string_color2);
	unsigned int string_color3 = GetColor(255, 100, 100);
	DrawString(0, 67, sstr3.str().data(), string_color3);
	sstr2.str("");
	sstr3.str("");
}
//==========================================================================//
//	@end		   				ループ設定									//
//==========================================================================//

//==========================================================================//
//	@start		   			グループ作成用関数								//
//==========================================================================//
std::vector< TetraElementD* > UseBlockObjectDouble::Create_Group_Candidate() {
	//グループの候補の作成
	std::vector< TetraElementD* > temp_tetras;
	for (auto _t : tetras) {
		if (!tetra_map[_t]) {
			_t->Clear_Same_Group();
			if (temp_tetras.size() == 0) { temp_tetras.push_back(_t); }
		}
	}

	for (unsigned int i = 0; i < temp_tetras.size(); ++i) {
		for (auto _t : temp_tetras[i]->Get_Side_Tetra()) {
			if (!tetra_map[_t] && !temp_tetras[i]->Is_Same_Group(_t)) {
				temp_tetras[i]->Set_Same_Group(_t, true);
				_t->Set_Same_Group(temp_tetras[i], true);
				temp_tetras.push_back(_t);
			}
		}
	}

	std::sort(temp_tetras.begin(), temp_tetras.end());
	temp_tetras.erase(std::unique(temp_tetras.begin(), temp_tetras.end()), temp_tetras.end());
	return temp_tetras;
}

// グループを分割するための判断用関数
bool UseBlockObjectDouble::GroupdivideNML(TetraElementD *t, int a, int b, int c) {
	for (auto _p : t->Get_Particle()) {
		Eigen::Vector3d pk = _p->Get_Grid();
		if (pk[0] > (a + 1) * ((xsize - 1) / dividedbyn) * sidelength) {
			return 0;
		}
		if (pk[1] >  (b + 1) * ((ysize - 1)  / dividedbym) * sidelength) {
			return 0;
		}
		if (pk[2] >  (c + 1) * ((zsize - 1) / dividedbyl) * sidelength) {
			return 0;
		}
	}
	return 1;
}

//==========================================================================//
//	@end		   			グループ作成用関数								//
//==========================================================================//