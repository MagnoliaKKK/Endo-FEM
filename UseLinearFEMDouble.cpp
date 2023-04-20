//===========================================================================//
//@author IKUTO ENDO
//@brief UseLinearFEMDouble.hの実装
//===========================================================================//
#include "UseLinearFEMDouble.h"
//オブジェクトをN(x),M(y),L(z)で分割する
//===========================================================================//
UseLinearFEMDouble::UseLinearFEMDouble(std::vector<ParticleD*> p, ObjectData data) //コンストラクタ
	: ObjectD(p, data)//親クラスはObjectクラス
{
	std::cout << "Create Object Using Local Stiffness Matrix with Group" << std::endl;
	Init();			//初期設定
}
UseLinearFEMDouble::~UseLinearFEMDouble() { //デコンストラクタ
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
void UseLinearFEMDouble::Init() {
	Delaunay_Triangulation();  //ドロネー三角形分割をする
	Create_Groups();		   //グループを作る
	data_name = "Use_LinearLocalFEM";//オブジェクトの名前をUse_GroupOneTetraとする
}
//==========================================================================//
//	@start		   			  グループ作成								    //
//==========================================================================//
void UseLinearFEMDouble::Create_Groups() {
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
	//--------------------------N個目までの作成-----------------------------//
	Create_Group(remain, 0);//remain=tetras(すべての四面体要素)
	remain = this->Create_Group_Candidate();
	if (remain.empty()) {
		std::cout << "All complete : " << std::endl;
	}
	else {
		std::cout << "ERROR(Create Group)" << std::endl;
	}
	//----------------------------------------------------------------------//
	std::cout << "success create All-Group. All-Group size : " << groups.size() << std::endl;
	//グループ分け完了

	//共有節点行列の作成
	//linearなのでしないはず
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
		std::cout << "オブジェクトの共有節点行列 is " << std::endl;
		std::cout << Node_Sharing_Matrix.block(0, 0, 3 * Share_particlenum, 3 * Sum_particlenum) << std::endl;
		std::cout << "オブジェクトの共有数 is " << Share_particlenum << std::endl;
		std::cout << "オブジェクトの共有節点 is " << std::endl;
		for (int i = 0; i < Share_particlenum; i++) {
			std::cout << Share_particle_id[i] << std::endl;
		}
		for (auto _g : groups) {
			AllParticlenum += _g->particle_num;
		}
		std::cout << "オブジェクトの全数 is " << AllParticlenum << std::endl;
	}
	else {
		std::cout << "You are  Linear!!" << std::endl;
	}

	//各種パラメータの設定(質量やKマトリックス)
	//質量の設定
	std::cout << particles.size() << std::endl;
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
	double tempM = 0.0;
	for (auto _g : groups) {
		tempM = 0.0;
		//モデル全体の質量から計算したM行列の作成
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			_g->M_Matrix_C.block(3 * pi, 3 * pi, 3, 3) = _g->particles[pi]->Get_Mass() * Eigen::Matrix3d::Identity(3, 3);
			tempM += _g->particles[pi]->Get_Mass();
		}
		//std::cout << _g->M_Matrix_C << std::endl;
		//本当はここでtempMで割ったやつをSUM_Mにしたほうがいいかもしれない
		_g->Create_SUM_M_Matrix();
		//_g->Create_SUM_M_Matrix(tempM);
		_g->Set_Group_Mass(tempM);
		//std::cout << _g->Get_Group_Mass() << std::endl;
		_g->Create_Center_Grid();
		_g->Create_Local_Stiffness_Matrix();
		_g->Create_Damping_Matrix();

		_g->Create_Information();
		//_g->Mode_Analytic();
		_g->OP_python_MC();
		_g->OP_python_Stiff();
		_g->OP_python_Damping();
	}
}
//==========================================================================//
//	@end		   			  グループ作成								    //
//==========================================================================//


// 物体の中心座標
Eigen::Vector3d UseLinearFEMDouble::Get_Center_Grid(std::vector<ParticleD*> p) {
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	for (auto _p : p) {
		center += _p->Get_Grid();
	}
	return center / p.size();
}
// グループの体積を計算して、グループの質量、剛性行列を作成する
void UseLinearFEMDouble::Create_Group(std::vector<TetraElementD*> tetra_set, int tetra_group_id) {
	if (tetra_set.size() == 0) { return; }
	double volume = 0;
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
void UseLinearFEMDouble::Update() {
	for (auto _g : groups) {
		//弾性力以外の力による位置更新を計算
		//std::cout << "Hello STEP" << countup << std::endl;
		//std::cout << std::fixed << std::setprecision(10) << particles[3]->Get_Grid() << std::endl;
		//初期速度代入
		//particles[3]->Update_Velocity(Eigen::Vector3d(0.0, 0.1, 0.0));

		mtCEPos.startMyTimer();
		_g->Calc_Exp_Pos2();
		mtCEPos.endMyTimer();

		//回転を計算
		mtUpRotate.startMyTimer();
		_g->Update_Rotate();
		mtUpRotate.endMyTimer();

	}
	mtCconstr.startMyTimer();
	for (auto _g : groups) {
		_g->iterativeVector = Eigen::VectorXd::Zero(3 * _g->particles.size());
		_g->Calc_Jacobi_Matrix_iteration();
		//差分のほうの定数値を計算
		_g->Calc_Constant_term_iteration2();
	}
	Solve_Constraints10(PBD_LOOP);
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
	sstr << "Rotate is " << std::setprecision(5) << mtUpRotate.getDt() << ",EXFis " << std::setprecision(5) << mtCEPos.getDt() 
		<< ",Jacobi is " << std::setprecision(5) << mtCP_3.getDt()
		<< ",FEM is " << std::setprecision(5) << mtCP_2.getDt() << ",Const" << std::setprecision(5) << mtCconstr.getDt() << std::endl;
	//     回転の計算時間,弾性力以外の力による位置更新の計算時間(ms)
	//	   一回目の有限要素法による位置更新の計算時間,制約条件による位置の修正にかかる時間
	DrawString(0, 15, sstr.str().data(), string_color);
	sstr.str("");

	//マウスの座標を出力(外力が働く)
	int Mouse = GetMouseInput();
	if (Mouse & MOUSE_INPUT_LEFT) {
		int x, y;
		GetMousePoint(&x, &y);
		sstr2 << "MousePoint(" << x << ", " << y << ")" << std::endl;
		sstr3 << "Force is (" << Outofforce[0] << "," << Outofforce[1] << "," << Outofforce[3] << ")" << "N" << std::endl;
	}
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
std::vector< TetraElementD* > UseLinearFEMDouble::Create_Group_Candidate() {
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

//==========================================================================//
//	@end		   			グループ作成用関数								//
//==========================================================================//