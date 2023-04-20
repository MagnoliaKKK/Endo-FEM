//===========================================================================//
//@author IKUTO ENDO
//@brief UseGroupOneTetraObject.hの実装
//===========================================================================//
#include "UseGroupOneTetraObjectDouble.h"
//オブジェクトをN(x),M(y),L(z)で分割する
//===========================================================================//
UseGroupOneTetraObjectDouble::UseGroupOneTetraObjectDouble(std::vector<ParticleD*> p, ObjectData data) //コンストラクタ
	: ObjectD(p, data)//親クラスはObjectクラス
{
	std::cout << "Create Object Using Local Stiffness Matrix with Group" << std::endl;
	Init();			//初期設定
}
UseGroupOneTetraObjectDouble::~UseGroupOneTetraObjectDouble() { //デコンストラクタ
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
void UseGroupOneTetraObjectDouble::Init() {
	Delaunay_Triangulation();  //ドロネー三角形分割をする
	Create_Groups();		   //グループを作る
	data_name = "Use_GroupOneTetra";//オブジェクトの名前をUse_GroupOneTetraとする
}
//==========================================================================//
//	@start		   			  グループ作成								    //
//==========================================================================//
void UseGroupOneTetraObjectDouble::Create_Groups() {
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
	Create_Group(tetras, 0);
	remain = this->Create_Group_Candidate();
	if (remain.empty()) {
		std::cout << "All complete : " << std::endl;
	}
	//----------------------------------------------------------------------//
	std::cout << "success create All-Group. All-Group size : " << groups.size() << std::endl;
}
//==========================================================================//
//	@end		   			  グループ作成								    //
//==========================================================================//


// 物体の中心座標
Eigen::Vector3d UseGroupOneTetraObjectDouble::Get_Center_Grid(std::vector<ParticleD*> p) {
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	for (auto _p : p) {
		center += _p->Get_Grid();
	}
	return center / p.size();
}
// グループの体積を計算して、グループの質量、剛性行列を作成する
void UseGroupOneTetraObjectDouble::Create_Group(std::vector<TetraElementD*> tetra_set, int tetra_group_id) {
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
void UseGroupOneTetraObjectDouble::Update() {
	for (auto _g : groups) {
		//弾性力以外の力による位置更新を計算
		//std::cout << "Hello STEP" << countup << std::endl;
		//std::cout << std::fixed << std::setprecision(10) << particles[3]->Get_Grid() << std::endl;
		//初期速度代入
		//particles[3]->Update_Velocity(Eigen::Vector3d(0.0, 0.1, 0.0));
		
		mtCEPos.startMyTimer();
		_g->Calc_Exp_Pos();
		mtCEPos.endMyTimer();
		//std::cout<<"Hello EXP"<<countup<<std::endl;
		//std::cout << std::fixed << std::setprecision(10) <<particles[3]->Get_Exp_Pos() <<std::endl;

		
		//回転を計算
		mtUpRotate.startMyTimer();
		_g->Update_Rotate();
		mtUpRotate.endMyTimer();
	
		//拘束力の中身が0かどうか確認
		//std::cout << "bindforce 0 is " << std::endl;
		//_g->Write_bind_force();
	}
	mtCconstr.startMyTimer();
	for (auto _g : groups) {
		_g->iterativeVector = Eigen::VectorXd::Zero(3 * _g->particles.size());
		_g->Calc_Jacobi_Matrix_iteration();
		//差分のほうの定数値(solve9)
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
std::vector< TetraElementD* > UseGroupOneTetraObjectDouble::Create_Group_Candidate() {
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
bool UseGroupOneTetraObjectDouble::GroupdivideNML(TetraElementD *t, int a, int b, int c) {
	for (auto _p : t->Get_Particle()) {
		Eigen::Vector3d pk = _p->Get_Grid();
		if (pk[0] < 100 + (dividedbyn - a)*(xsize - 1) * sidelength / dividedbyn) {
			return 0;
		}
		if (pk[1] < 200 + (dividedbym - b)*(ysize - 1) * sidelength / dividedbym) {
			return 0;
		}
		if (pk[2] < 0 + (dividedbyl - c)*(zsize - 1) * sidelength / dividedbyl) {
			return 0;
		}
	}
	return 1;
}

//==========================================================================//
//	@end		   			グループ作成用関数								//
//==========================================================================//