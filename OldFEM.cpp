//===========================================================================//
//@author IKUTO ENDO
//@brief OldFEM.hの実装
//===========================================================================//
#include "OldFEM.h"

//===========================================================================//
OldFEM::OldFEM(std::vector<Particle*> p, ObjectData data)
	: Object(p, data)
{
	std::cout << "Create Object Using Local Stiffness Matrix with Group" << std::endl;
	Init();
}
OldFEM::~OldFEM() {
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
void OldFEM::Init() {
	Delaunay_Triangulation();
	Create_Groups();
	data_name = "OldFEM";
	mtUpRotate.setid(4);
	mtCEPos.setid(5);
	mtCFEM.setid(6);
	mtCconstr.setid(7);
}

void OldFEM::Create_Groups() {
	//----------------------------------------------------------------------//
	//-------------------------TetraSet 初期設定----------------------------//
	int tetra_group_id = 0;

	for (auto _t : tetras) {
		_t->Create_Faces();
		tetra_map.insert(std::map<TetraElement*, bool>::value_type(_t, false));
	}
	for (int t1 = 0; t1 < tetras.size() - 1; ++t1) {
		for (int t2 = t1 + 1; t2 < tetras.size(); ++t2) {
			if (tetras[t1]->has_Common_Face(tetras[t2])) {
				tetras[t1]->Push_Side_Tetra(tetras[t2]);
				tetras[t2]->Push_Side_Tetra(tetras[t1]);
			}
		}
	}

	Create_Group(tetras, 0);
	std::cout << "success create All-Group. All-Group size : " << groups.size() << std::endl;
}




Eigen::Vector3f OldFEM::Get_Center_Grid(std::vector<Particle*> p) {
	Eigen::Vector3f center = Eigen::Vector3f::Zero();
	for (auto _p : p) {
		center += _p->Get_Grid();
	}
	return center / p.size();
}

void OldFEM::Create_Group(std::vector<TetraElement*> tetra_set, int tetra_group_id) {
	if (tetra_set.size() == 0) { return; }
	double volume = 0;
	for (auto _t : tetra_set) {
		volume += _t->Get_Volume();
		tetra_map[_t] = true;
	}
	TetraGroup* g = new TetraGroup(tetra_set, data, data.density * volume, tetra_group_id);
	groups.push_back(g);
	std::vector<Particle*> temp_p = g->Get_Particle();
	for (int i = 0; i < temp_p.size(); i++) {
		//belong_group[temp_p[i]].push_back(g);
		temp_p[i]->p_belong_TetraGroup_ids.push_back(tetra_group_id);
	}
}
//==========================================================================//
//	@end		   				初期設定									//
//==========================================================================//
//==========================================================================//
//	@start		   				ループ設定									//
//==========================================================================//
void OldFEM::Update() {
	for (auto _g : groups) {
		Create_Groups();
		mtCFEM.startMyTimer();
		_g->Calc_FEM();
		mtCFEM.endMyTimer();
	}
	mtCconstr.startMyTimer();
	Solve_Constraints(PBD_LOOP);
	mtCconstr.endMyTimer();

	std::ostringstream sstr;
	unsigned int string_color = GetColor(255, 255, 255);
	sstr << "Rotate is " << mtUpRotate.getDt() << ",EXFis " << mtCEPos.getDt()
		<< ",FEM is " << mtCFEM.getDt() << ",Const" << mtCconstr.getDt();

	DrawString(0, 15, sstr.str().data(), string_color);
	sstr.str("");

}
//==========================================================================//
//	@end		   				ループ設定									//
//==========================================================================//
//==========================================================================//
//	@start		   			グループ作成用関数								//
//==========================================================================//
std::vector< TetraElement* > OldFEM::Create_Group_Candidate() {
	//グループの候補の作成
	std::vector< TetraElement* > temp_tetras;
	for (auto _t : tetras) {
		if (!tetra_map[_t]) {
			_t->Clear_Same_Group();
			if (temp_tetras.size() == 0) { temp_tetras.push_back(_t); }
		}
	}

	for (int i = 0; i < temp_tetras.size(); ++i) {
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

bool OldFEM::GroupdivideNML(TetraElement *t, int a, int b, int c) {
	for (auto _p : t->Get_Particle()) {
		Eigen::Vector3f pk = _p->Get_Grid();
		if (pk[0] < 100 + (dividedbyn - a)*(13 - 1) * 40 / dividedbyn) {
			return 0;
		}
		if (pk[1] < 200 + (dividedbym - b)*(5 - 1) * 40 / dividedbym) {
			return 0;
		}
		if (pk[2] < 0 + (dividedbyl - c)*(5 - 1) * 40 / dividedbyl) {
			return 0;
		}
	}
	return 1;
}

//==========================================================================//
//	@end		   			グループ作成用関数								//
//==========================================================================//