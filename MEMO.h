//依然使ってたLocalStiffnessObject.hとLocalStiffnessObject.cpp

//LocalStiffnessObject.h
//===========================================================================//
//@author KatsuyaKikuchi
//@brief 局所合成行列を用いた手法でシミュレーションを行うObject
//===========================================================================//
//#ifndef _LOCALSTIFFNESSOBJECT
//#define _LOCALSTIFFNESSOBJECT
//
//#include "Object.h"
//class LocalStiffnessObject : public Object{
//public:
//	LocalStiffnessObject(std::vector<Particle*> p, ObjectData data, unsigned int flag);
//	~LocalStiffnessObject();
//
//	void Update();
//
//private:
//	void Init();
//	void Not_Using_Group();
//	void Create_Groups();
//	void Big_Group();
//	void Box_Group();
//	void Linear();
//
//
//	void Ortho_Skew_Method(std::vector<TetraElement*> tetra_set);
//	void Create_Group(std::vector<TetraElement*> tetra_set);
//	std::vector<Triangle*> Get_Outside_Face(std::vector<TetraElement*> tetra_set);
//	std::vector< TetraElement* > Get_Max_Group(std::vector<TetraElement*> tetra_set);
//	double Calc_Threshhold(std::vector<TetraElement*> tetra_set);
//	double Calc_Ortho_Skew(Eigen::Vector3f center, Triangle* face);
//	double Calc_Ortho_Skew(std::vector<TetraElement*> tetra_set, std::vector<Triangle*> outside_face, std::vector<int> &max_indexs);
//	Eigen::Vector3f Get_Center_Grid(std::vector<Particle*> p);
//	std::map<TetraElement*, bool> tetra_map;//TetraElementがすでにグループに使用されているかどうかのmap
//	std::vector< TetraElement* > Create_Group_Candidate();
//	std::vector< Particle* > Get_Outside_Particles(std::vector<TetraElement*> temp_tetras);
//
//
//	//debug
//	unsigned int flag;
//
//};
//#endif
//LocalStiffnessObject.cpp

//===========================================================================//
//@author KatsuyaKikuchi
//@brief LocalStiffnessObject.hの実装
//===========================================================================//
//#include "LocalStiffnessObject.h"
//
////===========================================================================//
//LocalStiffnessObject::LocalStiffnessObject(std::vector<Particle*> p, ObjectData data, unsigned int flag)
//	: Object(p, data), flag(flag)
//{
//	std::cout << "Create Object Using Local Stiffness Matrix" << std::endl;
//	Init();
//}
//LocalStiffnessObject::~LocalStiffnessObject(){
//	for (auto _g : groups){
//		delete _g;
//	}
//	tetras.clear();
//}
////===========================================================================//
////==========================================================================//
////	@start		   				初期設定									//
////==========================================================================//
//void LocalStiffnessObject::Init(){
//	Delaunay_Triangulation();
//	switch (flag){
//	case 0:
//		Not_Using_Group();
//		data_name = "No_Group";
//		break;
//	case 1:
//		Create_Groups();
//		data_name = "Use_Group";
//		break;
//	case 2:
//		Big_Group();
//		data_name = "Big_Group";
//		break;
//	case 3:
//		Box_Group();
//		data_name = "Box_Group";
//		break;
//	case 5:
//		Linear();
//		data_name = "Linear";
//		break;
//	default:
//		break;
//	}
//
//
//	std::vector<Triangle*> faces[6];
//
//	for (auto _f : outside_faces){
//		for (int i = 0; i < 6; i++){
//			if (faces[i].size() == 0){
//				faces[i].push_back(_f);
//				break;
//			}
//			if ((faces[i])[0]->Get_Area_Vec().dot(_f->Get_Area_Vec()) / ((faces[i])[0]->Get_Area_Vec().norm()*(_f->Get_Area_Vec()).norm()) > 0.90){
//				faces[i].push_back(_f);
//				break;
//			}
//		}
//	}
//
//
//	Eigen::Vector3f x, y, z;
//	x = Eigen::Vector3f(1, 0, 0);
//	y = Eigen::Vector3f(0, 1, 0);
//	z = Eigen::Vector3f(0, 0, 1);
//	for (auto _f : faces){
//		if (_f.size() == 18 && _f[0]->Get_Area_Vec().dot(x) > 0.9){
//			draw_faces[1] = _f;
//			continue;
//		}
//
//		if (_f[0]->Get_Area_Vec().dot(-z) > 0.9){
//			draw_faces[2] = _f;
//			continue;
//		}
//		if (_f[0]->Get_Area_Vec().dot(-y) > 0.9){
//			draw_faces[0] = _f;
//			continue;
//		}
//
//	}
//}
//
//void LocalStiffnessObject::Create_Groups(){
//	//----------------------------------------------------------------------//
//	//-------------------------TetraSet 初期設定----------------------------//
//	for (auto _t : tetras){
//		_t->Create_Faces();
//		tetra_map.insert(std::map<TetraElement*, bool>::value_type(_t, false));
//	}
//	for (int t1 = 0; t1 < tetras.size() - 1; ++t1){
//		for (int t2 = t1 + 1; t2 < tetras.size(); ++t2){
//			if (tetras[t1]->has_Common_Face(tetras[t2])){
//				tetras[t1]->Push_Side_Tetra(tetras[t2]);
//				tetras[t2]->Push_Side_Tetra(tetras[t1]);
//			}
//		}
//	}
//	//----------------------------------------------------------------------//
//	//-------------------------最外面Groupの作成----------------------------//
//	std::vector< TetraElement* > temp = this->Create_Group_Candidate();
//	std::vector<Particle*> outside_particles = this->Get_Outside_Particles(temp);
//	this->outside_faces = this->Get_Outside_Face(temp);
//
//	//最外面に配置されているparticleを所持しているTetraElementは、そのままグループとして利用する。
//	for (auto _p : outside_particles){
//		for (auto _t : tetras){
//			if (_t->has_Particle(_p) && !tetra_map[_t]){
//				_t->Clear_Same_Group();
//				//				tetra_map[_t] = true;
//				std::vector<TetraElement*> temp_t;
//				temp_t.push_back(_t);
//				Create_Group(temp_t);
//			}
//		}
//	}
//	std::cout << "success create Outside-Group. Outside-Group size : " << groups.size() << std::endl;
//	//----------------------------------------------------------------------//
//	//--------------------------内部Groupの作成-----------------------------//
//	std::vector< TetraElement* > next_t = this->Create_Group_Candidate();
//	//	Ortho_Skew_Method(next_t);
//	Create_Group(next_t);
//
//
//	std::cout << "success create All-Group. All-Group size : " << groups.size() << std::endl;
//}
//
//void LocalStiffnessObject::Ortho_Skew_Method(std::vector<TetraElement*> tetra_set){
//	if (tetra_set.size() == 0) {
//		std::cout << "Finish" << std::endl;
//		return;
//	}
//	std::vector<TetraElement*> copy_tetras = tetra_set;
//
//	while (true){
//		if (copy_tetras.size() == 0) {
//			std::cout << "BUG!!! : tetras size is 0" << std::endl;
//			return;
//		}
//
//		std::vector<int> max_indexs;
//		std::vector<Triangle*> outside_face = this->Get_Outside_Face(copy_tetras);
//		double ortho_skew = Calc_Ortho_Skew(copy_tetras, outside_face, max_indexs);
//		double threshhold = Calc_Threshhold(copy_tetras);
//
//		if (ortho_skew <= threshhold){
//			Create_Group(copy_tetras);
//			break;
//		}
//		else{
//			//copy_tetrasの修正
//			if (max_indexs.size() == copy_tetras.size()){
//				//最小のグループでも形が悪い
//				for (auto _t : copy_tetras){
//					std::vector<TetraElement*> temp;
//					temp.push_back(_t);
//					TetraGroup* g = new TetraGroup(temp, data, data.density * _t->Get_Volume());
//					groups.push_back(g);
//					std::vector<Particle*> temp_p = _t->Get_Particle();
//					for (auto _p : temp_p){
//						belong_group[_p].push_back(g);
//					}
//				}
//				break;
//			}
//			else{
//
//				for (auto _mi : max_indexs){
//					auto it = copy_tetras.begin();
//					while (it != copy_tetras.end()) {
//						if ((*it)->has_Face(outside_face[_mi])) {
//							for (auto _t : (*it)->Get_Side_Tetra()){
//								_t->Set_Same_Group((*it), false);
//							}
//							it = copy_tetras.erase(it);
//							break;
//						}
//						else { ++it; }
//					}
//				}
//				copy_tetras = this->Get_Max_Group(copy_tetras);
//			}
//		}
//	}
//
//	Ortho_Skew_Method(this->Create_Group_Candidate());
//}
//
//
//Eigen::Vector3f LocalStiffnessObject::Get_Center_Grid(std::vector<Particle*> p){
//	Eigen::Vector3f center = Eigen::Vector3f::Zero();
//	for (auto _p : p){
//		center += _p->Get_Grid();
//	}
//	return center / p.size();
//}
//
//void LocalStiffnessObject::Create_Group(std::vector<TetraElement*> tetra_set){
//	if (tetra_set.size() == 0){ return; }
//	double volume = 0;
//	for (auto _t : tetra_set){
//		volume += _t->Get_Volume();
//		tetra_map[_t] = true;
//	}
//	TetraGroup* g = new TetraGroup(tetra_set, data, data.density * volume);
//	groups.push_back(g);
//	std::vector<Particle*> temp_p = g->Get_Particle();
//	for (int i = 0; i < temp_p.size(); i++){
//		belong_group[temp_p[i]].push_back(g);
//	}
//}
//
//void LocalStiffnessObject::Not_Using_Group(){
//	//----------------------------------------------------------------------//
//	//-------------------------TetraSet 初期設定----------------------------//
//	for (auto _t : tetras){
//		_t->Create_Faces();
//		tetra_map.insert(std::map<TetraElement*, bool>::value_type(_t, false));
//	}
//	for (int t1 = 0; t1 < tetras.size() - 1; ++t1){
//		for (int t2 = t1 + 1; t2 < tetras.size(); ++t2){
//			if (tetras[t1]->has_Common_Face(tetras[t2])){
//				tetras[t1]->Push_Side_Tetra(tetras[t2]);
//				tetras[t2]->Push_Side_Tetra(tetras[t1]);
//			}
//		}
//	}
//	//----------------------------------------------------------------------//
//	//-------------------------最外面Groupの作成----------------------------//
//	std::vector< TetraElement* > temp = this->Create_Group_Candidate();
//	std::vector<Particle*> outside_particles = this->Get_Outside_Particles(temp);
//	this->outside_faces = this->Get_Outside_Face(temp);
//	for (auto _t : tetras){
//		std::vector<TetraElement*> temp;
//		temp.push_back(_t);
//		TetraGroup* g = new TetraGroup(temp, data, data.density * _t->Get_Volume());
//		groups.push_back(g);
//		std::vector<Particle*> temp_p = _t->Get_Particle();
//		for (int i = 0; i < 4; i++){
//			belong_group[temp_p[i]].push_back(g);
//		}
//	}
//}
//void LocalStiffnessObject::Big_Group(){
//	//----------------------------------------------------------------------//
//	//-------------------------TetraSet 初期設定----------------------------//
//	for (auto _t : tetras){
//		_t->Create_Faces();
//		tetra_map.insert(std::map<TetraElement*, bool>::value_type(_t, false));
//	}
//	for (int t1 = 0; t1 < tetras.size() - 1; ++t1){
//		for (int t2 = t1 + 1; t2 < tetras.size(); ++t2){
//			if (tetras[t1]->has_Common_Face(tetras[t2])){
//				tetras[t1]->Push_Side_Tetra(tetras[t2]);
//				tetras[t2]->Push_Side_Tetra(tetras[t1]);
//			}
//		}
//	}
//	//----------------------------------------------------------------------//
//	//-------------------------最外面Groupの作成----------------------------//
//	std::vector< TetraElement* > temp = this->Create_Group_Candidate();
//	std::vector<Particle*> outside_particles = this->Get_Outside_Particles(temp);
//	this->outside_faces = this->Get_Outside_Face(temp);
//
//	//最外面に配置されているparticleを所持しているTetraElementは、そのままグループとして利用する。
//	for (auto _p : outside_particles){
//		for (auto _t : tetras){
//			if (_t->has_Particle(_p) && !tetra_map[_t]){
//				_t->Clear_Same_Group();
//				//				tetra_map[_t] = true;
//				std::vector<TetraElement*> temp_t;
//				temp_t.push_back(_t);
//				Create_Group(temp_t);
//			}
//		}
//	}
//	std::cout << "success create Outside-Group. Outside-Group size : " << groups.size() << std::endl;
//	//----------------------------------------------------------------------//
//	//--------------------------内部Groupの作成-----------------------------//
//	std::vector< TetraElement* > next_t = this->Create_Group_Candidate();
//	Create_Group(next_t);
//	std::cout << "success create All-Group. All-Group size : " << groups.size() << std::endl;
//}
//void LocalStiffnessObject::Linear(){
//	//----------------------------------------------------------------------//
//	//-------------------------TetraSet 初期設定----------------------------//
//	for (auto _t : tetras){
//		_t->Create_Faces();
//		tetra_map.insert(std::map<TetraElement*, bool>::value_type(_t, false));
//	}
//	for (int t1 = 0; t1 < tetras.size() - 1; ++t1){
//		for (int t2 = t1 + 1; t2 < tetras.size(); ++t2){
//			if (tetras[t1]->has_Common_Face(tetras[t2])){
//				tetras[t1]->Push_Side_Tetra(tetras[t2]);
//				tetras[t2]->Push_Side_Tetra(tetras[t1]);
//			}
//		}
//	}
//	//----------------------------------------------------------------------//
//	//-------------------------最外面Groupの作成----------------------------//
//	std::vector< TetraElement* > temp = this->Create_Group_Candidate();
//	std::vector<Particle*> outside_particles = this->Get_Outside_Particles(temp);
//	this->outside_faces = this->Get_Outside_Face(temp);
//
//	Create_Group(tetras);
//	std::cout << "success create All-Group. All-Group size : " << groups.size() << std::endl;
//}
//void LocalStiffnessObject::Box_Group(){
//	//----------------------------------------------------------------------//
//	//-------------------------TetraSet 初期設定----------------------------//
//	for (auto _t : tetras){
//		_t->Create_Faces();
//		tetra_map.insert(std::map<TetraElement*, bool>::value_type(_t, false));
//	}
//	for (int t1 = 0; t1 < tetras.size() - 1; ++t1){
//		for (int t2 = t1 + 1; t2 < tetras.size(); ++t2){
//			if (tetras[t1]->has_Common_Face(tetras[t2])){
//				tetras[t1]->Push_Side_Tetra(tetras[t2]);
//				tetras[t2]->Push_Side_Tetra(tetras[t1]);
//			}
//		}
//	}
//	//----------------------------------------------------------------------//
//	//-------------------------最外面Groupの作成----------------------------//
//	std::vector< TetraElement* > temp = this->Create_Group_Candidate();
//	std::vector<Particle*> outside_particles = this->Get_Outside_Particles(temp);
//	this->outside_faces = this->Get_Outside_Face(temp);
//
//
//	int x_size = particles.size() / 25;
//	std::cout << x_size << std::endl;
//	std::vector< TetraElement* > temp_t;
//	int group_size = 1;//tetras.size() / 3;
//	for (unsigned int i = 0; i < tetras.size(); ++i){
//		if (i % group_size == 0){ temp_t.clear(); }
//		temp_t.push_back(tetras[i]);
//		if (i % group_size == group_size - 1 || i == tetras.size() - 1){ Create_Group(temp_t); }
//	}
//	std::cout << "success create All-Group. All-Group size : " << groups.size() << std::endl;
//}
////==========================================================================//
////	@end		   				初期設定									//
////==========================================================================//
////==========================================================================//
////	@start		   				ループ設定									//
////==========================================================================//
//void LocalStiffnessObject::Update(){
//	for (auto _g : groups){
//		_g->Update_Rotate();
//		_g->Calc_Exp_Pos();
//		_g->Calc_FEM();
//	}
//
//	Solve_Constraints(PBD_LOOP);
//}
////==========================================================================//
////	@end		   				ループ設定									//
////==========================================================================//
////==========================================================================//
////	@start		   			グループ作成用関数								//
////==========================================================================//
//std::vector< TetraElement* > LocalStiffnessObject::Create_Group_Candidate(){
//	//グループの候補の作成
//	std::vector< TetraElement* > temp_tetras;
//	for (auto _t : tetras){
//		if (!tetra_map[_t]){
//			_t->Clear_Same_Group();
//			if (temp_tetras.size() == 0) { temp_tetras.push_back(_t); }
//		}
//	}
//
//	for (int i = 0; i < temp_tetras.size(); ++i){
//		for (auto _t : temp_tetras[i]->Get_Side_Tetra()){
//			if (!tetra_map[_t] && !temp_tetras[i]->Is_Same_Group(_t)){
//				temp_tetras[i]->Set_Same_Group(_t, true);
//				_t->Set_Same_Group(temp_tetras[i], true);
//				temp_tetras.push_back(_t);
//			}
//		}
//	}
//
//	std::sort(temp_tetras.begin(), temp_tetras.end());
//	temp_tetras.erase(std::unique(temp_tetras.begin(), temp_tetras.end()), temp_tetras.end());
//	return temp_tetras;
//}
//std::vector< Particle* > LocalStiffnessObject::Get_Outside_Particles(std::vector<TetraElement*> temp_tetras){
//	//最外面のParticleを取得
//	//XXX::Same Groupを更新してからじゃないと使えない
//	std::vector< Particle* > outside_particles;
//
//	for (auto _t : temp_tetras){
//		std::vector<Triangle*> temp_f = _t->Get_Free_Face();
//		for (auto _f : temp_f){
//			std::vector<Particle*> temp_p = _f->Get_Particle();
//			for (auto _p : temp_p){
//				outside_particles.push_back(_p);
//			}
//		}
//	}
//
//	std::sort(outside_particles.begin(), outside_particles.end());
//	outside_particles.erase(std::unique(outside_particles.begin(), outside_particles.end()), outside_particles.end());
//
//	return outside_particles;
//}
//std::vector<Triangle*> LocalStiffnessObject::Get_Outside_Face(std::vector<TetraElement*> tetra_set){
//	//Same Groupを更新してからじゃないと使えない
//	std::vector< Triangle* > outside_triangle;
//	for (auto _t : tetra_set){
//		std::vector< Triangle* > temp_e = _t->Get_Free_Face();
//		for (auto _e : temp_e){
//			outside_triangle.push_back(_e);
//		}
//	}
//
//	std::sort(outside_triangle.begin(), outside_triangle.end());
//	outside_triangle.erase(std::unique(outside_triangle.begin(), outside_triangle.end()), outside_triangle.end());
//	return outside_triangle;
//}
//std::vector< TetraElement* > LocalStiffnessObject::Get_Max_Group(std::vector<TetraElement*> tetra_set){
//	std::vector<TetraElement*> max_tetras;
//	std::map<TetraElement*, bool> counted;
//	for (auto _t : tetra_set) { counted[_t] = false; }
//
//	for (auto _t : tetra_set){
//		if (!counted[_t]){
//			std::vector<TetraElement*> temp_tetras;
//			temp_tetras.push_back(_t);
//			counted[_t] = true;
//			for (size_t i = 0; i < temp_tetras.size(); ++i){
//				for (auto _side_t : temp_tetras[i]->Get_Side_Tetra()){
//					if (!counted[_side_t] && temp_tetras[i]->Is_Same_Group(_side_t)){
//						temp_tetras.push_back(_side_t);
//						counted[_side_t] = true;
//					}
//				}
//			}
//			std::sort(temp_tetras.begin(), temp_tetras.end());
//			temp_tetras.erase(std::unique(temp_tetras.begin(), temp_tetras.end()), temp_tetras.end());
//			if (max_tetras.size() < temp_tetras.size()) max_tetras = temp_tetras;
//		}
//	}
//	return max_tetras;
//}
//
//double LocalStiffnessObject::Calc_Threshhold(std::vector<TetraElement*> tetra_set){
//	double threshhold = 0;
//	for (auto _t : tetra_set){
//		Eigen::Vector3f t_center = Get_Center_Grid(_t->Get_Particle());
//		std::vector<Triangle*> outside_face = _t->Get_Face();
//		double ortho_skew = 0;
//		for (auto _f : outside_face){
//			double temp = Calc_Ortho_Skew(t_center, _f);
//			if (ortho_skew < temp){ ortho_skew = temp; }
//		}
//		threshhold += ortho_skew;
//	}
//	return threshhold / tetra_set.size();
//}
//
//double LocalStiffnessObject::Calc_Ortho_Skew(Eigen::Vector3f center, Triangle* face){
//	Eigen::Vector3f s_vec = face->Get_Area_Vec();
//	Eigen::Vector3f centroid_vec = face->Get_Center_Grid() - center;
//	return 1.0 - s_vec.dot(centroid_vec) / (s_vec.norm()*centroid_vec.norm());
//}
//double LocalStiffnessObject::Calc_Ortho_Skew(std::vector<TetraElement*> tetra_set, std::vector<Triangle*> outside_face, std::vector<int> &max_indexs){
//	std::vector<Particle*> outside_p = this->Get_Outside_Particles(tetra_set);
//
//	Eigen::Vector3f center = this->Get_Center_Grid(outside_p);
//	double ortho_skew = 0.0;
//	for (int i = 0; i < outside_face.size(); ++i){
//		double temp = Calc_Ortho_Skew(center, outside_face[i]);
//
//		if (temp == ortho_skew) max_indexs.push_back(i);
//		if (temp > ortho_skew){
//			ortho_skew = temp;
//			max_indexs.clear();
//			max_indexs.push_back(i);
//		}
//	}
//	return ortho_skew;
//}
//==========================================================================//
//	@end		   			グループ作成用関数								//
//==========================================================================//