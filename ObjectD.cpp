//===========================================================================//
//@author KatsuyaKikuchi
//@brief Implementation of the ObjectD.h
//===========================================================================//

#include "ObjectD.h"
#include "Eigen/SVD"

//===========================================================================//
ObjectD::ObjectD(std::vector<ParticleD*> p, ObjectData data)
	: particles(p), data(data)//変数の初期化
	 ,Outofforce(Eigen::Vector3d::Zero())
{
	mtUpRotate.setid(4), mtCEPos.setid(5), mtCFEM.setid(6), mtCconstr.setid(7), mtCP_1.setid(8), mtCP_2.setid(9), mtCP_3.setid(10), mtEnergyConstraint.setid(11);
}	//stopwatchのidをセット
ObjectD::~ObjectD() {
}
//===========================================================================//

void ObjectD::Draw()const {
	for (auto _g : groups) {
		_g->Draw();
	}
}

void ObjectD::Solve_Constraints13(unsigned int loop) {
	for (auto _g : groups)
	{
		_g->ReSet_Fbind_Pos();
		//_g->LHS0();
	}
	int aa = 1;
	for (int i = 0; i < loop; i++)
	{
		
		//#pragma omp parallel for //作用最大 135 to 30
		for (int j = 0; j < static_cast<int>(groups.size()); j++) {
			auto _g = groups[j];

			_g->RHS0();

			_g->CalcDeltax();

			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID Deltax" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				_p->Set_Deltax_In_Model(CalcMeanPos(_p));
				_p->Set_Velocity_In_Model(Calc_Mean_Vel(_p));
			}
			else {
				_p->Set_Deltax_In_Model(CalcMeanPos(_p));
				_p->Set_Velocity_In_Model(Calc_Mean_Vel(_p));
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD 263" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		//#pragma omp parallel for //没明显作用
		for (int i = 0; i < static_cast<int>(groups.size()); i++) {
			auto _g = groups[i];
			_g->Update_Fbind_Pos8();
			
		}
		aa++;
	
	}
	for (auto _g : groups) {
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			//速度をいれてみた
			_g->GroupVelVector.block(3 * pi, 0, 3, 1) = (_g->PrimeVector.block(3 * pi, 0, 3, 1) + _g->DeltaxNew.block(3 * pi, 0, 3, 1) -_g->GroupGridVector.block(3 * pi, 0, 3, 1)) / TIME_STEP;
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID Calc_Exp6" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
			_g->GroupGridVector.block(3 * pi, 0, 3, 1) = _g->PrimeVector.block(3 * pi, 0, 3, 1) + _g->DeltaxNew.block(3 * pi, 0, 3, 1);
			if (_g->particles[pi]->Is_Fixed()) {
				_g->GroupVelVector.block(3 * pi, 0, 3, 1) = Eigen::Vector3d::Zero();
				_g->GroupGridVector.block(3 * pi, 0, 3, 1) = _g->InitialVector.block(3 * pi, 0, 3, 1);
			}
		}
		
	}
	
	
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Posi_set" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);
}
double ObjectD::Get_V() {
	double v = 0;
	for (auto _g : groups) {
		v += _g->Get_Volume();
	}
	return v;
}
double ObjectD::Get_M() {
	double m = 0.0;
	for (auto _p : particles) {
		m += _p->p_mass;
	}
	return m;
}
int ObjectD::Get_ezolg() {
	return ezolg;
}
Eigen::Vector3d ObjectD::Calc_New_Exp_Pos(ParticleD* p) {
	//if (belong_group[p].size() == 0){ return p->Get_Grid(); }
	if (p->p_belong_TetraGroup_ids.size() == 0) { return p->Get_Grid(); }
#ifdef _DEBUGMODE2 
	double C = 0;
	Eigen::Vector3f nabla_C = Eigen::Vector3f::Zero();

	for (size_t g1 = 0; g1 < belong_group[p].size() - 1; ++g1) {
		for (size_t g2 = g1 + 1; g2 < belong_group[p].size(); ++g2) {
			double temp = (belong_group[p][g1]->Get_Exp_Pos(p) - belong_group[p][g2]->Get_Exp_Pos(p)).norm();
			if (0 == temp) { continue; }
			C += temp;
			nabla_C += (belong_group[p][g1]->Get_Exp_Pos(p) - belong_group[p][g2]->Get_Exp_Pos(p)) / temp;
		}
	}
	Eigen::Vector3f exp_pos;
	if (0 == C) { exp_pos = belong_group[p][0]->Get_Exp_Pos(p); }
	else { exp_pos = (-C / (nabla_C.norm()*nabla_C.norm()))*nabla_C; }

	return exp_pos;
#else
	/*Eigen::Vector3f delta_p = Eigen::Vector3f::Zero();
	for (auto _g : belong_group[p]){
	delta_p += _g->Get_Exp_Pos(p->p_id);
	}
	return delta_p / belong_group[p].size();*/
	Eigen::Vector3d delta_p = Eigen::Vector3d::Zero();
	Eigen::Vector3d delta_k = Eigen::Vector3d::Zero();
	for (auto _g : p->p_belong_TetraGroup_ids) {
		TetraGroupD* tg = groups[_g];
		/*delta_p += tg->Get_Exp_Pos(p->p_id);*/
		delta_p += tg->Get_Exp_Pos(p->p_id);
		delta_k += tg->Get_Exp_Pos(p->p_id) - p->Get_Initial_Pos();
	}
	/*return delta_p - (p->p_belong_TetraGroup_ids.size()-1)*p->Get_Prime_Pos();*/
	/*return   delta_p/ p->p_belong_TetraGroup_ids.size() + delta_k;*/
	if (p->p_belong_TetraGroup_ids.size() == 1) {
		return delta_p;
	}
	else {
		/*std::cout << p->p_belong_TetraGroup_ids.size() << std::endl;*/
		return   delta_k + p->Get_Initial_Pos();
	}
	//return p->Get_Exp_Pos();
#endif

}
//グループごとの頂点に対して平均をとる(制約条件)
Eigen::Vector3d ObjectD::Calc_New_Exp_Pos_Mean(ParticleD* p) {
	if (p->p_belong_TetraGroup_ids.size() == 0) { return p->Get_Grid(); }
	Eigen::Vector3d delta_p = Eigen::Vector3d::Zero();
	for (auto _g : p->p_belong_TetraGroup_ids) {
		TetraGroupD* tg = groups[_g];
		delta_p += tg->Get_Exp_In_Group(p->p_id);
	}
	if (p->p_belong_TetraGroup_ids.size() == 1) {
		return delta_p;
	}
	else {
		delta_p = delta_p / p->p_belong_TetraGroup_ids.size();
		return   delta_p;
	}
}
//グループの節点の座標を物体の座標に変換する(制約条件)
Eigen::Vector3d ObjectD::Set_New_Exp_Pos(ParticleD* p) {
	if (p->p_belong_TetraGroup_ids.size() == 0) { return p->Get_Grid(); }
	Eigen::Vector3d delta_p = Eigen::Vector3d::Zero();
	for (auto _g : p->p_belong_TetraGroup_ids) {
		TetraGroupD* tg = groups[_g];
		delta_p += tg->Get_Exp_Pos(p->p_id);
		//std::cout << "tg" << tg->tetra_group_id << std::endl;
		//std::cout << "p_id" << p->p_id << std::endl;
	}
	if (p->p_belong_TetraGroup_ids.size() == 1) {
		return delta_p;
	}
	else {
		return   delta_p / p->p_belong_TetraGroup_ids.size();
	}
}
//グループの節点の座標を物体の座標に変換する(制約条件)
Eigen::Vector3d ObjectD::Set_New_Exp_Pos(ParticleD* p,Eigen::VectorXd Delta) {
	if (p->p_belong_TetraGroup_ids.size() == 0) { return p->Get_Grid(); }
	Eigen::Vector3d delta_p = Eigen::Vector3d::Zero();
	int Delnum = 0;
	if (p->p_belong_TetraGroup_ids.size() == 1) {
		return  groups[p->p_belong_TetraGroup_ids[0]]->Get_Exp_Pos(p->p_id);
	}
	//std::cout << "Delta" << std::endl;
	//std::cout << Delta << std::endl;
	//std::cout << "p_id is " << p->p_id << std::endl;
	//std::cout << "p's group is " << p->p_belong_TetraGroup_ids[0] << std::endl;
	while (p->p_id != Share_particle_id[Delnum]) {
		Delnum++;
	}
	//std::cout << "Delnum" << Delnum << std::endl;
	if (p->p_belong_TetraGroup_ids.size() > 1) {
		for (unsigned int i = 0; i < groups[p->p_belong_TetraGroup_ids[0]]->particle_num;i++) {
			//std::cout << "particle " << groups[p->p_belong_TetraGroup_ids[0]]->particles[i]->p_id<< std::endl;
			//std::cout << groups[p->p_belong_TetraGroup_ids[0]]->particles[i]->Get_Exp_Pos() << std::endl;
		}
		//std::cout << "particle EXP" << std::endl;
		//std::cout << groups[p->p_belong_TetraGroup_ids[0]]->Get_Exp_Pos(p->p_id) << std::endl;
		//std::cout << "particle Delta" << std::endl;
		//std::cout << Delta.block(3 * Delnum, 0, 3, 1) << std::endl;
		//std::cout << "particle " << std::endl;
		//std::cout << groups[p->p_belong_TetraGroup_ids[0]]->Get_Exp_Pos(p->p_id) + Delta.block(3 * Delnum, 0, 3, 1) << std::endl;
		for (unsigned int i = 0; i < groups[p->p_belong_TetraGroup_ids[1]]->particle_num; i++) {
			//std::cout << "particle " << groups[p->p_belong_TetraGroup_ids[1]]->particles[i]->p_id << std::endl;
			//std::cout << groups[p->p_belong_TetraGroup_ids[1]]->particles[i]->Get_Exp_Pos() << std::endl;
		}
		//std::cout << "particle EXP" << std::endl;
		//std::cout << groups[p->p_belong_TetraGroup_ids[1]]->Get_Exp_Pos(p->p_id) << std::endl;
		//std::cout << "particle Delta" << std::endl;
		//std::cout << Delta.block(3 * Delnum, 0, 3, 1) << std::endl;
		//std::cout << "particle " << std::endl;
		//std::cout << groups[p->p_belong_TetraGroup_ids[0]]->Get_Exp_Pos(p->p_id) + Delta.block(3 * Delnum, 0, 3, 1) << std::endl;
		return   groups[p->p_belong_TetraGroup_ids[0]]->Get_Exp_Pos(p->p_id) + Delta.block(3 * Delnum,0,3,1);
	}
	else{
		std::cout << "ERORR" << std::endl;
		return   groups[p->p_belong_TetraGroup_ids[0]]->Get_Exp_Pos(p->p_id);
	}
}
//グループごとの頂点に対して平均をとる(制約条件)
Eigen::Vector3d ObjectD::Calc_New_Delatax_Mean(ParticleD* p) {
	if (p->p_belong_TetraGroup_ids.size() == 0) { std::cout << "ERROR486" << std::endl; return Eigen::Vector3d::Zero(); }
	Eigen::Vector3d delta_p = Eigen::Vector3d::Zero();
	for (auto _g : p->p_belong_TetraGroup_ids) {
		TetraGroupD* tg = groups[_g];
		delta_p += tg->Get_Deltax_In_Group(p->p_id);
		// = tg->Get_Deltax_In_Group(p->p_id);

	}
	//std::cout << p->p_id << "is" << std::setprecision(10) << delta_p << std::endl;
	//所属するグループが1つ
	if (p->p_belong_TetraGroup_ids.size() == 1) {
		//std::cout << p->p_id <<"is"<<std::setprecision(10) << delta_p << std::endl;
		return delta_p;
	}
	//所属するグループが複数
	else {
		delta_p = delta_p / p->p_belong_TetraGroup_ids.size();
		//std::cout << p->p_id << "is" << std::setprecision(10) << delta_p <<std::endl;
		return delta_p;
	}
}

Eigen::Vector3d ObjectD::CalcMeanPos(ParticleD* p) {
	if (p->p_belong_TetraGroup_ids.size() == 0) { std::cout << "ERROR486" << std::endl; return Eigen::Vector3d::Zero(); }
	Eigen::Vector3d delta_p = Eigen::Vector3d::Zero();
	for (auto _g : p->p_belong_TetraGroup_ids) {
		TetraGroupD* tg = groups[_g];
		delta_p += tg->Get_DeltaxNew(p->p_id);
		// = tg->Get_Deltax_In_Group(p->p_id);

	}
	//std::cout << p->p_id << "is" << std::setprecision(10) << delta_p << std::endl;
	//所属するグループが1つ
	if (p->p_belong_TetraGroup_ids.size() == 1) {
		//std::cout << p->p_id <<"is"<<std::setprecision(10) << delta_p << std::endl;
		return delta_p;
	}
	//所属するグループが複数
	else {
		delta_p = delta_p / p->p_belong_TetraGroup_ids.size();
		//std::cout << p->p_id << "is" << std::setprecision(10) << delta_p <<std::endl;
		return delta_p;
	}
}
Eigen::Vector3d ObjectD::Calc_Mean_Vel(ParticleD* p) {
	if (p->p_belong_TetraGroup_ids.size() == 0) { std::cout << "ERROR486" << std::endl; return Eigen::Vector3d::Zero(); }
	Eigen::Vector3d delta_p = Eigen::Vector3d::Zero();
	for (auto _g : p->p_belong_TetraGroup_ids) {
		TetraGroupD* tg = groups[_g];
		delta_p += tg->Get_Velocity(p->p_id);
	}
	//std::cout << p->p_id << "is" << std::setprecision(10) << delta_p << std::endl;
	//所属するグループが1つ
	if (p->p_belong_TetraGroup_ids.size() == 1) {
		//std::cout << p->p_id <<"is"<<std::setprecision(10) << delta_p << std::endl;
		return delta_p;
	}
	//所属するグループが複数
	else {
		delta_p = delta_p / p->p_belong_TetraGroup_ids.size();
		//std::cout << p->p_id << "is" << std::setprecision(10) << delta_p <<std::endl;
		return delta_p;
	}
}
//TimeStepのはじめ毎に外力を0にする
void ObjectD::Timestep_Init() {
	Eigen::Vector3d zero_vec = Eigen::Vector3d::Zero();
	for (auto _p : particles) {
		_p->Set_Force(zero_vec);
	}
}

void ObjectD::Set_Force(Eigen::Vector3d grid) {
	//particles[particles.size() - 1]->Set_Force(grid - particles[particles.size() - 1]->Get_Grid());
	particles[particles.size() - 1]->Set_Force(grid);
	this->Outofforce = particles[particles.size() - 1]->Get_Force();
	MyDrawLine2(particles[particles.size() - 1]->Get_Grid(), grid);
}
void ObjectD::Volume_consevation(unsigned int loop) {
	//初期化
	for (auto _g : groups) {
		for (auto _e : _g->elements) {
			_e->Kappa = 0.0;
		}
	}
	for (unsigned int i = 0; i < loop; ++i) {
		for (auto _g:groups) {
			for (auto _e: _g->elements) {
				_e->Calc_Conservation(_g->particles, _g->m_In_Group);
			}
		}
	}
}

//===========================================================================//
//	@start		   				四面体要素作成							     //
//===========================================================================//

//===========================================================================//
//							  ドロネー三角形分割							 //
//===========================================================================//
void ObjectD::Delaunay_Triangulation() {
	//------------------------------------------------------------------//
	//オブジェクトを囲む巨大四辺形を作成
	//------------------------------------------------------------------//
	Eigen::Vector3d p[4];
	Eigen::Vector3d max_grid, min_grid;
	for (unsigned int i = 0; i < 3; i++) {
		max_grid[i] = FLT_MIN;
		min_grid[i] = FLT_MAX;
	}
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Eigen::Vector3d pos = (*it)->Get_IM_Grid();
		for (unsigned int i = 0; i < 3; i++) {
			if (max_grid[i] < pos[i]) { max_grid[i] = pos[i]; }
			if (min_grid[i] > pos[i]) { min_grid[i] = pos[i]; }
		}
	}
	Eigen::Vector3d center_grid = 0.5*(max_grid + min_grid);
	double dx = center_grid[0] - min_grid[0];
	double dy = center_grid[1] - min_grid[1];
	double dz = center_grid[2] - min_grid[2];
	double radius = sqrt(dx*dx + dy*dy + dz*dz) + 1.0;

	// 4つの頂点(すべての点を内包する球に外接する正四面体)
	p[0].x() = center_grid[0];
	p[0].y() = center_grid[1] + 3.0 * radius;
	p[0].z() = center_grid[2];

	p[1].x() = center_grid[0] - 2.0 * sqrt(2.0) * radius;
	p[1].y() = center_grid[1] - radius;
	p[1].z() = center_grid[2];

	p[2].x() = center_grid[0] + sqrt(2.0) * radius;
	p[2].y() = center_grid[1] - radius;
	p[2].z() = center_grid[2] + sqrt(6.0) * radius;

	p[3].x() = center_grid[0] + sqrt(2.0) * radius;
	p[3].y() = center_grid[1] - radius;
	p[3].z() = center_grid[2] - sqrt(6.0) * radius;

	std::vector<ParticleD*> huge_vertex;
	for (int i = 0; i < 4; i++) {
		ParticleD *temp_p = new ParticleD(p[i]);
		huge_vertex.push_back(temp_p);
	}
	for (auto _p : huge_vertex) {
		_p->Set_IM_Grid(2);
	}
	TetraElementD huge_tetra = TetraElementD(huge_vertex);
	//------------------------------------------------------------------//
	//------------------------------------------------------------------//

	std::set < TetraElementD > tetra_set;
	tetra_set.insert(huge_tetra);
	for (auto _p : particles) {
		Eigen::Vector3d pos = _p->Get_IM_Grid();
		//std::cout << pos << std::endl;
		std::map<TetraElementD, bool> tetra_map;

		for (auto tIter = tetra_set.begin(); tIter != tetra_set.end();) {
			//std::cout << tetra_set.size() << std::endl;
			TetraElementD tSet = *tIter;
			CircleD c = tSet.Get_Circum_Circle();
			//std::cout << c.center << std::endl;
			Eigen::Vector3d x = c.center - pos;
			long double len = x.norm();

			if (len < c.radius + 0.1) {//外接球の内部にパーティクルが存在する
				for (int i = 0; i < tSet.Get_Vertex_Num(); ++i) {
					std::vector<ParticleD*> vertex;
					vertex = tSet.Get_Particle();
					//sort(vertex.begin(), vertex.end());

					vertex[i] = _p;
					TetraElementD new_tetra = TetraElementD(vertex);

					std::map<TetraElementD, bool>::iterator it = tetra_map.find(new_tetra);
					if (it != tetra_map.end() && it->second) {
						tetra_map.erase(it);
						tetra_map.insert(std::map<TetraElementD, bool>::value_type(new_tetra, false));
					}
					else {
						tetra_map.insert(std::map<TetraElementD, bool>::value_type(new_tetra, true));
					}
				}

				tetra_set.erase(tIter++);
			}
			else ++tIter;
		}

		for (auto iter = tetra_map.begin(); iter != tetra_map.end(); ++iter) {
			if (iter->second) {
				tetra_set.insert(iter->first);
			}
		}
	}

	//最初に作ったオブジェクトを囲む巨大四面体に関係する四面体を削除
	for (auto tIter = tetra_set.begin(); tIter != tetra_set.end();) {
		if (huge_tetra.has_Common_Points(*tIter)) tetra_set.erase(tIter++);
		else ++tIter;
	}

	//巨大四面体作成に使ったParticleのメモリを開放
	for (auto _p : huge_vertex) {
		delete _p;
	}
	huge_vertex.clear();

	for (auto tetra : tetra_set) {
		std::vector<ParticleD*> temp_p = tetra.Get_Particle();

		TetraElementD* t = new TetraElementD(temp_p);
		tetras.push_back(t);
	}

	//四面体要素数を表示
	std::cout << "success create TetraSet. Tetra Set size : " << tetras.size() << std::endl;
	std::cout << std::endl;
	for (auto _t : tetras) {
		for (auto _p : _t->Get_Particle()) {
			std::cout << _p->p_id<<"particle, ";
		}
		std::cout << std::endl;
	}
}
void ObjectD::Triprism_Triangulation() {
	std::vector<ParticleD*> temp_p;
	TetraElementD* temp_t;
	for (int i = 0; i < Tri_prismnum; i++) {
		//1つ目の四面体要素
		temp_p.push_back(particles[0 + 3*i]);
		temp_p.push_back(particles[1 + 3*i]);
		temp_p.push_back(particles[2 + 3*i]);
		temp_p.push_back(particles[4 + 3*i]);
		temp_t = new TetraElementD(temp_p);
		tetras.push_back(temp_t);
		temp_p.clear();

		//2つ目の四面体要素
		temp_p.push_back(particles[2 + 3*i]);
		temp_p.push_back(particles[3 + 3*i]);
		temp_p.push_back(particles[4 + 3*i]);
		temp_p.push_back(particles[5 + 3*i]);
		temp_t = new TetraElementD(temp_p);
		tetras.push_back(temp_t);
		temp_p.clear();

		//3つ目の四面体要素
		temp_p.push_back(particles[0 + 3*i]);
		temp_p.push_back(particles[2 + 3*i]);
		temp_p.push_back(particles[3 + 3*i]);
		temp_p.push_back(particles[4 + 3*i]);
		temp_t = new TetraElementD(temp_p);
		tetras.push_back(temp_t);
		temp_p.clear();
	}
	//四面体要素数を表示
	std::cout << "success create TetraSet. Tetra Set size : " << tetras.size() << std::endl;
	std::cout << std::endl;
	for (auto _t : tetras) {
		for (auto _p : _t->Get_Particle()) {
			std::cout << _p->p_id<<"particle, ";
		}
		std::cout << std::endl;
	}
}
//===========================================================================//
//	@end		   				四面体要素作成							     //
//===========================================================================//

const std::string& ObjectD::Get_Name()const { //オブジェクトの名前を取得
	return data_name;
}
void ObjectD::CalcMassMatrix() {
	M_MatrixBody = Eigen::MatrixXd::Zero(3 * particles.size(), 3 * particles.size());
	for (auto _e : tetras) {
		for (auto p1_it = particles.begin(); p1_it != particles.end(); ++p1_it) {
			size_t p1_index = std::distance(particles.begin(), p1_it);
			for (auto p2_it = particles.begin(); p2_it != particles.end(); ++p2_it) {
				size_t p2_index = std::distance(particles.begin(), p2_it);

				M_MatrixBody.block(3 * p1_index, 3 * p2_index, 3, 3) =
					M_MatrixBody.block(3 * p1_index, 3 * p2_index, 3, 3) + _e->Get_M_Submatrix(*p1_it, *p2_it);
			}
		}
	}

}
void ObjectD::Initialization()
{
	InitialPos = Eigen::VectorXd::Zero(3 * particles.size());
	InitialVel = Eigen::VectorXd::Zero(3 * particles.size());
	for (int i = 0; i < particles.size(); i++)
	{
		InitialPos.segment(3 * i, 3) = particles[i]->Get_Initial_Pos();
		InitialVel.segment(3 * i, 3).setZero();
		//x_corrected.segment(3 * i, 3).setZero();
	}
}
void ObjectD::CalcPrePos() {
	
		//初期化
	f_Local = Eigen::VectorXd::Zero(3 * particles.size());
	x_Local = Eigen::VectorXd::Zero(3 * particles.size());
	v_Local = Eigen::VectorXd::Zero(3 * particles.size());
	for (unsigned int pi = 0; pi < particles.size(); pi++) {
		f_Local.block(3 * pi, 0, 3, 1) = Eigen::Vector3d(0.0, Gravity, 0.0);
		x_Local.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Grid();
		v_Local.block(3 * pi, 0, 3, 1) = particles[pi]->Get_Vel();

	}
	
	for (unsigned int pi = 0; pi < particles.size(); pi++) {
		auto aaa = M_MatrixBody(3 * pi, 3 * pi);
		v_Local.segment(3 * pi, 3) = v_Local.segment(3 * pi, 3) + f_Local.segment(3 * pi, 3) * TIME_STEP / M_MatrixBody(3 * pi, 3 * pi);
		x_Local.segment(3 * pi, 3) = x_Local.segment(3 * pi, 3) + v_Local.segment(3 * pi, 3) * TIME_STEP;
		particles[pi]->Set_Exp_Pos(x_Local.segment(3 * pi, 3));
	}

}

void ObjectD::Assemble_EnergyGradGlobal() {
	EnergyGradGlobal = Eigen::VectorXd::Zero(3* particles.size());
	EnergyGradMatrix = Eigen::MatrixXd::Zero(3 , particles.size());//组装好的能量梯度  
	//for (auto _e : tetras) {
	//	for (auto p_it = particles.begin(); p_it != particles.end(); ++p_it) {
	//		size_t p_index = std::distance(particles.begin(), p_it);

	//		// 从当前元素获取此粒子的能量梯度。
	//		Eigen::Vector3d local_gradient = _e->Get_SubEnergyGrad(*p_it);

	//		// 将局部梯度加到全局梯度的适当位置。
	//		temp.col(p_index) += local_gradient;

	//	}
	//}
	for (auto e : tetras) {
		for (int i = 0; i < 4; i++) {
			EnergyGradMatrix.col(e->Get_Particle()[i]->p_id) += e->EnergyGrad.block(0, i, 3, 1);	
		}
	}
	for (int i = 0; i < 3 * particles.size(); ++i) {
		int row = i % 3;
		int col = i / 3;
		EnergyGradGlobal(i) = EnergyGradMatrix(row, col);
	}
}
void ObjectD::CreateEnergyBody() {
	for (auto _e : tetras) {
		_e->CreateEnegyDensity();
		_e->CreatePotentialEnergy();
		EnergyBody += _e->PotentialEnergy;
	}
}
void ObjectD::CreateLagrangeMulti() {
	Bottom = 0;
	int c = 0;
	for (int i = 0; i < particles.size(); i++) {
		auto tmpa = (1 / M_MatrixBody(3 * i, 3 * i));
		auto tmpb = EnergyGradGlobal.segment(3 * i, 3).squaredNorm();
		Bottom = Bottom + (1 / M_MatrixBody(3 * i, 3 * i)) * EnergyGradGlobal.segment(3 * i, 3).squaredNorm();
		c += 1;
	}
	LagrangeMulti = (-1) * EnergyBody / Bottom;
}
void ObjectD::UpdatePos() {
	x_corrected = Eigen::VectorXd::Zero(3 * particles.size());
	Deltax = Eigen::VectorXd::Zero(3 * particles.size());

	for (int i = 0; i < particles.size(); i++) {
		Deltax.segment(3 * i, 3) = (1 / M_MatrixBody(3 * i, 3 * i)) * EnergyGradGlobal.segment(3 * i, 3) * LagrangeMulti;
		if (!(particles[i]->Is_Fixed()))
		{
			x_corrected.segment(3 * i, 3) = x_Local.segment(3 * i, 3)+Deltax.segment(3 * i, 3);

			particles[i]->Update_Grid(x_corrected.segment(3 * i, 3));
		}
		else {
			x_corrected.segment(3 * i, 3) = InitialPos.segment(3 * i, 3);
			//particles[i]->Update_Grid(x_corrected.segment(3 * i, 3));
		}
		
		//particles[i]->Update_Grid(x_corrected.segment<3>(i));
		//该写velocity了
	}
	
	//x_corrected = x_Local + Deltax;

}

void ObjectD::UpdateVel() {
	for (int i = 0; i < particles.size(); i++) {
		if (!(particles[i]->Is_Fixed()))
		{
			particles[i]->Set_Velocity_In_Model((x_corrected.segment(3 * i, 3) - particles[i]->Get_Grid()) / TIME_STEP);
		}
		else {
			particles[i]->Set_Velocity_In_Model(Eigen::Vector3d::Zero(3));
		}
		
	}
}
Eigen::Vector3d ObjectD::Get_Grid_In_Object(int pid) {
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	for (unsigned int pi = 0; pi < particles.size(); pi++) {
		if (particles[pi]->p_id == pid) {
			temp = x_corrected.segment(3 * pi, 3);
			return temp;
		}
	}
	std::cout << "ERROR1177 NOT FOUND" << std::endl;
	return temp;
}

void ObjectD::PBDCalculation() {
	//CalcMassMatrix();
	CalcPrePos();
	mtEnergyConstraint.startMyTimer();
	for (int i = 0; i < 1; i++) {
		for (auto _e : tetras) {
			_e->CreateDs();
			_e->CreateDefTensor();
			_e->CreateStrain();
			_e->CreateStress(data.young, data.poisson);
			_e->CreateEnegyDensity();
			_e->CreatePKFirstStress();
			_e->CreateEnergyGrad();
		}
		Assemble_EnergyGradGlobal();
		CreateEnergyBody();
		CreateLagrangeMulti();
		UpdatePos(); //e->Get_Particle()[0]->p_id
	}
	mtEnergyConstraint.endMyTimer();
	//std::cout << std::setprecision(4) << mtEnergyConstraint.getDt() << std::endl;
	
	UpdateVel();
	
		
}
void ObjectD::UpdateOldFEM() {
	for (auto g : groups) {
		g->OldFEM();
	}
}