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
	mtUpRotate.setid(4), mtCEPos.setid(5), mtCFEM.setid(6), mtCconstr.setid(7), mtCP_1.setid(8), mtCP_2.setid(9), mtCP_3.setid(10);
}	//stopwatchのidをセット
ObjectD::~ObjectD() {
}
//===========================================================================//

void ObjectD::Draw()const {
	for (auto _g : groups) {
		_g->Draw();
	}
}
//反復法
//LU反復
void ObjectD::Solve_Constraints8(unsigned int loop) {
	ezolg = 0;
	convite = 0.0;
	mtCP_1.startMyTimer();
	//一つ前のexpを0(固定点以外)に初期化する
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Set_ExpAgo_Pos(Eigen::Vector3d::Zero());
		}
		if (fetestexcept(FE_INVALID)) {
			std::cout << "FE_INVALID PBD" << std::endl;
		}
		feclearexcept(FE_ALL_EXCEPT);
	}
	for (unsigned int i = 0; i < loop; ++i) {
		//std::cout << "To say Goodbye" << ezolg << std::endl;
		double delta = 0;
		mtCP_2.startMyTimer();
		//ひずみの部分の計算
		for (auto _g : groups) {
			//_g->Update_Rotate();			
			if (!useCRS) {
				_g->Calc_iterative_FEM_Fbind_pivot();
			}
			else {
				//まだできていない
				_g->Calc_CRSFEM();
			}
		}
		mtCP_2.endMyTimer();
		//位置合わせの部分
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				_p->Set_Exp_Pos(Calc_New_Exp_Pos_Mean(_p));
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		//拘束力の更新部分,
		for (auto _g : groups) {
			//EXPはグループで同じ値（afterなので位置合わせした後の値）
			//X_in_groupはグループで違う値のはず
			/*std::cout << "After shared Groval " << "groupnumbe is " << _g->tetra_group_id << std::endl;
			std::cout << std::setprecision(6) << _g->Get_Exp_Pos(3) << std::endl;
			std::cout << "After shared Local " << "groupnumbe is " << _g->tetra_group_id << std::endl;
			std::cout << std::setprecision(6) << _g->Get_X_In_Group(3) << std::endl;*/

			_g->Update_Fbind_Pos();
			//収束用の値を加算する
			//convite = _g->Add_convergence_iteration(convite);
			//particleで計算するようにしよ
		}
		//一つ前kのexp(expago)を現在k+1のexpに更新
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				convite = _p->Add_convergence_iteration(convite);
				_p->Set_ExpAgo_Pos(_p->Get_Exp_Pos());
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		std::cout << convite << ",";
		ezolg++;
	}
	//拘束力のリセット
	for (auto _g : groups) {
		_g->ReSet_Fbind_Pos();
		std::cout << "Did Bind Reset? "<<loop <<"times"<< std::endl;
		_g->Write_bind_force();
	}
	mtCP_1.endMyTimer();
	//制約条件により更新した位置を代入する
	mtCP_3.startMyTimer();
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Update(_p->Get_Exp_Pos());
		}
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Posi_set" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);
	mtCP_3.endMyTimer();
}
//反復法
//LU反復
//変位で計算
void ObjectD::Solve_Constraints9(unsigned int loop) {
	ezolg = 0;
	convite = 0.0;
	//std::cout << "np.array([";
	mtCP_1.startMyTimer();
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Set_ExpAgo_Pos(_p->Get_Exp_Pos());
		}
		if (fetestexcept(FE_INVALID)) {
			std::cout << "FE_INVALID PBD" << std::endl;
		}
		feclearexcept(FE_ALL_EXCEPT);
	}
	for (unsigned int i = 0; i < loop; ++i) {
		//std::cout << "To say Goodbye" << ezolg << std::endl;
		double delta = 0;
		convite = 0.0;
		mtCP_2.startMyTimer();
		//ひずみの部分の計算
		for (auto _g : groups) {
			//_g->Update_Rotate();			
			if (!useCRS) {
				_g->Calc_iterative_FEM_Fbind_pivot2();
			}
			else {
				//まだできていない
				//_g->Calc_CRSFEM();

				//ヤコビの実験に使用
				_g->Calc_iterative_FEM_Fbind_Jacobi2();
				//係数行列の性質でjacobi時間ステップをうまく調整しないとできない
			}
		}
		mtCP_2.endMyTimer();
		//位置合わせの部分
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				_p->Set_Exp_Pos(Calc_New_Exp_Pos_Mean(_p));
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		//拘束力の更新部分
		for (auto _g : groups) {
			//EXPはグループで同じ値（afterなので位置合わせした後の値）
			//X_in_groupはグループで違う値のはず
			/*std::cout << "After shared Groval " << "groupnumbe is " << _g->tetra_group_id << std::endl;
			std::cout << std::setprecision(6) << _g->Get_Exp_Pos(3) << std::endl;
			std::cout << "After shared Local " << "groupnumbe is " << _g->tetra_group_id << std::endl;
			std::cout << std::setprecision(6) << _g->Get_X_In_Group(3) << std::endl;*/

			_g->Update_Fbind_Pos();
		}
		//一つ前kのexp(expago)を現在k+1のexpに更新
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				convite = _p->Add_convergence_iteration(convite);
				_p->Set_ExpAgo_Pos(_p->Get_Exp_Pos());
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		//std::cout << std::setprecision(5)<< convite << ",";
		ezolg++;
	}
	//std::cout<<std::endl;
	//拘束力のリセット
	for (auto _g : groups) {
		_g->ReSet_Fbind_Pos();
		//拘束力をリセットできたか出力
		//std::cout << "Did Bind Reset? " << loop << "times" << std::endl;
		//_g->Write_bind_force();
	}
	mtCP_1.endMyTimer();
	//制約条件により更新した位置を代入する
	mtCP_3.startMyTimer();
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Update(_p->Get_Exp_Pos());
		}
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Posi_set" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);
	mtCP_3.endMyTimer();
}
//debug用
void ObjectD::Solve_Constraints10(unsigned int loop) {
	ezolg = 0;
	convite = 0.0;
	MicroSecondTimer mtOldIni;
	MicroSecondTimer mtShare;
	MicroSecondTimer mtUpBind;
	MicroSecondTimer mtUpAgo;
	MicroSecondTimer mtReBind;
	MicroSecondTimer mtUpPos;
	MicroSecondTimer mtGMRES;
	MicroSecondTimer mtLOOP;
	mtOldIni.setid(12);
	mtShare.setid(13);
	mtUpBind.setid(14);
	mtUpAgo.setid(15);
	mtReBind.setid(16);
	mtUpPos.setid(17);
	mtGMRES.setid(18);
	mtLOOP.setid(19);
	//std::cout << std::endl;
	//std::cout << "Rerate" << std::endl;
	mtCP_1.startMyTimer();
	//ノルムを計算するための値(1つ前の解)の初期化
	/*
	mtOldIni.startMyTimer();
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Set_ExpAgo_Pos(_p->Get_Prime_Pos());
			_p->Set_DeltaxAgo_In_Model(Eigen::Vector3d::Zero());
		}
		if (fetestexcept(FE_INVALID)) {
			std::cout << "FE_INVALID PBD 238" << std::endl;
		}
		feclearexcept(FE_ALL_EXCEPT);
	}
	mtOldIni.endMyTimer();
	*/
	//std::cout << "Initial Ago is " << mtOldIni.getDt() << std::endl;
	//反復法により解を求める
	mtLOOP.startMyTimer();
	for (unsigned int i = 0; i < loop; ++i) {
		convite = 0.0;
		//弾性力の部分の線形方程式を計算
		mtGMRES.startMyTimer();
		for (auto _g : groups) {	
			mtCP_2.startMyTimer();
			_g->Calc_iterative_LocalFEM();
			mtCP_2.endMyTimer();
		}
		mtGMRES.endMyTimer();
		//std::cout << "Calc GMRES is " << mtGMRES.getDt() << std::endl;
		//位置合わせの部分(Delta)
		mtShare.startMyTimer();
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				_p->Set_Deltax_In_Model(Calc_New_Delatax_Mean(_p));
			}
			else {
				_p->Set_Deltax_In_Model(Eigen::Vector3d::Zero());
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD 263" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		for (auto _p : particles) {
			//std::cout << _p->p_id<<" Deltax_In_Model"<< std::endl;
			//std::cout << _p->Get_Deltax_In_Model() << std::endl;
		}
		mtShare.endMyTimer();
		//std::cout << "Sahre potitionis " << mtShare.getDt() << std::endl;
		//拘束力の更新部分
		mtUpBind.startMyTimer();
		for (auto _g : groups) {
			_g->Update_Fbind_Pos2();
		}
		mtUpBind.endMyTimer();
		//std::cout << "Upadate Bind is " << mtUpBind.getDt() << std::endl;
		//一つ前kのDeltaxModel(expago)を現在k+1のDeltamodelに更新
		/*
		mtUpAgo.startMyTimer();
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				convite = _p->Add_convergence_iteration2(convite);
				_p->Set_DeltaxAgo_In_Model(_p->Get_Deltax_In_Model());
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD 284" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		*/
		//std::cout << "Update Ago Pos is " << mtUpAgo.getDt() << std::endl;
		mtUpAgo.endMyTimer();
		//一回の更新でどれだけ変化したか出力する
		//std::cout << std::setprecision(5)<< convite << ",";
		//ConbiteGMRES[i] += convite;
		ezolg++;
	}
	mtLOOP.endMyTimer();
	//std::cout << "Loop is " << mtLOOP.getDt() << std::endl;
	/*
	GMREScount++;
	if (GMREScount == 2500) {
		ConbiteGMRES = ConbiteGMRES / 2500;
		std::cout << "Conbite is " << std::setprecision(10) << ConbiteGMRES << std::endl;
	}
	*/
	mtReBind.startMyTimer();
	//拘束力のリセット
	/*
	for (auto _g : groups) {
		_g->ReSet_Fbind_Pos();
	}
	*/
	//::cout << "Reset Bind is " << mtReBind.getDt() << std::endl;
	mtReBind.endMyTimer();
	mtCP_1.endMyTimer();
	//制約条件により更新した位置を代入する
	mtUpPos.startMyTimer();
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Update(_p->Get_Deltax_In_Model() + _p->Get_Prime_Pos());
		}
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Posi_set" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);
	mtUpPos.endMyTimer();
	//std::cout << "Update Pos is " << mtUpPos.getDt() << std::endl;
}
//debug用
void ObjectD::Solve_Constraints10_LU(unsigned int loop) {
	ezolg = 0;
	convite = 0.0;
	MicroSecondTimer mtOldIni;
	MicroSecondTimer mtShare;
	MicroSecondTimer mtUpBind;
	MicroSecondTimer mtUpAgo;
	MicroSecondTimer mtReBind;
	MicroSecondTimer mtUpPos;
	MicroSecondTimer mtGMRES;
	MicroSecondTimer mtLOOP;
	mtOldIni.setid(12);
	mtShare.setid(13);
	mtUpBind.setid(14);
	mtUpAgo.setid(15);
	mtReBind.setid(16);
	mtUpPos.setid(17);
	mtGMRES.setid(18);
	mtLOOP.setid(19);
	//反復法により解を求める
	mtLOOP.startMyTimer();
	Eigen::MatrixXd SuperA = Eigen::MatrixXd::Zero(3 * AllParticlenum + 3 * Share_particlenum, 3 * AllParticlenum);
	unsigned int tempposu;
	tempposu = 0;
	for (unsigned int gi = 0; gi < groups.size(); gi++) {
		SuperA.block(tempposu,tempposu, 3 * groups[gi]->particle_num, 3 * groups[gi]->particle_num) = groups[gi]->Jacobi_Matrix;
		tempposu += 3 * groups[gi]->particle_num;
	}
	//std::cout << "SuperA" << std::endl;
	//std::cout << SuperA << std::endl;
	if (groups.size() > 1) {
		SuperA.block(tempposu, 0, 3 * Share_particlenum, 3 * AllParticlenum) = Node_Sharing_Matrix;
	}
	Eigen::MatrixXd HyperA = Eigen::MatrixXd::Zero(3 * AllParticlenum , 3 * AllParticlenum);
	HyperA = SuperA.transpose() * SuperA;
	Eigen::SparseMatrix<double> HyperA_Sparse;
	HyperA_Sparse = HyperA.sparseView();
	Eigen::VectorXd Superb = Eigen::VectorXd::Zero(3 * AllParticlenum + 3 * Share_particlenum);
	tempposu = 0;
	for (unsigned int gi = 0; gi < groups.size(); gi++) {
		Superb.block(tempposu, 0, 3 * groups[gi]->particle_num, 1) = groups[gi]->Constant_term_iteration;
		tempposu += 3 * groups[gi]->particle_num;
	}
	Superb = SuperA.transpose() * Superb;
	Eigen::VectorXd vector_u = Eigen::VectorXd::Zero(3 * AllParticlenum);
	//FEM計算
	//GMRESなどの反復法を用いて反復法(Local)を解く
	if (useGMRES) {
		//GMRES

		//前処理なし
		if (!usePreIte) {//执行这里
			Eigen::GMRES< Eigen::SparseMatrix<double>> gmresFEM;
			gmresFEM.setMaxIterations(outerGMRES);//外部反復の設定
			gmresFEM.set_restart(innerGMRES);//内部反復の設定
			gmresFEM.compute(HyperA_Sparse);// Compute 

			//初期値0
			//vector_u = gmresFEM.solve(Constant_term_iteration + bind_force_iterative);

			//初期値は一つ前の値から計算したものをいれる
			vector_u = gmresFEM.solve(Superb);
			//std::cout << "GMRES前なし：#iterations：" << gmresFEM.iterations() << "、推定エラー：" << gmresFEM.error() << std::endl;
		}
		else {
			//前処理あり

			Eigen::GMRES< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > gmresFEM;
			gmresFEM.preconditioner().setFillfactor(7);  //Get the reference of the preconditioner and set properties
			//gmresFEM.setTolerance(10e-3);//許容値の設定
			gmresFEM.setMaxIterations(outerGMRES);//外部反復の設定
			gmresFEM.set_restart(innerGMRES);//内部反復の設定
			gmresFEM.compute(HyperA_Sparse);// Compute the ILUT factorization

			//初期値0
			//vector_u = gmresFEM2.solve(Constant_term_iteration + bind_force_iterative);

			//初期値は一つ前の値から計算したものをいれる
			vector_u = gmresFEM.solve(Superb);
			//std::cout << "GMRES前あり：#iterations：" << gmresFEM2.iterations() << "、推定エラー：" << gmresFEM2.error() << std::endl;
		}
	}
	//LU分解を用いて反復法(Local)を解く
	else {
		Eigen::FullPivLU<Eigen::MatrixXd> lu(HyperA);
		vector_u = lu.solve(Superb);
	}
	//Groupの解に代入
	tempposu = 0;
	for (unsigned int gi = 0; gi < groups.size(); gi++) {
		groups[gi]->Set_Deltax_In_Group(vector_u.block(tempposu, 0, 3 * groups[gi]->particle_num, 1));
		tempposu += 3 * groups[gi]->particle_num;
	}
	//位置修正
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Set_Deltax_In_Model(Calc_New_Delatax_Mean(_p));
		}
		else {
			_p->Set_Deltax_In_Model(Eigen::Vector3d::Zero());
		}
		if (fetestexcept(FE_INVALID)) {
			std::cout << "FE_INVALID PBD 263" << std::endl;
		}
		feclearexcept(FE_ALL_EXCEPT);
	}
	//制約条件により更新した位置を代入する
	mtUpPos.startMyTimer();
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Update(_p->Get_Deltax_In_Model() + _p->Get_Prime_Pos());
		}
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Posi_set" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);
	mtUpPos.endMyTimer();
	//std::cout << "Update Pos is " << mtUpPos.getDt() << std::endl;
}
//debug用
//グループごとに予測位置が違う
void ObjectD::Solve_Constraints11(unsigned int loop) {
	ezolg = 0;
	convite = 0.0;
	MicroSecondTimer mtOldIni;
	MicroSecondTimer mtShare;
	MicroSecondTimer mtUpBind;
	MicroSecondTimer mtUpAgo;
	MicroSecondTimer mtReBind;
	MicroSecondTimer mtUpPos;
	MicroSecondTimer mtGMRES;
	MicroSecondTimer mtLOOP;
	mtOldIni.setid(12);
	mtShare.setid(13);
	mtUpBind.setid(14);
	mtUpAgo.setid(15);
	mtReBind.setid(16);
	mtUpPos.setid(17);
	mtGMRES.setid(18);
	mtLOOP.setid(19);
	//std::cout << std::endl;
	//std::cout << "Rerate" << std::endl;
	mtCP_1.startMyTimer();
	//ノルムを計算するための値(1つ前の解)の初期化
	/*
	mtOldIni.startMyTimer();
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Set_ExpAgo_Pos(_p->Get_Prime_Pos());
			_p->Set_DeltaxAgo_In_Model(Eigen::Vector3d::Zero());
		}
		if (fetestexcept(FE_INVALID)) {
			std::cout << "FE_INVALID PBD 238" << std::endl;
		}
		feclearexcept(FE_ALL_EXCEPT);
	}
	mtOldIni.endMyTimer();
	*/
	//std::cout << "Initial Ago is " << mtOldIni.getDt() << std::endl;
	//反復法により解を求める
	mtLOOP.startMyTimer();
	for (unsigned int i = 0; i < loop; ++i) {
		convite = 0.0;
		//弾性力の部分の線形方程式を計算
		mtGMRES.startMyTimer();
		for (auto _g : groups) {
			mtCP_2.startMyTimer();
			//_g->Calc_iterative_LocalFEM();
			_g->Calc_GMRES_FEM();
			mtCP_2.endMyTimer();
		}
		mtGMRES.endMyTimer();
		//std::cout << "Calc GMRES is " << mtGMRES.getDt() << std::endl;
		//位置合わせの部分(Delta)
		mtShare.startMyTimer();
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				_p->Set_Deltax_In_Model(Calc_New_Delatax_Mean(_p));
			}
			else {
				_p->Set_Deltax_In_Model(Eigen::Vector3d::Zero());
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD 263" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		for (auto _p : particles) {
			//std::cout << _p->p_id<<" Deltax_In_Model"<< std::endl;
			//std::cout << _p->Get_Deltax_In_Model() << std::endl;
		}
		//Deltaを含んだやつ
		
		std::cout<<"Delta"<<std::endl;
		for (auto _g : groups) {
			_g->Calc_Convergence2();
		}
		
		mtShare.endMyTimer();
		//std::cout << "Sahre potitionis " << mtShare.getDt() << std::endl;
		//拘束力の更新部分
		mtUpBind.startMyTimer();
		for (auto _g : groups) {
			_g->Update_Fbind_Pos5();
		}
		mtUpBind.endMyTimer();
		//std::cout << "Upadate Bind is " << mtUpBind.getDt() << std::endl;
		//一つ前kのDeltaxModel(expago)を現在k+1のDeltamodelに更新
		/*
		mtUpAgo.startMyTimer();
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				convite = _p->Add_convergence_iteration2(convite);
				_p->Set_DeltaxAgo_In_Model(_p->Get_Deltax_In_Model());
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD 284" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		*/
		//std::cout << "Update Ago Pos is " << mtUpAgo.getDt() << std::endl;
		mtUpAgo.endMyTimer();
		//一回の更新でどれだけ変化したか出力する
		//std::cout << std::setprecision(5)<< convite << ",";
		//ConbiteGMRES[i] += convite;
		ezolg++;
	}
	mtLOOP.endMyTimer();
	//std::cout << "Loop is " << mtLOOP.getDt() << std::endl;
	/*
	GMREScount++;
	if (GMREScount == 2500) {
		ConbiteGMRES = ConbiteGMRES / 2500;
		std::cout << "Conbite is " << std::setprecision(10) << ConbiteGMRES << std::endl;
	}
	*/
	mtReBind.startMyTimer();
	//拘束力のリセット
	for (auto _g : groups) {
		_g->ReSet_Fbind_Pos();
	}
	//::cout << "Reset Bind is " << mtReBind.getDt() << std::endl;
	mtReBind.endMyTimer();
	mtCP_1.endMyTimer();
	//制約条件により更新した位置を代入する
	mtUpPos.startMyTimer();
	//Delta
	/*
	std::cout << "Delta" << std::endl;
	for (auto _g : groups) {
		_g->Calc_Convergence2();
	}
	*/
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Update(_p->Get_Deltax_In_Model() + _p->Get_Exp_Pos());
			//_p->Set_Exp_Pos(_p->Get_Deltax_In_Model() + _p->Get_Exp_Pos());
		}
		else {
			_p->Set_Exp_Pos(_p->Get_Grid());
		}
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Posi_set" << std::endl;
	}
	//体積保存制約
	//Volume_consevation(2);
	/*
	for (auto _p : particles) {
		if (!(_p->Is_Fixed())) {
			_p->Update(_p->Get_Exp_Pos());
		}
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Posi_set" << std::endl;
	}
	feclearexcept(FE_ALL_EXCEPT);
	*/
	mtUpPos.endMyTimer();
	//std::cout << "Update Pos is " << mtUpPos.getDt() << std::endl;
}
//debug用
//グループごとに予測位置が違う
void ObjectD::Solve_Constraints12(unsigned int loop) {
	ezolg = 0;
	convite = 0.0;
	MicroSecondTimer mtOldIni;
	MicroSecondTimer mtShare;
	MicroSecondTimer mtUpBind;
	MicroSecondTimer mtUpAgo;
	MicroSecondTimer mtReBind;
	MicroSecondTimer mtUpPos;
	MicroSecondTimer mtGMRES;
	MicroSecondTimer mtLOOP;
	mtOldIni.setid(12);
	mtShare.setid(13);
	mtUpBind.setid(14);
	mtUpAgo.setid(15);
	mtReBind.setid(16);
	mtUpPos.setid(17);
	mtGMRES.setid(18);
	mtLOOP.setid(19);
	//std::cout << std::endl;
	//std::cout << "Rerate" << std::endl;
	mtCP_1.startMyTimer();
	//std::cout << "Initial Ago is " << mtOldIni.getDt() << std::endl;
	//反復法により解を求める
	mtLOOP.startMyTimer();
	for (unsigned int i = 0; i < loop; ++i) {
		convite = 0.0;
		//弾性力の部分の線形方程式を計算
		mtGMRES.startMyTimer();
		for (auto _g : groups) {
			mtCP_2.startMyTimer();
			_g->Calc_GMRES_FEM();
			mtCP_2.endMyTimer();
		}
		mtGMRES.endMyTimer();
		//std::cout<<"663"<<std::endl;
		//std::cout << "Calc GMRES is " << mtGMRES.getDt() << std::endl;
		//位置合わせの部分(Delta)
		mtShare.startMyTimer();
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				_p->Set_Deltax_In_Model(Calc_New_Delatax_Mean(_p));
			}
			else {
				_p->Set_Deltax_In_Model(Calc_New_Delatax_Mean(_p));
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD 263" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		//Deltaを含んだやつ
		/*
		std::cout << "Delta" << std::endl;
		for (auto _g : groups) {
			_g->Calc_Convergence2();
		}
		*/
		mtShare.endMyTimer();
		//std::cout << "Sahre potitionis " << mtShare.getDt() << std::endl;
		//拘束力の更新部分
		mtUpBind.startMyTimer();
		for (auto _g : groups) {
			//Eigen::VectorXd GMRES_Bind = Eigen::VectorXd(3 * _g->particle_num);
			//GMRES_Bind = _g->bind_force_iterative;
			_g->Update_Fbind_Pos6();
			//GMRES_Bind = GMRES_Bind - _g->bind_force_iterative;
			//ConbiteGMRES[i] += GMRES_Bind.squaredNorm();
			
			//if (GMREScount == 90 && _g->tetra_group_id==0) {
			/*
			if (i== 19 && _g->tetra_group_id == 0) {
				//GMRES_Bind = _g->bind_force_iterative;
				//Eigen::Vector3d Distance = Eigen::Vector3d::Zero();
				ConbiteGMRES[GMREScount] = (_g->Calc_Distance()).squaredNorm();
				std::cout << std::setprecision(10) << ConbiteGMRES[GMREScount] << ",";
				//ConbiteGMRES[i] = (_g->Calc_Distance()).squaredNorm();
				//std::cout << std::setprecision(10) << ConbiteGMRES[i] << std::endl;
				//ConbiteGMRES[i] = (GMRES_Bind.squaredNorm())/10000;
			}
			*/
		}
		mtUpBind.endMyTimer();
		//std::cout << "Update Ago Pos is " << mtUpAgo.getDt() << std::endl;
		mtUpAgo.endMyTimer();
		//一回の更新でどれだけ変化したか出力する
		//std::cout << std::setprecision(5)<< convite << ",";
		//ConbiteGMRES[i] += convite;
		ezolg++;
	}
	mtLOOP.endMyTimer();
	//std::cout << "Loop is " << mtLOOP.getDt() << std::endl;
	
	GMREScount++;
	/*
	if (GMREScount == 101) {
		std::cout << "Conbite is " << std::setprecision(10) << ConbiteGMRES << std::endl;
	}
	*/
	mtReBind.startMyTimer();
	//拘束力のリセット
	for (auto _g : groups) {
		_g->ReSet_Fbind_Pos();
	}
	//::cout << "Reset Bind is " << mtReBind.getDt() << std::endl;
	mtReBind.endMyTimer();
	mtCP_1.endMyTimer();
	//制約条件により更新した位置を代入する
	mtUpPos.startMyTimer();
	//Delta
	/*
	std::cout << "Delta" << std::endl;
	for (auto _g : groups) {
		_g->Calc_Convergence2();
	}
	*/
	for (auto _g: groups) {
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			//速度をいれてみた
			_g->GroupVelVector.block(3 * pi, 0, 3, 1) = (_g->PrimeVector.block(3 * pi, 0, 3, 1) + _g->Deltax_In_Group.block(3 * pi, 0, 3, 1) - _g->GroupGridVector.block(3 * pi, 0, 3, 1)) / TIME_STEP;
			_g->GroupGridVector.block(3 * pi, 0, 3, 1) = _g->PrimeVector.block(3 * pi, 0, 3, 1) + _g->Deltax_In_Group.block(3 * pi, 0, 3, 1);
			if (_g->particles[pi]->Is_Fixed()) {
				_g->GroupVelVector.block(3 * pi, 0, 3, 1) = Eigen::Vector3d::Zero();
				_g->GroupGridVector.block(3 * pi, 0, 3, 1) = _g->InitialVector.block(3 * pi, 0, 3, 1);
			}
		}
		std::ofstream outputfile("Deltax_In_Group.txt", std::ios_base::app);
		outputfile << "Deltax_In_Group" << _g->tetra_group_id << " is " << std::endl;
		outputfile << std::setprecision(3) << _g->Deltax_In_Group << std::endl;
		outputfile.close();
	}
	if (fetestexcept(FE_INVALID)) {
		std::cout << "FE_INVALID Posi_set" << std::endl;
	}
	mtUpPos.endMyTimer();
	//std::cout << "Update Pos is " << mtUpPos.getDt() << std::endl;
}

void ObjectD::Solve_Constraints13(unsigned int loop) {
	for (auto _g : groups)
	{
		_g->ReSet_Fbind_Pos();
		
	}

	for (int i = 0; i < loop; i++)
	{
		for (auto _g : groups) {
			_g->LHS0();
			_g->RHS0();

			//Eigen::GMRES<Eigen::SparseMatrix<double>> solver;
			////solver.preconditioner().setF
			//solver.setMaxIterations(outerGMRES);
			//solver.set_restart(innerGMRES);
			////solver.setTolerance(1e-10);
			//solver.compute(_g->Jacobi_Matrix_Sparse);
			//_g->DeltaxNew = solver.solve(_g->Constant_term_iteration);
			

			_g->CalcDeltax();
			//std::ofstream outputfile("deltaxx.txt", std::ios_base::app);
		    //outputfile << "deltax" << _g->tetra_group_id << " is " << std::endl;
			//outputfile << std::setprecision(3) << _g->DeltaxNew << std::endl;
			//outputfile.close();
			

			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID Deltax" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		for (auto _p : particles) {
			if (!(_p->Is_Fixed())) {
				_p->Set_Deltax_In_Model(CalcMeanPos(_p));
			}
			else {
				_p->Set_Deltax_In_Model(CalcMeanPos(_p));
			}
			if (fetestexcept(FE_INVALID)) {
				std::cout << "FE_INVALID PBD 263" << std::endl;
			}
			feclearexcept(FE_ALL_EXCEPT);
		}
		for (auto _g : groups) {
			//Eigen::VectorXd GMRES_Bind = Eigen::VectorXd(3 * _g->particle_num);
			//GMRES_Bind = _g->bind_force_iterative;
			_g->Update_Fbind_Pos8();
		}
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