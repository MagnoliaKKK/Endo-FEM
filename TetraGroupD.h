//-------------------------------------------------------------//
//@author	KatsuyaKikuchi
//@brief	四面体要素のグループ
//-------------------------------------------------------------//

#ifndef _TETRAGROUPD
#define _TETRAGROUPD

#include "TetraElementD.h"
#include <algorithm>

class TetraGroupD {
public:
	TetraGroupD(std::vector< TetraElementD* > elements, ObjectData data, double mass, int tetra_Group_id);
	~TetraGroupD();

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	void Init();
	void BuildMatrix();
	const std::vector< ParticleD* > particles;//groupが持っているparticle
	const unsigned int particle_num;	     //particle数
	void TetraGroupD::Set_Size_para(int particle_num);					 //各変数の要素数を決定する
	void TetraGroupD::Set_Size_para2(std::vector< ParticleD* > particles); //各変数の要素数を決定する

																		 //剛性行列を作成するための変数、関数
	ObjectData data;			      //材料パラメータ
	Eigen::MatrixXd R_Matrix;		  // グループの回転行列
	Eigen::MatrixXd trans_R_Matrix;	  // グループの回転行列の逆行列
	Eigen::MatrixXd M_Matrix;		  // グループの質量行列
	Eigen::MatrixXd M_Matrix_C;		  // グループの質量行列(モデル全体を考慮したもの)
	Eigen::MatrixXd inv_M_Matrix;	  // グループの質量行列の逆行列
	Eigen::MatrixXd stiffness_matrix; // グループの剛性行列
	Eigen::MatrixXd Damping_Matrix;	  // グループの減衰行列(*Timestep)
	Eigen::MatrixXd MassDamInv_Matrix;// グループの質量行列+減衰行列の逆行列
	Eigen::SparseMatrix<double> MassDamInvSparse;//(M+C')^-1 sparse
	Eigen::SparseMatrix<double> StiffnessSparse;//stiffness sparse 
	Eigen::SparseMatrix<double> DammingT_Matrix_Sparse; // 予測位置計算用行列((mass+damping)^-1 * mass * timestep)(Sparse)

	Eigen::SparseMatrix<double> Damm_Matrix_Sparse;//Jacobi計算用(mass+damping)
	Eigen::SparseMatrix<double> StiffnessTT_Matrix_Sparse;//Jacobi計算用(stiffness * Timestep * Timestep)
	Eigen::SparseMatrix<double> Rn_Matrix_Sparse;//Jacobi計算用(rotation_Matrix)
	Eigen::SparseMatrix<double> Rn_MatrixTR_Sparse;//Jacobi計算用(rotation_Matrix.T)
	Eigen::SparseMatrix<double> MassCondi_Sparse;//Jacobi計算用(I-(sumn/groupmass))

	Eigen::MatrixXd Diag_M_Matrix;	  // グループの対角質量行列
	Eigen::MatrixXd SUM_M_Matrix;	  // グループの帯状質量行列
	Eigen::SparseMatrix<double> SUM_M;//SUM_M sparse
	Eigen::MatrixXd Jacobi_Matrix; // 反復法の係数行列
	Eigen::MatrixXd Jacobi_Matrix_Inv;
	Eigen::SparseMatrix<double> Jacobi_Matrix_Sparse; // 反復法の係数行列(Sparse)
	Eigen::MatrixXd DiagFEM_Matrix_iteration; // 反復法のガウスザイデル用の対角化行列
	Eigen::MatrixXd F_FEM_Matrix_iteration;   // 反復法のガウスザイデル用の対角化じゃないところの行列
	Eigen::VectorXd Constant_term_iteration;  // 反復法のガウスザイデル用の定数
    
	Eigen::VectorXd PrimeVector;  // constanttermを計算するときに使うprimeのベクトル(予測位置ベクトル)
	Eigen::VectorXd OrigineVector;  // constanttermを計算するときに使うorigine_localのベクトル
	Eigen::VectorXd InitialVector;  // constanttermを計算するときに使う初期座標のベクトル(OldFEM用)
	Eigen::VectorXd GroupGridVector;  // グループにおける位置
	Eigen::VectorXd GroupVelVector;  // グループにおける速度

	Eigen::VectorXd iterativeVector; //反復法のGMRES用の初期値

	void Create_Stiffness_Matrix();
	void Create_M_Matrix();
	void Create_Diag_M_Matrix(); 
	void Create_SUM_M_Matrix();
	void Create_Center_Grid();
	void Create_Local_Stiffness_Matrix();
	void Create_Damping_Matrix();
	void Create_Information();

	double f_damping, v_damping;
	double mass;				// グループの質量(体積x密度)
	Eigen::VectorXd gravity;	// 重力加速度ベクトル(=m/s2, 3Nx1)
	Eigen::Vector3d origin_center_grid;//t=0でのグループの質量中心グローバルベクトル(重心ベクトル)
	Eigen::Vector3d center_grid;	   //t!=0でのグループの質量中心グローバルベクトル(重心ベクトル)
	std::vector<Eigen::Vector3d> origin_center_distance;//各particleにおける初めの重心からの距離(global)(変化しない)
	std::vector<Eigen::Vector3d> TetraGroupD::Get_origin_center_distance();//上の変数を取得
	std::vector<Eigen::Vector3d> center_distance;		//各particleにおける現在の重心からの距離(global)
	std::vector<Eigen::Vector3d> origin_local_grid;		//各particleのローカル座標初期位置
	std::vector<Eigen::Vector3d> TetraGroupD::Get_origin_local_grid();	  //上の変数を取得
	Eigen::VectorXd TetraGroupD::Get_bind_force();//拘束力を取得
	void TetraGroupD::Write_bind_force();//拘束力を取得

	//拘束力の変更
	void ReSet_Fbind_Pos(); //各節点における拘束力ベクトルを0にリセットする
	void Update_Fbind_Pos(); //各節点における拘束力ベクトルを更新する
	void Update_Fbind_Pos2(); //各節点における拘束力ベクトルを更新する(差を利用)
	void Update_Fbind_Pos3(); //各節点における拘束力ベクトルを更新する(差を利用.予測位置はグループで異なる)
	void Update_Fbind_Pos4();//予測位置の計算(対応する節点のみに力を加える)
	void Update_Fbind_Pos5();//予測位置が同じとき
	void Update_Fbind_Pos6();//予測位置が異なる(節点でも同じ計算)
	void Update_Fbind_Pos7();//予測位置が同じとき(節点は別の計算)
	void Update_Fbind_Pos8();
	Eigen::Vector3d Calc_Distance();
	double Add_convergence_iteration(double convite);//収束するかどうか値(前回と今回のDeltax)を計算する

	//FEMで計算するための変数
	std::vector<Eigen::Vector3d> prime_pos_g;	    //事前計算した予測位置(弾性力を考慮しない)
	std::vector<Eigen::Vector3d> initial_pos_g;		//事前計算した予測位置にFEMを適用しだした位置(1ループで同じ値)
	std::vector<Eigen::Vector3d> exp_pos_g;			//制約条件により更新した位置

	//idから各位置を取得する
	Eigen::Vector3d Get_Prime_Pos(int pid);
	Eigen::Vector3d Get_Initial_Pos(int pid);
	Eigen::Vector3d Get_Exp_Pos(int pid);
	Eigen::Vector3d Get_X_In_Group(int pid);
	Eigen::Vector3d Get_Deltax_In_Group(int pid);
	Eigen::Vector3d Get_DeltaxNew(int pid);
	Eigen::Vector3d Get_Velocity(int pid);
	Eigen::Vector3d Get_Exp_In_Group(int pid);
	Eigen::Vector3d Get_Grid_In_Group(int pid);
	Eigen::Vector3d Get_Vel_In_Group(int pid);

	//idからグループのインデックスを返す
	unsigned int Get_Group_Index(int pid);

	//idからグループだけで計算した質量を取得する
	//idから全体で考えた質量を取得する
	double Get_GMass_In_Group(int pid);
	double Get_Mass_In_Group(int pid);

	//idからexp_posを変数のベクトル値に変更する
	void Update_CoExp_Pos(int pid, Eigen::Vector3d a);

	//CRSで行列を扱うための変数、関数
	int tetra_group_id;											//tetraグループのid
	Eigen::MatrixXd inverted_stiffness_matrix;

	//rotated stiffmatrixes are matrixes from stiffmatrixes * rotate matrix, and been linearly stored, only the half triangle matrix is stored
	std::vector<std::vector<int>> stiffmatrix_valued_list;		//ある点に隣接する点のリスト
	std::vector<std::vector<int>> stiffmatrix_valued_list_sym;  //上のリストの対称となるリスト(対角は含まない)
	std::vector<int> valued_sym_rotated_indexes;				//particleの位置を記録し、剛性行列を計算する

	std::vector<Eigen::Matrix3d*> rotated_stiffmatrixs;			//隣接リストから求めた剛性行列のリスト

	void Find_Edges();					  //各particleに隣接する点を記録する
	void Create_Local_Stiffmatrix_Array();//行列を値が存在するところのみ格納する



    //シミュレーションを行うための変数、関数
	void Update_Rotate();				    //回転行列の更新
	void Create_Rotate_Matrix();		    //回転行列を極分解によって計算する
	void Create_Rotate_Matrix_APD();		//APDによって回転行列を極分解によって計算する
	
	Eigen::Vector3d clamp2(Eigen::Vector3d x, double y, double z);
	Eigen::Quaterniond Cay2(Eigen::Vector3d a);
	Eigen::Quaterniond Cay(Eigen::Vector3d a);
	Eigen::Quaterniond Exp(Eigen::Vector3d a);
	Eigen::Quaterniond Exp2(Eigen::Vector3d a);
	Eigen::Matrix3d rotate_matrix;	//回転行列3*3
	Eigen::MatrixXd rotate_matrix3N;//回転行列3n*3n
	Eigen::Quaterniond quaternion;			//回転行列のクォータニオン
	Eigen::Matrix3d rotate_matrix_trans;	//回転行列の逆行列
	std::vector<Eigen::Vector3d*> local_xs; //ローカル座標における変位ベクトル

	void Calc_Exp_Pos();//予測位置の計算
	void Calc_Exp_Pos2();//予測位置の計算(damperあり)
	void Calc_Exp_Pos3();//予測位置の計算(damperをいれる+グループごとの予測位置をもつ)
	void Calc_Exp_Pos_Group();//グループごとに位置を計算

	void Calc_CRSFEM();

	void Calc_FEM();

	void Calc_implicit_FEM_FF7();
	Eigen::VectorXd Calc_implicit_FEM_FF9();
	void Set_implicit_FEM_FF9(Eigen::VectorXd Delta, int num);
	void Set_implicit_FEM_FF12(Eigen::VectorXd Delta, int num, Eigen::MatrixXd Share,int AllParticlenum,std::vector<int> Sharepacle);//うしろはベクターほうがいいかも
	Eigen::MatrixXd Calc_Zidane_Tribal();
	Eigen::MatrixXd Calc_VIVI_Ornitier(int a);

	void Calc_iterative_FEM_Fbind();//反復法のFEM部分の計算
	void Calc_iterative_FEM_Fbind_Jacobi();//反復法のFEM部分の計算,ヤコビ反復
	void Calc_iterative_FEM_Fbind_Gauss();//反復法のFEM部分の計算, ガウスザイデル反復
	void Calc_iterative_FEM_Fbind_pivot();//反復法のFEM部分の計算,LU反復
	void Calc_iterative_FEM_Fbind_pivot2();//反復法のFEM部分の計算,LU反復(差分の値)
	void Calc_iterative_FEM_Fbind_Jacobi2(); //反復法のFEM部分の計算, ヤコビ反復(差分の値)

	void Calc_iterative_LocalFEM();	//debug用,反復法のFEM部分の計算(差分,LUとGMRES)
	void Calc_GMRES_Pre();//debug用(前処理)
	void Calc_GMRES_FEM();//debug用(前処理済)
	void CalcDeltax();

	void Calc_Jacobi_Matrix_iteration();//反復法、ヤコビで使う行列の更新(1ステップで一回)
	void Calc_Jacobi_Matrix_iteration_Sparse();//反復法、ヤコビで使う行列の更新(1ステップで一回)(Sparse)
	void Calc_Constant_term_iteration();//反復法、ヤコビで使う定数ベクトルの更新(1ステップで一回)
	void Calc_Constant_term_iteration2();//反復法、ヤコビで使う定数ベクトルの更新(1ステップで一回)(差分の値)
	void Calc_Constant_term_iteration_Sparse();//反復法、ヤコビで使う定数ベクトルの更新(1ステップで一回)(差分の値)(Sparse)
	void LHS();
	void RHS();
	void LHS0();
	void RHS0();
	void Calc_Jacobi_Matrix_iteration_Old();//OldFEMの係数行列作成
	void Calc_Constant_term_iteration_Old ();//OldFEMの定数値作成

	void Draw()const;

	double Get_Volume()const;
	void Set_Mass(double whole_mass);
	void Set_Gravity();

	void Set_Group_Mass(double a);
	double Get_Group_Mass();

	double Get_Volume();

	void Set_Deltax_In_Group(Eigen::VectorXd a);

	std::vector<ParticleD*> Get_Particle()const;


	const double  TIME_STEP2 = TIME_STEP * TIME_STEP;
	const double  TIME_STEP23 = TIME_STEP * TIME_STEP * f_damping;
	std::vector< TetraElementD* > elements;


	int key = 1;

	std::vector< ParticleD* > Create_Particles(std::vector< TetraElementD* > elements);

	//python用の出力
	void OP_python_M();
	void OP_python_MC();
	void OP_python_R(Eigen::MatrixXd rotate_matrix3N);
	void OP_python_SumM();
	void OP_python_I(Eigen::MatrixXd Ident);
	void OP_python_Stiff();
	void OP_python_Damping();
	void OP_python_Jacobi1(Eigen::MatrixXd rotate_matrix3N, Eigen::MatrixXd Ident);
	void OP_python_Jacobi();
	void OP_python_ConstantDelta();
	void OP_python_File();
	//その他のdebug用の関数
	void OP_eigenvalue(Eigen::MatrixXd Matrix);
	void OP_diag_advantage();
	void OP_diag_advantage2();
	void OP_Symetric(Eigen::MatrixXd rotate_matrix3N);
	void OP_CommutativeSR(Eigen::MatrixXd Matrix);
	void OP_OtherMatrix(Eigen::MatrixXd rotate_matrix3N);

	void OP_debug_iterative1(Eigen::VectorXd v);

	void Mode_Analytic();

	void Create_Rotate_Matrix_APD_Debug(Eigen::Matrix3d temp);//APDによって回転行列を極分解によって計算する(計測用)
	void InsertAPD(Eigen::Matrix3d p,Eigen::Matrix3d q, int c);//APDのMAXMINを代入(反復回数ごとに)(svdのR,APDのR,反復回数)
	Eigen::Vector3d axlAPD(Eigen::Matrix3d a);//行列aの軸ベクトルを計算する
	Eigen::VectorXd MaxRotationVector;//APDの精度を図るための反復回数ごとの値を格納する(最大値)
	Eigen::VectorXd MinRotationVector;//APDの精度を図るための反復回数ごとの値を格納する(最小値)
	int APDcount;
	double APDtime;//PDにかかる時間をはかる
	double APDSparsetime;//PDのsparse化にかかる時間をはかる
	double Apqtime;//Apqにかかる時間をはかる

	Eigen::VectorXd m_In_Group;//グループにおける質量ベクトル（グループによって異なる値をもつ）
	void Calc_Convergence();//拘束力によって位置が収束したかどうか確認する
	void Calc_Convergence2();//拘束力によって位置が収束したかどうか確認する

	Eigen::VectorXd Deltax_In_Group;//グループにおける解ベクトル（FEMで計算）（グループによって異なる値をもつ）
	Eigen::VectorXd Deltax_CoFEM;//グループにおける解ベクトル（ひずみの部分）（グループによって異なる値をもつ）
	Eigen::VectorXd Deltax_Bind;//グループにおける解ベクトル（拘束力の部分）（グループによって異なる値をもつ）
	Eigen::VectorXd bind_force_iterative;//各節点における拘束力ベクトル（変化する）（ループ0時点ではいつも0）
	Eigen::VectorXd DeltaxNew;
	//Debug用
	void Update_Rotate2();

private:
	//CRSを使わない実装の時の変数
	Eigen::VectorXd pLocal; //計算する時点でのparticleの位置(最初は弾性力を除いて更新した位置)
	Eigen::VectorXd f_inLocal;
	Eigen::VectorXd xgLocal;
	Eigen::VectorXd centerLocal;
	//事前計算のための変数

	Eigen::VectorXd f_Local;	//Nニュートン
	Eigen::VectorXd x_Local;
	Eigen::VectorXd v_Local;

	
	Eigen::VectorXd x_In_Group;//グループにおける予測位置ベクトル（FEMで計算）（グループによって異なる値をもつ）
	

	double Group_Mass;//グループの質量（節点の質量の和）

	double _mass_Radio;

	Eigen::GMRES< Eigen::SparseMatrix<double>> gmresFEM_Pre;
	Eigen::GMRES< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > gmresFEM_Pre2;
	Eigen::GMRES< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > gmresFEM_Pre3;
};

#endif