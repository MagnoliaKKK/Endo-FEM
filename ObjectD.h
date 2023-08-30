//===========================================================================//
//@author KatsuyaKikuchi
//@brief シミュレーションを行うモデルを扱う仮想クラス
//===========================================================================//
#ifndef _OBJECTD
#define _OBJECTD

#include "TetraGroupD.h"
#include "TetraElementD.h"

class ObjectD {
public:
	ObjectD(std::vector<ParticleD*> p, ObjectData data);
	virtual ~ObjectD();

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	virtual void Update() = 0;	//位置更新
	virtual void Draw()const;	//描画

	void Timestep_Init();	//各タイムステップで外力を0にする
	double Get_V();			//オブジェクトの体積を取得
	double Get_M();
	int ezolg;				//制約条件による更新回数を記録する
	int Get_ezolg();		//制約条件による更新回数を取得
	double convite;			//反復においてどれだけ値が変化しているかを記録する(前回と今の佐野ノルムの2乗)
	const std::string& Get_Name()const;//オブジェクトの名前を取得

	void Set_Force(Eigen::Vector3d grid);//最後の頂点(一番右下の頂点)にマウスのポインタ分だけ力をかける
	Eigen::Vector3d Outofforce;
	MicroSecondTimer mtUpRotate; // 回転の計算時間
	MicroSecondTimer mtCEPos;	 // 弾性力以外の力による位置更新の計算時間
	MicroSecondTimer mtCFEM;	 // 一回目の有限要素法による位置更新の計算時間
	MicroSecondTimer mtCconstr;	 // 制約条件による位置の修正にかかる時間
	MicroSecondTimer mtCP_1;     //省略法の反復計算にかかる時間
	MicroSecondTimer mtCP_2;     //省略法の局所剛性行列のFEMによる制約での計算時間
	MicroSecondTimer mtCP_3;     //省略法の位置の更新にかかる時間

	std::vector<TetraGroupD*> groups;
	std::vector<ParticleD*> particles; //頂点粒子
	std::vector<TetraElementD*> tetras;//四面体要素
									  //std::map < Particle*, std::vector<TetraGroup*> > belong_group;
	ObjectData data;				  //材料パラメータ

	int Sum_particlenum;                //グループの節点数の合計
	Eigen::MatrixXd Node_Sharing_Matrix;//Node_Sharing_Matrixの作成
	int Share_particlenum;                //共有する節点の組み合わせの総和=l
	std::vector<int>Share_particle_id;  //共有する節点の各座標
	int AllParticlenum;					//節点の総和=l
	Eigen::MatrixXd M_MatrixBody;		//質量行列
	//Debiug用
	Eigen::VectorXd ConbiteGMRES;//GMRESの収束をみるために20反復を値として格納
	int GMREScount;//出力するためにSTep数を記録
	
	Eigen::VectorXd f_Local;	//Nニュートン
	Eigen::VectorXd x_Local;
	Eigen::VectorXd v_Local;

protected:

	virtual void Init() = 0;			//初期化(オーバーライドする)
	virtual void Create_Groups() = 0;	//グループの生成(オーバーライドする)
	void Delaunay_Triangulation();		//ドロネー三角形分割
	void Triprism_Triangulation();		//三角柱を3個に分割
	std::string data_name;				//オブジェクトの名前

	void Solve_Constraints(unsigned int loop);	  //制約条件による更新スキーム
	void Solve_Constraints2(unsigned int loop);	  //制約条件による更新スキーム(位置制約を行列の中にいれる)
	void Solve_Constraints3(unsigned int loop);	  //制約条件による更新スキーム(位置制約を行列の中にいれる,正方)
	void Solve_Constraints4(unsigned int loop);	  //制約条件による更新スキーム(位置制約を行列の中にいれる,正方,ガウスザイデル)
	void Solve_Constraints5(unsigned int loop);	  //制約条件による更新スキーム(位置制約を行列の中にいれる,正方,三角柱)
	void Solve_Constraints6(unsigned int loop);	  //反復法による更新スキーム
	void Solve_Constraints7(unsigned int loop);	  //反復法による更新スキーム(ガウスザイデル,ヤコビ)
	void Solve_Constraints8(unsigned int loop);	  //反復法による更新スキーム(係数行列fulpivot)
	void Solve_Constraints9(unsigned int loop);	  //反復法による更新スキーム(係数行列fulpivot)(変位)
	void Solve_Constraints10(unsigned int loop);  //debug用
	void Solve_Constraints10_LU(unsigned int loop);  //debug用
	void Solve_Constraints11(unsigned int loop);  //debug用(予測位置がグループごとに異なる)
	void Solve_Constraints12(unsigned int loop); 
	void Solve_Constraints13(unsigned int loop);//new method

	void Volume_consevation(unsigned int loop);
	void CalcMassMatrix();
	void CalcPrePos();//New prediction position
	
	Eigen::Vector3d Calc_New_Exp_Pos(ParticleD* p);//位置修正(差を考える)
	Eigen::Vector3d Calc_New_Exp_Pos_Mean(ParticleD* p);//位置修正(現在は足して平均をとる)
	Eigen::Vector3d Calc_Mean_Vel(ParticleD* p);//平均速度
	Eigen::Vector3d Calc_New_Delatax_Mean(ParticleD* p);//解の修正(平均をとる)
	Eigen::Vector3d CalcMeanPos(ParticleD* p);
	Eigen::Vector3d Set_New_Exp_Pos(ParticleD* p);//グループの位置をオブジェクト位置にのSet
	Eigen::Vector3d Set_New_Exp_Pos(ParticleD* p,Eigen::VectorXd Delta);//グループの位置をオブジェクト位置にのSet.Deltaつき
	
};

#endif
