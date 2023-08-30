//==========================================================================//
//@author	KatsuyaKikuchi
//@brief	有限要素法四面体要素
//==========================================================================//

#ifndef _TETRAELEMENTD
#define _TETRAELEMENTD

#include "TriangleD.h"


class TetraElementD {
public:
	//Eigenをメンバ変数に使うためのマクロ定義
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	TetraElementD(std::vector<ParticleD*> particles);
	~TetraElementD();

	void Draw()const;
	bool has_Common_Points(TetraElementD t)const;

	std::vector<ParticleD*> Get_Particle()const;
	int Get_Vertex_Num()const;
	CircleD Get_Circum_Circle();

	double Get_Volume();//四面体の体積の取得
	double Calc_Volume(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, Eigen::Vector3d d);
	double Get_Exp_Volume();
						//剛性行列の作成
	void Create_Stiffness_Matrix(const Eigen::Vector3d& origin, const double& young, const double& poisson);
	void Create_Stiffness_Matrix2(const Eigen::Vector3d& origin, const double& young, const double& poisson);
	void Create_M_Matrix(double& density);	// 質量行列の作成	
	void Create_Faces();
	Eigen::Matrix<double, 6, 12> Create_B_Martix(const Eigen::Vector3d& origin);//(6x12)
	Eigen::Matrix<double, 6, 6> Create_D_Martix(const double& young, const double& poisson);//(6x6)
	void CreateDm();//create Dm
	void CreateDs();//create Ds
	void CreateDefTensor();//create deformation tensor
	void CreateStrain();//create strain tensor
	void CreateStress(const double& young, const double& poisson);//create PK stress tensor
	void CreateEnegyDensity();//create energy density
	void CreatePKFirstStress();//create PK first stress tensor
	void CreatePotentialEnergy();//create potential energy
	Eigen::Vector3d Get_Center();
	Eigen::Matrix3d Get_K_Submatrix(ParticleD* p1, ParticleD* p2);
	int Get_K_Submatrix_Edge(ParticleD* p1, ParticleD* p2);
	Eigen::Matrix3d Get_M_Submatrix(ParticleD* p1, ParticleD* p2);
	Eigen::MatrixXd TetraElementD::Get_M_Matrix();
	int Get_M_Submatrix_Edge(ParticleD* p1, ParticleD* p2);
	bool Is_Adequate()const;


	//==========================================================================//
	const std::vector<TriangleD*>& Get_Face()const;
	bool has_Common_Face(const TetraElementD *t)const;
	bool Push_Side_Tetra(TetraElementD *t);
	bool has_Face(const TriangleD *t)const;
	void Set_Same_Group(TetraElementD* t, bool b);
	void Clear_Same_Group();
	std::vector<TetraElementD*> Get_Side_Tetra();
	bool Is_Same_Group(TetraElementD* t);
	std::vector<TriangleD*> Get_Free_Face();
	bool has_Particle(const ParticleD *p)const;
	//==========================================================================//
	bool operator ==(const TetraElementD& t) const;
	bool operator<(const TetraElementD& t) const;

	double Kappa;
	void Calc_Conservation(std::vector<ParticleD*> Gp,Eigen::VectorXd m);
	Eigen::MatrixXd Ds;
	Eigen::MatrixXd Dm;
	Eigen::MatrixXd DefTensor;
	Eigen::MatrixXd Strain;
	Eigen::MatrixXd Stress;
	double EnergyDensity;
	Eigen::MatrixXd PKFirstStress;
	double PotentialEnergy;
private:
	std::vector<ParticleD*> particles;				  // 持っているparticle(4個)
	std::vector<TriangleD*> faces;					  // 持っている面(4面)
	std::vector<TetraElementD*> side_tetras;			 // 隣接する四面体の要素群
	std::map< TetraElementD*, TriangleD* > common_face; // 隣接する四面体の共通面
	std::map< TetraElementD*, bool > same_group;
	unsigned int VERTEX_NUM;
	Eigen::MatrixXd local_stiffness_matrix;			  // 四面体の剛性行列(12x12)
	Eigen::MatrixXd m_matrix;						  // 四面体の質量行列(対角行列)(12x12)
	double Ini_volume;//静止状態の体積
};

	

#endif
