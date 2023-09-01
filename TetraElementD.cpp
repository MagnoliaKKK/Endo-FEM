//==========================================================================//
//@author	KatsuyaKikuchi
//@brief	TetraElementD.hの実装
//==========================================================================//

#include "TetraElementD.h"

//==========================================================================//
TetraElementD::TetraElementD(std::vector<ParticleD*> particles)
	:particles(particles),//変数の初期化, 四面体要素の頂点(4こ), 頂点数
	VERTEX_NUM(4)
{
}
TetraElementD::~TetraElementD() {
	//particlesに格納されているポインタは外部で定義。
	//TetraElementD内では解放しないように注意。
	for (auto _f : faces) {
		delete _f;
	}
	faces.clear();
}
//==========================================================================//
//==========================================================================//
//	@start		   				初期設定									//
//==========================================================================//
void TetraElementD:: Create_Faces() {
	for (unsigned int i = 0; i < VERTEX_NUM; ++i) {
		std::vector<ParticleD*> copy_p = particles;
		copy_p.erase(copy_p.begin() + i);
		TriangleD* f = new TriangleD(copy_p);
		faces.push_back(f);
	}
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	for (auto _p : particles) {
		center += _p->Get_Grid();
	}
	center /= double(particles.size());
	for (auto _f : faces) {
		_f->Create_Area_Vec(_f->Get_Center_Grid() - center);
	}
}
//==========================================================================//
//	@end		   				初期設定									//
//==========================================================================//


//==========================================================================//
//	@start		   				ループ設定									//
//==========================================================================//
void TetraElementD::Draw()const {
	Eigen::Vector3d line = particles[particles.size() - 1]->Get_Grid();

	int color = WHITE;

	MyDrawLine3(particles[0]->Get_Grid(), particles[2]->Get_Grid(), color);
	MyDrawLine3(particles[0]->Get_Grid(), particles[1]->Get_Grid(), color);
	MyDrawLine3(particles[1]->Get_Grid(), particles[2]->Get_Grid(), color);
	MyDrawLine3(particles[0]->Get_Grid(), particles[3]->Get_Grid(), color);
	MyDrawLine3(particles[2]->Get_Grid(), particles[3]->Get_Grid(), color);
	MyDrawLine3(particles[1]->Get_Grid(), particles[3]->Get_Grid(), color);
	for (auto it = particles.begin(), end = particles.end(); it != end; it++) {
		(*it)->Draw();
	}
}
//==========================================================================//
//	@end		   				ループ設定									//
//==========================================================================//
//===========================================================================//
//	@start		   			ドロネー分割用関数							     //
//===========================================================================//
bool TetraElementD::has_Common_Points(TetraElementD t)const {
	for (unsigned int i = 0; i < VERTEX_NUM; ++i) {
		for (unsigned int j = 0; j < VERTEX_NUM; ++j) {
			if (particles[i] == t.particles[j]) return true;
		}
	}
	return false;
}
CircleD TetraElementD::Get_Circum_Circle() {
	CircleD c;
	c.radius = 0.0;
	c.center = Eigen::Vector3d::Zero();

	std::vector<Eigen::Vector3d> v(VERTEX_NUM);
	for (unsigned int i = 0; i < VERTEX_NUM; ++i) {
		v[i] = (this->particles[i]->Get_IM_Grid());
	}

	Eigen::Matrix3d A;
	Eigen::Vector3d b;

	for (unsigned int i = 0; i < 3; ++i) {
		b[i] = 0.5*(v[i + 1].x()*v[i + 1].x() - v[0].x()*v[0].x()
			+ v[i + 1].y()*v[i + 1].y() - v[0].y()*v[0].y()
			+ v[i + 1].z()*v[i + 1].z() - v[0].z()*v[0].z());
		for (unsigned int j = 0; j < 3; ++j) {
			A(i, j) = v[i + 1][j] - v[0][j];
		}
	}

	Eigen::Vector3d x = A.fullPivLu().solve(b);
	c.center = x;
	double dx = x[0] - v[0][0];
	double dy = x[1] - v[0][1];
	double dz = x[2] - v[0][2];
	c.radius = sqrt(dx*dx + dy*dy + dz*dz);
	return c;
}
//===========================================================================//
//	@end		   			ドロネー分割用関数							     //
//===========================================================================//


//===========================================================================//
//	@start		   			局所剛性行列の作成							     //
//===========================================================================//

//B_Matrixの作成
Eigen::Matrix<double, 6, 12> TetraElementD::Create_B_Martix(const Eigen::Vector3d& origin) {
	double V = Get_Volume();
	Eigen::Matrix<double, 6, 12> B_Matrix = Eigen::Matrix<double, 6, 12>::Zero();
	Eigen::Vector3d p1, p2, p3, p4;
	p1 = particles[0]->Get_Grid() - origin;
	p2 = particles[1]->Get_Grid() - origin;
	p3 = particles[2]->Get_Grid() - origin;
	p4 = particles[3]->Get_Grid() - origin;

	double N_x[4], N_y[4], N_z[4];
	N_x[0] = (-p3.y()*p4.z() + p3.z()*p4.y() + p2.y()*p4.z()
		- p2.y()*p3.z() - p2.z()*p4.y() + p2.z()*p3.y());
	N_y[0] = (p3.x()*p4.z() - p3.z()*p4.x() - p2.x()*p4.z()
		+ p2.x()*p3.z() + p2.z()*p4.x() - p2.z()*p3.x());
	N_z[0] = (-p3.x()*p4.y() + p3.y()*p4.x() + p2.x()*p4.y()
		- p2.x()*p3.y() - p2.y()*p4.x() + p2.y()*p3.x());

	N_x[1] = (p3.y()*p4.z() - p3.z()*p4.y() - p1.y()*p4.z()
		+ p1.y()*p3.z() + p1.z()*p4.y() - p1.z()*p3.y());
	N_y[1] = (-p3.x()*p4.z() + p3.z()*p4.x() + p1.x()*p4.z()
		- p1.x()*p3.z() - p1.z()*p4.x() + p1.z()*p3.x());
	N_z[1] = (p3.x()*p4.y() - p3.y()*p4.x() - p1.x()*p4.y()
		+ p1.x()*p3.y() + p1.y()*p4.x() - p1.y()*p3.x());

	N_x[2] = (-p2.y()*p4.z() + p2.z()*p4.y() + p1.y()*p4.z()
		- p1.y()*p2.z() - p1.z()*p4.y() + p1.z()*p2.y());
	N_y[2] = (p2.x()*p4.z() - p2.z()*p4.x() - p1.x()*p4.z()
		+ p1.x()*p2.z() + p1.z()*p4.x() - p1.z()*p2.x());
	N_z[2] = (-p2.x()*p4.y() + p2.y()*p4.x() + p1.x()*p4.y()
		- p1.x()*p2.y() - p1.y()*p4.x() + p1.y()*p2.x());

	N_x[3] = (p2.y()*p3.z() - p2.z()*p3.y() - p1.y()*p3.z()
		+ p1.y()*p2.z() + p1.z()*p3.y() - p1.z()*p2.y());
	N_y[3] = (-p2.x()*p3.z() + p2.z()*p3.x() + p1.x()*p3.z()
		- p1.x()*p2.z() - p1.z()*p3.x() + p1.z()*p2.x());
	N_z[3] = (p2.x()*p3.y() - p2.y()*p3.x() - p1.x()*p3.y()
		+ p1.x()*p2.y() + p1.y()*p3.x() - p1.y()*p2.x());


	for (unsigned int i = 0; i < 4; i++) {
		B_Matrix(0, 3 * i) = N_x[i];
		B_Matrix(1, 3 * i + 1) = N_y[i];
		B_Matrix(2, 3 * i + 2) = N_z[i];
		B_Matrix(3, 3 * i) = N_y[i]; B_Matrix(3, 3 * i + 1) = N_x[i];
		B_Matrix(4, 3 * i) = N_z[i]; B_Matrix(4, 3 * i + 2) = N_x[i];
		B_Matrix(5, 3 * i + 1) = N_z[i]; B_Matrix(5, 3 * i + 2) = N_y[i];
	}
	B_Matrix = B_Matrix / (6.0 * V);//ここのVはKを作るときに打ち消すことができる
	return B_Matrix;
}

//D_Matrixの作成
Eigen::Matrix<double, 6, 6> TetraElementD::Create_D_Martix(const double& young, const double& poisson) {
	Eigen::MatrixXd D_Matrix = Eigen::MatrixXd::Zero(6, 6);
	double temp = young / ((1 + poisson)* (1 - 2 * poisson));
	double d1 = temp*(1 - poisson);
	double d2 = temp*poisson;
	double G = (1 - 2 * poisson)*temp / 2;

	D_Matrix(0, 0) = d1; D_Matrix(1, 1) = d1; D_Matrix(2, 2) = d1;
	D_Matrix(1, 0) = d2; D_Matrix(2, 0) = d2;
	D_Matrix(0, 1) = d2; D_Matrix(2, 1) = d2;
	D_Matrix(0, 2) = d2; D_Matrix(1, 2) = d2;
	D_Matrix(3, 3) = G; D_Matrix(4, 4) = G; D_Matrix(5, 5) = G;

	//Dの中身をチェック
	/*std::cout << "D" << std::endl;
	std::cout << std::setprecision(10) << D_Matrix << std::endl;
	for (int qii = 0; qii< 6; qii++) {
		for (int qij = 0; qij < 6; qij++) {
			std::cout << (D_Matrix).block(qii, qij, 1, 1) - ((D_Matrix).transpose()).block(qii, qij, 1, 1) << ",";
		}
		std::cout << "F" << std::endl;
	}*/
	//Dは確かに対称行列

	return D_Matrix;
}

//合成行列(K_Matrix)の作成
void TetraElementD::Create_Stiffness_Matrix(const Eigen::Vector3d& origin, const double& young, const double& poisson) {
	local_stiffness_matrix = Eigen::MatrixXd::Zero(12, 12);
	//std::cout << "local0" << std::endl;
	//std::cout << std::setprecision(10) << local_stiffness_matrix << std::endl;
	Eigen::MatrixXd D = Create_D_Martix(young, poisson);
	Eigen::MatrixXd B = Create_B_Martix(origin);
	double V = Get_Volume();

	//std::cout << "B" << std::endl;
	//std::cout << std::setprecision(10) << B << std::endl;

	local_stiffness_matrix = B.transpose()*D*B*V;//ここのVはBと打ち消すことができる

	std::cout << "local" << std::endl;
	std::cout << std::setprecision(10) << local_stiffness_matrix << std::endl;
	std::cout << "localtt" << std::endl;
	std::cout << std::setprecision(10) << local_stiffness_matrix- local_stiffness_matrix.transpose() << std::endl;
}
//異なる方法で作った剛性行列(originはグループの重心)
void TetraElementD::Create_Stiffness_Matrix2(const Eigen::Vector3d& origin, const double& young, const double& poisson) {
	//定数
	double lambda = (young * poisson) /((1 + poisson)*(1 - 2 * poisson)) ;
	double G = young / (2*(1 + poisson));

	Eigen::MatrixXd XX = Eigen::MatrixXd::Zero(4, 4);
	Eigen::MatrixXd YY = Eigen::MatrixXd::Zero(4, 4);
	Eigen::MatrixXd ZZ = Eigen::MatrixXd::Zero(4, 4);
	Eigen::MatrixXd XY = Eigen::MatrixXd::Zero(4, 4);
	Eigen::MatrixXd XZ = Eigen::MatrixXd::Zero(4, 4);
	Eigen::MatrixXd YX = Eigen::MatrixXd::Zero(4, 4);
	Eigen::MatrixXd YZ = Eigen::MatrixXd::Zero(4, 4);
	Eigen::MatrixXd ZX = Eigen::MatrixXd::Zero(4, 4);
	Eigen::MatrixXd ZY = Eigen::MatrixXd::Zero(4, 4);

	Eigen::Vector3d p1, p2, p3, p4;
	p1 = particles[0]->Get_Grid() - origin;
	p2 = particles[1]->Get_Grid() - origin;
	p3 = particles[2]->Get_Grid() - origin;
	p4 = particles[3]->Get_Grid() - origin;

	double N_x[4], N_y[4], N_z[4];
	N_x[0] = (-p3.y()*p4.z() + p3.z()*p4.y() + p2.y()*p4.z()
		- p2.y()*p3.z() - p2.z()*p4.y() + p2.z()*p3.y());
	N_y[0] = (p3.x()*p4.z() - p3.z()*p4.x() - p2.x()*p4.z()
		+ p2.x()*p3.z() + p2.z()*p4.x() - p2.z()*p3.x());
	N_z[0] = (-p3.x()*p4.y() + p3.y()*p4.x() + p2.x()*p4.y()
		- p2.x()*p3.y() - p2.y()*p4.x() + p2.y()*p3.x());

	N_x[1] = (p3.y()*p4.z() - p3.z()*p4.y() - p1.y()*p4.z()
		+ p1.y()*p3.z() + p1.z()*p4.y() - p1.z()*p3.y());
	N_y[1] = (-p3.x()*p4.z() + p3.z()*p4.x() + p1.x()*p4.z()
		- p1.x()*p3.z() - p1.z()*p4.x() + p1.z()*p3.x());
	N_z[1] = (p3.x()*p4.y() - p3.y()*p4.x() - p1.x()*p4.y()
		+ p1.x()*p3.y() + p1.y()*p4.x() - p1.y()*p3.x());

	N_x[2] = (-p2.y()*p4.z() + p2.z()*p4.y() + p1.y()*p4.z()
		- p1.y()*p2.z() - p1.z()*p4.y() + p1.z()*p2.y());
	N_y[2] = (p2.x()*p4.z() - p2.z()*p4.x() - p1.x()*p4.z()
		+ p1.x()*p2.z() + p1.z()*p4.x() - p1.z()*p2.x());
	N_z[2] = (-p2.x()*p4.y() + p2.y()*p4.x() + p1.x()*p4.y()
		- p1.x()*p2.y() - p1.y()*p4.x() + p1.y()*p2.x());

	N_x[3] = (p2.y()*p3.z() - p2.z()*p3.y() - p1.y()*p3.z()
		+ p1.y()*p2.z() + p1.z()*p3.y() - p1.z()*p2.y());
	N_y[3] = (-p2.x()*p3.z() + p2.z()*p3.x() + p1.x()*p3.z()
		- p1.x()*p2.z() - p1.z()*p3.x() + p1.z()*p2.x());
	N_z[3] = (p2.x()*p3.y() - p2.y()*p3.x() - p1.x()*p3.y()
		+ p1.x()*p2.y() + p1.y()*p3.x() - p1.y()*p2.x());
	//本当は上の値は6Vでわらないといけないけど、あとで積分してVをかけるからあとで計算
	for (int ki = 0; ki < 4;ki++) {
		for (int kj = 0; kj < 4; kj++) {
			XX(ki, kj) = N_x[ki] * N_x[kj];
			YY(ki, kj) = N_y[ki] * N_y[kj];
			ZZ(ki, kj) = N_z[ki] * N_z[kj];

			XY(ki, kj) = N_x[ki] * N_y[kj];
			XZ(ki, kj) = N_x[ki] * N_z[kj];
			
			YX(ki, kj) = N_y[ki] * N_x[kj];
			YZ(ki, kj) = N_y[ki] * N_z[kj];

			ZX(ki, kj) = N_z[ki] * N_x[kj]; 
			ZY(ki, kj) = N_z[ki] * N_y[kj];
		}
	}

	Eigen::MatrixXd ONE = Eigen::MatrixXd::Zero(12, 12);
	Eigen::MatrixXd TWO = Eigen::MatrixXd::Zero(12, 12);
	Eigen::MatrixXd THREE = Eigen::MatrixXd::Zero(12, 12);
	Eigen::MatrixXd FOUR = Eigen::MatrixXd::Zero(12, 12);

	//一つ目の行列作成
	ONE.block(0, 0, 4, 4) = XX;
	                            ONE.block(4, 4, 4, 4) = YY;
	                                                        ONE.block(8, 8, 4, 4) = ZZ;
	ONE = 2 * G * ONE;
	//二つ目の行列作成
	TWO.block(0, 0, 4, 4) = XX; TWO.block(0, 4, 4, 4) = XY; TWO.block(0, 8, 4, 4) = XZ;
	TWO.block(4, 0, 4, 4) = YX; TWO.block(4, 4, 4, 4) = YY; TWO.block(4, 8, 4, 4) = YZ;
	TWO.block(8, 0, 4, 4) = ZX; TWO.block(8, 4, 4, 4) = ZY; TWO.block(8, 8, 4, 4) = ZZ;
	TWO = lambda * TWO;

	//三つ目の行列作成
	THREE.block(0, 0, 4, 4) = ZZ; THREE.block(0, 4, 4, 4) = YX;
	THREE.block(4, 0, 4, 4) = XY; THREE.block(4, 4, 4, 4) = XX; THREE.block(4, 8, 4, 4) = ZY;
							      THREE.block(8, 4, 4, 4) = YZ; THREE.block(8, 8, 4, 4) = YY;
	THREE = G * THREE;

	//四つ目の行列作成
	FOUR.block(0, 0, 4, 4) = YY;							  FOUR.block(0, 8, 4, 4) = ZX;
								 FOUR.block(4, 4, 4, 4) = ZZ; 
	FOUR.block(8, 0, 4, 4) = XZ;							  FOUR.block(8, 8, 4, 4) = XX;
	FOUR = G * FOUR;

	//四つの行列の合成
	local_stiffness_matrix = Eigen::MatrixXd::Zero(12, 12);
	local_stiffness_matrix = ONE + TWO + THREE + FOUR;
	//微分した時に6Vで割ってなくて(2乗する)、そのあと積分してVをかけるから、総じて36Vで割る
	local_stiffness_matrix = local_stiffness_matrix / (36.0 * this->Get_Volume());

	Ini_volume = this->Get_Volume();
	//出力
	//std::cout << "local" << std::endl;
	//std::cout << local_stiffness_matrix << std::endl;
	//std::cout << "localtt" << std::endl;
	//std::cout << std::setprecision(10) << local_stiffness_matrix - local_stiffness_matrix.transpose() << std::endl;
	
	//pythonの出力
	/*
	std::cout << "local_stiffness_matrix is"  << std::endl;
	std::cout << "np.array([";
	for (unsigned int pi = 0; pi < 3 * 4 - 1; pi++) {
		//値
		std::cout << "[";
		for (unsigned int pj = 0; pj < 3 * 4 - 1; pj++) {
			//値
			std::cout << local_stiffness_matrix.block(pi, pj, 1, 1) << ",";
		}
		std::cout << local_stiffness_matrix.block(pi, 3 * 4 - 1, 1, 1) << "]," << std::endl;
	}
	std::cout << "[";
	for (unsigned int pj = 0; pj < 3 * 4 - 1; pj++) {
		//値
		std::cout << local_stiffness_matrix.block(3 * 4 - 1, pj, 1, 1) << ",";
	}
	std::cout << local_stiffness_matrix.block(3 * 4 - 1, 3 * 4 - 1, 1, 1) << "]])" << std::endl;
	*/
	
}
void TetraElementD::CreateDm() {
	Eigen::Vector3d p1, p2, p3, p4;
	Dm = Eigen::Matrix<double, 3, 3>::Zero();
	p1 = particles[0]->Get_Initial_Pos();
	p2 = particles[1]->Get_Initial_Pos();
	p3 = particles[2]->Get_Initial_Pos();
	p4 = particles[3]->Get_Initial_Pos();
	
	Dm << p1[0] - p4[0], p2[0] - p4[0], p3[0] - p4[0],
		p1[1] - p4[1], p2[1] - p4[1], p3[1] - p4[1],
		p1[2] - p4[2], p2[2] - p4[2], p3[2] - p4[2];
	
}
void TetraElementD::CreateDs() {
	Ds = Eigen::Matrix<double, 3, 3>::Zero();
	Eigen::Vector3d p1, p2, p3, p4;
	p1 = particles[0]->Get_Grid();
	p2 = particles[1]->Get_Grid();
	p3 = particles[2]->Get_Grid();
	p4 = particles[3]->Get_Grid();
	Eigen::Matrix<double, 3, 3> Ds;
	Ds << p1[0] - p4[0], p2[0] - p4[0], p3[0] - p4[0],
		p1[1] - p4[1], p2[1] - p4[1], p3[1] - p4[1],
		p1[2] - p4[2], p2[2] - p4[2], p3[2] - p4[2];
	
}
void TetraElementD::CreateDefTensor() {
	DefTensor = Ds * Dm.inverse();
}
void TetraElementD::CreateStrain() {
	Strain = 0.5 * (DefTensor.transpose() * DefTensor - Eigen::Matrix<double, 3, 3>::Identity());
}
void TetraElementD::CreateStress(const double& young, const double& poisson) {
	//Create stress tensor, D_Matrix is the strain-displacement matrix,strain is calculated in the function CreateStrain
//stress = D_Matrix * strain
	Eigen::MatrixXd D = Create_D_Martix(young, poisson);

	Eigen::VectorXd E_voigt(6);
	E_voigt << Strain(0, 0), Strain(1, 1), Strain(2, 2), 2 * Strain(1, 2), 2 * Strain(0, 2), 2 * Strain(0, 1);
	Eigen::VectorXd S_voigt = D * E_voigt;

	// Convert S from Voigt notation back to 3x3 form
	
	Stress = Eigen::MatrixXd::Zero(3, 3);
	Stress(0,0) = S_voigt(0);
	Stress(1,1) = S_voigt(1);
	Stress(2,2) = S_voigt(2);
	Stress(1,2) = S_voigt(3);
	Stress(0,2) = S_voigt(4);
	Stress(0,1) = S_voigt(5);	
	Stress(2,1) = S_voigt(3);
	Stress(2,0) = S_voigt(4);
	Stress(1,0) = S_voigt(5);

}
void TetraElementD::CreateEnegyDensity() {
	EnergyDensity = 0.5 * (Strain.transpose() * Stress).trace();
}
void TetraElementD::CreatePKFirstStress() {
	PKFirstStress = DefTensor * Stress;
}
void TetraElementD::CreatePotentialEnergy() {
	PotentialEnergy = EnergyDensity * Ini_volume;
}
void TetraElementD::CreateEnergyGrad() {
	EnergyGrad = Eigen::Matrix<double, 3, 4>::Zero();
	//EnergyGrad = Ini_volume * PKFirstStress * Dm.inverse().transpose()
	
	Eigen::Matrix<double, 3, 3> tempResult = Ini_volume * PKFirstStress * Dm.inverse().transpose();
	// 设置前三列
	EnergyGrad.block<3, 3>(0, 0) = tempResult;
	// 计算前三列的和
	Eigen::Matrix<double, 3, 1> sumOfColumns = tempResult.rowwise().sum();
	// 设置第四列为前三列的和的负数
	EnergyGrad.col(3) = -sumOfColumns;
}

Eigen::Vector3d TetraElementD::Get_SubEnergyGrad(ParticleD* p) {
	Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
	// Find the index of particle p in the particles list.
	size_t p_index = -1;
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		if (p == *it) {
			p_index = std::distance(particles.begin(), it);
			break;
		}
	}

	if (-1 == p_index) {
		return gradient; // Returns zero vector.
	}
	// Extract energy gradient for particle p from EnergyGrad matrix.
	gradient = EnergyGrad.col(p_index);

	return gradient;
}

//===========================================================================//
//	@end		   			局所剛性行列の作成							     //
//===========================================================================//


//===========================================================================//
//	@start		   				質量行列の作成							     //
//===========================================================================//
void TetraElementD::Create_M_Matrix(double& density) {
	//質量行列が対角成分のみのとき
	if (mdiag == TRUE) {
		//集中質量マトリクス
		double mass = density * this->Get_Volume();//四面体の質量
		m_matrix = (mass / 4.0) * Eigen::MatrixXd::Identity(3 * particles.size(), 3 * particles.size());
	}
	//対角成分のみでないとき（まだ実装できていない）
	else {
		double V = Get_Volume();
		double mass = V*density;
		Eigen::MatrixXd n_matrix = Eigen::MatrixXd::Zero(3, 3 * VERTEX_NUM);
		for (unsigned int i = 0; i < VERTEX_NUM; ++i) {
			n_matrix.block(0, i * 3, 3, 3) = Eigen::Matrix3d::Identity();
		}
		m_matrix = mass / 4 * n_matrix.transpose()*n_matrix;
		//std::cout << " M Matrix Of Group" << m_matrix << std::endl;
	}
}
//===========================================================================//
//	@end		   				質量行列の作成							     //
//===========================================================================//


//===========================================================================//
//	@start		   				グループ作成用関数						     //
//===========================================================================//
bool TetraElementD::has_Common_Face(const TetraElementD *t)const {
	int common_vertex_num = 0;
	for (unsigned int i = 0; i < VERTEX_NUM; ++i) {
		for (unsigned int j = 0; j < VERTEX_NUM; ++j) {
			if (particles[i] == t->particles[j]) ++common_vertex_num;
		}
	}
	if (3 == common_vertex_num) return true;
	return false;
}
bool TetraElementD::Push_Side_Tetra(TetraElementD *t) {
	for (auto _f : faces) {
		if (t->has_Face(_f)) {
			this->side_tetras.push_back(t);
			common_face[t] = _f;
			return true;
		}
	}
	return false;
}
bool TetraElementD::has_Face(const TriangleD *t)const {
	std::vector<ParticleD*> tri_p = t->Get_Particle();
	for (auto _p : tri_p) {
		if (!this->has_Particle(_p)) { return false; }
	}
	return true;
}
void TetraElementD::Clear_Same_Group() {
	for (auto _t : side_tetras) {
		same_group[_t] = false;
	}
}
void TetraElementD::Set_Same_Group(TetraElementD* t, bool b) {
	same_group[t] = b;
}
std::vector<TriangleD*> TetraElementD::Get_Free_Face() {
	std::vector<TriangleD*> copy_face = faces;

	for (auto side_t : side_tetras) {
		if (same_group[side_t]) {
			for (auto it = copy_face.begin(); it != copy_face.end();) {
				if (*it == common_face[side_t]) { it = copy_face.erase(it); }
				else { ++it; }
			}
		}
	}
	return copy_face;
}
bool TetraElementD::has_Particle(const ParticleD *p)const {
	for (auto _p : this->particles) {
		if (p == _p) return true;
	}
	return false;
}
const std::vector<TriangleD*>& TetraElementD::Get_Face()const {
	return faces;
}
std::vector<TetraElementD*> TetraElementD::Get_Side_Tetra() {
	return side_tetras;
}
bool TetraElementD::Is_Same_Group(TetraElementD* t) {
	return same_group[t];
}
//===========================================================================//
//	@end		   				グループ作成用関数						     //
//===========================================================================//


//===========================================================================//
//	@brief	グループの全体行列作成のため部分行列を取得
//===========================================================================//
Eigen::Matrix3d TetraElementD::Get_K_Submatrix(ParticleD* p1, ParticleD* p2) {
	Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
	size_t p1_index = -1, p2_index = -1;
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		size_t index = std::distance(particles.begin(), it);
		if (p1 == *it) { p1_index = index; }
		if (p2 == *it) { p2_index = index; }
	}
	if (-1 == p1_index || -1 == p2_index) {
		return m;
	}

	m = local_stiffness_matrix.block<3, 3>(3 * p1_index, 3 * p2_index);
	return m;
}
//隣接行列を作るための関数
int TetraElementD::Get_K_Submatrix_Edge(ParticleD* p1, ParticleD* p2) {
	size_t p1_index = -1, p2_index = -1;
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		size_t index = std::distance(particles.begin(), it);
		if (p1 == *it) { p1_index = index; }
		if (p2 == *it) { p2_index = index; }
	}
	if (-1 == p1_index || -1 == p2_index) {
		return 0;
	}
	return 1;
}
Eigen::Matrix3d TetraElementD::Get_M_Submatrix(ParticleD* p1, ParticleD* p2) {
	Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
	size_t p1_index = -1, p2_index = -1;
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		size_t index = std::distance(particles.begin(), it);
		if (p1 == *it) { p1_index = index; }
		if (p2 == *it) { p2_index = index; }
	}
	if (-1 == p1_index || -1 == p2_index) { return m; }

	m = m_matrix.block<3, 3>(3 * p1_index, 3 * p2_index);
	return m;
}
//隣接行列を作るための関数(たぶんKと一緒)
int TetraElementD::Get_M_Submatrix_Edge(ParticleD* p1, ParticleD* p2) {
	size_t p1_index = -1, p2_index = -1;
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		size_t index = std::distance(particles.begin(), it);
		if (p1 == *it) { p1_index = index; }
		if (p2 == *it) { p2_index = index; }
	}
	if (-1 == p1_index || -1 == p2_index) {
		return 0;
	}
	return 1;
}

//===========================================================================//
std::vector<ParticleD*> TetraElementD::Get_Particle()const {
	return this->particles;
}
int TetraElementD::Get_Vertex_Num()const {
	return VERTEX_NUM;
}
bool TetraElementD::Is_Adequate()const {
	if (particles.size() == VERTEX_NUM)return true;
	return false;
}
//四面体の体積の取得
double TetraElementD::Get_Volume() {
	Eigen::Matrix<double, 4, 4> A;
	A = Eigen::Matrix4d::Ones();
	Eigen::Vector3d l_p1, l_p2, l_p3, l_p4;
	l_p1 = particles[0]->Get_Grid();
	l_p2 = particles[1]->Get_Grid();
	l_p3 = particles[2]->Get_Grid();
	l_p4 = particles[3]->Get_Grid();

	A(0, 1) = l_p1.x(); A(0, 2) = l_p1.y(); A(0, 3) = l_p1.z();
	A(1, 1) = l_p2.x(); A(1, 2) = l_p2.y(); A(1, 3) = l_p2.z();
	A(2, 1) = l_p3.x(); A(2, 2) = l_p3.y(); A(2, 3) = l_p3.z();
	A(3, 1) = l_p4.x(); A(3, 2) = l_p4.y(); A(3, 3) = l_p4.z();

	return fabs(A.determinant() / 6.0);
}
double TetraElementD::Calc_Volume(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, Eigen::Vector3d d) {
	Eigen::Matrix<double, 4, 4> A;
	A = Eigen::Matrix4d::Ones();
	Eigen::Vector3d l_p1, l_p2, l_p3, l_p4;
	l_p1 = a;
	l_p2 = b;
	l_p3 = c;
	l_p4 = d;

	A(0, 1) = l_p1.x(); A(0, 2) = l_p1.y(); A(0, 3) = l_p1.z();
	A(1, 1) = l_p2.x(); A(1, 2) = l_p2.y(); A(1, 3) = l_p2.z();
	A(2, 1) = l_p3.x(); A(2, 2) = l_p3.y(); A(2, 3) = l_p3.z();
	A(3, 1) = l_p4.x(); A(3, 2) = l_p4.y(); A(3, 3) = l_p4.z();

	return fabs(A.determinant() / 6.0);
}
double TetraElementD::Get_Exp_Volume() {
	Eigen::Matrix<double, 4, 4> A;
	A = Eigen::Matrix4d::Ones();
	Eigen::Vector3d l_p1, l_p2, l_p3, l_p4;
	l_p1 = particles[0]->Get_Exp_Pos();
	l_p2 = particles[1]->Get_Exp_Pos();
	l_p3 = particles[2]->Get_Exp_Pos();
	l_p4 = particles[3]->Get_Exp_Pos();
	A(0, 1) = l_p1.x(); A(0, 2) = l_p1.y(); A(0, 3) = l_p1.z();
	A(1, 1) = l_p2.x(); A(1, 2) = l_p2.y(); A(1, 3) = l_p2.z();
	A(2, 1) = l_p3.x(); A(2, 2) = l_p3.y(); A(2, 3) = l_p3.z();
	A(3, 1) = l_p4.x(); A(3, 2) = l_p4.y(); A(3, 3) = l_p4.z();

	return fabs(A.determinant() / 6.0);
}
Eigen::Vector3d TetraElementD::Get_Center() {
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	for (auto _p : particles) {
		center += _p->Get_Grid();
	}
	return center / VERTEX_NUM;
}
//四面体の体積の取得
Eigen::MatrixXd TetraElementD::Get_M_Matrix() {
	return m_matrix;
}
//===========================================================================//
//===========================================================================//
bool TetraElementD::operator < (const TetraElementD& t) const {
	std::vector<ParticleD*> p1 = this->particles, p2 = t.particles;
	std::sort(p1.begin(), p1.end(), [](const ParticleD* p1, const ParticleD* p2) { return p1->Get_Grid() < p2->Get_Grid(); });
	std::sort(p2.begin(), p2.end(), [](const ParticleD* p1, const ParticleD* p2) { return p1->Get_Grid() < p2->Get_Grid(); });

	for (unsigned int i = 0; i < VERTEX_NUM - 1; ++i) {
		if (p1[i]->Get_Grid() != p2[i]->Get_Grid()) return p1[i]->Get_Grid() < p2[i]->Get_Grid();
	}
	return p1[VERTEX_NUM - 1]->Get_Grid() < p2[VERTEX_NUM - 1]->Get_Grid();
}

bool TetraElementD::operator ==(const TetraElementD& t) const {
	std::vector<ParticleD*> p1 = this->particles, p2 = t.particles;

	std::sort(p1.begin(), p1.end(), [](const ParticleD* p1, const ParticleD* p2) { return p1->Get_Grid() < p2->Get_Grid(); });
	std::sort(p2.begin(), p2.end(), [](const ParticleD* p1, const ParticleD* p2) { return p1->Get_Grid() < p2->Get_Grid(); });
	for (unsigned int i = 0; i < VERTEX_NUM; ++i) {
		if (p1[i]->Get_Grid() != p2[i]->Get_Grid()) return false;
	}
	return true;
}
//==========================================================================//
void TetraElementD::Calc_Conservation(std::vector<ParticleD*> Gp,Eigen::VectorXd Gm) {
	double lambda = 1.0e+10;//材料パラメータ
	//四点の獲得
	Eigen::Vector3d l_p0, l_p1, l_p2, l_p3;
	l_p0 = particles[0]->Get_Exp_Pos();
	l_p1 = particles[1]->Get_Exp_Pos();
	l_p2 = particles[2]->Get_Exp_Pos();
	l_p3 = particles[3]->Get_Exp_Pos();
	Eigen::VectorXd m = Eigen::VectorXd::Zero(4);
	for (unsigned int pi = 0; pi < Gp.size();pi++) {
		if (particles[0]->p_id == Gp[pi]->p_id) {
			m[0] = Gm[pi];
		}
		else if (particles[1]->p_id == Gp[pi]->p_id) {
			m[1] = Gm[pi];
		}
		else if (particles[2]->p_id == Gp[pi]->p_id) {
			m[2] = Gm[pi];
		}
		else if (particles[3]->p_id == Gp[pi]->p_id) {
			m[3] = Gm[pi];
		}
	}
	if ((m[0] * m[1] * m[2] * m[3]) <= 0.0) {
		std::cout << "ERROR M 598" << std::endl;
	}
	double DKappa = 0.0;
	double C = 0.0;
	double alpha = 0.0;
	//制約条件Ctの計算
	C = ((l_p1 - l_p0).cross(l_p2 - l_p0)).dot(l_p3- l_p0);
	C = (C / (6 * Ini_volume)) - 1.0;
	//コンプライアンス係数アルファの計算
	alpha = 1.0/(lambda * Ini_volume * TIME_STEP * TIME_STEP);
	alpha = 0.0001;
	//勾配計算
	Eigen::MatrixXd Ct = Eigen::MatrixXd::Zero(3,4);
	Ct.block(0, 1, 3, 1) = ((l_p2 - l_p0).cross(l_p3 - l_p0)) / 6.0 / this->Get_Exp_Volume();
	Ct.block(0, 2, 3, 1) = ((l_p3 - l_p0).cross(l_p1 - l_p0)) / 6.0 / this->Get_Exp_Volume();
	Ct.block(0, 3, 3, 1) = ((l_p1 - l_p0).cross(l_p2 - l_p0)) / 6.0 / this->Get_Exp_Volume();
	Ct.block(0, 0, 3, 1) = -1 * (Ct.block(0, 1, 3, 1) + Ct.block(0, 2, 3, 1) + Ct.block(0, 3, 3, 1));
	//Delta Kappaの計算
	double Frac = 0.0;
	Frac = (Ct.block(0, 1, 3, 1).squaredNorm() / m[1]) + (Ct.block(0, 2, 3, 1).squaredNorm() / m[2]) + (Ct.block(0, 3, 3, 1).squaredNorm() / m[3]) + (Ct.block(0, 0, 3, 1).squaredNorm() / m[0]);
	Frac = Frac + alpha;

	DKappa = C - alpha * Kappa;
	DKappa = DKappa / Frac;
	//Kappaの更新
	Kappa = Kappa + DKappa;
	for (unsigned int pi = 0; pi < 4; pi++) {
		if (!(particles[pi]->Is_Fixed())) {
			particles[pi]->Set_Exp_Pos(particles[pi]->Get_Exp_Pos() + (DKappa * Ct.block(0, pi, 3, 1))/m[pi]);
		}
	}
}