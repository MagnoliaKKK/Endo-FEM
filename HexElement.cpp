//==========================================================================//
//@author	IkutoEndo
//@brief	HexElement.hの実装
//==========================================================================//
#include"HexElement.h"

//==========================================================================//
TetraElement::TetraElement(std::vector<Particle*> particles)
	:particles(particles),
	VERTEX_NUM(4)
{
}
TetraElement::~TetraElement() {
	//particlesに格納されているポインタは外部で定義。
	//TetraElement内では解放しないように注意。
	for (auto _f : faces) {
		delete _f;
	}
	faces.clear();
}
//==========================================================================//
//==========================================================================//
//	@start		   				初期設定									//
//==========================================================================//
void TetraElement::Create_Faces() {
	for (unsigned int i = 0; i < VERTEX_NUM; ++i) {
		std::vector<Particle*> copy_p = particles;
		copy_p.erase(copy_p.begin() + i);
		Triangle* f = new Triangle(copy_p);
		faces.push_back(f);
	}
	Eigen::Vector3f center = Eigen::Vector3f::Zero();
	for (auto _p : particles) {
		center += _p->Get_Grid();
	}
	center /= particles.size();
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
void TetraElement::Draw()const {
	Eigen::Vector3f line = particles[particles.size() - 1]->Get_Grid();

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
bool TetraElement::has_Common_Points(TetraElement t)const {
	for (int i = 0; i < VERTEX_NUM; ++i) {
		for (int j = 0; j < VERTEX_NUM; ++j) {
			if (particles[i] == t.particles[j]) return true;
		}
	}
	return false;
}
Circle TetraElement::Get_Circum_Circle() {
	Circle c;
	c.radius = 0;
	c.center = Eigen::Vector3f::Zero();

	std::vector<Eigen::Vector3f> v(VERTEX_NUM);
	for (unsigned int i = 0; i < VERTEX_NUM; ++i) {
		v[i] = this->particles[i]->Get_Grid();
	}

	Eigen::Matrix3f A;
	Eigen::Vector3f b;

	for (unsigned int i = 0; i < 3; ++i) {
		b[i] = 0.5*(v[i + 1].x()*v[i + 1].x() - v[0].x()*v[0].x()
			+ v[i + 1].y()*v[i + 1].y() - v[0].y()*v[0].y()
			+ v[i + 1].z()*v[i + 1].z() - v[0].z()*v[0].z());
		for (unsigned int j = 0; j < 3; ++j) {
			A(i, j) = v[i + 1][j] - v[0][j];
		}
	}

	Eigen::Vector3f x = A.fullPivLu().solve(b);
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
Eigen::Matrix<float, 6, 12> TetraElement::Create_B_Martix(const Eigen::Vector3f& origin) {
	double V = Get_Volume();
	Eigen::Matrix<float, 6, 12> B_Matrix = Eigen::Matrix<float, 6, 12>::Zero();
	Eigen::Vector3f p1, p2, p3, p4;
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
	B_Matrix /= 6 * V;
	return B_Matrix;
}

Eigen::Matrix<float, 6, 6> TetraElement::Create_D_Martix(const double& young, const double& poisson) {
	Eigen::MatrixXf D_Matrix = Eigen::MatrixXf::Zero(6, 6);
	double temp = young / ((1 + poisson)* (1 - 2 * poisson));
	double d1 = temp*(1 - poisson);
	double d2 = temp*poisson;
	double G = (1 - 2 * poisson)*temp / 2;

	D_Matrix(0, 0) = d1; D_Matrix(1, 1) = d1; D_Matrix(2, 2) = d1;
	D_Matrix(1, 0) = d2; D_Matrix(2, 0) = d2;
	D_Matrix(0, 1) = d2; D_Matrix(2, 1) = d2;
	D_Matrix(0, 2) = d2; D_Matrix(1, 2) = d2;
	D_Matrix(3, 3) = G; D_Matrix(4, 4) = G; D_Matrix(5, 5) = G;

	return D_Matrix;
}
void TetraElement::Create_Stiffness_Matrix(const Eigen::Vector3f& origin, const double& young, const double& poisson) {
	Eigen::MatrixXf D = Create_D_Martix(young, poisson);
	Eigen::MatrixXf B = Create_B_Martix(origin);
	double V = Get_Volume();

	local_stiffness_matrix = B.transpose()*D*B*V;
}
//===========================================================================//
//	@end		   			局所剛性行列の作成							     //
//===========================================================================//
//===========================================================================//
//	@start		   				質量行列の作成							     //
//===========================================================================//
void TetraElement::Create_M_Matrix(double& density) {
#ifdef _DEBUGMODE
	//集中質量マトリクス
	double mass = density*this->Get_Volume();
	m_matrix = (mass / particles.size())*Eigen::MatrixXf::Identity(3 * particles.size(), 3 * particles.size());
#else
	double V = Get_Volume();
	double mass = V*density;
	Eigen::MatrixXf n_matrix = Eigen::MatrixXf::Zero(3, 3 * VERTEX_NUM);
	for (int i = 0; i < VERTEX_NUM; ++i) {
		n_matrix.block(0, i * 3, 3, 3) = Eigen::Matrix3f::Identity();
	}

	m_matrix = mass * n_matrix.transpose()*n_matrix;

#endif
}
//===========================================================================//
//	@end		   				質量行列の作成							     //
//===========================================================================//
//===========================================================================//
//	@start		   				グループ作成用関数						     //
//===========================================================================//
bool TetraElement::has_Common_Face(const TetraElement *t)const {
	int common_vertex_num = 0;
	for (int i = 0; i < VERTEX_NUM; ++i) {
		for (int j = 0; j < VERTEX_NUM; ++j) {
			if (particles[i] == t->particles[j]) ++common_vertex_num;
		}
	}
	if (3 == common_vertex_num) return true;
	return false;
}
bool TetraElement::Push_Side_Tetra(TetraElement *t) {
	for (auto _f : faces) {
		if (t->has_Face(_f)) {
			this->side_tetras.push_back(t);
			common_face[t] = _f;
			return true;
		}
	}
	return false;
}
bool TetraElement::has_Face(const Triangle *t)const {
	std::vector<Particle*> tri_p = t->Get_Particle();
	for (auto _p : tri_p) {
		if (!this->has_Particle(_p)) { return false; }
	}
	return true;
}
void TetraElement::Clear_Same_Group() {
	for (auto _t : side_tetras) {
		same_group[_t] = false;
	}
}
void TetraElement::Set_Same_Group(TetraElement* t, bool b) {
	same_group[t] = b;
}
std::vector<Triangle*> TetraElement::Get_Free_Face() {
	std::vector<Triangle*> copy_face = faces;

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
bool TetraElement::has_Particle(const Particle *p)const {
	for (auto _p : this->particles) {
		if (p == _p) return true;
	}
	return false;
}
const std::vector<Triangle*>& TetraElement::Get_Face()const {
	return faces;
}
std::vector<TetraElement*> TetraElement::Get_Side_Tetra() {
	return side_tetras;
}
bool TetraElement::Is_Same_Group(TetraElement* t) {
	return same_group[t];
}
//===========================================================================//
//	@end		   				グループ作成用関数						     //
//===========================================================================//
//===========================================================================//
//	@brief	グループの全体行列作成のため部分行列を取得
//===========================================================================//
Eigen::Matrix3f TetraElement::Get_K_Submatrix(Particle* p1, Particle* p2) {
	Eigen::Matrix3f m = Eigen::Matrix3f::Zero();
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
Eigen::Matrix3f TetraElement::Get_M_Submatrix(Particle* p1, Particle* p2) {
	Eigen::Matrix3f m = Eigen::Matrix3f::Zero();
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

//===========================================================================//
std::vector<Particle*> TetraElement::Get_Particle()const {
	return this->particles;
}
const int& TetraElement::Get_Vertex_Num()const {
	return VERTEX_NUM;
}
bool TetraElement::Is_Adequate()const {
	if (particles.size() == VERTEX_NUM)return true;
	return false;
}
double TetraElement::Get_Volume() {
	Eigen::Matrix<float, 4, 4> A;
	A = Eigen::Matrix4f::Ones();
	Eigen::Vector3f l_p1, l_p2, l_p3, l_p4;
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
//===========================================================================//
//===========================================================================//
bool TetraElement::operator < (const TetraElement& t) const {
	std::vector<Particle*> p1 = this->particles, p2 = t.particles;
	std::sort(p1.begin(), p1.end(), [](const Particle* p1, const Particle* p2) { return p1->Get_Grid() < p2->Get_Grid(); });
	std::sort(p2.begin(), p2.end(), [](const Particle* p1, const Particle* p2) { return p1->Get_Grid() < p2->Get_Grid(); });

	for (int i = 0; i < VERTEX_NUM - 1; ++i) {
		if (p1[i]->Get_Grid() != p2[i]->Get_Grid()) return p1[i]->Get_Grid() < p2[i]->Get_Grid();
	}
	return p1[VERTEX_NUM - 1]->Get_Grid() < p2[VERTEX_NUM - 1]->Get_Grid();
}

bool TetraElement::operator ==(const TetraElement& t) const {
	std::vector<Particle*> p1 = this->particles, p2 = t.particles;

	std::sort(p1.begin(), p1.end(), [](const Particle* p1, const Particle* p2) { return p1->Get_Grid() < p2->Get_Grid(); });
	std::sort(p2.begin(), p2.end(), [](const Particle* p1, const Particle* p2) { return p1->Get_Grid() < p2->Get_Grid(); });
	for (int i = 0; i < VERTEX_NUM; ++i) {
		if (p1[i]->Get_Grid() != p2[i]->Get_Grid()) return false;
	}
	return true;
}
//==========================================================================//

Eigen::Vector3f TetraElement::Get_Center() {
	Eigen::Vector3f center = Eigen::Vector3f::Zero();
	for (auto _p : particles) {
		center += _p->Get_Grid();
	}
	return center / VERTEX_NUM;
}