//===========================================================================//
//@author	KatsuyaKikuchi
//@brief	Triangle.hの実装
//===========================================================================//
#include "TriangleD.h"

//===========================================================================//
TriangleD::TriangleD(std::vector<ParticleD*> vertices)
	:VERTEX_NUM(3),
	particles(vertices)
{
}
TriangleD::~TriangleD() {}
//===========================================================================//
//==========================================================================//
//	@start		   				ループ設定									//
//==========================================================================//
void TriangleD::Draw(int color)const {
	VECTOR p[3];
	for (unsigned int i = 0; i < 3; ++i) {
		Eigen::Vector3d grid = particles[i]->Get_Grid();
		p[i] = VGet(float(grid.x()), 480 - float(grid.y()), float(grid.z()));
	}
	Eigen::Vector3d normal = this->area_vec.normalized();
	DrawTriangle3D(p[0], p[1], p[2], color, true);
}
//==========================================================================//
//	@end		   				ループ設定									//
//==========================================================================//
bool TriangleD::operator ==(const TriangleD& t) const {
	std::vector<ParticleD*> v1 = this->particles, v2 = t.particles;

	std::sort(v1.begin(), v1.end(), [](const ParticleD* v1, const ParticleD* v2) { return v1->Get_Grid() < v2->Get_Grid(); });
	std::sort(v2.begin(), v2.end(), [](const ParticleD* v1, const ParticleD* v2) { return v1->Get_Grid() < v2->Get_Grid(); });
	for (int i = 0; i < VERTEX_NUM; ++i) {
		if (v1[i]->Get_Grid() != v2[i]->Get_Grid()) return false;
	}
	return true;
}


std::vector<ParticleD*> TriangleD::Get_Particle()const {
	return particles;
}
Eigen::Vector3d TriangleD::Get_Center_Grid()const {
	return (particles[0]->Get_Grid() + particles[1]->Get_Grid() + particles[2]->Get_Grid()) / VERTEX_NUM;
}
Eigen::Vector3d TriangleD::Get_Area_Vec()const {
	return area_vec;
}

void TriangleD::Create_Area_Vec(Eigen::Vector3d centroid_vec) {
	area_vec = ((particles[2]->Get_Grid() - particles[0]->Get_Grid()).cross(particles[1]->Get_Grid() - particles[0]->Get_Grid())).normalized();

	if (area_vec.dot(centroid_vec) == 0) {
		area_vec = centroid_vec.normalized()*area_vec.norm();
	}
	else if (area_vec.dot(centroid_vec) < 0) {
		area_vec *= -1;
	}
}
