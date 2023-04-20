//===========================================================================//
//@author	KatsuyaKikuchi
//@brief	四面体要素の面として使用する三角形要素
//===========================================================================//
#ifndef _TRIANGLED
#define _TRIANGLED

#include "ParticleD.h"

class TriangleD {
public:
	TriangleD(std::vector<ParticleD*> vertices);
	~TriangleD();

	void Draw(int color)const;

	Eigen::Vector3d Get_Center_Grid()const;
	Eigen::Vector3d Get_Area_Vec()const;
	std::vector<ParticleD*> Get_Particle()const;

	void Create_Area_Vec(Eigen::Vector3d centroid_vec);

	bool operator ==(const TriangleD& t)const;
private:
	const int VERTEX_NUM;
	std::vector<ParticleD*> particles;
	Eigen::Vector3d area_vec;
};
#endif
