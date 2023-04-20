#pragma once
//===========================================================================//
//@author Ikuto Endo
//@brief 合成行列を用いた手法でシミュレーションを行うObject
//===========================================================================//
#ifndef _OldFEM
#define _OldFEM

#include "Object.h"
class OldFEM : public Object {
public:
	OldFEM(std::vector<Particle*> p, ObjectData data);
	~OldFEM();

	void Update();

private:
	void Init();
	void Create_Groups();

	void Create_Group(std::vector<TetraElement*> tetra_set, int tetra_group_id);
	Eigen::Vector3f Get_Center_Grid(std::vector<Particle*> p);


	std::map<TetraElement*, bool> tetra_map;//TetraElementがすでにグループに使用されているかどうかのmap
	std::vector< TetraElement* > Create_Group_Candidate();
	bool OldFEM::GroupdivideNML(TetraElement *t, int a, int b, int c);

	//debug
	unsigned int flag;

};
#endif#pragma once
#pragma once
#pragma once
