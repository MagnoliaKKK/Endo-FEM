//===========================================================================//
//@author Ikuto Endo
//@brief 局所合成行列を用いた手法でシミュレーションを行うObject(四面体一つをgroupとして使用)
//===========================================================================//
#ifndef _UseGroupOneTetraObjectDouble
#define _UseGroupOneTetraObjectDouble

#include "ObjectD.h"
class UseGroupOneTetraObjectDouble : public ObjectD {//親クラスはObjectクラス
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	UseGroupOneTetraObjectDouble(std::vector<ParticleD*> p, ObjectData data);//コンストラクタ
	~UseGroupOneTetraObjectDouble();										 //デコンストラクタ

	void Update();												 //更新作業

private:
	void Init();												 //初期化
	void Create_Groups();										 //グループ作成

	void Create_Group(std::vector<TetraElementD*> tetra_set, int tetra_group_id);
	Eigen::Vector3d Get_Center_Grid(std::vector<ParticleD*> p);
	std::map<TetraElementD*, bool> tetra_map;//TetraElementがすでにグループに使用されているかどうかのmap
	std::vector< TetraElementD* > Create_Group_Candidate();
	bool UseGroupOneTetraObjectDouble::GroupdivideNML(TetraElementD *t, int a, int b, int c);

	//debug
	unsigned int flag;

};
#endif#pragma once
#pragma once
#pragma once

#pragma once
