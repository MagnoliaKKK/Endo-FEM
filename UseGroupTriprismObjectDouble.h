#pragma once
//===========================================================================//
//@author Ikuto Endo
//@brief 局所合成行列を用いた手法でシミュレーションを行うObject(三角柱をgroupとして使用)
//===========================================================================//
#ifndef _UseGroupTriprismObjectDouble
#define _UseGroupTriprismObjectDouble

#include "ObjectD.h"
class UseGroupTriprismObjectDouble : public ObjectD {//親クラスはObjectクラス
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	 UseGroupTriprismObjectDouble(std::vector<ParticleD*> p, ObjectData data);//コンストラクタ
	~UseGroupTriprismObjectDouble();										 //デコンストラクタ

	void Update();												 //更新作業

private:
	void Init();												 //初期化
	void Create_Groups();										 //グループ作成

	void Create_Group(std::vector<TetraElementD*> tetra_set, int tetra_group_id);
	Eigen::Vector3d Get_Center_Grid(std::vector<ParticleD*> p);
	std::map<TetraElementD*, bool> tetra_map;//TetraElementがすでにグループに使用されているかどうかのmap
	std::vector< TetraElementD* > Create_Group_Candidate();
	bool UseGroupTriprismObjectDouble::GroupdivideNML(TetraElementD *t, int a, int b, int c);
	bool UseGroupTriprismObjectDouble::BlockPrismDivBy(TetraElementD *t, int a);

	//debug
	unsigned int flag;

};
#endif#pragma once
#pragma once
#pragma once

#pragma once

