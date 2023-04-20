#pragma once
//===========================================================================//
//@author Ikuto Endo
//@brief従来のFEMでシミュレーションを行うObject(モデルを1つのグループにする)
//===========================================================================//
#ifndef _UseOldFEMDouble
#define _UseOldFEMDouble

#include "ObjectD.h"
class UseOldFEMDouble : public ObjectD {//親クラスはObjectクラス
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	UseOldFEMDouble(std::vector<ParticleD*> p, ObjectData data);//コンストラクタ
	~UseOldFEMDouble();										    //デコンストラクタ

	void Update();												 //更新作業

private:
	void Init();												 //初期化
	void Create_Groups();										 //グループ作成

	void Create_Group(std::vector<TetraElementD*> tetra_set, int tetra_group_id);
	Eigen::Vector3d Get_Center_Grid(std::vector<ParticleD*> p);
	std::map<TetraElementD*, bool> tetra_map;//TetraElementがすでにグループに使用されているかどうかのmap
	std::vector< TetraElementD* > Create_Group_Candidate();

	//debug
	unsigned int flag;

};
#endif#pragma once
