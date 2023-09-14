#include "Base.h"


#include "UseBlockObjectDouble.h"

#include "UseLinearFEMDouble.h"

#include "InputKey.h"
#include <windows.h>

#include <stdio.h>

#include <stdlib.h>

#include <fenv.h>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

//#pragma fenv_access (on)
#pragma comment(lib, "winmm.lib")
//カメラの回転スピード
static const float ROTATE_SPEED = DX_PI_F / 90;

std::vector<ParticleD*> Create_ParticlesD(Eigen::Vector3d origin, ObjectSize size_data);
std::vector<ParticleD*> Create_Particles_OneTetraD(Eigen::Vector3d origin, ObjectSize size_data);
std::vector<ParticleD*> Create_Particles_TwoTetraD(Eigen::Vector3d origin, ObjectSize size_data);
std::vector<ParticleD*> Create_Particles_Tri_prismD(Eigen::Vector3d origin, ObjectSize size_data);
void Draw_Mesh();
void Draw_Rotation(ObjectD* obj,float SinParam, float CosParam, float CameraVAngle, float CameraHAngle,double cameraZoom);
void Draw_Group_Grid(ObjectD* obj, float SinParam, float CosParam, float CameraVAngle, float CameraHAngle, double cameraZoom);
Eigen::Vector3d Calc_Draw_Grid(Eigen::Vector3d a, float SinParam, float CosParam, float CameraVAngle, float CameraHAngle, double cameraZoom);
void DrawRotationNew(ObjectD* obj, float SinParam, float CosParam, float CameraVAngle, float CameraHAngle, double cameraZoom);
void Draw_Group_Grid_New(ObjectD* obj, float SinParam, float CosParam, float CameraVAngle, float CameraHAngle, double cameraZoom);

// (x,y)の点を(mx,my)を中心にang角回転する
void rotate(float *x, float *y, const float ang, const float mx, const float my) {
	const float ox = *x - mx, oy = *y - my;
	*x = ox * cos(ang) + oy * sin(ang);
	*y = -ox * sin(ang) + oy * cos(ang);
	*x += mx;
	*y += my;
}
double MOUSE_RANGE;
double TIME_STEP;  //秒, second(リアルタイム時間なら0.03ぐらい)
double THRESHHOLD; // 1.0e+8;
int PBD_LOOP;      //制約条件における更新式(int)
double config;	  //修正値和による閾値の設定

float Gravity;  //重力加速度(m/s2)
bool useCRS;    //CRS(1)か昔の実装(0)か
bool mdiag;    //質量行列が対角(1)か否か(0)
int mSys; //質量行列がモデル全体の平均か(0)グループ別か(1)グループ和か(2)
bool loopite;    //更新数を決め打ち(1)か否か(0)
bool rollcamera;   //カメラが回転する(1)か否か(0)
bool fixedion;   //節点を固定する(1)か否か(0)
bool whichmethodused;//差分で反復する(1)かしないか(0)
bool useGMRES;//反復計算を非定常反復(GMRESなど)で解くか(1)かLUか(0)
int innerGMRES;//GMRESの内部反復の数(restart)
int outerGMRES;//GMRESの外部反復の数
bool usePreIte;//GMRESなどの反復法で前処理をするか(1)しないか(0)
bool useSparse; //行列をSparseで計算するか(1)しないか(0)
bool useUpRotate;//回転行列を更新するか(1)しないか(0)
bool useAPD;//回転行列をAPDで更新するか(1)しないか(0)
bool useUpBind;//拘束力を更新するか(1)しないか(0)
bool UseIniForce;//最初に力をかけるか(1)いなか(0)
int howmanytetras; //物体の形、箱型(3以上)、四面体一つ(1)、四面体二つ(2)
double M_damping;	  //バネダンパ係数(M,alpha)
double K_damping;	  //バネダンパ係数(K,Beta)
int dividedbyn;//x方向をn分割
int dividedbym;//y方向をm分割
int dividedbyl;//z方向をl分割　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　
double F_bind_coeff;//拘束力の係数
double F_bind_damping;//拘束力のダンパー係数

int Tri_prismnum;//三角柱の長さ
int GroupDividedby;//ブロックごとに切り分ける　

unsigned int xsize;	  //x方向のサイズ
unsigned int ysize;	  //y方向のサイズ
unsigned int zsize;	  //z方向のサイズ
double sidelength;	//一辺の長さ
double cipher;//スケーリング用（描画）
 // カメラの回転速度
#define CAMERA_ANGLE_SPEED		2.0
 // カメラの注視点の高さ
#define CAMERA_LOOK_AT_HEIGHT	400.0f
 // カメラと注視点の距離
#define CAMERA_LOOK_AT_DISTANCE	2150.0f
 // ラインを描く範囲
#define LINE_AREA_SIZE	10000.0f

// ラインの数
#define LINE_NUM 500

// Dataから数値を読む(Read numbers from Data.txt)
void BasicInformation() {

	std::ifstream ifs("Data.txt");
	std::string line;

	if (ifs.fail()) {
		std::cerr << "File Open Error" << std::endl;
	}
	while (getline(ifs, line)) {
		if (line[0] == '#')
			continue;
		if (line.find('=') == std::string::npos)
			continue;
		std::stringstream ss(line);
		std::string name;
		ss >> name;
		ss.ignore(line.size(), '=');
		std::cout << name << " = ";
		float temp;
		ss >> temp;
		std::cout << temp << ",";
		if (name == "MOUSE_RANGE")
			MOUSE_RANGE = temp;
		else if (name == "TIME_STEP")
			TIME_STEP = temp;
		else if (name == "THRESHHOLD")
			THRESHHOLD = temp;
		else if (name == "PBD_LOOP")
			PBD_LOOP = int(temp);
		else if (name == "config")
			config = temp;
		else if (name == "Gravity")
			Gravity = temp;
		else if (name == "useCRS")
			if (temp)  useCRS = TRUE;
			else  useCRS = FALSE;
		else if (name == "mdiag")
			if (temp)  mdiag = TRUE;
			else  mdiag = FALSE;
		else if (name == "mSys")
			mSys = int(temp);
		else if (name == "loopite")
			if (temp)  loopite = TRUE;
			else  loopite = FALSE;
		else if (name == "rollcamera")
			if (temp)  rollcamera = TRUE;
			else  rollcamera = FALSE;
		else if (name == "fixedion")
			if (temp)  fixedion = TRUE;
			else  fixedion = FALSE;
		else if (name == "howmanytetras")
			howmanytetras = int(temp);
		else if (name == "M_damping")
			M_damping = temp;
		else if (name == "K_damping")
			K_damping = temp;
		else if (name == "F_bind_coeff")
			F_bind_coeff = temp;
		else if (name == "F_bind_damping")
			F_bind_damping = temp;
		else if (name == "whichmethodused")
			whichmethodused = temp;
		else if (name == "useGMRES")
			useGMRES = temp;
		else if (name == "innerGMRES")
			innerGMRES = int(temp);
		else if (name == "outerGMRES")
			outerGMRES = int(temp);
		else if (name == "usePreIte")
			usePreIte = temp;
		else if (name == "useSparse")
			useSparse = temp;
		else if (name == "useUpRotate")
			useUpRotate = temp;
		else if (name == "useAPD")
			useAPD = temp;
		else if (name == "useUpBind")
			useUpBind = temp;
		else if (name == "UseIniForce")
			UseIniForce = temp;
		else if (name == "xsize")
			xsize = int(temp);
		else if (name == "ysize")
			ysize = int(temp);
		else if (name == "zsize")
			zsize = int(temp);
		else if (name == "sidelength")
			sidelength = double(temp);
		else if (name == "dividedbyn")
			dividedbyn = int(temp);
		else if (name == "dividedbym")
			dividedbym = int(temp);
		else if (name == "dividedbyl")
			dividedbyl = int(temp);
		else if (name == "Tri_prismnum")
			Tri_prismnum = int(temp);
		else if (name == "GroupDividedby")
			GroupDividedby = int(temp);
	}
	std::cout << std::endl;
}
//===========================================================================//
//===========================================================================//
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {
	ChangeWindowMode(TRUE), DxLib_Init(), SetDrawScreen(DX_SCREEN_BACK);
	//ウィンドウモードを非全画面にし、DXライブラリの初期化、裏画面設定

	if (!::AttachConsole(ATTACH_PARENT_PROCESS)) {
		::AllocConsole();           // WEBアプリケーションなので、標準入出力のためにコマンドプロンプトを表示する
									// コマンドプロンプトからの起動じゃない時は新たに開く
	}
	freopen("CON", "r", stdin);     // 標準入力の割り当て
	freopen("CON", "w", stdout);    // 標準出力の割り当て

	// Dataから数値を読む
	BasicInformation();				

	//	X+が右、Y+が下、Z+が手前
	std::vector<ObjectD*> obj;									//シミュレーションで生成するオブジェクト群
	ObjectSize size_data = { xsize, ysize, zsize ,sidelength };	//	モデルの大きさ(x,y,z方向), 1辺の長さ
	ObjectData almi = { 2.7e+03, 6.9e+10, 0.3, size_data };	    //	アルミの（密度、ヤング率、ポアソン比）
	ObjectData almin = { 2.7e+03, 1.0e+7, 0.49, size_data };	//	ゴムの（密度、ヤング率、ポアソン比）
	ObjectData gum = { 0.91e+03, 1.0e+06, 0.49, size_data };	//	ゴムの（密度、ヤング率、ポアソン比）
	ObjectData gum2 = { 0.91e+03, 1.0e+07, 0.49, size_data };	//	堅いゴムの（密度、ヤング率、ポアソン比）
	ObjectData gum3 = { 0.91e+03, 1.0e+08, 0.49, size_data };	//	堅いゴムの（密度、ヤング率、ポアソン比）
	ObjectData orihar = { 2.7e-03, 6.9e+12, 0.3, size_data };	//	幻想物体オリハルコンの（密度、ヤング率、ポアソン比）
    
    //スケールが変わっても、分割と描画ができるように係数を求める
	if (sidelength >= 40.0) {
		cipher = 40.0 / sidelength;
	}
	else if (sidelength>1) {
		cipher = 40.0 / sidelength;
	}
	else {
		cipher = 40.0 / sidelength;
	}
	//描画の見た目を小さくする係数cameraZoomの設定
	//Setting the coefficient cameraZoom to reduce the appearance of the drawing.
	double cameraZoom = 0.0;
	if (sidelength == 4.0) {
		if (xsize == 3) {
			cameraZoom = 3 * 40.0 / sidelength;
		}
		else if(xsize == 5){
			cameraZoom = 0.5 * 40.0 / sidelength;
		}
		else if (xsize == 7) {
			cameraZoom = 1.0 * 40.0 / sidelength;
		}
		else if (xsize == 9) {
			cameraZoom = 0.5 * 40.0 / sidelength;
		}
		else if (xsize == 11) {
			cameraZoom = 0.5 * 40.0 / sidelength;
		}
		else if (xsize == 13) {
			cameraZoom = 0.25 * 40.0 / sidelength;
		}
	}
	else if (sidelength >= 40.0) {
		if (xsize==3) {
			cameraZoom = 3 * 40.0 / sidelength;
		}
		else {
			cameraZoom = 0.5 * 40.0 / sidelength;
		}
	}
	else {
		cameraZoom = 5 * 40.0 / sidelength;
	}
	

	// オブジェクトの頂点群
	std::vector<ParticleD*> particles; 
	//オブジェクトのインスタンス
	ObjectD* o;						  
	
	particles = Create_ParticlesD(Eigen::Vector3d(0.0, 0.0, 0.0), size_data);
	//Vector3dが原点で、右下奥に直方体ができる。
	//オブジェクトのインスタンスを生成する
	//o = new UseBlockObjectDouble(particles, gum);  // モデルを従来のFEMでシミュレーションする
	o = new UseBlockObjectDouble(particles, gum3);

	obj.push_back(o);// 生成したオブジェクトをシミュレーションで使うオブジェクト群にpushする
	o->CalcMassMatrix();

	//DXライブラリで描画する色(白),随時設定する
	unsigned int string_color = GetColor(255, 255, 255);

	//事前計算の終了
	//End of precomputation
	std::cout << "Let's start!! " << std::endl;

	//===========================================================================//
	//描画の設定（Drawing Settings）
	//===========================================================================//
	float  CameraHAngle;
	float  CameraVAngle;
	float  SinParam;
	float  CosParam;
	VECTOR Position = VGet(0.0f, 0.0f, 0.0f);
	// カメラの向きを初期化
	CameraHAngle = -0.0f;
	CameraVAngle = 0.0f;

	SetCameraNearFar(0.1f, 1000.0f);//奥行0.1～1000までをカメラの描画範囲とする

	mt.setid(0);			//１ステップを測るstopwatchのidをセット
	mtUpdate.setid(1);		//１ステップ中の位置更新の計算時間を測るstopwatchのidをセット
	mtDraw.setid(3);		//描写する時間を測るstopwatchのidをセット
	mt.startMyTimer();		//１ステップを測るstopwatchをスタート

	std::ostringstream sstr;	//ウィンドウで出力する文字列1列目
	std::ostringstream sstr3;	//ウィンドウで出力する文字列3列目
	std::ostringstream sstr4;	//ウィンドウで出力する文字列4列目
	std::ostringstream sstr5;	//ウィンドウで出力する文字列5列目
	// 描画先を裏画面にする
	//SetDrawScreen(DX_SCREEN_BACK);

	// 背景の色を灰色にする
	//SetBackgroundColor(128, 128, 128);

	//===========================================================================//
	//実行時処理(run-time processing)
	//===========================================================================//

	while (ScreenFlip() == 0 && ProcessMessage() == 0 && ClearDrawScreen() == 0 && KeyBoard::gpUpdateKey() == 0) {
		//↑裏画面を表画面に反映   ↑ﾒｯｾｰｼﾞ処理			   ↑画面をｸﾘｱ               ↑keyが押されていない

		mt.endMyTimer();	//１ステップを測るstopwatchを終了
		mt.startMyTimer();	//１ステップを測るstopwatchをスタート
		SetCameraPositionAndTarget_UpVecY(VGet(320.0f, 150.0f, 0.0f), VGet(0, 0.0f, 0));
		int Mouse = GetMouseInput();	//マウスのボタンの入力状態
		
		// ZCSXキーでカメラの操作(回転)
		if (CheckHitKey(KEY_INPUT_C) == 1) {
			CameraHAngle += CAMERA_ANGLE_SPEED;
			if (CameraHAngle >= 360.0f) {
				CameraHAngle -= 360.0f;
			}
		}
		if (CheckHitKey(KEY_INPUT_Z) == 1) {
			CameraHAngle -= CAMERA_ANGLE_SPEED;
			if (CameraHAngle <= -360.0f) {
				CameraHAngle += 360.0f;
			}
		}
		if (CheckHitKey(KEY_INPUT_S) == 1) {
			CameraVAngle += CAMERA_ANGLE_SPEED;
			if (CameraVAngle >= 360.0f) {
				CameraVAngle -= 360.0f;
			}
		}
		if (CheckHitKey(KEY_INPUT_X) == 1) {
			CameraVAngle -= CAMERA_ANGLE_SPEED;
			if (CameraVAngle <= -360.0f) {
				CameraVAngle += 360.0f;
			}
		}
		if (CheckHitKey(KEY_INPUT_R) == 1) {
			CameraVAngle = 0.0f;
			CameraHAngle = -80.0f;
		}
		if (CheckHitKey(KEY_INPUT_T) == 1) {
			CameraVAngle += 90.0f;
			if (CameraVAngle >= 360.0f) {
				CameraVAngle -= 360.0f;
			}
		}
		//===========================================================================//
		// カメラの位置と向きを設定(Set the position and orientation of the camera.)
		//===========================================================================//
		VECTOR TempPosition1;
		VECTOR TempPosition2;
		VECTOR CameraPosition;
		VECTOR CameraLookAtPosition;

		// 注視点はキャラクターモデルの座標から CAMERA_LOOK_AT_HEIGHT 分だけ高い位置
		CameraLookAtPosition = Position;
		CameraLookAtPosition.y += CAMERA_LOOK_AT_HEIGHT;

		// カメラの位置はカメラの水平角度と垂直角度から算出

		// 最初に垂直角度を反映した位置を算出
		SinParam = sin(CameraVAngle / 180.0f * DX_PI_F);
		CosParam = cos(CameraVAngle / 180.0f * DX_PI_F);
		TempPosition1.x = 0.0f;
		TempPosition1.y = SinParam * CAMERA_LOOK_AT_DISTANCE;
		TempPosition1.z = -CosParam * CAMERA_LOOK_AT_DISTANCE;

		// 次に水平角度を反映した位置を算出
		SinParam = sin(CameraHAngle / 180.0f * DX_PI_F);
		CosParam = cos(CameraHAngle / 180.0f * DX_PI_F);
		TempPosition2.x = CosParam * TempPosition1.x - SinParam * TempPosition1.z;
		TempPosition2.y = TempPosition1.y;
		TempPosition2.z = SinParam * TempPosition1.x + CosParam * TempPosition1.z;

		// 算出した座標に注視点の位置を加算したものがカメラの位置
		CameraPosition = VAdd(TempPosition2, CameraLookAtPosition);

		// カメラの設定に反映する
		SetCameraPositionAndTarget_UpVecY(CameraPosition, CameraLookAtPosition);
		//SetCameraPositionAndTarget_UpVecY();



		//===========================================================================//
		// 実行時計算(run-time calculation)
		//===========================================================================//
		for (unsigned int i = 0; i < obj.size(); ++i) {	//現在,オブジェクトの数は一つ

			obj[i]->Timestep_Init();//1ステップの外力をリセット
			

			if (Mouse & MOUSE_INPUT_LEFT) {	//マウスの左ボタンが押されたら
				int x, y;
				GetMousePoint(&x, &y);		//マウスの座標を取得
			    //obj[i]->Set_Force(Eigen::Vector3d(double(x), double(y), 0.0));// 外力を一番端の点にかける
			}
			if (UseIniForce) {
				if (countup < 20000) {
					//obj[i]->Set_Force(Eigen::Vector3d(0.0, 1.0e+7, 0.0));// 外力を一番端の点にかける
				}
			}
			
			mtUpdate.startMyTimer();	//１ステップ中の位置更新の計算時間を測るstopwatchをスタート
			if (countup>50 && countup<100000) {
				
				//obj[i]->Update();	
			

				obj[i]->PBDCalculation();
			}

			mtUpdate.endMyTimer();		//１ステップ中の位置更新の計算時間を測るstopwatchを終了

			countup++;						 // 経過したステップ数を記録する
			if (countup<50) {				 // 50ステップのときの総計算時間を測定
				TIMERof50 += mtUpdate.getDt();// 1ステップ中の位置更新の計算時間を取得し加算
			}

			//実験用の変数
			int countup2 = 0;
			if (countup>50) {
				countup2 = countup;
			}
			if (countup > 100000) {
				countup2 = 100000;
			}
			//実験用の変数（終了）

			sstr3 << std::fixed;
			unsigned int string_color2 = GetColor(0, 255, 0); // 3列目の文字列の色を緑色にする
			sstr3 << "TimeStep is " << std::setprecision(5) << countup2 << ", Mean Of TIME(" << (TIME_STEP) << "s )is " << std::setprecision(5) << TIMERof50 / 50 << ", PBDLOOP" << obj[i]->ezolg;
			//3列目  ↑経過したタイムステップ数	     ↑50ステップにおける計算時間の平均値(s)      ↑各ステップのPBDにおける位置更新の回数
			DrawString(0, 32, sstr3.str().data(), string_color2);// 3列目を描画する
			sstr3.str("");				// 3列目のバッファをクリア
			mtDraw.startMyTimer();		// 描写する時間を測るstopwatchをスタート

			//座標系の表示
			Draw_Rotation(obj[0],SinParam,CosParam, CameraVAngle,CameraHAngle, cameraZoom);
			//グループごとに節点を描画
			
			Draw_Group_Grid(obj[0], SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			// 描写する時間を測るstopwatchを終了
			//DrawRotationNew(obj[0], SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			//Draw_Group_Grid_New(obj[0], SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			mtDraw.endMyTimer();		
			
		}
		sstr << std::fixed;
		sstr << std::setprecision(5) << "DeltaT=" << std::setprecision(5) << mt.getDt() << ",Update is " << std::setprecision(5) << mtUpdate.getDt() << ",Drawing is " << std::setprecision(5) << mtDraw.getDt()
			<< ",Others is " << std::setprecision(5) << mt.getDt() - mtUpdate.getDt() - mtDraw.getDt();
		//1列目	↑浮動小数点は4まで表示	↑１ステップにかかる時間　　　↑１ステップ中の位置更新の計算時間	　　↑描写する時間
		DrawString(0, 0, sstr.str().data(), string_color);// 1列目を描画する
		sstr.str("");	// 1列目のバッファをクリア
		//動く節点3の座標を追跡する(四面体要素1,2つのとき用)
		sstr4 << std::fixed;
		//Eigen::Vector3d Grige = Eigen::Vector3d::Zero();
		//Grige = obj[0]->groups[obj[0]->groups.size() - 1]->Get_Grid_In_Group(61);
		//Eigen::Vector3d Grige2 = Eigen::Vector3d::Zero();
		//Grige2 = obj[0]->groups[obj[0]->groups.size() - 1]->Get_Vel_In_Group(61);
		//sstr4 << Grige;
		//unsigned int string_color2 = GetColor(0, 255, 0); // 4列目の文字列の色を緑色にする
		////DrawString(0, 80, sstr4.str().data(), string_color2);// 4列目を描画する
		//sstr4.str("");				// 4列目のバッファをクリア
		////動く節点3の速度を追跡する
		//sstr5 << std::fixed;
		//sstr5 << Grige2;
		unsigned int string_color3 = GetColor(0, 255, 255); // 5列目の文字列の色を緑色にする
		//DrawString(100, 80, sstr5.str().data(), string_color3);// 5列目を描画する
		sstr5.str("");				// 5列目のバッファをクリア
		std::ofstream file("frame_times.txt");
	
		// 位置関係が分かるように地面にラインを描画する
		//Draw_Mesh();
	}
	::FreeConsole();	// コマンドプロンプトを閉じる
	DxLib_End();		// DXライブラリ終了処理
	return 0;
}

//===========================================================================//
//size_dataから頂点の集合を作成(Create a set of vertices from size_data)
//===========================================================================//

//Nodes for a box model
std::vector<ParticleD*> Create_ParticlesD(Eigen::Vector3d origin, ObjectSize size_data) {
	std::vector<ParticleD*> particles;//オブジェクトの頂点群
	unsigned int num = size_data.x_vertex_num * size_data.y_vertex_num * size_data.z_vertex_num; //頂点の数は各軸の要素の積 xyz
	std::vector<ParticleD*> p(num);//頂点の集合
	std::vector<Eigen::Vector3d> phys;//頂点の各座標
	for (unsigned int xi = 0; xi < size_data.x_vertex_num; ++xi) {
		for (unsigned int yi = 0; yi < size_data.y_vertex_num; ++yi) {
			for (unsigned int zi = 0; zi < size_data.z_vertex_num; ++zi) {
				phys.push_back(Eigen::Vector3d(origin.x() + xi * size_data.size, origin.y() + yi * size_data.size, origin.z() + zi * size_data.size));
			}
		}
	}
	//各頂点にidと座標を入れる
	for (unsigned int i = 0; i < num; i++) {
		p[i] = new ParticleD(phys[i]);
		p[i]->p_id = i;
		particles.push_back(p[i]);
		//if (i >= size_data.x_vertex_num * size_data.y_vertex_num * size_data.z_vertex_num - size_data.y_vertex_num * size_data.z_vertex_num) {
		//	p[i]->Set_Fixed(fixedion);// 一番端の頂点を固定する
		//}
		if (i < size_data.y_vertex_num * size_data.z_vertex_num) {
			p[i]->Set_Fixed(fixedion);// 一番端の頂点を固定する
		}
	}
	std::cout << "success create particle" << ":" << num << std::endl;//頂点の作成に成功したことを出力
	return particles;
}
//Nodes for one Tetra model
std::vector<ParticleD*> Create_Particles_OneTetraD(Eigen::Vector3d origin, ObjectSize size_data) {
	std::vector<ParticleD*> particles;//オブジェクトの頂点群
	unsigned int num = 4; //頂点の数は4で決め打ち
	std::vector<ParticleD*> p(num);//頂点の集合
	std::vector<Eigen::Vector3d> phys;//頂点の各座標
	phys.push_back(Eigen::Vector3d(origin.x(), origin.y(), origin.z()));
	phys.push_back(Eigen::Vector3d(origin.x() + size_data.size, origin.y(), origin.z()));
	phys.push_back(Eigen::Vector3d(origin.x() + size_data.size * 0.5, origin.y() + size_data.size * sqrt(3.0)*0.5, origin.z()));
	phys.push_back(Eigen::Vector3d(origin.x() + size_data.size * 0.5, origin.y() + size_data.size * sqrt(3.0)*0.5 / 3.0, origin.z() + size_data.size * sqrt(6.0) / 3.0));
	//各頂点にidと座標を入れる
	for (unsigned int i = 0; i < num; i++) {
		p[i] = new ParticleD(phys[i]);
		p[i]->p_id = i;
		particles.push_back(p[i]);
		if (i != 3) {
			p[i]->Set_Fixed(fixedion);// 一番端の頂点を固定する
		}
		std::cout << i << "particle is " << std::endl;
		std::cout << p[i]->Get_Grid() << std::endl;
	}
	//初期速度代入
	//particles[3]->Update_Velocity(Eigen::Vector3d(0.0, 1.0, 0.0));	
	std::cout << "success create particle" << ":" << num << std::endl;//頂点の作成に成功したことを出力
	return particles;
}
//Nodes for two Tetra model
std::vector<ParticleD*> Create_Particles_TwoTetraD(Eigen::Vector3d origin, ObjectSize size_data) {
	std::vector<ParticleD*> particles;//オブジェクトの頂点群
	unsigned int num = 5; //頂点の数は5で決め打ち
	std::vector<ParticleD*> p(num);//頂点の集合
	std::vector<Eigen::Vector3d> phys;//頂点の各座標
	phys.push_back(Eigen::Vector3d(origin.x(), origin.y(), origin.z()));
	phys.push_back(Eigen::Vector3d(origin.x() + size_data.size, origin.y(), origin.z()));
	phys.push_back(Eigen::Vector3d(origin.x() + size_data.size * 0.5, origin.y() + size_data.size * sqrt(3.0)*0.5, origin.z()));
	phys.push_back(Eigen::Vector3d(origin.x() + size_data.size * 0.5, origin.y() + size_data.size * sqrt(3.0)*0.5 / 3.0, origin.z() + size_data.size * sqrt(6.0) / 3.0));
	phys.push_back(Eigen::Vector3d(origin.x() + size_data.size * 0.5, origin.y() - size_data.size * sqrt(3.0)*0.5 / 3.0, origin.z() + size_data.size * sqrt(6.0) / 3.0));
	//各頂点にidと座標を入れる
	for (unsigned int i = 0; i < num; i++) {
		p[i] = new ParticleD(phys[i]);
		p[i]->p_id = i;
		particles.push_back(p[i]);
		if (i <= 2) {
			p[i]->Set_Fixed(fixedion);// 一番端の頂点を固定する
		}
		//if (i==4) {
			//p[i]->Set_Fixed(fixedion);
		//}
		std::cout << i << "particle is " << std::endl;
		std::cout << p[i]->Get_Grid() << std::endl;
	}
	//初期速度代入
	//particles[3]->Update_Velocity(Eigen::Vector3d(0.0, 1.0, 0.0));	
	std::cout << "success create particle of Two Tetra" << ":" << num << std::endl;//頂点の作成に成功したことを出力
	return particles;
}
//Nodes for Triangular columns model
std::vector<ParticleD*> Create_Particles_Tri_prismD(Eigen::Vector3d origin, ObjectSize size_data) {
	std::vector<ParticleD*> particles;//オブジェクトの頂点群
	unsigned int num = 3*(Tri_prismnum + 1); //頂点の数
	//std::cout << "Tri_particle = " << num << std::endl;
	std::vector<ParticleD*> p(num);//頂点の集合
	std::vector<Eigen::Vector3d> phys;//頂点の各座標
    for (int i = 0; i < Tri_prismnum + 1;i++) {
		phys.push_back(Eigen::Vector3d(origin.x(), origin.y(), origin.z() + size_data.size*i));
		phys.push_back(Eigen::Vector3d(origin.x() + size_data.size * sqrt(3.0)*0.5, origin.y() + size_data.size*0.5, origin.z() + size_data.size*i));
		phys.push_back(Eigen::Vector3d(origin.x() , origin.y() + size_data.size, origin.z() + size_data.size*i));
	}
	//各頂点にidと座標を入れる
	for (unsigned int i = 0; i < num; i++) {
		p[i] = new ParticleD(phys[i]);
		p[i]->p_id = i;
		particles.push_back(p[i]);
		if (i <= 2) {
			p[i]->Set_Fixed(fixedion);// 一番端の頂点3つを固定する
		}
		std::cout << i << "particle is " << std::endl;
		std::cout << p[i]->Get_Grid() << std::endl;
	}
	//初期速度代入
	//particles[14]->Update_Velocity(Eigen::Vector3d(0.0, 100.0, 0.0));	
	std::cout << "success create particle of Tri prism" << ":" << num << std::endl;//頂点の作成に成功したことを出力
	return particles;
}


//===========================================================================//
//節点の描画用の関数(Functions for drawing nodes.)
//===========================================================================//

//Draw_Mesh
void Draw_Mesh() {
	{
		int i = 0;
		VECTOR Pos1;
		VECTOR Pos2;

		SetUseZBufferFlag(TRUE);

		Pos1 = VGet(-LINE_AREA_SIZE / 2.0f, 0.0f, -LINE_AREA_SIZE / 2.0f);
		Pos2 = VGet(-LINE_AREA_SIZE / 2.0f, 0.0f, LINE_AREA_SIZE / 2.0f);
		for (i = 0; i <= LINE_NUM; i++)
		{
			DrawLine3D(Pos1, Pos2, GetColor(255, 255, 255));
			Pos1.x += LINE_AREA_SIZE / LINE_NUM;
			Pos2.x += LINE_AREA_SIZE / LINE_NUM;
		}

		Pos1 = VGet(-LINE_AREA_SIZE / 2.0f, 0.0f, -LINE_AREA_SIZE / 2.0f);
		Pos2 = VGet(LINE_AREA_SIZE / 2.0f, 0.0f, -LINE_AREA_SIZE / 2.0f);
		for (i = 0; i < LINE_NUM; i++)
		{
			DrawLine3D(Pos1, Pos2, GetColor(255, 255, 255));
			Pos1.z += LINE_AREA_SIZE / LINE_NUM;
			Pos2.z += LINE_AREA_SIZE / LINE_NUM;
		}

		SetUseZBufferFlag(FALSE);
	}
}
void DrawRotationNew(ObjectD* obj, float SinParam, float CosParam, float CameraVAngle, float CameraHAngle, double cameraZoom) {
	Eigen::VectorXd CentorVector = Eigen::Vector3d::Zero();
	Eigen::VectorXd PosVector = Eigen::VectorXd::Zero(3 * obj->particles.size());
	//Calculate Centor of mass of the object and save it in CentorVector
	for (int i = 0; i < obj->particles.size(); i++)
	{
		CentorVector[0] += obj->M_MatrixBody(3 * i, 3 * i) * obj->x_corrected(3 * i);
		CentorVector[1] += obj->M_MatrixBody(3 * i + 1, 3 * i + 1) * obj->x_corrected(3 * i + 1);
		CentorVector[2] += obj->M_MatrixBody(3 * i + 2, 3 * i + 2) * obj->x_corrected(3 * i + 2);
	}
	Eigen::Vector3d Draw_centor = Eigen::Vector3d::Zero();
	Draw_centor = Calc_Draw_Grid(CentorVector, SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
	DrawCircle(int(Draw_centor.x()), int(Draw_centor.y()), 3, RED, TRUE);
	std::ostringstream sstrC;
	sstrC << std::setprecision(3) << "  " << CentorVector[0] << "," << CentorVector[1] << "," << CentorVector[2];
	//DrawString(0, 49 + 16 * _g->tetra_group_id, sstrC.str().data(), RED);
	sstrC.str("");
	//座標系を出力
	double Range = 1.0 * cameraZoom;
	Eigen::Vector3d Xcard(Range, 0, 0);
	Eigen::Vector3d Ycard(0, Range, 0);
	Eigen::Vector3d Zcard(0, 0, Range);
	Xcard = CentorVector + Xcard;
	Ycard = CentorVector + Ycard;
	Zcard = CentorVector + Zcard;
	Eigen::Vector3d Draw_Xcard = Eigen::Vector3d::Zero();
	Eigen::Vector3d Draw_Ycard = Eigen::Vector3d::Zero();
	Eigen::Vector3d Draw_Zcard = Eigen::Vector3d::Zero();
	//回転したあとの座標系を計算


	
	Xcard = CentorVector + (Xcard - CentorVector);
	Ycard = CentorVector + (Ycard - CentorVector);
	Zcard = CentorVector + (Zcard - CentorVector);
	
	//描画用の座標に変換
	Draw_Xcard = Calc_Draw_Grid(Xcard, SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
	Draw_Ycard = Calc_Draw_Grid(Ycard, SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
	Draw_Zcard = Calc_Draw_Grid(Zcard, SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);

	//座標系を出力
	int color = GetColor(0, 255, 0);
	MyDrawLine3(Draw_centor, Draw_Xcard, RED);
	MyDrawLine3(Draw_centor, Draw_Ycard, GREEN);
	MyDrawLine3(Draw_centor, Draw_Zcard, BLUE);


	
}
//Draw_coordinate system of each Groups
void Draw_Rotation(ObjectD* obj,float SinParam, float CosParam, float CameraVAngle, float CameraHAngle,double cameraZoom) {
	int parti = 0;
	//まずは予測位置の重心を出力
	for (auto _g : obj->groups) {
		Eigen::Vector3d centor = Eigen::Vector3d::Zero();
		Eigen::VectorXd PosVector = Eigen::VectorXd::Zero(3 * _g->particle_num);
		//重心をSUMから求める
		
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			PosVector.block(3 * pi, 0, 3, 1) = _g->GroupGridVector.block(3 * pi, 0, 3, 1);
		}
		centor = _g->SUM_M_Matrix.block(0,0,3, 3 * _g->particle_num) * PosVector;
		
		if (countup == 100000) {
			std::cout <<_g->tetra_group_id << " is centor"<<std::setprecision(10) << centor << std::endl;
		}
		
		//重心を現在の座標から足して求める
		/*
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			//centor += _g->particles[pi]->Get_Mass() * _g->particles[pi]->Get_Initial_Pos();
			centor += _g->particles[pi]->Get_Mass() * _g->particles[pi]->Get_Grid();
		}
		centor = centor / _g->Get_Group_Mass();
		*/
		/*
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			centor += _g->particles[pi]->Get_Grid();
		}
		centor = centor / _g->particle_num;
		*/
		//これはできてそう
		/*
		for (unsigned int pi = 0; pi < _g->particle_num; pi++) {
			centor += _g->particles[pi]->Get_Initial_Pos();
		}
		centor = centor / _g->particle_num;
		*/


		//centor = centor + Eigen::Vector3d(100,100,100);
		//重心描画
		Eigen::Vector3d Draw_centor = Eigen::Vector3d::Zero();
		Draw_centor = Calc_Draw_Grid(centor,SinParam,CosParam ,CameraVAngle,CameraHAngle,cameraZoom);
		DrawCircle(int(Draw_centor.x()), int(Draw_centor.y()), 3, RED, TRUE);

		std::ostringstream sstrC;
		sstrC <<std::setprecision(3)<<"  " <<centor.x() <<","<< centor.y() << "," << centor.z();
		DrawString(0, 49+16 * _g->tetra_group_id, sstrC.str().data(), RED);
		sstrC.str("");
		//座標系を出力
		double Range = 1.0 * cameraZoom;
		Eigen::Vector3d Xcard(Range,0,0);
		Eigen::Vector3d Ycard(0,Range,0);
		Eigen::Vector3d Zcard(0,0,Range);
		Xcard = centor + Xcard;
		Ycard = centor + Ycard;
		Zcard = centor + Zcard;
		Eigen::Vector3d Draw_Xcard = Eigen::Vector3d::Zero();
		Eigen::Vector3d Draw_Ycard = Eigen::Vector3d::Zero();
		Eigen::Vector3d Draw_Zcard = Eigen::Vector3d::Zero();
		//回転したあとの座標系を計算
		
		Xcard = centor + _g->rotate_matrix * (Xcard - centor);
		Ycard = centor + _g->rotate_matrix * (Ycard - centor);
		Zcard = centor + _g->rotate_matrix * (Zcard - centor);
		/*
		Xcard = centor + (Xcard - centor);
		Ycard = centor + (Ycard - centor);
		Zcard = centor + (Zcard - centor);
		*/
		//描画用の座標に変換
		Draw_Xcard = Calc_Draw_Grid(Xcard, SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
		Draw_Ycard = Calc_Draw_Grid(Ycard, SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
		Draw_Zcard = Calc_Draw_Grid(Zcard, SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);

		//座標系を出力
		int color = GetColor(0, 255, 0);
		MyDrawLine3(Draw_centor, Draw_Xcard, RED);
		MyDrawLine3(Draw_centor, Draw_Ycard, GREEN);
		MyDrawLine3(Draw_centor, Draw_Zcard, BLUE);
	}
	/*
	for (auto _t : obj->tetras) {
		std::vector<ParticleD*> it = _t->Get_Particle();
		int color = WHITE;
		MyDrawLine3(it[0]->Get_Draw_Grid(), it[2]->Get_Draw_Grid(), color);
		MyDrawLine3(it[0]->Get_Draw_Grid(), it[1]->Get_Draw_Grid(), color);
		MyDrawLine3(it[1]->Get_Draw_Grid(), it[2]->Get_Draw_Grid(), color);
		MyDrawLine3(it[0]->Get_Draw_Grid(), it[3]->Get_Draw_Grid(), color);
		MyDrawLine3(it[2]->Get_Draw_Grid(), it[3]->Get_Draw_Grid(), color);
		MyDrawLine3(it[1]->Get_Draw_Grid(), it[3]->Get_Draw_Grid(), color);
	}
	*/
}
void Draw_Group_Grid_New(ObjectD* obj, float SinParam, float CosParam, float CameraVAngle, float CameraHAngle, double cameraZoom) {
	double Volume_i = 0.0;
	for (auto e : obj->tetras) {
		Eigen::Vector3d Draw_particle0 = Eigen::Vector3d::Zero();
		Eigen::Vector3d Draw_particle1 = Eigen::Vector3d::Zero();
		Eigen::Vector3d Draw_particle2 = Eigen::Vector3d::Zero();
		Eigen::Vector3d Draw_particle3 = Eigen::Vector3d::Zero();

		//Volume_i += _e->Calc_Volume(_g->Get_Grid_In_Group((_e->Get_Particle())[0]->p_id), _g->Get_Grid_In_Group((_e->Get_Particle())[1]->p_id), _g->Get_Grid_In_Group((_e->Get_Particle())[2]->p_id), _g->Get_Grid_In_Group((_e->Get_Particle())[3]->p_id));
		Draw_particle0 = Calc_Draw_Grid(obj->Get_Grid_In_Object(e->Get_Particle()[0]->p_id), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
		Draw_particle1 = Calc_Draw_Grid(obj->Get_Grid_In_Object(e->Get_Particle()[1]->p_id), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
		Draw_particle2 = Calc_Draw_Grid(obj->Get_Grid_In_Object(e->Get_Particle()[2]->p_id), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
		Draw_particle3 = Calc_Draw_Grid(obj->Get_Grid_In_Object(e->Get_Particle()[3]->p_id), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);

		DrawCircle(int(Draw_particle0.x()), int(Draw_particle0.y()), 3, WHITE, TRUE);
		DrawCircle(int(Draw_particle1.x()), int(Draw_particle1.y()), 3, WHITE, TRUE);
		DrawCircle(int(Draw_particle2.x()), int(Draw_particle2.y()), 3, WHITE, TRUE);
		DrawCircle(int(Draw_particle3.x()), int(Draw_particle3.y()), 3, WHITE, TRUE);

		MyDrawLine3(Draw_particle0, Draw_particle1, WHITE);
		MyDrawLine3(Draw_particle0, Draw_particle2, WHITE);
		MyDrawLine3(Draw_particle0, Draw_particle3, WHITE);
		MyDrawLine3(Draw_particle1, Draw_particle2, WHITE);
		MyDrawLine3(Draw_particle1, Draw_particle3, WHITE);
		MyDrawLine3(Draw_particle2, Draw_particle3, WHITE);
	}
}
//Draw node position  of each Groups 
void Draw_Group_Grid(ObjectD* obj, float SinParam, float CosParam, float CameraVAngle, float CameraHAngle, double cameraZoom) {
	//グループごとの座標を出力
	double Volume_i = 0.0;
	Eigen::Vector3d Draw_particle0 = Eigen::Vector3d::Zero();
	Eigen::Vector3d Draw_particle1 = Eigen::Vector3d::Zero();
	Eigen::Vector3d Draw_particle2 = Eigen::Vector3d::Zero();
	Eigen::Vector3d Draw_particle3 = Eigen::Vector3d::Zero();
	for (auto _g : obj->groups) {
		for (auto _e : _g->elements) {
			

			//Volume_i += _e->Calc_Volume(_g->Get_Grid_In_Group((_e->Get_Particle())[0]->p_id), _g->Get_Grid_In_Group((_e->Get_Particle())[1]->p_id), _g->Get_Grid_In_Group((_e->Get_Particle())[2]->p_id), _g->Get_Grid_In_Group((_e->Get_Particle())[3]->p_id));
			Draw_particle0 = Calc_Draw_Grid(_g->Get_Grid_In_Group((_e->Get_Particle())[0]->p_id), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			Draw_particle1 = Calc_Draw_Grid(_g->Get_Grid_In_Group((_e->Get_Particle())[1]->p_id), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			Draw_particle2 = Calc_Draw_Grid(_g->Get_Grid_In_Group((_e->Get_Particle())[2]->p_id), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			Draw_particle3 = Calc_Draw_Grid(_g->Get_Grid_In_Group((_e->Get_Particle())[3]->p_id), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);

			DrawCircle(int(Draw_particle0.x()), int(Draw_particle0.y()), 3, WHITE, TRUE);
			DrawCircle(int(Draw_particle1.x()), int(Draw_particle1.y()), 3, WHITE, TRUE);
			DrawCircle(int(Draw_particle2.x()), int(Draw_particle2.y()), 3, WHITE, TRUE);
			DrawCircle(int(Draw_particle3.x()), int(Draw_particle3.y()), 3, WHITE, TRUE);

			MyDrawLine3(Draw_particle0, Draw_particle1, WHITE);
			MyDrawLine3(Draw_particle0, Draw_particle2, WHITE);
			MyDrawLine3(Draw_particle0, Draw_particle3, WHITE);
			MyDrawLine3(Draw_particle1, Draw_particle2, WHITE);
			MyDrawLine3(Draw_particle1, Draw_particle3, WHITE);
			MyDrawLine3(Draw_particle2, Draw_particle3, WHITE);

		}

		/*
			Eigen::Vector3d Draw_particle0 = Eigen::Vector3d::Zero();
			Eigen::Vector3d Draw_particle1 = Eigen::Vector3d::Zero();
			Eigen::Vector3d Draw_particle2 = Eigen::Vector3d::Zero();
			Eigen::Vector3d Draw_particle3 = Eigen::Vector3d::Zero();
			//Draw_particle0 = Calc_Draw_Grid(_g->particles[0]->Get_Grid(), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			//Draw_particle1 = Calc_Draw_Grid(_g->particles[1]->Get_Grid(), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			//Draw_particle2 = Calc_Draw_Grid(_g->particles[2]->Get_Grid(), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			//Draw_particle3 = Calc_Draw_Grid(_g->particles[3]->Get_Grid(), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			Draw_particle0 = Calc_Draw_Grid(_g->GroupGridVector.block(0,0,3,1), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			Draw_particle1 = Calc_Draw_Grid(_g->GroupGridVector.block(3, 0, 3, 1), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			Draw_particle2 = Calc_Draw_Grid(_g->GroupGridVector.block(6, 0, 3, 1), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);
			Draw_particle3 = Calc_Draw_Grid(_g->GroupGridVector.block(9, 0, 3, 1), SinParam, CosParam, CameraVAngle, CameraHAngle, cameraZoom);

			DrawCircle(int(Draw_particle0.x()), int(Draw_particle0.y()), 3, RED, TRUE);
			DrawCircle(int(Draw_particle1.x()), int(Draw_particle1.y()), 3, RED, TRUE);
			DrawCircle(int(Draw_particle2.x()), int(Draw_particle2.y()), 3, RED, TRUE);
			DrawCircle(int(Draw_particle3.x()), int(Draw_particle3.y()), 3, RED, TRUE);

			//線を描画
			MyDrawLine3(Draw_particle0, Draw_particle1, WHITE);
			MyDrawLine3(Draw_particle0, Draw_particle2, WHITE);
			MyDrawLine3(Draw_particle0, Draw_particle3, WHITE);
			MyDrawLine3(Draw_particle1, Draw_particle2, WHITE);
			MyDrawLine3(Draw_particle1, Draw_particle3, WHITE);
			MyDrawLine3(Draw_particle2, Draw_particle3, WHITE);
		*/
	}
	//std::cout << Volume_i << ",";
}
//Calc drawing node position  of each Groups using Camera parameters
Eigen::Vector3d Calc_Draw_Grid(Eigen::Vector3d a, float SinParam, float CosParam, float CameraVAngle, float CameraHAngle, double cameraZoom) {
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	double gridx1;
	double gridy1;
	double gridz1;
	double gridx2;
	double gridy2;
	double gridz2;
	SinParam = sin(CameraVAngle / 180.0f * DX_PI_F);
	CosParam = cos(CameraVAngle / 180.0f * DX_PI_F);
	gridx1 = a.x() * cameraZoom;
	gridy1 = a.y() * cameraZoom * CosParam - a.z() * cameraZoom * SinParam;
	gridz1 = a.y() * cameraZoom * SinParam + a.z() * cameraZoom * CosParam;
	// カメラの角度に合わせて移動ベクトルを回転してから加算
	SinParam = sin(CameraHAngle / 180.0f * DX_PI_F);
	CosParam = cos(CameraHAngle / 180.0f * DX_PI_F);
	gridx2 = gridx1 * CosParam - gridz1 * SinParam;
	gridy2 = gridy1;
	gridz2 = gridx1 * SinParam + gridz1 * CosParam;
	temp[0] = gridx2 + 100.0;
	temp[1] = gridy2 + 200.0;
	temp[2] = gridz2;

	return temp;
}
