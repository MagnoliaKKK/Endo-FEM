#ifndef _MYDEBUG
#define _MYDEBUG

#define _DEBUGMODE
#include "Base.h"
#define _CREATE_GROUP

#pragma warning(disable : 4996)
static void PrintV(int x, int y, const Eigen::Vector2d& v){
	DrawFormatString(x, y, WHITE, "[%9.3f,%9.3f]", v(0), v(1));
}
static void PrintM(int x, int y, const Eigen::Matrix2d& v){
	DrawFormatString(x, y, WHITE, "%9.3f %9.3f \n%9.3f %9.3f", v(0), v(1));
}
static void MyDrawLine(const Eigen::Vector3d& o, const Eigen::Vector3d& v) {
	DrawLine((int)o[0], (int)o[1], (int)o[0] + (int)v[0], (int)o[1] + (int)v[1], WHITE);
}
static void MyDrawLine2(const Eigen::Vector3d& o, const Eigen::Vector3d& v) {
	DrawLine((int)o[0], (int)o[1], (int)v[0], (int)v[1], GREEN);
}
static void MyDrawLine3(const Eigen::Vector3d& o, const Eigen::Vector3d& v, int color) {
	DrawLine((int)o[0], (int)o[1], (int)v[0], (int)v[1], color);
}
static void MyDrawLine31(const Eigen::Vector3d& o, const Eigen::Vector3d& v, int color) {
	DrawLine((int)o[0], (int)o[2], (int)v[0], (int)v[2], color);
}
static void PrintV(const Eigen::VectorXd& v) {
	std::cout << v << std::endl << std::endl;
}
static void PrintM(const Eigen::MatrixXd& m) {
	std::cout << m << std::endl << std::endl;
}

static void MyDrawLine(const Eigen::Vector3f& o, const Eigen::Vector3f& v){
	DrawLine((int)o[0], (int)o[1], (int)o[0] + (int)v[0], (int)o[1] + (int)v[1], WHITE);
}
static void MyDrawLine2(const Eigen::Vector3f& o, const Eigen::Vector3f& v){
	DrawLine((int)o[0], (int)o[1], (int)v[0], (int)v[1],GREEN);
}
static void MyDrawLine3(const Eigen::Vector3f& o, const Eigen::Vector3f& v, int color){
	DrawLine((int)o[0], (int)o[1], (int)v[0], (int)v[1], color);
}
static void MyDrawLine31(const Eigen::Vector3f& o, const Eigen::Vector3f& v, int color) {
	DrawLine((int)o[0], (int)o[2], (int)v[0], (int)v[2], color);
}
static void PrintV(const Eigen::VectorXf& v){
	std::cout << v << std::endl << std::endl;
}
static void PrintM(const Eigen::MatrixXf& m){
	std::cout << m << std::endl << std::endl;
}

template<typename T>
static void P(T t){
	std::cout << t << std::endl;
}
#endif