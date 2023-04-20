#ifndef _BASE
#define _BASE

//#define PY_SSIZE_T_CLEAN
//#include <python.h>

#include <vector>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <map>
#include <set>
#include <fstream>
#include <windows.h>
#include <iomanip>
#include <fenv.h>

#include "DxLib.h"
#include "Eigen/Dense"
#include "Eigen/LU"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "StopwatchTimer.h"
#include "Eigen/SVD"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

//typedef unsigned int size_t;

//Dxlib 色関係定数
const unsigned int WHITE = GetColor(255, 255, 255);
const unsigned int BLACK = GetColor(0, 0, 0);
const unsigned int RED = GetColor(255, 0, 0);
const unsigned int GREEN = GetColor(0, 255, 0);
const unsigned int BLUE = GetColor(0, 0, 255);

//シミュレーション定数
extern double MOUSE_RANGE;
extern double TIME_STEP;  //秒, second(リアルタイム時間なら0.03ぐらい)
extern double THRESHHOLD; // 1.0e+8;
extern int PBD_LOOP;      //制約条件における更新式(int)
extern double config;	  //修正値和による閾値の設定

//物理量定数
const  double PI = 3.14159265358979265358979;
extern float Gravity;  //重力加速度(m/s2)
extern bool UseIniForce;//最初に力をかけるか(1)いなか(0)
extern bool useCRS;    //CRS(1)か昔の実装(0)か
extern bool mdiag;    //質量行列が対角(1)か否か(0)
extern int mSys; //質量行列がモデル全体の平均か(0)グループ別か(1)グループ和か(2)
extern bool loopite;    //更新数を決め打ち(1)か否か(0)
extern bool rollcamera;    //カメラが回転する(1)か否か(0)
extern bool fixedion;   //節点を固定する(1)か否か(0)
extern bool whichmethodused;//省略法(1)か反復法か(0)
extern bool useGMRES;//反復計算を非定常反復(GMRESなど)で解くか(1)かLUか(0)
extern bool useGMRES;//反復計算を非定常反復(GMRESなど)で解くか(1)かLUか(0)
extern int innerGMRES;//GMRESの内部反復の数(restart)
extern int outerGMRES;//GMRESの外部反復の数
extern bool usePreIte;//GMRESなどの反復法で前処理をするか(1)しないか(0)
extern bool useSparse; //行列をSparseで計算するか(1)しないか(0)
extern bool useUpRotate;//回転行列を更新するか(1)しないか(0)
extern bool useAPD;//回転行列をAPDで更新するか(1)しないか(0)
extern bool useUpBind;//拘束力を更新するか(1)しないか(0)

extern int howmanytetras; //物体の形、箱型(3以上)、四面体一つ(1)、四面体二つ(2)
extern double M_damping;	  //バネダンパ係数
extern double K_damping;	  //バネダンパ係数
extern double F_bind_coeff;	  //反復法における拘束力の係数
extern double F_bind_damping; //反復法における拘束力のダンパー係数
extern int dividedbyn;//x方向をn分割
extern int dividedbym;//y方向をm分割
extern int dividedbyl;//z方向をl分割
extern int Tri_prismnum;//三角柱の長さ
extern int GroupDividedby;//ブロックごとに切り分ける

extern int unsigned xsize;	  //x方向のサイズ
extern int unsigned ysize;	  //y方向のサイズ
extern int unsigned zsize;	  //z方向のサイズ
extern int immortallength;//一辺の長さ(仮)
extern double sidelength; //一辺の長さ
extern double cipher; //一辺の長さ

static int countup;
static double TIMERof50;

struct Circle{
	Eigen::Vector3f center;  // 中心座標
	long double radius;		 // 半径
};

struct CircleD {
	Eigen::Vector3d center;  // 中心座標
	long double radius;		 // 半径
};

struct ObjectSize{
	unsigned int x_vertex_num ;
	unsigned int y_vertex_num ;
	unsigned int z_vertex_num ;
	double size ;
};

struct ObjectData{
	double density; //密度
	double young;	//ヤング率
	double poisson; //ポアソン比
	ObjectSize os;
};

const static bool operator<(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2){
	return v1.x() != v2.x() ? v1.x() < v2.x() :
		v1.y() != v2.y() ? v1.y() < v2.y() :
		v1.z() < v2.z();
}
const static bool operator<(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
	return v1.x() != v2.x() ? v1.x() < v2.x() :
		v1.y() != v2.y() ? v1.y() < v2.y() :
		v1.z() < v2.z();
}

/*
//前処理してるかんすう
class MatrixReplacement;
using Eigen::SparseMatrix;
namespace Eigen {
	namespace internal {
		// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
		template<>
		struct traits<MatrixReplacement> : public Eigen::internal::traits<Eigen::SparseMatrix<double> >
		{};
	}
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement:public Eigen::EigenBase<MatrixReplacement> {
public:
	// Required typedefs, constants, and method:
	typedef double Scalar;
	typedef double RealScalar;
	typedef int StorageIndex;
	enum {
		ColsAtCompileTime = Eigen::Dynamic,
		MaxColsAtCompileTime = Eigen::Dynamic,
		IsRowMajor = false
	};

	Index rows() const { return mp_mat->rows(); }
	Index cols() const { return mp_mat->cols(); }

	template<typename Rhs>
	Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
		return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
	}

	// Custom API:
	MatrixReplacement() : mp_mat(0) {}

	void attachMyMatrix(const Eigen::SparseMatrix<double> &mat) {
		mp_mat = &mat;
	}
	const Eigen::SparseMatrix<double> my_matrix() const { return *mp_mat; }

private:
	const Eigen::SparseMatrix<double> *mp_mat;
};

// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
	namespace internal {

		template<typename Rhs>
		struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
			: generic_product_impl_base<MatrixReplacement, Rhs, generic_product_impl<MatrixReplacement, Rhs> >
		{
			typedef typename Product<MatrixReplacement, Rhs>::Scalar Scalar;

			template<typename Dest>
			static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
			{
				// This method should implement "dst += alpha * lhs * rhs" inplace,
				// however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
				assert(alpha == Scalar(1) && "scaling is not implemented");
				EIGEN_ONLY_USED_FOR_DEBUG(alpha);

				// Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
				// but let's do something fancier (and less efficient):
				for (Index i = 0; i<lhs.cols(); ++i)
					dst += rhs(i) * lhs.my_matrix().col(i);
			}
		};
	}
}
*/
/*
template<unsigned int NI, unsigned int NJ, typename Blocks>
class MatrixReplacement :
	public Eigen::EigenBase< MatrixReplacement< NI, NJ, Blocks > >
{
public:
	// types
	typedef first_block_type::Scalar Scalar;
	typedef Scalar                   RealScalar;
	typedef size_t                   Index;
	typedef int                      StorageIndex;
	typedef unspecified              InnerIterator;

	enum @0 { ColsAtCompileTime = = Eigen::Dynamic,
		RowsAtCompileTime = = Eigen::Dynamic,
		MaxColsAtCompileTime = = Eigen::Dynamic,
		MaxRowsAtCompileTime = = Eigen::Dynamic, IsRowMajor = = false };

	// construct/copy/destruct
	MatrixReplacement(const Blocks &);
	MatrixReplacement(Blocks &&);

	// public member functions
	Index rows() const;
	Index cols() const;
	Index innerSize() const;
	Index outerSize() const;
	void resize(Index, Index);
	Scalar coeff(const Index, const Index) const;
	template<typename Rhs>
	Eigen::Product< MatrixReplacement, Rhs, Eigen::AliasFreeProduct >
		operator*(const Eigen::MatrixBase< Rhs > &) const;
	template<unsigned int I, unsigned int J>
	const std::tuple_element< I *NJ + J, Blocks >::type & get_kernel() const;
	const std::tuple_element< 0, Blocks >::type & get_first_kernel() const;
	template<unsigned int I, unsigned int J>
	std::tuple_element< I *NJ + J, Blocks >::type & get_kernel();
	std::tuple_element< 0, Blocks >::type & get_first_kernel();
	template<typename Derived>
	void assemble(Eigen::DenseBase< Derived > &) const;
	template<int _Options, typename _StorageIndex>
	void assemble(Eigen::SparseMatrix< Scalar, _Options, _StorageIndex > &);
	template<std::size_t... I> Index rows_impl(unspecified) const;
	template<std::size_t... J> Index cols_impl(unspecified) const;
	template<int I> Index start_col() const;
	template<int I> Index size_col() const;
	template<int I> Index start_row() const;
	template<int I> Index size_row() const;
	template<typename block_type>
	Scalar coeff_impl_block(const Index, const Index, const block_type &) const;
	template<std::size_t... I>
	Scalar coeff_impl(const Index, const Index, unspecified) const;
	template<typename Block>
	void assemble_block_impl(const size_t, const size_t,
		std::vector< Eigen::Triplet< Scalar >> &,
		const Block &) const;
	template<typename Block, typename Derived>
	void assemble_block_impl(const Eigen::MatrixBase< Derived > &,
		const Block &) const;
	template<std::size_t... I>
	void assemble_impl(std::vector< Eigen::Triplet< Scalar >> &, unspecified) const;
	template<typename Derived, std::size_t... I>
	void assemble_impl(Eigen::DenseBase< Derived > &, unspecified) const;

	// public data members
	Blocks m_blocks;
};
*/
#endif