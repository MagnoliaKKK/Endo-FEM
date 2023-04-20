#ifndef UT_CLAPACK_H
#define UT_CLAPACK_H

//#include <Springhead.h>


#define BOOST_NUMERIC_BINDINGS_USE_CLAPACK
#pragma warning(push)
#pragma warning(disable:4267)
#pragma warning(disable:4005)
#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/bindings/lapack/driver/sygv.hpp>
#include <boost/numeric/bindings/lapack/driver/sygvx.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <boost/numeric/bindings/lapack/driver/gels.hpp>
#include <boost/numeric/bindings/lapack/driver/gelsd.hpp>
#include <boost/numeric/bindings/lapack/driver/gesdd.hpp>
#include <boost/numeric/bindings/noop.hpp>
#include <boost/numeric/bindings/ublas/banded.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#pragma warning(pop)

#ifdef TRACE		// Trace
#if (_MSC_VER == 1500)    // Visual Studio 2008
#ifdef _WIN64 
# pragma comment(lib, "LIBF2C14.0x64.lib")
# pragma comment(lib, "BLAS14.0x64.lib")
# pragma comment(lib, "CLAPACK14.0x64.lib")
#else
# pragma comment(lib, "LIBF2C14.0Win32.lib")
# pragma comment(lib, "BLAS14.0Win32.lib")
# pragma comment(lib, "CLAPACK14.0Win32.lib")
#endif
#else			    // after Visual Studio 2010
#ifdef _WIN64 
# pragma comment(lib, "LIBF2C14.0Tx64.lib")
# pragma comment(lib, "BLAS14.0Tx64.lib")
# pragma comment(lib, "CLAPACK14.0Tx64.lib")
#else
# pragma comment(lib, "LIBF2C14.0TWin32.lib")
# pragma comment(lib, "BLAS14.0TWin32.lib")
# pragma comment(lib, "CLAPACK14.0TWin32.lib")
#endif
#endif
#else /* TRACE */
#ifdef _DEBUG
#ifdef _DLL		// Debug (former DebugDll)
#ifdef _WIN64 
# pragma comment(lib, "LIBF2C14.0MDx64.lib")
# pragma comment(lib, "BLAS14.0MDx64.lib")
# pragma comment(lib, "CLAPACK14.0MDx64.lib")
#else
# pragma comment(lib, "LIBF2C14.0MDWin32.lib")
# pragma comment(lib, "BLAS14.0MDWin32.lib")
# pragma comment(lib, "CLAPACK14.0MDWin32.lib")
#endif
#else		// (former Debug)
#ifdef _WIN64 
# pragma comment(lib, "LIBF2C14.0Dx64.lib")
# pragma comment(lib, "BLAS14.0Dx64.lib")
# pragma comment(lib, "CLAPACK14.0Dx64.lib")
#else
# pragma comment(lib, "LIBF2CD.lib")
# pragma comment(lib, "BLASD.lib")
# pragma comment(lib, "LAPACKD.lib")
#endif
#endif
#else /* _DEBUG */
#ifdef _DLL		// Release (former ReleaseDll)
#ifdef _WIN64 
# pragma comment(lib, "LIBF2C14.0Mx64.lib")
# pragma comment(lib, "BLAS14.0Mx64.lib")
# pragma comment(lib, "CLAPACK14.0Mx64.lib")
#else
# pragma comment(lib, "LIBF2C14.0MWin32.lib")
# pragma comment(lib, "BLAS14.0MWin32.lib")
# pragma comment(lib, "CLAPACK14.0MWin32.lib")
#endif
#else		// (former Release)
#ifdef _WIN64 
# pragma comment(lib, "LIBF2C14.0x64.lib")
# pragma comment(lib, "BLAS14.0x64.lib")
# pragma comment(lib, "CLAPACK14.0x64.lib")
#else
# pragma comment(lib, "LIBF2C.lib")
# pragma comment(lib, "BLAS.lib")
# pragma comment(lib, "LAPACK.lib")
#endif
#endif
#endif /* _DEBUG */
#endif /* TRACE */

/*
/ リンクするlibファイルは_cdeclで呼び出し
/ _fastcallでリンクしたい場合はpringhead2\src\boost\numeric\bindings\lapack\lapack.h
/ の関数をすべて_cdecl呼び出しにすること
*/

#define ublas boost::numeric::ublas
#define bindings boost::numeric::bindings
#define lapack bindings::lapack


#endif
#pragma once
