#pragama comment(lib,"LIBF2C14.0Win32.lib")

exturn "C"{
	void gbsv(const fortran_int_t n, const fortran_int_t kl,
		const fortran_int_t ku, const fortran_int_t nrhs, float* ab,
		const fortran_int_t ldab, fortran_int_t* ipiv, float* b,
		const fortran_int_t ldb);
}