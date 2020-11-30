#include <cstdint>
#include <cstddef>
using std::int32_t;
using std::int64_t;
using std::size_t;


// mps float
// mpd fdouble
// mpi int32_t
// mpl int64_t
// nm mphistab.o | grep hmpwrt
// https://docs.oracle.com/cd/E19422-01/819-3685/11_cfort.html
// http://www.mathcs.emory.edu/~cheung/Courses/561/Syllabus/5-Fortran/param-passing.html

#if __GNUC__ > 7
typedef size_t fortran_charlen_t;
#else
typedef int fortran_charlen_t;
#endif

extern "C" {
  // __bar_MOD_foo
  extern int32_t __mpmod_MOD_ictest;
  extern int32_t __mpmod_MOD_lunmon;
  extern int32_t __mpmod_MOD_lunlog;
  extern int32_t __mpmod_MOD_lvllog;
  extern int32_t __mpmod_MOD_icheck;
  extern int32_t __mpmod_MOD_mprint;
  extern int32_t __mpmod_MOD_memdbg;
  extern int32_t __mpmod_MOD_ncache;
  extern int32_t __mpmod_MOD_mthrd;
  extern int32_t __mpmod_MOD_mthrdr;
  extern int32_t __mpmod_MOD_nrecpr;
  extern int32_t __mpmod_MOD_nrecp2;
  extern float __mpmod_MOD_chicut;
  extern float __mpmod_MOD_chirem;
  extern int32_t __mpmod_MOD_lhuber;
  extern float __mpmod_MOD_dwcut;
  extern int32_t __mpdalc_MOD_printflagalloc;
  void loop1_();
  void loop2_();
  void mvopen_(int32_t* unit, const char* name, fortran_charlen_t);
  void peend_(int32_t* n, const char* str, fortran_charlen_t);
  void filetc_();
  void filetx_();
  void hmplun_(int32_t* unit);
  void gmplun_(int32_t* unit);
  void xloopn_();
  void mstart_(const char* str, fortran_charlen_t);
  void hmprnt_(int32_t* ih);
  void hmpwrt_(int32_t* ih);
  void gmpwrt_(int32_t* gh);

  // extern int32_t* __mpmod_MOD_readBufferInfo; // dynamic varialbe is not able to exprot
  // void __mpdalc_MOD_mpallociarr(int32_t* array, int64_t* rows, int64_t* cols,  const char* str, fortran_charlen_t);
  void preparereadbufferinfo_(int64_t* rows, int64_t* cols); // workaround
}
