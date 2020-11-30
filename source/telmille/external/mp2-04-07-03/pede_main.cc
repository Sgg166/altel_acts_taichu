#include <string>
#include <cstdio>
#include <cstdint>
#include <cinttypes>

#include "pede.h"

using std::int32_t;
using std::int64_t;
using std::size_t;


int main(int argc, char** args){

  __mpmod_MOD_lunmon=0; // !     millepede monitoring file
  __mpmod_MOD_lunlog=8; // !     millepede.log file

  __mpmod_MOD_lvllog=1;
  std::string logfile_name("millepede.log");
  mvopen_(&__mpmod_MOD_lunlog, logfile_name.c_str(), logfile_name.size());
  std::string running_str("Still running or crashed");
  int32_t error_code_n1 =  -1;
  peend_(&error_code_n1,running_str.c_str(), running_str.size());

  // !     read command line and text files
  //TODO: filetc has internal  CALL getarg(i,text) and open streering file
  filetc_();   // command line and steering file analysis
  filetx_();   // read text files

  if(__mpmod_MOD_icheck > 0){
    throw;//TODO: not yet implenmented
  }
  __mpmod_MOD_lvllog=__mpmod_MOD_mprint; // export print level

  if(__mpmod_MOD_memdbg > 0){
    __mpdalc_MOD_printflagalloc=1; // debug memory management
  }

  if(__mpmod_MOD_ncache < 0){
    __mpmod_MOD_ncache=25000000*__mpmod_MOD_mthrd;  // default cache size (100 MB per thread)
  }


  int64_t rows = 6;
  int64_t cols =__mpmod_MOD_mthrdr;
  preparereadbufferinfo_(&rows, &cols);

  int32_t lunhis = 7;
  std::string hisfile_name("millepede.his");
  mvopen_(&lunhis, hisfile_name.c_str(), hisfile_name.size());
  hmplun_(&lunhis);// unit for histograms
  gmplun_(&lunhis);//  unit for xy data

  int32_t lundeb = 1;
  if(__mpmod_MOD_nrecpr != 0 || __mpmod_MOD_nrecp2 != 0){
    std::string debugfile_name("mpdebug.txt");
    mvopen_(&lundeb, debugfile_name.c_str(), debugfile_name.size());
  }

  loop1_();
  loop2_();
  if(__mpmod_MOD_chicut != 0.0){
    std::fprintf(stdout, "Chi square cut equiv 3 st.dev applied ...");
    std::fprintf(stdout, "in  first iteration with factor %f",__mpmod_MOD_chicut);
    std::fprintf(stdout, "in second iteration with factor %f",__mpmod_MOD_chirem);
    std::fprintf(stdout, "(reduced by sqrt in next iterations)");
  }

  if(__mpmod_MOD_lhuber != 0){
   std::fprintf(stdout, "Down-weighting of outliers in %d iterations" , __mpmod_MOD_lhuber);
   std::fprintf(stdout, "Cut on downweight fraction %f",__mpmod_MOD_dwcut);
  }

  std::string iter_str("Iteration");
  mstart_(iter_str.c_str(), iter_str.size());   // Solution module starting
  xloopn_(); // all methods;
  //TODO: remove prtglo_ at the end of xloopn_ , millepede.res;


/*
  if(nloopn > 2  && nhistp != 0){ //      ! last iteration
    hmprnt_(3);  //! scaled residual of single measurement (with global deriv.)
    hmprnt_(12); //! scaled residual of single measurement (no global deriv.)
    hmprnt_(4);  //! chi^2/Ndf
  }

  if(nloopn > 2){
    hmpwrt_(3);
    hmpwrt_(12);
    hmpwrt_(4);
    gmpwrt_(4); //! location, dispersion (res.) as a function of record nr
    if(nloopn <= lfitnp){
      hmpwrt_(13);
      hmpwrt_(14);
      gmpwrt_(5);
    }
  }

  if(nhistp != 0){
    gmprnt_(1);
    gmprnt_(2);
  }

  gmpwrt_(1);             //! output of xy data
  gmpwrt_(2);             //! output of xy data
*/
}
