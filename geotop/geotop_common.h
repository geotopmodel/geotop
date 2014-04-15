#ifndef __GEOTOP_COMMON__
#define __GEOTOP_COMMON__

/*----------   Global variables  ------------*/

#include "times.h"
#include <string>
#include <vector>

namespace geotop
{

    namespace common
    {

        class Variables
        {
        public:
            static std::string WORKING_DIRECTORY ;

            static  std::vector<std::string> hpnt ;
            static std::vector<std::string> hbsn ;
            static std::vector<std::string> hsnw ;
            static std::vector<std::string> hglc ;
            static std::vector<std::string> hsl ;

            static const double TZ;

            static TInit *UV;

            static std::string logfile;
            static std::vector<std::string> files;

            static long Nl ;
            static long Nr ;
            static long Nc ;

            static double t_meteo ;
            static double t_energy ;
            static double t_water ;
            static double t_sub ;
            static double t_sup ;
            static double t_blowingsnow ;
            static double t_out ;

            static double **odpnt ;
            static double **odp ;

            static long *opnt ;
            static long nopnt ;
            static short *ipnt ;
            static short *ibsn ;

            static double *odbsn  ;
            static double *odb  ;
            static long *obsn  ;
            static long nobsn  ;

            static long *osnw  ;
            static long nosnw  ;

            static long *oglc  ;
            static long noglc  ;

            static long *osl  ;
            static long nosl  ;

            static FILE *ffbas  ;
            static FILE *ffpoint  ;
            static FILE *ffT  ;
            static FILE *ffTav  ;
            static FILE *ffpsi  ;
            static FILE *ffpsitot  ;
            static FILE *ffliq  ;
            static FILE *ffliqav  ;
            static FILE *ffice  ;
            static FILE *fficeav  ;
            static FILE *ffsnowT  ;
            static FILE *ffsnowl  ;
            static FILE *ffsnow  ;
            static FILE *ffsnowi  ;
            static FILE *ffsnowd  ;
            static FILE *ffglac  ;

            static long i_sim  ;
            static long i_run  ;
            static long i_sim0 ;
            static long i_run0  ;

            //char *SuccessfulRunFile, *FailedRunFile;

            static std::string SuccessfulRunFile ;
            static std::string FailedRunFile ;

            static time_t start_time ;
            static double elapsed_time  ;
            static double elapsed_time_start  ;
            static double cum_time  ;
            static double max_time  ;
        } ;

    } // end namespace common
} // end namespace geotop

#endif //__GEOTOP_COMMON__
