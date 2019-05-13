#ifndef _GEOTOP_METEOIO_H
#define _GEOTOP_METEOIO_H

#include <meteoio/MeteoIO.h>
#include <string>

class MeteoioWrapper {
  
 public:
  std::string cfgfile;

  void which_cfgfile(){
    std::cout << "cfgfile = " << cfgfile << std::endl;
  }

    MeteoioWrapper(const std::string _cfgfile):
            cfgfile{_cfgfile} {}

 // MeteoioWrapper(){};
};

#endif // _GEOTOP_METEOIO_H
