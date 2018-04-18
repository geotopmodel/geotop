#include "turtle.h"
#include <iostream>
/*-----------------------------------------------------------------------*/
void t_error(const char *error_text)
/* Error handling */
{

  std::cerr << "\nError::Run time error\n"
            <<"Error:: " << error_text << "\n"
            << "........exiting...\n" << std::endl;
  exit(1);
}

