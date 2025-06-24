//
// Auxiliary common functions
//
// Author: Jose Javier Gonzalez Ortiz, 6/6/2017
//

#ifndef UTILS_H
#define UTILS_H

#include <cstring>

void stripcrlf(char *tok);

void stripcrlf(char *tok)
{
   int l = strlen(tok);
   if (l > 0 && tok[l - 1] == '\n') {
      if (l > 1 && tok[l - 2] == '\r')
         tok[l - 2] = '\0';
      else
         tok[l - 1] = '\0';
   }
}

#endif // UTILS_H
