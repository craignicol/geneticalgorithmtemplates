/***************************************************************************
                          testhappyfaces.cpp  -  
			  for solving the happy faces puzzle
                             -------------------
    begin                : Tue Mar 30 2004
    copyright            : (C) 2004 by Craig Nicol
    email                : craig.nicol@gmail.com
 ***************************************************************************/
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "GAhappyfaces.h"

int main(int argc, char **argv)
{
  HappyFaces test;
  int rseed;

  const int solution650[] = { 9,23,20,22,11,17,7,18,11,14,5,7,6,2,8,21,12,13,23,1,2,16,24,0,20,12,6,5,5,13,13,1,15,7,1,3,17,4,17,4,12,22,4,22,20,13,10,20,8,17};
  const int solution850[] = { 12,9,1,18,14,6,8,7,17,7,24,13,5,16,0,22,21,21,17,18,1,10,10,11,15,13,1,13,4,11,23,21,2,1,24,2,4,2,7,12,7,11,0,5,20,22,7,17,19,14};
  const int solution20[] = { 8,26,40,46,15,35,15,48,25,16,12,37,17,1,41,17,45,7,49,38,18,10,35,22,33,14,6,4,2,15,10,13,19,32,41,45,46,36,20,21,22,0,33,25,47,11,11,7,3,8 };

  const int * solution = solution650;

  for (int idx = 0; (idx < 50) && (!test.all_happy()); ++idx) {
    std::cout << "============== Move " << idx << " ===" << std::endl;
    std::cout << test.to_string() << std::endl;

    test.make_move(*(solution++));  
  }

  std::cout << "========== Final Board ===" << std::endl;
  std::cout << test.to_string() << std::endl;  

  return 0;
}
