/***************************************************************************
                          GAhappyfaces.h  -  
			  for solving the happy faces puzzle
                             -------------------
    begin                : Wed Mar 03 2004
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

/***************************************************************************
 * The happy faces puzzle is played on a 5x5 grid. Each square can be in   *
 * one of two states (happy or sad, but we will use 0 and 1 here).         *
 * On each move, the chosen square and it's immediate neighbours above,    *
 * below, left and right are all switched into their opposite state.       *
 * The game starts with 25 sad faces. The game ends when the grid contains *
 * 25 happy faces. It is possible to complete the game in 25 moves.        *
 ***************************************************************************/

#include "GAhappyfaces.h"

grid_t HappyFaces::moves[25] = {
  // All possible moves
  // xxxxxxx xxxxx xxxxx xxxxx xxxxx xxxxx
  
  // 0000000 11000 10000 00000 00000 00000
  0x01880000,
  // 0000000 11100 01000 00000 00000 00000
  0x01C40000,
  // 0000000 01110 00100 00000 00000 00000
  0x00E20000,
  // 0000000 00111 00010 00000 00000 00000
  0x00710000,
  // 0000000 00011 00001 00000 00000 00000
  0x00308000,
  // 0000000 10000 11000 10000 00000 00000
  0x010C4000,
  // 0000000 01000 11100 01000 00000 00000
  0x008E2000,
  // 0000000 00100 01110 00100 00000 00000
  0x00471000,
  // 0000000 00010 00111 00010 00000 00000
  0x00238800,
  // 0000000 00001 00011 00001 00000 00000
  0x00118400,
  // 0000000 00000 10000 11000 10000 00000
  0x00086200,
  // 0000000 00000 01000 11100 01000 00000
  0x00047100,
  // 0000000 00000 00100 01110 00100 00000
  0x00023880,
  // 0000000 00000 00010 00111 00010 00000
  0x00011C40,
  // 0000000 00000 00001 00011 00001 00000
  0x00008C20,
  // 0000000 00000 00000 10000 11000 10000
  0x00004310,
  // 0000000 00000 00000 01000 11100 01000
  0x00002388,
  // 0000000 00000 00000 00100 01110 00100
  0x000011C4,
  // 0000000 00000 00000 00010 00111 00010
  0x000008E2,
  // 0000000 00000 00000 00001 00011 00001
  0x00000461,
  // 0000000 00000 00000 00000 10000 11000
  0x00000218,
  // 0000000 00000 00000 00000 01000 11100
  0x0000011C,
  // 0000000 00000 00000 00000 00100 01110
  0x0000008E,
  // 0000000 00000 00000 00000 00010 00111
  0x00000047,
  // 0000000 00000 00000 00000 00001 00011
  0x00000023
};

int HappyFaces::make_move(unsigned int x, unsigned int y) 
{ 
  if((x>4)&&(y>4)) 
    return -1; 

  int square = 5*x+y; 
  _current ^= moves[square]; 
  return 0;
};

int HappyFaces::make_move(unsigned int square) 
{ 
  if(square>=25) 
    return -1; 

  _current ^= moves[square]; 
  return 0;
};

std::string HappyFaces::to_string(grid_t &grid)
{
  std::string out;
  grid_t mask = 1;

  for(int i=0; i < 25; i++, mask <<= 1) {
    out = out + ((mask & grid)?'O':'-');
    if((i+1)%5 == 0)
      out = out + '\n';
  }

  return out;
}

int HappyFaces::count_happy() {
  if (all_sad())
    return 0;

  if (all_happy())
    return 25;

  int count = 0;
  grid_t mask = 1;

  for(int i = 0; i < 25; ++i, mask <<= 1) {
    if(mask & _current) {
      ++count;
    }
  }
  return count;
}

// GA fitness function
double HappyFaces::operator()(int chrom[50])
{
  int moves = 0;

  _current = 0;

  for(;(moves<50) && !all_happy(); ++moves) {
    make_move(chrom[moves]);
  }
  
  return count_happy() + (50-moves)*25;
}
