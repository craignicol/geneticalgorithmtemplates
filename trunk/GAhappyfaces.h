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

#include <string>

#include "GApopulation.h"

#ifndef _GAHAPPYFACES_H
#define _GAHAPPYFACES_H

typedef int grid_t; // Where grid contains at least 25 boolean values

// class HappyFaces { // Non-GA version
class HappyFaces : public mg_GA::fitness_base<int,50> { // GA version with a 50-move sequence

 public:
  static grid_t moves[25]; // List of possible moves

  int make_move(unsigned int x, unsigned int y);
  int make_move(unsigned int square);
  std::string to_string(grid_t &grid);
  std::string to_string() {return to_string(_current); };

  int count_happy();
  int count_sad() {return 25-count_happy();};
  bool all_happy() {return (_current == (1<<25)-1); };
  bool all_sad() {return (_current == 0); };

  HappyFaces(grid_t start = 0) : _current(start) {};

  // GA fitness function
  double operator()(int chrom[50]);

 private:
  grid_t _current;

};

#endif // defined _GAHAPPYFACES_H
