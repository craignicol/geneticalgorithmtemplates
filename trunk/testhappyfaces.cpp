/***************************************************************************
                          testhappyfaces.cpp  -  
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

#include <stdlib.h>
#include <stdio.h>

#include "GAhappyfaces.h"
//#include "GApopulation.h"

int main(int argc, char **argv)
{
  HappyFaces test;
  int rseed;

  for (int i=0; i < 25; ++i)
    std::cout << i << ":\n" << test.to_string(test.moves[i]);

  std::cout << "==============\n";
  std::cout << test.to_string() << "==" << test.count_happy() << std::endl;;
  
  test.make_move(2,3);
  
  std::cout << "==============\n";
  std::cout << test.to_string() << "==" << test.count_happy() << std::endl;;

  std::cout << "========================\nGA test\n========================\n";

  double mutationrate = 0.7; // was 0, 0.7
  double chrommutationrate = 0.6; // was 0.8, 0.6
  mg_GA::population<mg_GA::chrom_cached_lazy<int, 50, 25> > happy_pop(mutationrate,chrommutationrate,0.1,1000,mg_GA::PC_RANK,mg_GA::SM_RANK);
  happy_pop.setverbose(1,50);

  std::cout << "Enter Seed Value: ";
  std::cin >> rseed;
  srand(rseed);

  HappyFaces happy_func;
  mg_GA::twopointx<int, 50> happy_xover;
  mg_GA::random_mutate<int, 50> mut_func;
   happy_pop.initialise((mg_GA::fitness_base<int,50> *)&happy_func, (mg_GA::crossover_base<int,50> *)&happy_xover, &mut_func);
  std::cout << "Range before: " << happy_pop.min() << " - " << happy_pop.max() << std::endl;
  happy_pop.run_once();
  std::cout << "Range after: " << happy_pop.min() << " - " << happy_pop.max() << std::endl;
  happy_pop.run(1000);
  std::cout << "Range final: " << happy_pop.min() << " - " << happy_pop.max() << std::endl;
  std::cout << "Chromosome range: " << happy_pop.min_chrom().showchrom() << " < " << happy_pop.max_chrom().showchrom() << std::endl;

  return EXIT_SUCCESS;
};
