/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Fri Sep 27 21:26:34 BST 2002
    copyright            : (C) 2002 by Craig Nicol
    email                : mad_goldfish@yahoo.co.uk
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
 *                                                                         *
 *   This set of templates is provided as source code only. If you have    *
 *    compilation queries then please use the e-mail address above, but    *
 *    support for this source is not guaranteed.                           *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "GApopulation.h"
#include "GAfitness.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <iomanip> // For setw/setfill

using namespace std;

int difftimeb(timeb *late, timeb *early, timeb *result)
{
  if (!((late->timezone == early->timezone) && (late->dstflag == early->dstflag))) {
    return -1;
  } else {
    if (late->millitm > early->millitm) {
      result->time = late->time - early->time;
      result->millitm = late->millitm - early->millitm;
    } else {
      result->time = late->time - early->time - 1;
      result->millitm = 1000 + late->millitm - early->millitm;
    }
    result->timezone = late->timezone;
    result->dstflag = late->dstflag;
    return 0;
  }
  return -2; // Should never return
}

template<typename T>
inline void print_time(const T text, const timeb& time, std::ostream& out = cout)
{
  out << text << time.time << "." << setw(3) << setfill('0') << time.millitm << endl;
}

/****** This does not work. setw and setfill not found... */
ostream& operator<<(ostream& stream, timeb& time)
{
  stream << time.time << "." << setw(3) << setfill('0') << time.millitm;
  return stream;
}

template<typename T, int Sz>
void testpop(int nruns, double popmrate, double chrommrate, double crate, long psize, mg_GA::population_control_t pc, mg_GA::selection_method_t sm, int npops = 1, double popcrate = 0, int verbose = 0, int step_size = 10)
{
  mg_GA::dj1<T, Sz> pop_func;
  mg_GA::twopointx<T, Sz> cross_func;
  mg_GA::random_mutate<T, Sz> mut_func;
   
  mg_GA::population<mg_GA::chrom_cached_lazy<T, Sz> > pop_test(popmrate, chrommrate, crate, psize, pc, sm, npops, popcrate, verbose, step_size);
  pop_test.setverbose(verbose);
  std::cout << "Populating world..." << std::endl;
  pop_test.initialise(&pop_func, &cross_func, &mut_func);
  std::cout << "Range before: " << pop_test.min() << " - " << pop_test.max() << std::endl;
  if (verbose > 1)
    std::cout << pop_test.showfitness() << std::endl;
  pop_test.run_once();
  std::cout << "Range after: " << pop_test.min() << " - " << pop_test.max() << std::endl;
  if (verbose > 1) 
    std::cout << pop_test.showfitness() << std::endl;

  timeb now,then,change;
  ftime(&then);
  //std::cout << "Time before: " << then.time << "." << setw(3) << setfill('0') << then.millitm << endl;
  //print_time("Time before: ", then);
  cout << "Time before: " << then << endl;
  pop_test.run(nruns);
  ftime(&now);
  //  std::cout << "Time after: " << now.time << "." << now.millitm << endl;
  print_time("Time after: ", now);
  if(difftimeb(&now,&then,&change) != 0) {
    std::cout << "Incompatable times." << endl;
  } else {
    // std::cout << "Time diff: " << change.time << "." << change.millitm << endl;
    print_time("Time diff: ", change);
  }
  std::cout << "Range final: " << pop_test.min() << " - " << pop_test.max() << std::endl;
  std::cout << "Chromosome range: " << pop_test.min_chrom().showchrom() << " - " << pop_test.max_chrom().showchrom() << std::endl;
  if (verbose > 1)
    std::cout << pop_test.showfitness() << std::endl;
}

void test_chrom()
{
  mg_GA::maxones<bool, 32> fz;
  mg_GA::maxones<double, 16> f16;
  mg_GA::allpointx<bool, 32> apx;
  mg_GA::random_mutate<bool, 32> rb3;

  mg_GA::chromosome<bool,32,2> testz(&fz,&apx,&rb3,0.1);
  //  testz.set_fitfunc(&fz);
  std::cout << "Hello, " << testz.tostring() << " = " << testz.showchrom() << " = " << testz.f() << std::endl;
  std::cout << mg_GA::bool2double(testz.getchrom(),32,-5,7) << '\t' << mg_GA::dj1<bool, 32>()(testz.getchrom()) << std::endl;

  for (int idx = 0; idx < 64; ++idx) {
    bool data[6] = {idx & 32, idx & 16, idx & 8, idx & 4, idx & 2, idx & 1};
    bool ans[6];
    mg_GA::twopointx<bool,6>()(data,data,ans);
    std::cout << idx << ": " << data[0] << data[1] << data[2] << data[3] << data[4] << data[5] 
	      << " = (" << mg_GA::bool2long(data,6) << ", " 
	      << mg_GA::bool2double(data,6,-5.12,5.12) << ", ["
 	      << mg_GA::bool2double(data,2,-5.12,5.12) << ", " 
	      << mg_GA::bool2double(data+2,2,-5.12,5.12) << ", " 
	      << mg_GA::bool2double(data+4,2,-5.12,5.12) << "], " 
	      << mg_GA::dj1<bool, 6>()(data) << " ! "
              << mg_GA::dj1<bool, 6>()(ans) << ")" << std::endl;
  }

  mg_GA::onepointx<double, 16> apx_d16;
  mg_GA::random_mutate<double, 16> rm_d16;
   
  mg_GA::chromosome<double, 16, 100> test1(&f16,&apx_d16,&rm_d16,0.1);
  mg_GA::chromosome<double, 16, 100> test2(&f16,&apx_d16,&rm_d16,0.1);

  std::cout << "Hello, " << test1.tostring() << " = " << test1.showchrom() << " = " << test1.f() << std::endl;
  std::cout << "Hello, " << test2.tostring() << " = " << test2.showchrom() << " = " << test2.f() << std::endl;

  mg_GA::chromosome<double, 16, 100> test3 = test1 + test2;
  std::cout << "test1 + test2 = " << test3.showchrom() << " = " << test3.f() << std::endl;
  ~test1;
  ~test2;
  ~test3;
  std::cout << "~test1 = " << test1.showchrom() << " = " <<  test1.f() << std::endl;
  std::cout << "~test2 = " << test2.showchrom() << " = " <<  test2.f() << std::endl;
  std::cout << "~test3 = " << test3.showchrom() << " = " <<  test3.f() << std::endl;
}

int main(int argc, char *argv[])
{
  double popmrate = 0.1;
  double chrommrate = 0.01;
  double crossrate = 0.2;
  int popsize = 10000;
  mg_GA::population_control_t pop_control = mg_GA::PC_REPLACE;
  int verbose = 0;
  int generations = 100;

  get_type<int>();
  get_type<double>();

  if (argc > 1) {
    srand(atoi(argv[1]));
  }
  if (argc > 2) {
    popmrate = atof(argv[2]);
  }
  if (argc > 3) {
    chrommrate = atof(argv[3]);
  }
  if (argc > 4) {
    crossrate = atof(argv[4]);
  }
  if (argc > 5) {
    popsize = atoi(argv[5]);
  }
  if (argc > 6) {
    if (atoi(argv[6]) != 0) {
      pop_control = mg_GA::PC_RANK;
    }
  }
  if (argc > 7) {
    verbose = atoi(argv[7]);
  }

  cout << "Args: " << endl;
  cout << "Seed: " << argv[1] << ", popmrate: " << popmrate
     << ", chrommrate: " << chrommrate << ", crossrate: " << crossrate
     << ", popsize: " << popsize << ", pop_control: " << pop_control
     << ", verbose: " << verbose << endl;
   
  // test_chrom();
  // cout << "--------------------------------------" << endl;
  cout << "USER2 (random) test:" << endl;
  testpop<bool,36>(generations, popmrate, chrommrate, crossrate, popsize, pop_control, mg_GA::SM_USER2, 1, 0, verbose);
  cout << "USER0 test:" << endl;
  testpop<bool,36>(generations, popmrate, chrommrate, crossrate, popsize, pop_control, mg_GA::SM_USER0, 1, 0, verbose);
  cout << "USER1 test:" << endl;
  testpop<bool,36>(generations, popmrate, chrommrate, crossrate, popsize, pop_control, mg_GA::SM_USER1, 1, 0, verbose);
  cout << "SUS test:" << endl;
  testpop<bool,36>(generations, popmrate, chrommrate, crossrate, popsize, pop_control, mg_GA::SM_SUS, 1, 0, verbose);

  return EXIT_SUCCESS;
}
