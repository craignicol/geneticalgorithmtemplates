/***************************************************************************
                          GAfitness.h  -  description
                             -------------------
    begin                : Wed Feb 04 2004
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

#ifndef _GAFITNESS_H
#define _GAFITNESS_H

#include "GAchromosome.h"

namespace mg_GA {

/**
  *@author Craig Nicol
  */

#ifdef MG_GA_PGA
  // Use PGA offsets for functions where global maxima occurs at the origin
  const double offset = 0.053;
#else
  const double offset = 0;
#endif


  /***********************************************************************
   *            AUXILIARY FUNCTIONS                                      *
   ***********************************************************************/
  template<typename T>
    T sq(T x) {
    return x*x;
  }

  unsigned long bool2long(bool * barray, unsigned int size)
    {
      bool * firstbool = barray;
      bool * finalbool = barray + size;
      unsigned long result = 0;
      unsigned long mult = 1;

      if (size > sizeof(unsigned long)*8) {
	// Only convert the LSB
	firstbool += (size - sizeof(unsigned long)*8);
      }

      for (--finalbool; firstbool <= finalbool; --finalbool) {
	if (*finalbool)
	  result += mult;
	mult *= 2;
      }

      return result;
    };

  // Returns the double represented by the first 'size' bools of barray
  // Such that if barray is the unsigned representation of x
  // bool2double(x=0) == min, bool2double(x=2**size) == max
  double bool2double(bool * barray, unsigned int size, double min, double max) {
    unsigned long rawint = bool2long(barray, size);
    double range = max - min;
    double drange = 0;

    if (size >= sizeof(unsigned long)*8)
      drange = ULONG_MAX;
    else
      drange = (1<<size)-1;

    return ((rawint*range)/drange) + min; 
  }

  /***********************************************************************
   *            FITNESS FUNCTIONS                                        *
   ***********************************************************************/

  template<class _Ctype, int _Csize, _Ctype N>
  struct max_n : public fitness_base<_Ctype,_Csize>{
    double operator()(_Ctype chrom[_Csize]) {
      return count(chrom, chrom+_Csize, N);    
    };
  } ;
  
  // De Jong's first function
  template<class _Ctype, int _Csize>
  struct dj1 : public fitness_base<_Ctype,_Csize>{
    double operator()(_Ctype chrom[_Csize]) {
      const double rangemax = 5.12;
      const double rangemin = -rangemax;
      const int var_size = _Csize/3;

      double n1 = bool2double(chrom, var_size, rangemin, rangemax) + offset;
      double n2 = bool2double(chrom + var_size, var_size, rangemin, rangemax) + offset;
      double n3 = bool2double(chrom + 2*var_size, var_size, rangemin, rangemax) + offset;

      return 100 - sq(n1) - sq(n2) - sq(n3);    
    };
  } ;

  // De Jong's second function
  template<class _Ctype, int _Csize>
  struct dj2 : public fitness_base<_Ctype,_Csize>{
    double operator()(_Ctype chrom[_Csize]) {
      const double rangemax = 2.048;
      const double rangemin = -rangemax;
      const int var_size = _Csize/2;

      double n0 = bool2double(chrom, var_size, rangemin, rangemax) + offset;
      double n1 = bool2double(chrom + var_size, var_size, rangemin, rangemax) + offset;

      return 1000 - 100*sq(sq(n0) - n1) - sq(1-n1);    
    };
  } ;

  // De Jong's third function
  template<class _Ctype, int _Csize>
  struct dj3 : public fitness_base<_Ctype,_Csize>{
    double operator()(_Ctype chrom[_Csize]) {
      const double rangemax = 5.12;
      const double rangemin = -rangemax;
      const int var_size = _Csize/5;

      double n1 = bool2double(chrom, var_size, rangemin, rangemax) + offset;
      double n2 = bool2double(chrom + var_size, var_size, rangemin, rangemax) + offset;
      double n3 = bool2double(chrom + 2*var_size, var_size, rangemin, rangemax) + offset;
      double n4 = bool2double(chrom + 3*var_size, var_size, rangemin, rangemax) + offset;
      double n5 = bool2double(chrom + 4*var_size, var_size, rangemin, rangemax) + offset;

      return 25 - floor(n1) - floor(n2) - floor(n3) - floor(n4) - floor(n5);    
    };
  } ;

  // Binary f6 (Schaffer et. al. 
  // - Proceedings of the Third International Conference on Genetic Algorithms, 1989
  template<class _Ctype, int _Csize>
  struct bf6 : public fitness_base<_Ctype,_Csize>{
    double operator()(_Ctype chrom[_Csize]) {
      const double rangemax = 100;
      const double rangemin = -rangemax;
      const int var_size = _Csize/2;

      double n1 = bool2double(chrom, var_size, rangemin, rangemax) + offset;
      double n2 = bool2double(chrom + var_size, var_size, rangemin, rangemax) + offset;

      return sq(sin(sqrt(sq(n1)+sq(n2))))/(1+0.001*sq(sq(n1)+sq(n2)));    
    };
  } ;

  // Himmelblau function
  template<class _Ctype, int _Csize>
  struct himm : public fitness_base<_Ctype,_Csize>{
    double operator()(_Ctype chrom[_Csize]) {
      const double rangemax = 6;
      const double rangemin = -rangemax;
      const int var_size = _Csize/2;

      double n1 = bool2double(chrom, var_size, rangemin, rangemax) + offset;
      double n2 = bool2double(chrom + var_size, var_size, rangemin, rangemax) + offset;

      return fabs(200-sq(sq(n1)+n2-11)-sq(n1+sq(n2)-7));    
    };
  } ;

}

#endif // defined _GAFITNESS_H_
