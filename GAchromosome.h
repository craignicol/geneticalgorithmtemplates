/***************************************************************************
                          GAchromosome.h  -  description
                             -------------------
    begin                : Fri Sep 27 2002
    copyright            : (C) 2002 by Craig Nicol
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
 * CURRENT STATUS                                                          *
 *   # chromosome is fully templated, but the sep string for printing a    *
 *     chromosome can lead to a segfault.                                  *
 *   # population_bitfield is templated against chromosome and             *
 *     its children.                                                       *
 *   # population_bitfield has a decent stochastic selector.               *
 *   # population_bitfield needs a decent wander selector.                 *
 ***************************************************************************/

#ifndef __GACHROMOSOME_H__
#define __GACHROMOSOME_H__
	
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <strstream>
#include <string>
#include <algorithm>
#include <type_traits>
#include <limits.h>
#include <typeinfo>

template<typename T> std::string is_true(std::string truestr, std::string falsestr, std::__false_type) { return falsestr; };
template<typename T> std::string is_true(std::string truestr, std::string falsestr, std::__true_type) { return truestr; };

template<typename T>
std::string get_type() {
  // std::strstream out;
  // out << is_true<T>("int", "float", typename _Is_integer<T>::_Integral()) << sizeof(T)*8 << '\0';
  return typeid(T).name(); // Only returns one character
  // return out.str();
}

namespace mg_GA {

  /***********************************************************************
   *            SAMPLE FITNESS FUNCTIONS                                 *
   ***********************************************************************/

  template<class _Ctype, int _Csize>
  struct fitness_base {
    virtual double operator()(_Ctype chrom[_Csize]) = 0;
  } ;

  template<class _Ctype, int _Csize>
  struct maxones : public fitness_base<_Ctype,_Csize>{
    double operator()(_Ctype chrom[_Csize]) {
      return count(chrom, chrom+_Csize, 1);    
    };
  } ;

  template<class _Ctype, int _Csize>
  struct crossover_base {
    virtual _Ctype * operator()(_Ctype chrom_one[_Csize], _Ctype chrom_two[_Csize], _Ctype result[_Csize]) = 0;
  } ;

  template<class _Ctype, int _Csize>
  struct allpointx: public crossover_base<_Ctype,_Csize> {
    _Ctype * operator()(_Ctype chrom_one[_Csize], _Ctype chrom_two[_Csize], _Ctype result[_Csize]) {

      for(int i=0; i<_Csize; i++) 
	result[i] = rand()%2?chrom_one[i]:chrom_two[i];
      return result;
    }
  } ;

  template<class _Ctype, int _Csize>
  struct onepointx: public crossover_base<_Ctype,_Csize> {
    _Ctype * operator()(_Ctype chrom_one[_Csize], _Ctype chrom_two[_Csize], _Ctype result[_Csize]) {

      _Ctype* first;
      _Ctype* second;

      int point = rand()%_Csize;
      int start_chrom = rand()%2;
      if(start_chrom == 1) {
	first = chrom_one;
	second = chrom_two;
      } else {
	first = chrom_two;
	second = chrom_one;
      }

      // would memcpy be more efficient here?
      for(int i=0; i<point; ++i) 
	result[i] = first[i];
      for(int i=point; i<_Csize; ++i)
	result[i] = second[i];

      return result;
    }
  } ;

  template<class _Ctype, int _Csize>
  struct twopointx: public crossover_base<_Ctype,_Csize> {
    _Ctype * operator()(_Ctype chrom_one[_Csize], _Ctype chrom_two[_Csize], _Ctype result[_Csize]) {

      _Ctype* first;
      _Ctype* second;

      int point1 = rand()%_Csize;
      int point2 = rand()%_Csize;

      if(point1>point2) {
	int temp = point1;
	point1 = point2;
	point2 = temp;
      }

       //assert(point1 < _Csize);
       //assert(point2 < _Csize);

      int start_chrom = rand()%2;
      if(start_chrom == 1) {
	first = chrom_one;
	second = chrom_two;
      } else {
	first = chrom_two;
	second = chrom_one;
      }

      // would memcpy be more efficient here?
      for(int i=0; i<point1; ++i) 
	result[i] = first[i];
      for(int i=point1; i<point2; ++i)
	result[i] = second[i];
      for(int i=point2; i<_Csize; ++i)
	result[i] = first[i];

      return result;
    }
  } ;

  template<class _Ctype, int _Csize>
  struct mutate_base {
//     mutate_base() { std::cout << "<" << typeid(_Ctype) << ", " << _Csize << ">" << std::endl; };
    virtual _Ctype * operator()(_Ctype chrom[_Csize], float mr) = 0;
  } ;

  template<class _Ctype, int _Csize>
  struct random_mutate: public mutate_base<_Ctype, _Csize> {
//     random_mutate() { std::cout << "r<" << typeid(_Ctype) << ", " << _Csize << ">" << std::endl; };
     _Ctype * operator()(_Ctype chrom[_Csize], float mr) {
	for(int i=0; i<_Csize; i++) 
	  chrom[i] = ((rand() % INT_MAX)<(mr*INT_MAX))?(rand() % _Csize):chrom[i];
	return chrom;
     };
  };
   
  template<typename T>
  void arg_func(double ff(T * i) ) {} ;

  /***********************************************************************
   *                                                                     *
   *            CHROMOSOME TEMPLATE DEFINITION                           *
   *                                                                     *
   ***********************************************************************/

  template<class _CTYPE = bool, int _Csize = 32, int _Numvals = 2> 
    class chromosome {

    /***********************************************************************
     *            VARIABLE DEFINITIONS                                     *
     ***********************************************************************/
    
      protected:
      typedef _CTYPE _chrom_t[_Csize];
      _chrom_t _chromosome;
      double _mrate;
      char sep[2];
      fitness_base<_CTYPE,_Csize> * _fitfunc;
      mutate_base<_CTYPE, _Csize> * _mutatefunc;
      crossover_base<_CTYPE, _Csize> * _crossfunc;
	
      // Defaults
      // static allpointx<_CTYPE, _Csize> defaultcross;

    /***********************************************************************
     *            FUNCTION DEFINITIONS                                     *
     ***********************************************************************/
		
      private:
      //      void set_sep() {if ((_Numvals >= 10) || (_Numvals <= 0)) sep = '-'; else sep = 24; /* Non-printing CAN signal */} ;
      void set_sep() {sep[1] = '\0'; if ((_Numvals < 10) && (_Numvals > 0)) sep[0] = '\0'; else sep[0] = '-'; /* assert(_crossfunc!=NULL); */ } ;
      
      protected:
      chromosome() {} ; // Disable default constructor

      //      double calcfitness() {return count(_chromosome, _chromosome+_Csize, 1); };
      // Function below SEGFAULTs if result is returned directly when called
      // from chromosome < chromosome
      double calcfitness() { double result = (*_fitfunc)(_chromosome); return result; /* assert(_crossfunc!=NULL); */ };
      // double calcfitness() { std::cout << "("; double result = (*_fitfunc)(_chromosome); std::cout << ")"; return result; };
      void init_chrom() {for(int i=0; i<_Csize; i++) _chromosome[i] = rand() % _Numvals; /* assert(_crossfunc!=NULL); */ set_sep();};

      // merge_chrom returns the array as its third argument as c++ cannot have an
      // array as a return type
      //      void merge_chrom(_chrom_t first, _chrom_t second, _chrom_t result) {for(int i=0; i<_Csize; i++) result[i] = rand()%2?first[i]:second[i]; };
      void merge_chrom(_chrom_t first, _chrom_t second, _chrom_t result) { /* assert(_crossfunc!=NULL); */ (*_crossfunc)(first, second, result); /* assert(_crossfunc!=NULL); */ };
      chromosome(_chrom_t& chrom, fitness_base<_CTYPE,_Csize> * ff, crossover_base<_CTYPE,_Csize> * xf, mutate_base<_CTYPE,_Csize> * mf, double mrate = 0.01) : _mrate(mrate), _fitfunc(ff), _mutatefunc(mf), _crossfunc(xf) { /* assert(_crossfunc!=NULL); assert(_mutatefunc!=NULL); */ std::copy(chrom, chrom+_Csize, _chromosome); set_sep(); };
      void mutate() { (*_mutatefunc)(_chromosome, _mrate); }; // for(int i=0; i<_Csize; i++) _chromosome[i] = ((rand() % INT_MAX)<(_mrate*INT_MAX))?(rand() % _Numvals):_chromosome[i]; assert(_crossfunc!=NULL);};

      public:
      chromosome(fitness_base<_CTYPE,_Csize>* ff, crossover_base<_CTYPE,_Csize> * xf, mutate_base<_CTYPE,_Csize> * mf, double mrate = 0.01) : _mrate(mrate), _fitfunc(ff), _mutatefunc(mf), _crossfunc(xf) { /* assert(_crossfunc!=NULL); assert(_mutatefunc!=NULL); */ init_chrom(); };
      //      chromosome(double mrate = 0.01) : _mrate(mrate) {for(int i=0; i<_Csize; i++) _chromosome[i] = rand() % _Numvals; _fitfunc = new maxones<_CTYPE,_Csize>; set_sep(); };
      ~chromosome() {};
			
      void set_fitfunc(fitness_base<_CTYPE,_Csize>* ff) {_fitfunc = ff; /* assert(_crossfunc!=NULL); */ } ;

      chromosome operator+(chromosome second) { _chrom_t result; /* assert(_crossfunc!=NULL); */ merge_chrom(_chromosome,second._chromosome,result); return chromosome(result,_fitfunc,_crossfunc,_mutatefunc); };
      chromosome operator~() { mutate(); return *this; };
      bool operator[](size_t __pos) {return _chromosome[__pos]; };

      std::string tostring() { std::strstream out; out << '<' << get_type<_CTYPE>() << ", " << _Csize << ">" << " chromosome" << '\0'; /* assert(_crossfunc!=NULL); */ return out.str(); };
      std::string showchrom() { std::strstream out; for(int i=0; i<_Csize; i++) { out << _chromosome[i] << sep; }; out << '\0'; /* assert(_crossfunc!=NULL); */ return out.str(); };
      _CTYPE* getchrom() { return _chromosome; };
			
      //			int f() {return _chromosome.count(); }; // Fitness function. I expect this to be overwritten
      double f() {return calcfitness(); } ;
    };
  /***********************************************************************
   *            'FRIENDLY' FUNCTIONS                                     *
   ***********************************************************************/	

  template <class _CTYPE, int _Csize, int _Numvals>
      static bool operator>(chromosome<_CTYPE,_Csize,_Numvals> first, chromosome<_CTYPE,_Csize,_Numvals> second) { return first.f() > second.f(); };
  template <class _CTYPE, int _Csize, int _Numvals>
      static bool operator>=(chromosome<_CTYPE,_Csize,_Numvals> first, chromosome<_CTYPE,_Csize,_Numvals> second) { return first.f() >= second.f(); };
  template <class _CTYPE, int _Csize, int _Numvals>
      static bool operator<(chromosome<_CTYPE,_Csize,_Numvals> first, chromosome<_CTYPE,_Csize,_Numvals> second) { return first.f() < second.f(); };
  template <class _CTYPE, int _Csize, int _Numvals>
      static bool operator<=(chromosome<_CTYPE,_Csize,_Numvals> first, chromosome<_CTYPE,_Csize,_Numvals> second) { return first.f() <= second.f(); };


  // TODO: Need to check that the caching here actually stops calcfitness() being
  // called everytime...
  template <class _Ctype, int _Csize, int _Maxval=2>
  class chrom_cached_lazy : public chromosome<_Ctype, _Csize, _Maxval>
  {
  private:
    bool _isvalid;
    double _fitness;

  public:
    chrom_cached_lazy(fitness_base<_Ctype,_Csize>* ff, crossover_base<_Ctype,_Csize> * xf, mutate_base<_Ctype,_Csize> * mf, double mr = 0.01) : _isvalid(false) {this->_mrate = mr; this->_fitfunc = ff; this->_crossfunc = xf; this->_mutatefunc = mf; /* assert(_crossfunc!=NULL); assert(_mutatefunc!=NULL); */ this->init_chrom(); };

    double f() {if (!_isvalid) {_fitness = this->calcfitness(); _isvalid = true;} return _fitness; };
    chrom_cached_lazy operator+(chrom_cached_lazy second) {/* assert(_crossfunc!=NULL); */ chrom_cached_lazy newc(this->_fitfunc, this->_crossfunc, this->_mutatefunc, this->_mrate); merge_chrom(this->_chromosome, second._chromosome, newc._chromosome); return newc; };
    chrom_cached_lazy operator~() {this->mutate(); _isvalid = false; return *this;}
  };

  /*
  template <typename _Ctype=bool, int _Csize=32, int _Maxval=2>
  class chrom_findzeroes : public chromosome<_Ctype, _Csize, _Maxval> {
  protected:
    double calcfitness() {return count(this->_chromosome, _chromosome+_Csize, 0); };

  public:
    chrom_findzeroes(double mr = 0.01) {chromosome<_Ctype, _Csize, _Maxval>(); _mrate = mr;};
    chrom_findzeroes operator+(chrom_findzeroes second) { return chrom_findzeroes(chromosome<_Ctype,_Csize,_Maxval>(*this) + chromosome<_Ctype,_Csize,_Maxval>(second)); };
  };
  */


};

#endif
