/***************************************************************************
                          GAchromosome.h  -  description
                             -------------------
    begin                : Fri Sep 27 2002
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
 * CURRENT STATUS                                                          *
 *   # chromosome is fully templated, but the sep string for printing a    *
 *     chromosome can lead to a segfault.                                  *
 *   # chrom_findzeroes is a fully realised child class of chromosome      *
 *   # population_bitfield need to be templated against chromosome and     *
 *     its children.                                                       *
 *   # population_bitfield needs a decent stochastic selector.             *
 *   # population_bitfield needs a decent wander selector.                 *
 ***************************************************************************/
	
#include <stdlib.h>
#include <math.h>
#include <g++-3/bitset>
#include <g++-3/vector>
#include <g++-3/strstream>
#include <g++-3/algorithm>
#include <g++-3/type_traits.h>

#define MAX_INT 10000

template<typename T> std::string is_true(std::string truestr, std::string falsestr, __false_type) { return falsestr; };
template<typename T> std::string is_true(std::string truestr, std::string falsestr, __true_type) { return truestr; };

template<typename T>
std::string get_type() {
  std::strstream out;
  out << is_true<T>("int", "float", typename _Is_integer<T>::_Integral()) << sizeof(T)*8 << '\0';
  return out.str();
}

namespace mg_GA {

  /***********************************************************************
   *						NAMESPACE VARIABLE DEFINITIONS                           *
   ***********************************************************************/
	
  enum selection_method_list {
    SM_STOCHASTIC,
    SM_RANK,
    SM_WANDER,
    SM_USER0 = 1000,
    SM_USER1,
    SM_USER2,
    SM_USER3,
    SM_USER4,
    SM_USER5,
    SM_USER6,
    SM_USER7,
    SM_USER8,
    SM_USER9
  };
	
  typedef enum selection_method_list selection_method_t;
	
  enum population_control_list {
    PC_REPLACE,
    PC_RANK,
    PC_USER0 = 1000,
    PC_USER1,
    PC_USER2,
    PC_USER3,
    PC_USER4,
    PC_USER5,
    PC_USER6,
    PC_USER7,
    PC_USER8,
    PC_USER9
  };
	
  typedef enum population_control_list population_control_t;
	
	
	
  /***********************************************************************
   *                                                                     *
   *            CHROMOSOME DEFINITION                                    *
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
      std::string sep;
      chromosome(_chrom_t& chrom, double mrate = 0.01) : _chromosome(chrom), _mrate(mrate) {set_sep(); };
	
    /***********************************************************************
     *            FUNCTION DEFINITIONS                                     *
     ***********************************************************************/
		
      private:
      void set_sep() {if ((_Numvals >= 10) || (_Numvals <= 0)) sep = "-"; else sep = ""; } ;
      
      public:
      chromosome(double mrate = 0.01) : _mrate(mrate) {for(int i=0; i<_Csize; i++) _chromosome[i] = rand() % _Numvals; set_sep(); };
      ~chromosome() {};
			
      chromosome operator+(chromosome second) {_chrom_t result; for(int i=0; i<_Csize; i++) result[i] = rand()%2?_chromosome[i]:second[i]; return chromosome(result); };
      chromosome operator~() {for(int i=0; i<_Csize; i++) _chromosome[i] = ((rand() % INT_MAX)<(_mrate*INT_MAX))?(rand() % _Numvals):_chromosome[i]; };
      bool operator[](size_t __pos) {return _chromosome[__pos]; };

      std::string tostring() { std::strstream out; out << '<' << get_type<_CTYPE>() << ", " << _Csize << ">" << " chromosome" << '\0'; return out.str(); };

      // WARNING: This function can SEGFAULT if you use std::string sep in place of "-"...
      std::string showchrom() { std::strstream out; for(int i=0; i<_Csize; i++) { out << _chromosome[i] << "-"; }; out << '\0'; return out.str(); };
			
      //			int f() {return _chromosome.count(); }; // Fitness function. I expect this to be overwritten
      int f() {return count(_chromosome, _chromosome+_Csize, 1); };
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

  class chrom_findzeroes : public chromosome<bool, 32, 2> {
  public:
    chrom_findzeroes(double mr = 0.01) {chromosome<bool, 32, 2>(); _mrate = mr;};

    int f() {return count(_chromosome, _chromosome+32, 0); };
  };

  class chromosome_bitfield { // sample class upon which the template will be built

    /***********************************************************************
     *            VARIABLE DEFINITIONS                                     *
     ***********************************************************************/
	 		
  protected:
    // typedef std::bitset<32> chrom_t;
    typedef bool chrom_t[32];
    chrom_t _chromosome;
    double _mrate;
	
    /***********************************************************************
     *            FUNCTION DEFINITIONS                                     *
     ***********************************************************************/
		
  public:
    chromosome_bitfield(double mrate = 0.01) : _mrate(mrate) {for(int i=0; i<32; i++) _chromosome[i] = rand() % 2; };
    chromosome_bitfield(chrom_t& chrom, double mrate = 0.01) : _chromosome(chrom), _mrate(mrate) {};
    ~chromosome_bitfield() {};
			
    chromosome_bitfield operator+(chromosome_bitfield second) {chrom_t result; for(int i=0; i<32; i++) result[i] = rand()%2?_chromosome[i]:second[i]; return chromosome_bitfield(result); };
    chromosome_bitfield operator~() {for(int i=0; i<32; i++) _chromosome[i] = ((rand() % INT_MAX)<(_mrate*INT_MAX))?!_chromosome[i]:_chromosome[i]; };
    bool operator[](size_t __pos) {return _chromosome[__pos]; };

    std::string string() { return std::string("Is a chromosome"); }; // For debug purposes
    std::string showchrom() { std::strstream out; for(int i=0; i<32; i++) out << _chromosome[i]; out << '\0'; return out.str(); };
			
    //			int f() {return _chromosome.count(); }; // Fitness function. I expect this to be overwritten
    int f() {return count(_chromosome, _chromosome+32, true); };
  };

  /***********************************************************************
   *            'FRIENDLY' FUNCTIONS                                     *
   ***********************************************************************/	

  static bool operator>(chromosome_bitfield first, chromosome_bitfield second) { return first.f() > second.f(); };
  static bool operator>=(chromosome_bitfield first, chromosome_bitfield second) { return first.f() >= second.f(); };
  static bool operator<(chromosome_bitfield first, chromosome_bitfield second) { return first.f() < second.f(); };
  static bool operator<=(chromosome_bitfield first, chromosome_bitfield second) { return first.f() <= second.f(); };

  /***********************************************************************
   *                                                                     *
   *           POPULATION CONTROL                                        *
   *                                                                     *
   ***********************************************************************/
	
  class population_bitfield {
	
    /***********************************************************************
     *            VARIABLE DEFINITIONS                                     *
     ***********************************************************************/
	
  protected:
	 
    // Using vector here as we require a lot of random access
    // and not much resizing
    typedef chromosome_bitfield chrom_t;
    typedef vector<chrom_t> pop_t;
	
    pop_t _population;
    double _mutation_rate, _crossover_rate, _population_crossover_rate;
    long _population_size;
    int _npopulations; // Number of populations, = 1 for now
    selection_method_t _selection_method;
    population_control_t _population_control;
	 
    int _verb;

  private:
    int generation;
	
    /***********************************************************************
     *            FUNCTION DEFINITIONS                                     *
     ***********************************************************************/
	
  public:
    population_bitfield(double popmrate = 0.01, double chrommrate = 0.1,
			double crate = 0.2, long psize = 100,
			population_control_t pc = PC_REPLACE,
			selection_method_t sm = SM_STOCHASTIC,
			int npops = 1, double popcrate = 0, 
			int verbose = 0)
      {
	// check data is valid here
	 			
	_mutation_rate = popmrate;
	_crossover_rate = crate;
	_population_crossover_rate = popcrate;
	_population_size = psize;
	_npopulations = npops;
	_population_control = pc;
	_selection_method = sm;
	_verb = verbose;
	 					
	//initialise population here
	for (int i = 0; i < _population_size; i++) {
	  _population.push_back(chromosome_bitfield(chrommrate));
	}
						
	// Sort ascending
	std::sort(_population.begin(), _population.end());
	// std::reverse(_population.begin(), _population.end());
      };
	 				
    ~population_bitfield() {};
	
    void run_once() {
      switch (_population_control) {
      case PC_REPLACE: run_once_replace(); break; //this will be able to step through the populations
      case PC_RANK: run_once_rank(_population_size / 2); break;
      default: // raise an error here
	;
      }
      mutate();
      if (_selection_method != SM_WANDER) {
	// Sort ascending
	std::sort(_population.begin(), _population.end());
	// std::reverse(_population.begin(), _population.end());
      }
    };
	 	
    void run(int ngenerations) { // as above but with loops
      for (int i = 0; i < ngenerations; i++) {
	run_once();
      }
    };

    void setverbose(int verbose) {
      _verb = verbose;
    }

    std::string showfitness() {
      std::strstream out;
	   
      for (int i = 0; i < _population_size; i++) {
	out << _population[i].f() << ' ';
      }

      out << '\0';
      return out.str();
    }

    chrom_t first() {
      return _population[0];
    }

    chrom_t last() {
      return _population[_population_size - 1];
    }

    int min() {
      return min_element(_population.begin(), _population.end())->f();

      /*
      int result = 50;

      for (int i = 0; i < _population_size; i++) {
	if (_population[i].f() < result) {
	  result = _population[i].f();
	}
      }
      return result;
      */
    }

    int max() {
      return max_element(_population.begin(), _population.end())->f();

/*       int result = 0; */

/*       for (int i = 0; i < _population_size; i++) { */
/* 	if (_population[i].f() > result) { */
/* 	  result = _population[i].f(); */
/* 	} */
/*       } */
/*       return result; */
    }
	
  private:
    void mutate() {
      for(int i = 0; i < _population_size; i++) {
	if ((rand() % MAX_INT) < (_mutation_rate * MAX_INT)) {
	  ~_population[i];
	}
      }
    }

    void run_once_rank(int N) {			// one generation = (0 < N <= psize) crossovers
      chromosome_bitfield child;
      pop_t add_pop;
      int minf = this->min();
	 		
      for(int i = 0; i < N; i++) {
	child = select_chromosome() + select_chromosome();
	if (child.f() > minf ) {
	  add_pop.push_back(child);
	}
      }
      if (_verb > 1) {
	std::cerr << "Adding " << add_pop.size() << " children,,, ";
      }
      _population.erase(_population.begin(), _population.begin() + add_pop.size());
      _population.insert(_population.begin(), add_pop.begin(), add_pop.end());
    };
	 			
    void run_once_replace() { // one generation = psize crossovers
      pop_t newpop;
      chromosome_bitfield child;
	
      for(int i = 0; i < _population_size; ++i) {
	child = select_chromosome() + select_chromosome();
	if (_verb > 1) {
	  std::cerr << child.f() << '-';
	}
	newpop.push_back(child);
      }
			
      if (_verb > 1) {
	std::cerr << std::endl;
      }

      _population.clear();
      _population.insert(_population.begin(),newpop.begin(),newpop.end());
    };
	
    chromosome_bitfield select_chromosome() {
      // Assume _population is sorted s.t. largest value is at end
      switch (_selection_method) {
      case SM_RANK:
	// Prefer larger indices...
	return _population[floor(pow(rand() % int(pow(_population_size, 1.0/_crossover_rate)), _crossover_rate))];
      case SM_STOCHASTIC:
      case SM_WANDER:
      default:
	return _population[rand() % _population_size];
      }
    };
  };
	
};
