/***************************************************************************
                          GApopulation.h  -  description
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

#ifndef __GAPOPULATION_H__
#define __GAPOPULATION_H__

#include "GAchromosome.h"
#include <numeric>  // For accumulate
#include "MGalgo.h"
#include "GAsus.h" // New O(1) SUS algorithm
#include <iostream>

namespace mg_GA {
  /***********************************************************************
   *		NAMESPACE VARIABLE DEFINITIONS                           *
   ***********************************************************************/
	
  enum selection_method_list {
    SM_ROULETTE,
    SM_SUS,
    SM_RANK,
    SM_WANDER,
    SM_TOURNAMENT,
    SM_RANDOM,
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

  const bool SORT_ASCENDING = false;
  const bool SORT_DESCENDING = true;

  /***********************************************************************
   *                                                                     *
   *           POPULATION CONTROL                                        *
   *                                                                     *
   ***********************************************************************/

  template<class chrom_t = chromosome<bool, 32> >
  class population {
	
    /***********************************************************************
     *            VARIABLE DEFINITIONS                                     *
     ***********************************************************************/
	
  protected:
	 
    // Using vector here as we require a lot of random access
    // and not much resizing
    //    typedef chromosome<bool, 32> chrom_t;
    typedef std::vector<chrom_t> pop_t;

    pop_t _population;
    double _chrom_mutation_rate;
    double _mutation_rate, _crossover_rate, _population_crossover_rate;
    long _population_size;
    int _npopulations; // Number of populations, = 1 for now
    selection_method_t _selection_method;
    population_control_t _population_control;
    bool _autosort;
	 
    int _verb;
    int _verb_step;
    double _max_fitness;

    // select_chromosome cached data 
    double _sc_sum_fitness;
    std::vector<double> _sc_cumul_fitns;
    sus_search _sc_suss;
     
  private:
    int _generation;
	
    /***********************************************************************
     *            FUNCTION DEFINITIONS                                     *
     ***********************************************************************/
	
  public:
    template<class _Ctype, int _Csize>
    void initialise(fitness_base<_Ctype,_Csize> * fitfunc, crossover_base<_Ctype,_Csize> * crossfunc, mutate_base<_Ctype,_Csize> * mutfunc)
    {
	//initialise population here
	for (int i = 0; i < _population_size; i++) {
	  _population.push_back(chrom_t(fitfunc, crossfunc, mutfunc, _chrom_mutation_rate));
	}
      };
	 				
    ~population() { };
	
    void run_once() {
       prepare();
       switch (_population_control) {
	case PC_REPLACE: run_once_replace(); break; //this will be able to step through the populations
	case PC_RANK: run_once_rank(_population_size / 2); break;
	default: // raise an error here
	  ;
       }
       mutate();
    };
	 	
    void run(int ngenerations) { // as above but with loops
      if(_verb > 0) {
	std::cerr << "Running " << ngenerations << " generations." << std::endl;
	_max_fitness = this->max();
      }
      for (int i = 0; i < ngenerations; i++) {
	_generation = i;
	run_once();
	if (_verb > 0) {
	  if (this->max() > _max_fitness) {
	    _max_fitness = this->max();
	    std::cerr << "New fittest: " << max_chrom().showchrom() << " = " << _max_fitness << std::endl;
	  }
	  if (i%_verb_step == 0) {
	    std::cerr << "Completed " << i << " of " << ngenerations << " generations." << std::endl;
	  }
	}
      }
      _generation = -1; 
    };

    void setverbose(int verbose, int step_size = 10) {
      _verb = verbose;
      _verb_step = step_size; 
    }

    std::string showfitness() {
      std::strstream out;
	   
      for (typename pop_t::iterator it = _population.begin(); it != _population.end(); ++it) {
	out << it->f() << ' ';
      }

      out << '\0';
      return out.str();
    }

    void sort(bool dir=SORT_ASCENDING) {
	std::sort(_population.begin(), _population.end());
	if(dir == SORT_DESCENDING)
	  std::reverse(_population.begin(), _population.end());
    }

    chrom_t first() {
      return _population[0];
    }

    chrom_t last() {
      return _population[_population_size - 1];
    }

    chrom_t min_chrom() {
      return *std::min_element(_population.begin(), _population.end());
    }

    chrom_t max_chrom() {
      return *std::max_element(_population.begin(), _population.end());
    }
	
    double min() { return min_chrom().f(); };
    double max() { return max_chrom().f(); };

  private:
    bool rand_less(float select_rate) {
      return ((rand() % INT_MAX) < (select_rate * INT_MAX));
    }

    void mutate() {
      for(int i = 0; i < _population_size; i++) {
	if (rand_less(_mutation_rate)) {
	  ~_population[i];
	}
      }
    }

    void run_once_rank(int N) {			// one generation = (0 < N <= psize) crossovers
      static pop_t add_pop;
      add_pop.clear();
      add_pop.reserve(N);
      double minf = this->min();
	 		
      for(int i = 0; i < N; i++) {
	chrom_t child = select_chromosome() + select_chromosome();
	~child;
	if (child.f() > minf ) {
	  add_pop.push_back(child);
	}
      }
      if (_verb > 1) {
	std::cerr << "Adding " << add_pop.size() << " children..." << std::endl;
	std::cerr << "Original Range: " << this->min() << " < " << this->max() << std::endl;
	std::cerr << "Adding Range: " << std::min_element(add_pop.begin(), add_pop.end())->f() << " < " << std::max_element(add_pop.begin(), add_pop.end())->f() << std::endl;
      }
      if (!_autosort)
	sort(SORT_ASCENDING);
      copy(add_pop.begin(), add_pop.end(), _population.begin());

    };
	 			
    void run_once_replace() { // one generation = psize crossovers
      pop_t newpop;
      newpop.reserve(_population_size);

      for(int i = 0; i < _population_size; ++i) {
	chrom_t child = select_chromosome() + select_chromosome();
	if (_verb > 1) {
	  std::cerr << child.f() << ',';
	}
	newpop.push_back(child);
      }
			
      if (_verb > 1) {
	std::cerr << std::endl;
      }

      copy(newpop.begin(),newpop.end(),_population.begin());

    };

    struct add_fitness : public std::binary_function<double, chrom_t, double> {
       double operator()(double total, chrom_t C) { 
	  return total + C.f();
       };
    } _fitness_sum; 
     
    double total_fitness() { 
       return std::accumulate(_population.begin(), _population.end(), 0.0, _fitness_sum);  
    };
    
    std::vector<double> fitness_redistribution(std::vector<double> fitnesses) { 
       double cum_f = 0;
       int i = 0;
       for (typename pop_t::iterator pit = _population.begin();
	    pit != _population.end(); ++pit, ++i) {
	
	  cum_f += pit->f();
	  fitnesses[i] = cum_f;
       }
       return fitnesses;
    };

    std::vector<double> fitness_distribution() {
       std::vector<double> fitnesses(_population_size);
       return fitness_redistribution(fitnesses);
    };

    void prepare() {
      if (_autosort) {
	sort(SORT_ASCENDING);
      };
       
       switch(_selection_method) {
	case SM_ROULETTE:
	case SM_SUS:
	case SM_USER0:
	    {
	       _sc_sum_fitness = total_fitness();
	       _sc_cumul_fitns = fitness_redistribution(_sc_cumul_fitns);
	    } break;
	case SM_USER1:
	    {
	       _sc_sum_fitness = total_fitness();
	       _sc_cumul_fitns = fitness_redistribution(_sc_cumul_fitns);
	       _sc_suss.construct_data(_sc_cumul_fitns);
	    } break;
       default:
	 // Do nothing
	 break;
       };
    };
     
    chrom_t select_chromosome() {
 
      // Assume _population is sorted s.t. largest value is at end
      switch (_selection_method) {
      case SM_RANK:
	// Prefer larger indices...
	return _population[floor(pow(rand() % int(pow(_population_size, 1.0/_crossover_rate)), _crossover_rate))]; break;
      case SM_TOURNAMENT: {
	chrom_t first = _population[rand() % _population_size];
	chrom_t second = _population[rand() % _population_size];
	return ((first.f() >= second.f()) && !rand_less(_crossover_rate))?first:second;
      }
	break;
      case SM_ROULETTE:
      case SM_SUS: {
	 // Pick stuff - this is the slow way to do it
	 double target_f = (double(rand()) * _sc_sum_fitness)/RAND_MAX;
	 unsigned int i = 0;
	 for(typename std::vector<double>::iterator cum_it = _sc_cumul_fitns.begin();
	     (cum_it != _sc_cumul_fitns.end()) && (*cum_it < target_f); 
	     ++cum_it, ++i) {
	 }
	 if (i < _population.size())
	   return _population[i];
	 else
	   return _population[_population.size()-1];
      }
	 break;
	 
      case SM_USER0: // Test faster SM_SUS O(log n)
      {
	 unsigned int index;
	 double target_f = (double(rand()) * _sc_sum_fitness)/RAND_MAX;
	 //index = mgtl::binary_search(target_f, _sc_cumul_fitns);
	 index = mgtl::binary_search(target_f, _sc_cumul_fitns.begin(), _sc_cumul_fitns.end()-1);
	 if (index >= _population.size())
	   return _population[_population.size() - 1];
	 else
	   return _population[index];
      }
	 break;
      case SM_USER1: // Test O(1) SUS algorithm
	   {
	      unsigned int index;
	      double target_f = (double(rand()) * _sc_sum_fitness)/RAND_MAX;
	      index = _sc_suss.search_data(target_f);
	      if (index < _population.size())
		return _population[index];
	      else
		return _population[_population.size()-1];
	   }
	 break;
      case SM_WANDER:
      case SM_RANDOM:
      default:
	return _population[rand() % _population_size];
      }
    };
     
   public:
    population(double popmrate = 0.01, double chrommrate = 0.1,
			double crate = 0.2, long psize = 100,
			population_control_t pc = PC_REPLACE,
			selection_method_t sm = SM_SUS,
			int npops = 1, double popcrate = 0, 
	                int verbose = 0, int step_size = 10,
	                float sus_factor = 10.0) :
                        _sc_suss(sus_factor)
      {
	// check data is valid here
	 			
	_mutation_rate = popmrate;
	_chrom_mutation_rate = chrommrate;
	_crossover_rate = crate;
	_population_crossover_rate = popcrate;
	_population_size = psize;
	_npopulations = npops;
	_population_control = pc;
	_selection_method = sm;
	_verb = verbose;
	_verb_step = step_size;

	_population.reserve(_population_size);
	 
	_generation = -1;
	 
	switch (_selection_method) {
	 case SM_RANK:
	   _autosort = true;
	   break;
	 case SM_SUS:
	 case SM_ROULETTE:
	 case SM_USER0:
	 case SM_USER1:
	   _sc_cumul_fitns = fitness_distribution();
	 default:
	   _autosort = false;
	}
      }     
  };

  template<class chrom_t = chromosome<bool, 32>, int PSIZE = 1000 >
  class population_array : public population<chrom_t> {
    typedef chrom_t pop_t[PSIZE];
    chrom_t * _population_array;

  public:
    population_array(double popmrate = 0.01, double chrommrate = 0.1,
		        double crate = 0.2, int ignored_psize = PSIZE,
			population_control_t pc = PC_REPLACE,
			selection_method_t sm = SM_SUS,
			int npops = 1, double popcrate = 0, 
			int verbose = 0, int step_size = 10)
      {
	// check data is valid here
	// check data is valid here
	 			
	this->_mutation_rate = popmrate;
	this->_chrom_mutation_rate = chrommrate;
	this->_crossover_rate = crate;
	this->_population_crossover_rate = popcrate;
	this->_population_size = PSIZE;
	this->_npopulations = npops;
	this->_population_control = pc;
	this->_selection_method = sm;
	this->_verb = verbose;
	this->_verb_step = step_size;

	if (this->_selection_method == SM_RANK)
	  this->_autosort = true;
	else
	  this->_autosort = false;
      };
    
      template<class _Ctype, int _Csize>
      void initialise(fitness_base<_Ctype,_Csize> * fitfunc, crossover_base<_Ctype,_Csize> * crossfunc)
      {
	this->_population_array = new chrom_t(fitfunc, crossfunc, this->_chrom_mutation_rate);
      };
  };

};

#endif // __GAPOPULATION_H__
