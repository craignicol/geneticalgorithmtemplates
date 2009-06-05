// GAsus.h - Testing a possible O(1) implementation of SUS
// search. Based on hash table ideas to bound the range
// of elements in a binary search.

// See sus_select.idea for discussion of the concept here

#ifndef __GASUS_H__
#define __GASUS_H__

#include <map>
#include <vector>
#include <math.h>

namespace mg_GA {
   
 class sus_search {
    
    private:
    
    class sus_row {
     private:
       std::map<double,int> _row;
       
     public:
       int rowstart, rowend;
       void add(double freq, int idx) { _row.insert(std::pair<double,int>(idx,freq)); };
       int search(double target_freq) { return (_row.lower_bound(target_freq)->second); };
       sus_row(int _start = 0, int _end = 0) : rowstart(_start),rowend(_end) {};
    };
 
    typedef std::vector<sus_row> rowdata;
    typedef std::vector<double> fitdata;
    
    rowdata _fitnesses;
    float _scale;
    
  public:
    sus_search(double scaling = 1.0) : _scale(scaling) { _fitnesses.reserve(2000); };
    void reset() { _fitnesses.clear(); };
    void construct_data(fitdata& input)
      {
	 // Scale data
	 // Construct rows
 	 _fitnesses.clear();
	 // _fitnesses.reserve(2000);
	 sus_row blank(0,0);
	 
	 int n = 0;
	 int c = input.size(); // num of chromosomes
	 for (int i = 0; i < floor(input[c-1])/_scale; ++i) 
	   {
	      _fitnesses.push_back(blank);
	      // cout << "Adding row " << i << "...";
	      if((n==c) || (input[n+1]/_scale > i)) {
		 _fitnesses[i].rowstart = n;
		 _fitnesses[i].rowend = n;
	      } else {
		 _fitnesses[i].rowstart = n;
		 for(; (n <= c) && (input[n]/_scale <= i); ++n) {
		    _fitnesses[i].add(input[n]/_scale, n);
		 }
		 _fitnesses[i].rowend = n;
	      };
	      // cout << _fitnesses[i].rowend - _fitnesses[i].rowstart + 1 << " members." << endl;
	   };
      };
 
    int search_data(double target)
      {
	// Scale data here
 	// Search data
	int row = (int)(floor(target)/_scale);
	if(_fitnesses[row].rowstart == _fitnesses[row].rowend) {
	   return _fitnesses[row].rowstart;
	} else {
	   return _fitnesses[row].search(target/_scale);
	};
      }; 
  };
};

#endif // !defined __GASUS_H__
