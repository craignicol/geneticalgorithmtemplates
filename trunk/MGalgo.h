/************************************
 * MGalgo.h                         *
 * A collection of helpful          *
 * algorithms not included in the   *
 * stl.                             *
 ************************************/

#ifndef __MGalgo_h__
#define __MGalgo_h__

#include <vector>

namespace mgtl { // Mildly Generic Template Library
   
   // int binary_search(T target, vector<T> list [, start, end])
   // 
   // Inputs:
   //   target - the value to search for
   //   list - the sorted list to search in
   //   start - smallest index in list to use (default 0)
   //   end - largest index in list to use (default list.size()-1)
   // 
   // Outputs:
   //   index of the largest value in list no greater than target
   //   i.e. list[result] <= target < list[result+1]
   //   except where target < list[start] when result = start-1
   template <typename T>
     int binary_search_recurse(T target, std::vector<T> li, int start, int end) {
	if(target < li[start]) return start-1;
	if(start==end) return start;
	int middle = (start + end + 1)/2; // Add one to take the ceiling of the mean
	if(target < li[middle]) {
	   return binary_search_recurse(target, li, start, middle-1);
	} else {
	   return binary_search_recurse(target, li, middle, end);
	}
     };

   template <typename T>
     int binary_search_norecurse(T target, std::vector<T>& li, int start, int end) {
        int middle;
        while(start < end) {
	  middle = (start + end + 1)/2; // Add one to take the ceiling of the mean
	  if(target < li[middle]) {
	    end = middle-1;
	  } else {
	    start = middle;
	  }
	}
	return start;
   }
   
   template <typename T>
     inline int binary_search(T target, std::vector<T> li, bool recurse=false) {
        if(recurse)
	  return binary_search_recurse(target, li, 0, li.size()-1);
	else
	  return binary_search_norecurse(target, li, 0, li.size()-1);
     };

   // This version uses iterators to avoid copying a vector
    template <typename T>
   inline const int binary_search(const T target, typename std::vector<T>::iterator start, typename std::vector<T>::iterator end, typename std::vector<T>::iterator result) {
     typename std::vector<T>::iterator first = start;
     typename std::vector<T>::iterator last = end;
     unsigned long middle = ((unsigned long)(1 + last - first)>>1);
     unsigned int i = 0;
     while(first < last) {
//       middle >>= 4;
//       middle <<= 4; // make sure it's a multiple of 16

       // The mean should err on the larger side to guarantee
       // The next iteration shrinks the search
       // Ignore first version as it can give overflow...
       // middle = (T *)( ( ((unsigned long)first + (unsigned long)last) + 1) / 2);
       // This one is still too slow...
       // middle = first + ((unsigned long)(1 + last - first)>>1);
       if(target < first[middle]) {
	  last -= middle+1;
       } else {
	  first += middle;
       }
       middle = (middle+1) >> 1;	
       ++i;
     }
     result = first;
     return (int)(first - start);
   }
   template <typename T>
   inline const int binary_search(const T target, typename std::vector<T>::iterator start, typename std::vector<T>::iterator end) {
	return binary_search<T>(target, start, end, start);
	};


   // ToDo: Implement a binary_search using random_access_iterator
   // iter binary_search(T target, container<T> list, iter start, iter end)
   //    [With the priviso that start and end are valid iterators for list]
   // if(target < *start) return NULL;
   // if(start==end) return start;
   // iter middle = start + (1 + end - start)/2;
   // if (target < *middle) {
   //   return binary_search(target, list, start, --middle);
   // } else {
   //   return binary_search(target, list, middle, end);
   // }
   // [BUT, STL end() is one past the last element. Is this a problem?]
};

#endif // !defined __MGalgo_h__
