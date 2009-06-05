#include "GAsus.h"
#include <iostream>
#include <functional>

int main(int argc, char **argv)
{
   vector<float> cum_f(2000);
   int i;
   float j;
   float total = 0;
   
   cout << "Populating cum_f..." << endl;
   for (i = 0, j = 0.2; i < 2000; ++i, j=fmod(j*i,10.0f), total += j) {
      cum_f[i] = total;
   }
   
   cout << "Constructing test..." << endl;
   mg_GA::sus_search test;
   cout << "Constructing test data..." << endl;
   test.construct_data(cum_f);
   float suite[] = { 0.2, 4.3, 2.6, 10.67, 29.2 };
   int result;
   
   cout << "Searching test data..." << endl;
   for (int idx = 0; idx < (sizeof(suite)/sizeof(float)); ++idx) {
      result = test.search_data(suite[idx]);
      
      cout << "Target " << idx << ": " << suite[idx] << " found at " << result << "." << endl;
      cout << "\tResult--: " << (result==0?0.0:cum_f[result-1]) 
	<< ", Result: " << cum_f[result]
	<< ", Result++: " << (result==1999?2000000.0:cum_f[result+1]) << endl;
   };
   
   cout << "Completed test." << endl;
   
   return 0;
};
