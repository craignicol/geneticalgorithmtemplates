#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "MGalgo.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
   vector<int> sorted_list;
   
   sorted_list.push_back(1);
   sorted_list.push_back(4);
   sorted_list.push_back(8);
   sorted_list.push_back(12);
   sorted_list.push_back(15);
   sorted_list.push_back(22);
   sorted_list.push_back(25);
   sorted_list.push_back(30);
   
   std::cout << mgtl::binary_search(20, sorted_list) << ' ';
   std::cout << mgtl::binary_search(0, sorted_list) << ' ';
   std::cout << mgtl::binary_search(40, sorted_list) << ' ';
   std::cout << mgtl::binary_search(8, sorted_list) << ' ';
   std::cout << mgtl::binary_search(3, sorted_list) << std::endl;

   vector<int> incident_3;

   incident_3.push_back(10);
   incident_3.push_back(20);

   std::cout << mgtl::binary_search(5, incident_3) << std::endl;
   
   return EXIT_SUCCESS;
};
