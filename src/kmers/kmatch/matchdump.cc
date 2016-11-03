#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <cstdlib>
#include "matchresult.h"

#include "KMatch.h"
#include <sys/time.h>
#include <thread>
#include <functional>

int main(int argc, char ** argv){

  if (argc!=2) {
    std::cout<<"Usage: "<<argv[0]<<" kmatch_output_file"<<std::endl;
    return -1;
  }

  std::ifstream result_file(argv[1]);
  t_result r;

  while (!result_file.read((char *) &r,sizeof(r)).eof()){
    std::cout<<r.query_id<<" "<<r.target_id<<" "<<r.query_size<<" "<<r.qstart<<" "<<r.tstart<<" "<<r.len<<" "<<r.reverse<<" "<<r.prob<<" "<<r.ident<<std::endl;
  }
  result_file.close();
}
