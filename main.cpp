// 2D implicit FDM for steady state
// Neumann BC
// 03.2011. WW
//
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

#include "mat.h"
#include "geo.h"
#include "fdm.h"
using namespace std;
using namespace _FDM;

string file_name;
string file_path;



int main ( int argc, char *argv[] )
{
  const int max_size = 1028; 
  char str1[max_size];
  

  if(argc>1) 
     strcpy(str1,argv[1]);
  else 
  {
     cout<<"\tInput file name (without extension): ";
     scanf("%s%*[^\n]%*c",str1);
  }
  clock_t time_cpu;
  time_cpu = -clock();


  // Name and path of the input data file
  file_name = str1;
  basic_string <char>::size_type indexChWin, indexChLinux; 
  indexChWin = indexChLinux = 0;
  indexChWin = file_name.find_last_of('\\');
  indexChLinux = file_name.find_last_of('/');
  //
  if(indexChWin!=string::npos)
     file_path = file_name.substr(0,indexChWin)+"\\";
  else if(indexChLinux!=string::npos)
     file_path = file_name.substr(0,indexChLinux)+"/";


  GeoRead();
  FiniteDifference *fdm = new FiniteDifference();
  fdm->TimeSteping();
    
#ifdef TEST_OUT  ///TEST 
  string fname = file_name+"_dat.out";
  ofstream os(fname.c_str(), ios::trunc);
  fdm->Write(os);
  os.close();

  fname = file_name+"_geo.out";
  os.open(fname.c_str(), ios::trunc);
  WriteGeoData(os);
  os.close();

#endif


  cout<<"\n\tMemory usage: "<< HeapUsed()/1024<<"KB"<<endl;


  delete fdm;
  GeoReleaseMemory();

  time_cpu += clock();
  cout<<"\tCPU time elapsed: "  <<(double)time_cpu / CLOCKS_PER_SEC<<"s"<<endl;

  return 0;
}
