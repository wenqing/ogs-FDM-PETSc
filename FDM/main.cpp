// 2D implicit FDM for steady state
// Neumann BC
// 03.2011. WW
//
#include <iostream>
#include <fstream>
#include <string.h>
#include <string>
#include <time.h>

#ifdef USE_PETSC
#include "petscksp.h"
#endif


#include "MatProperty.h"
#include "Polyline.h"
#include "FiniteDifference.h"
using namespace std;
using namespace _FDM;




int main ( int argc, char *argv[] )
{


  const int max_size = 1028; 
  char str1[max_size];
  string file_name;
  string file_path;
  

  if(argc>1) 
     strcpy(str1,argv[1]);
  else 
  {
     cout<<"\tA 2-D FDM groundwater flow simulator (by WW@UFZ) "<<endl;
     cout<<"\tV2.0. 09.2011 "<<endl;
     cout<<"\tInput file name (without extension): ";
     scanf("%s%*[^\n]%*c",str1);
  }


#ifdef USE_PETSC
  int rank, r_size;
  PetscLogDouble v1,v2;
  char help[] = "FDM with PETSc \n";
  //PetscInitialize(argc, argv, help);
  PetscInitialize(&argc,&argv,(char *)0,help);
  PetscGetTime(&v1);

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &r_size);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Number of CPUs: %d, rank: %d\n", r_size, rank);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

#else
  clock_t time_cpu;
  time_cpu = -clock();
#endif




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


  FiniteDifference *fdm = new FiniteDifference(file_path, file_name);
#ifdef USE_PETSC
  fdm->set_MPI_rank_size(rank, r_size);
#endif


  fdm->Initialize();
  fdm->TimeSteping();

    
#ifdef TEST_OUT  ///TEST 
  string fname = file_name+"_dat.out";
  ofstream os(fname.c_str(), ios::trunc);
  fdm->Write(os);
  os.close();


#endif

    
#ifdef WIN
  cout<<"\n\tMemory usage: "<< AuxFunctions::HeapUsed()/1024<<"KB"<<endl;
#endif

  delete fdm;

#ifdef USE_PETSC
  PetscGetTime(&v2);
  PetscPrintf(PETSC_COMM_WORLD,"\t\n>>Total elapsed time:%f s\n",v2-v1);
#else
  time_cpu += clock();
  cout<<"\tCPU time elapsed: "  <<(double)time_cpu / CLOCKS_PER_SEC<<"s"<<endl;
#endif


#ifdef USE_PETSC

   PetscFinalize();
#endif


  return 0;
}
