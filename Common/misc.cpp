/*!
\fn string_To_lower
    
   Convert upper case to lower case 
*/
#include <iostream>
#include <sstream>
#include <malloc.h>
#include <cmath>

#include "misc.h"


using namespace std;

namespace AuxFunctions
{
 string string_To_lower(string strToConvert)
{
   for(unsigned int i=0;i<strToConvert.length();i++)
   {
      strToConvert[i] = tolower(strToConvert[i]);
   }
   return strToConvert;
}

/*!
    \brief Binary search a array
     The array must be sorted
     
     \param arr     an array 
     \param target  searched index
     \param start   the start index of the array 
     \parm  end     the end index of the array
     
     By WW. 03.2011
    
*/
long binarySearch(long *arr, long target, long start, long end) 
{
   long middle;
   while (start <= end) 
   {
      middle = (start + end) / 2;
      if (arr[middle] == target)
         return middle;
      else if (arr[middle] > target)
         end = middle - 1;
      else
         start = middle + 1;
   }
   return -1;
}




/*!
\fn Read_Block
    
   Read a block of data pairs 
   
   \param ins : file stream
   \param key : keyword
   \param key_value : value to keyword


   WW 04.2011
*/
void Read_Block(ifstream &ins, vector<string> &key, vector<float> & key_value)
{
     string aline;
     //ios::pos_type position;
     std::stringstream ss;

     int i, counter, size;

     size = (int)key.size();
     vector<bool> foundkey(size); 
     for(i=0; i<size; i++)
        foundkey[i] = false;

     counter = 0;
     while(!ins.eof())
     {       
 
        //position = ins.tellg();
        getline(ins, aline); 
        if(CheckComment(aline))
           continue;
        
        /* 
        if(aline.find("---")!=string::npos || aline.find("...")!=string::npos)
        {
           ins.seekg(position);
           return;
        }
        */ 

        for(i=0; i<size; i++)
        {
           if(foundkey[i])
             continue;
   
           if(aline.find(key[i])!=string::npos)
           {
              counter++;
              ss.str(aline);
              ss>> key[i];
			  key[i] = AuxFunctions::string_To_lower(key[i]);
              foundkey[i] = true;
              ss>> key_value[i];
              ss.clear();
              break; // Exit for(i=0; i<size; i++)          
           }           
           
        }
        
        if(counter == size)
          return;
                 
     } 
    
     /// Output message 
     for(i=0; i<size; i++)
     {
        if(!foundkey[i])
          cout<<"Keyword "<<key[i]<<" is not correctly defined"<<endl;
     }

      
}




/// returns used heap size in bytes or negative if heap is corrupted.
#ifdef WIN
long HeapUsed()
{
    _HEAPINFO info = { 0, 0, 0 };
    long used = 0;
    int rc;

    while ((rc=_heapwalk(&info)) == _HEAPOK)
    {
        if (info._useflag == _USEDENTRY)
            used += (long)info._size;
        }
    if (rc != _HEAPEND && rc != _HEAPEMPTY)
        used = (used?-used:-1);

    return used;
}
#endif

/*!
   \fn  bool CheckComment(string& string_line)

   Check whether string_line is a comment or contains comment 

   10.06.2011. WW 
*/
bool CheckComment(string& string_line)
{
   // The string is a comment
   if(string_line.c_str()[0]=='#')
       return true;

   // The string contains comment. Remove the comment
   string_line = string_line.substr(0, string_line.find("#"));  
   return false;      
}


/*!
     Computer area of a triangle
*/
double ComputeDetTri(const double *x1, const double *x2,
                                const double *x3)
{
    static double u[3], v[3], z[3];
    
    u[0] = x3[0] - x1[0];	
    u[1] = x3[1] - x1[1];
    u[2] = x3[2] - x1[2];

    v[0] = x2[0] - x1[0];	
    v[1] = x2[1] - x1[1];
    v[2] = x2[2] - x1[2];
 
    z[0] = u[1]*v[2] - u[2]*v[1];
    z[1] = u[2]*v[0] - u[0]*v[2];
    z[2] = u[0]*v[1] - u[1]*v[0];

    return 0.5*sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2] );   
} 


} //namespace

