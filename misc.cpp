/*!
\fn string_To_lower
    
   Convert upper case to lower case 
*/
#include <iostream>
#include <sstream>
#include "misc.h"
 string string_To_lower(string strToConvert)
{
   for(unsigned int i=0;i<strToConvert.length();i++)
   {
      strToConvert[i] = tolower(strToConvert[i]);
   }
   return strToConvert;
}; 


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
              key[i] = string_To_lower(key[i]);
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

