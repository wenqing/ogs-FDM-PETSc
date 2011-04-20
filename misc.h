#ifndef misc_INC
#define misc_INC

#include <string>
#include <vector>
#include <fstream>
#include <cctype>
typedef  double real;

using namespace std;

extern string file_name;
extern string file_path;

/*!
\fn string_To_lower
    
   Convert upper case to lower case
   WW 04.2011 
*/
extern string string_To_lower(string strToConvert);

/*!
    \brief Binary search a array
     The array must be sorted
     
     \param arr     an array 
     \param target  searched index
     \param start   the start index of the array 
     \parm  end     the end index of the array
     
     By WW. 03.2011
    
*/
extern long binarySearch(long *arr, long target, long start, long end); 

/*!
\fn Read_Block
    
   Read a block of data pairs 

   WW 04.2011
*/
extern void Read_Block(ifstream &ins, vector<string> &key, vector<float> & key_value);

#endif