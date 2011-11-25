#ifndef misc_INC
#define misc_INC

#include <string>
#include <vector>
#include <fstream>
#include <cctype>


typedef  double real;



/*!
\fn string_To_lower
    
   Convert upper case to lower case
   WW 04.2011 
*/

namespace AuxFunctions
{
std::string string_To_lower(std::string strToConvert);
/*!
    \brief Binary search a array
     The array must be sorted
     
     \param arr     an array 
     \param target  searched index
     \param start   the start index of the array 
     \parm  end     the end index of the array
     
     By WW. 03.2011
    
*/
long binarySearch(long *arr, long target, long start, long end); 

/*!
\fn Read_Block
    
   Read a block of data pairs 

   WW 04.2011
*/
void Read_Block(std::ifstream &ins, std::vector<std::string> &key, std::vector<float> & key_value);

/*!
\fn  DeleteArray(num *an_array)
    
  Release the memory of arrary allocated by using 'new'. 

  10.2010
*/
/*!
\fn  DeleteArray(num *an_array)
    
  Release the memory of arrary allocated by using 'new'. 

  10.2010
*/
template<class num> void  DeleteArray(num *an_array)
{
   if(an_array) delete [] an_array;
   an_array = NULL; 
}
template<class num> void  DeleteVector(std::vector<num*> &a_vec)
{
   
   while (a_vec.size()>0)
   {
      delete a_vec[(int)a_vec.size()-1];
      a_vec.pop_back();
   }
   
}

/// Check comment
bool CheckComment(std::string& string_line);

/// returns used heap size in bytes or negative if heap is corrupted.
long HeapUsed();

double ComputeDetTri(const double *x1, const double *x2,
                                const double *x3);


}

#endif

