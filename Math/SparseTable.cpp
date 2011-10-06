/*========================================================================
 GeoSys - class Matrix (Definition)
          class vec  
 Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation. 
 Function:   See the definition below
 programming:
  22/08/2004  WW  
==========================================================================*/

/// Matrix
#include "SparseTable.h"
#include <iomanip>
#include <float.h>
#include <cmath>
//

///PDE related classes:
#include "FiniteDifference.h"


namespace Math_Group{

////PDE related class:
//using _FDM::Point;


/*\!
********************************************************************
   Create sparse matrix table
   01/2006 WW
   08/2007 WW
   10/2007 WW
   03/2010 WW: CRS storage
********************************************************************
*/
SparseTable::SparseTable(FiniteDifference *fdm)
{
  
   long i=0, j=0;
   long **larraybuffer;
   larraybuffer = NULL;
   storage_type = CRS;
   symmetry = false;  
 
   //
   rows = (long)fdm->grid_point_in_use.size();

   size_entry_column = 0;
   diag_entry = new long[rows]; 

   row_index_mapping_n2o = NULL; 
   row_index_mapping_o2n = NULL;
 
   
   /// CRS storage
   /// num_column_entries saves vector ptr of CRS 
   num_column_entries = new long[rows+1];

   vector<long> A_index;
   long col_index;

   for(i=0; i<rows; i++)
   {
      num_column_entries[i] = (long)A_index.size();
 
      for(j=0; j<fdm->grid_point_in_use[i]->GetNumNeighborPoints(); j++)
      {
         col_index = fdm->grid_point_in_use[i]->GetNeighborIndex(j);
                        
         if(i == col_index)
            diag_entry[i] = (long)A_index.size();
         A_index.push_back(col_index);
      }
   }
      
   size_entry_column = (long)A_index.size(); 
   num_column_entries[rows] = size_entry_column;

   entry_column = new long[size_entry_column]; 
   for(i=0; i<size_entry_column; i++)
      entry_column[i] = A_index[i];

}
/*\!
********************************************************************
   Create sparse matrix table
   08/2007 WW
   10/2007 WW
********************************************************************/
void SparseTable::Write(ostream &os)
{
   long i, k, counter=0;

   os.width(10);
   os<<"Symmetry: "<<symmetry<<endl;
   os<<"\n*** Row index  "<<endl;
 
   if(storage_type == CRS)
   {
      os<<"\n*** Sparse entry  "<<endl;
      for (i = 0; i < rows; i++)
      {
          for (k = num_column_entries[i]; k < num_column_entries[i+1]; k++)
             os<<entry_column[k]+1<<" ";
          os<<endl; 
      }
   } 
   else if(storage_type == JDS)
   {
      for (i = 0; i < rows; i++)
        os<<row_index_mapping_n2o[i]+1<<endl;
      // 
      os<<"\n*** Sparse entry  "<<endl;
      for (k = 0; k < max_columns; k++)
      {
         os<<"--Column: "<<k+1<<endl;
         for (i = 0; i < num_column_entries[k]; i++)
         {          
            os<<entry_column[counter]+1<<endl;;
            counter++;
         } 
         os<<endl;        
      } 
   } 
}    

/*\!
********************************************************************
   Create sparse matrix table
   08/2007 WW
   10/2007 WW
********************************************************************/
SparseTable::~SparseTable()
{
  if(entry_column) delete [] entry_column;
  if(num_column_entries) delete [] num_column_entries; 
  if(row_index_mapping_n2o) delete [] row_index_mapping_n2o;    
  if(row_index_mapping_o2n) delete [] row_index_mapping_o2n;    
  if(diag_entry) delete [] diag_entry;
  entry_column = NULL;
  num_column_entries = NULL; 
  row_index_mapping_n2o = NULL;    
  row_index_mapping_o2n = NULL;    
  diag_entry = NULL;
}

}// Namespace

