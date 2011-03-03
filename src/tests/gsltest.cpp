/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

int main(int argc, char *argv[])
{
  const size_t NDIM = 3;
  double data[NDIM * NDIM];
  gsl_matrix_view m = gsl_matrix_view_array(data, NDIM, NDIM);
  gsl_vector *eval = gsl_vector_alloc (NDIM);
  gsl_matrix *evec = gsl_matrix_alloc (NDIM, NDIM);
  gsl_eigen_symmv_workspace * w = 
    gsl_eigen_symmv_alloc (NDIM);
  gsl_eigen_symmv (&m.matrix, eval, evec, w);
  gsl_eigen_symmv_free (w);
  //Sort Eigenvalues in ascending order
  gsl_eigen_symmv_sort (eval, evec, 
			GSL_EIGEN_SORT_VAL_ASC);
  
//  for (size_t i = 0; i < NDIM; i++)
//    {
//      retVal.EigenVal[i] = gsl_vector_get(eval, i) / range->size();
//
//      gsl_vector_view evec_i 
//	= gsl_matrix_column (evec, i);
//
//      //EigenVec Components
//      for (size_t j = 0; j < NDIM; j++)
//	retVal.EigenVec[i][j] = gsl_vector_get(&evec_i.vector, j);
//    }

  //Cleanup GSL
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
}

