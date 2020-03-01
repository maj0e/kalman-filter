#include "kalman.h"
#include <stdlib.h>

lah_Return kal_correct(kal_Filter *kal, const lah_mat *z)
{
    lah_Return result = lahReturnOk;
    /*lapack_int *ipiv;   */ 
    
    if (z == NULL || kal == NULL || z->nR != kal->zp->nR )
    {
        return lahReturnParameterError;
    }
    
    /* predict measurement and get Jacobi of Measurement: 
     * zp = h(x), Jh = dh(x)/dx */
    if (kal->h == NULL
        || lahReturnOk != kal->h(kal->zp, kal->x, 
                                 kal->Jh, kal->optData))
    {
        return lahReturnFPError;
    }
    
    /* cholesky factorization L of innovation covariance 
     * S = (Jh*P*Jh' + R) = L*L' */ 
    /* Impl: W = P*Jh', S = Jh*W, S = S + R */
    /* After this operation Ws holds Ws = P *Jh' */
    /*ipiv = calloc(z->nR * 2, sizeof(lapack_int));*/
    
    result += lah_matUpdate(kal->S, kal->P, kal->Jh, kal->R, kal->WView_U);
    /* assumes S is symmetric positive-definite for 
     * efficiency of cholesky factorization
     * Only the lower triangular matrix is written
     * */
    if(lahReturnOk != /* lah_LU(kal->S, 0, ipiv)) */ lah_chol(kal->S, 0))
    {  
        return result;
    }

    /* compute intermediate matrix U used for the updates */
    /* K = U / L -> U = (L \ (P *Jh')')'      */
    /* In implementation --> WView_U_trans = S \ WView_U_trans      */

    if( lahReturnOk != /*lah_solveLU(kal->WView_U_trans, kal->S, ipiv)) */
        lah_forwardSub(kal->WView_U_trans, kal->S))
    {  
        return result;
    }

    /* correct state vector
     * x = x + U * L \ (z - zp) 
     * Impl: x = x + U * (S \(z -zp)) */
    /* zw = z - zp */
    
    result += lah_matAdd(kal->zw, 1.0, z, -1.0, kal->zp);
    
    /* z = S \ zw */
    if(lahReturnOk != lah_forwardSub(kal->zw, kal->S))
    {  
        return result;
    }

    /* x = x + U * zw */
    result += lah_matMul(lahNorm, lahNorm, 1.0, 1.0, 
                         kal->x, kal->WView_U, kal->zw);
    /* correct covariance */
    /* P = Pp - U * U' */
    /* Use Ut_ gain to prevent some overhead */
    result += lah_matMul(lahNorm, lahNorm, 1.0, -1.0, 
                        kal->P, kal->WView_U, kal->WView_U_trans);

    if (result == lahReturnOk)
        return lahReturnOk;
    else
        return lahReturnMathError;

}
