#include "kalman.h"
#include <stdlib.h>
#include <string.h>

lah_Return kal_predict(kal_Filter *kal, const lah_mat *u)
{
    lah_Return result;
    /*lah_mat *temp = NULL; */
    
    if ( kal == NULL)
    {
        return lahReturnParameterError;
    }

    /* predict state and get Jacobi matrix (extended filter) 
     * x1 = f(x,u), Jf = df(x,u)/dx */
    if (kal->f == NULL 
        || lahReturnOk != kal->f(kal->xp, kal->x, u, 
                                 kal->Jf, kal->optData))
    {
        return lahReturnFPError;
    }
    
    /* Swap prediction with state */
    /*
    temp = kal->x;
    kal->x = kal->xp; 
    kal->xp = temp;
    */
    memcpy(kal->x->data, kal->xp->data, kal->xp->nR * sizeof(lah_value));
    
    /* predict covariance P = A*P*A' + Q */
    result = lah_matUpdate(kal->P, kal->P, kal->Jf, kal->Q, kal->WView_P);
    
    return result;
}


