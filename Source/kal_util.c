#include "kalman.h"
#include <math.h>
#include <stdlib.h>

/* allocation for the feature vector feat */
kal_Filter *kal_create(lah_index nState, lah_index nMeas,
                       kal_StateFun f, kal_MeasFun h,
                       void *optData)
{
    /* Allocate the struct */
    size_t max = 0;
    kal_Filter *kal = NULL;
    
    kal = (kal_Filter *) calloc(1, sizeof(kal_Filter));
    if (!kal) return (NULL) ;     /* out of memory */
    max = (size_t) LAH_MAX(nState, nMeas);
    /*--- Allocate matrices of the Kalman Filter ------------------*/
    /* State vector*/
    kal->x = lah_matAlloc(1, nState, 1);
    /* Predicted state vector*/
    kal->xp = lah_matAlloc(1, nState, 1);
    
    /* Pred. meas. vector*/
    kal->zp = lah_matAlloc(1, nMeas, 1);
    /* Workspace for meas. vector*/
    kal->zw = lah_matAlloc(1, nMeas, 1);
    
    /* Covariance Matrix*/
    kal->P = lah_matAlloc(nState, nState, 1);
    
    /* Matrices for updating filter and calculations */ 
    /* innovation covariance */
    kal->S = lah_matAlloc(nMeas, nMeas, 1);
#ifdef HAVE_LAPACK
    /*kal->ipiv = malloc(max * sizeof(lapack_int)); */
#else
    kal->ipiv = NULL;
#endif
    /* Update matrix*/
    
    /* workspace max(N,M)^2*/
    kal->W = lah_matAlloc(max, max, 1);
    
    kal->Q = lah_matAlloc(nState, nState, 1);
    kal->R = lah_matAlloc(nMeas, nMeas, 1);
    
    
    /* Prediction functions */
    /* Jacobi of state transition */
    kal->Jf = lah_matAlloc(nState, nState, 1); 
    /* Jacobi of measurement transition */
    kal->Jh = lah_matAlloc(nState, nMeas, 1);
    
    kal->f = f;
    kal->h = h;
    kal->optData = optData;
    
    
    kal->WView_P = kal_getWorkspaceView(kal, kal->P);
    kal->WView_U = lah_matAlloc(nMeas, nState, 0);
    kal->WView_U->data = kal->W->data;
    kal->WView_U_trans = lah_matTrans(kal->WView_U);
    
    /* Do some checks if everything worked*/
    if ( kal->x == NULL || kal->xp == NULL || kal->zp == NULL || kal->zw == NULL
        || kal->P == NULL || kal->Jf == NULL || kal->Jh == NULL || kal->S == NULL 
        || /*kal->ipiv == NULL ||*/ kal->W == NULL || kal->Q == NULL 
        || kal->R == NULL || kal->WView_P == NULL || kal->WView_U == NULL || kal->WView_U_trans == NULL )
    {
        return kal_free(kal);
    }
    else
        return kal;
}

kal_Filter *kal_free(kal_Filter *kal)
{
    if (!kal) return NULL;      /* do nothing if kal already NULL */
    
    /* Free all data of the Kalman-Filter 
     * According to C standard free(ptr) does nothing 
     * if ptr is already NULL
     * */
    lah_matFree(kal->x);
    lah_matFree(kal->xp);
    lah_matFree(kal->zp);
    lah_matFree(kal->zw);
    lah_matFree(kal->P);
   
    lah_matFree(kal->S);
    lah_matFree(kal->W);

    free(kal->WView_P);
    free(kal->WView_U);
    free(kal->WView_U_trans);
    
    lah_matFree(kal->Q);
    lah_matFree(kal->R);
    
    lah_matFree(kal->Jh);
    lah_matFree(kal->Jf);
    
    /*free(kal->ipiv); */
    
    /* free kal itself*/
    free(kal);
    return NULL; 
}

/* Returns a View on the workspace kal->W compatible to matrix A */
lah_mat *kal_getWorkspaceView(kal_Filter *kal, 
                              lah_mat *A)
{
    lah_mat *WView;
    if (A == NULL || kal == NULL || kal->W == NULL 
        || A->nC > kal->W->nC || A->nR > kal->W->nR)
    {
        return NULL;
    } 
    
    WView = lah_matAlloc(A->nC, A->nR, 0);

    WView->data = kal->W->data;
    return WView;
}
