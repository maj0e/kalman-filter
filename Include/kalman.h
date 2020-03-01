#ifndef _KALMAN_H_
#define _KALMAN_H_
#include <lah.h>

#ifdef __cplusplus
extern "C" {
#endif
/*****************************************************************************/
/*  Implementation of the simple/extended kalman filter                      */
/*****************************************************************************/

/*--- Definitions for the Kalman-Filter -------------------------------------*/

/* kal_StateFun
 * -------------
 * Function type to compute the predicted state of the Kalman-Filter
 * This function can also change the state prediction Matrix of the
 * filter. In which case the filter is an extended Kalman-Filter and
 * the Matrix Jf is the Jacobi Matrix describing the
 * linearization from current filter state.
 */
typedef lah_Return (*kal_StateFun)(lah_mat *xp, const lah_mat *x,
                                   const lah_mat *u, 
                                   lah_mat *Jf, void *optData);

/* kal_MeasFun
 * -------------
 * Measurement function for the (extended) Kalman filter.
 * Describes the function z = h(x), relating state to measurement.
 * This function can also change the measurement prediction Matrix of
 * the filter. In which case the filter is an extended Kalman-Filter 
 * and the Matrix Jh is the Jacobi Matrix describing the
 * linearization of the current measurement.
 */
typedef lah_Return (*kal_MeasFun)(lah_mat *zp, const lah_mat *x, 
                                  lah_mat *Jh, void *optData);

/* Data structure for the kalman filter */
typedef struct
{
    lah_mat *x;            /* state vector */
    lah_mat *xp;           /* predicted state vector */
    lah_mat *zp;           /* pred. meas. vector */
    lah_mat *zw;           /* workspace meas. vector */
    
    lah_mat *P;            /* covariance (nState x nState) */

    lah_mat *S;            /* innovation covariance */
    int *ipiv;       /* Permutations for LU factorization */
    lah_mat *W;            /* Workspace for operations */
    
    /* Noise Matrices */
    lah_mat *Q;
    lah_mat *R;
    
    /* Matrix Views wihtout own data*/
    lah_mat *WView_P;
    lah_mat *WView_U;
    lah_mat *WView_U_trans;

    /* Propagation */
    lah_mat *Jf;           /* Jacobi of state transition */
    lah_mat *Jh;           /* Jacobi of state transition */
    kal_StateFun f;        /* state prediction function */
	kal_MeasFun h;	       /* measurement prediction function */
	
    void *optData;         /* pointer to optional user data for 
                               prediction functions */

} kal_Filter;

/*--------------------------------------------------------------------------*/

/*--- Kalman-Filter functions ----------------------------------------------*/    
/* For details see "Dan Simon - Optimal State estimation"                   */
/* The implementation design is inspired by                                 */ 
/* https://github.com/dr-duplo/eekf/blob/master/                            */
/* because it can use the linear and extended kalman filter at once         */
/* This is achieved by defining an interface for the measurement/state      */
/* prediction function.                                                     */
/*--------------------------------------------------------------------------*/

/* kal_predict:
 * ------------
 * This function predicts the next filter state using the current state, 
 * the given action u and the given process covariance. 
 * It will call the function f for state transition. 
 * The dimensions of process noise Q should match the dimensions 
 * of covariance P (nState x nState).
 */
lah_Return kal_predict(kal_Filter *kal, const lah_mat *u);

/* kal_correct:
 * ------------
 * This function corrects the current filter state by using the given 
 * measurement z and the given measurement process noise covariance R.
 * The dimensions of measurement noise R should match the dimensions 
 * of covariance P (nMeas x nMeas).
 */
lah_Return kal_correct(kal_Filter *kal, const lah_mat *z);

/*--- Utility/Helper functions ---------------------------------------------*/

/* Allocates all data for the Kalman struct*/
kal_Filter *kal_create(lah_index nState, lah_index nMeas,
                       kal_StateFun f, kal_MeasFun h, 
                       void *optData);

/* Free all data used by the Kalman-Filter struct*/
kal_Filter *kal_free(kal_Filter *kal);

/* Returns a View on the workspace kal->W compatible to matrix A */
lah_mat *kal_getWorkspaceView(kal_Filter *kal, lah_mat *A);

/*--------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif /* _KALMAN_H_ */
