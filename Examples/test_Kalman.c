/*******************************************************************/
/**  Example to test the Kalman-Filter functions                  **/
/*******************************************************************/
#include "kalman.h"
#include "lah.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*the state prediction function: linear case for simplicity */

const lah_value accel   = 0.4;
const lah_value deltaT  = 0.1;
const lah_value sigma_x = 1.5;
const lah_value sigma_z = 10;

struct optData 
{
    lah_index inits;
    lah_index initm;
    lah_mat* B; /* Control matrix*/
};

lah_Return timeEvolve(lah_mat* xp, const lah_mat *x,
                      const lah_mat *u, lah_mat *A, 
                      struct optData *Data)
{
	/* Set the state transition Matrix A
     * if not done yet*/
    lah_mat *B = Data->B;
    lah_Return result;
    if (Data->inits == 0 )
    {
        A->data[0] = 1;
        A->data[A->incRow] = 0;
        A->data[A->incCol] = deltaT;
        A->data[A->incCol + A->incRow] = 1;
        
        B->data[0] = deltaT * deltaT /2;
        B->data[1] = deltaT;
        Data->inits = 1;
    }

    /* compute predicted state */
    result = lah_matMul(lahNorm, lahNorm, 0.0, 1.0 , xp, A, x);
    result += lah_matMul(lahNorm, lahNorm, 1.0, 1.0, xp, B, u);
        
    if (result != 0)
        return lahReturnMathError;
    else
        return lahReturnOk;
}

/*  measurement prediction function */
lah_Return measure(lah_mat *zp, const lah_mat *x,
                   lah_mat *H, struct optData *Data)
{
    /* Set the state transition Matrix A
     * if not done yet*/
    if (Data->initm == 0 )
    {
        H->data[0] = 1;
        H->data[H->incCol] = 0 ;

        Data->initm = 1;
    }

    /* compute the measurement from state x */
    zp->data[0] = x->data[0];

    return lahReturnOk;
}

int main ()
{
    lah_index i;
    lah_value vel = 0, pos = 0;
    lah_value noise;
    lah_mat *u = lah_matAlloc(1, 1, 1);
    lah_mat *z = lah_matAlloc(1, 1, 1);
    kal_Filter* kal;

    const lah_index nState = 2;
    const lah_index nMeas = 1;
    /*Initialize the struct for prediction */
    struct optData Data;
    Data.initm = 0;
    Data.inits = 0;
    Data.B = lah_matAlloc(1, 2, 1);
    kal = kal_create(nState, nMeas,
                           (kal_StateFun) timeEvolve, 
                           (kal_MeasFun) measure, (void*) &Data);
    
    /* Initialization */ 
    kal->P->data[0] = pow(sigma_x, 2) * pow(deltaT, 4) / 4;
    kal->P->data[kal->P->incCol] = pow(sigma_x, 2) * pow(deltaT, 3) / 2;
    kal->P->data[kal->P->incRow] = pow(sigma_x, 2) * pow(deltaT, 3) / 2; 
    kal->P->data[kal->P->incRow + kal->P->incCol]
        = pow(sigma_x, 2) * pow(deltaT, 2);
    /* input and process noise variables */
    u->data[0] = accel;
    kal->R->data[0] = sigma_z * sigma_z;
    LAH_ENTRY(kal->Q, 0, 0) = pow(sigma_x, 2) * pow(deltaT, 4) / 4;
    LAH_ENTRY(kal->Q, 0, 1) = pow(sigma_x, 2) * pow(deltaT, 3) / 2;
    LAH_ENTRY(kal->Q, 1, 0) = pow(sigma_x, 2) * pow(deltaT, 3) / 2;
    LAH_ENTRY(kal->Q, 1, 1) = pow(sigma_x, 2) * pow(deltaT, 2);

    /* initialize random number generator */
    srand(0);

    /* print header */
    printf("i pos vel realpos realvel P11 P12 P21 P22 z zp\n");
    /* loop over time and demonstrate Kalman filter */
    for (i = 0; i < 5000; i++)
    {
        /* virtual measurement */
        noise = lah_gaussNoise();
        pos = u->data[0] / 2.0 * pow(i * deltaT, 2.0);
        vel = u->data[0] * i * deltaT;
	z->data[0] = pos + noise * sigma_z;
        
        /* correct filter state */
        if (lahReturnOk != kal_correct(kal, z))
            printf("Error in kal_correct !!!\n");
        /* print out */
        printf("%d %g %g %g %g %g %g %g %g %g %g\n", i, 
                kal->x->data[0], kal->x->data[1], pos, vel,
                kal->P->data[0], kal->P->data[1],
                kal->P->data[2], kal->P->data[3],
                z->data[0]     , kal->zp->data[0]);

		u->data[0] = accel;
        /* predict the next filter state */
        if (lahReturnOk != kal_predict(kal, u))
            printf("Error in kal_predict !!!\n");
    }

    lah_matFree(u);
    lah_matFree(z);
    lah_matFree(Data.B);
    kal_free(kal);

    return 0;
}
