/* scale.c 

    From original code by Tom Steppe.
    "Well Tempered Linear Scales",  
    Computer Language, September, 1989.
    pp 49-65.
    see <https://archive.org/details/sim_computer-language_1989_6_index/mode/2up?q=steppe>

    This collection of functions provides 4 methods for selecting the
    scale values used to plot an arbitrary data-set.  The four functions
    each take the same arguments, and return values that are used to
    compute the values along a Y axis in a plot.
    
    The DixonKronmal algorithm returns a fixed number of intervals,
    and may have sub-optimal utilization of the data range.
    
    The Lewart's algorithm returns a less rigorously defined number of
    intervals, in return for improved utilization of the range.
    
    The MaxIntervals method returns a number of intervals nearest to
    the specified amount.
    
    ???
    The Internal Labelling algorithm provides high utilization rates
    at the expense of not using well tempered minimum and maximums.
    ???

    $Id: scale.c,v 1.4 2025/04/06 04:08:56 bomr Exp $

    $Log: scale.c,v $
    Revision 1.4  2025/04/06 04:08:56  bomr
    Minor comment typos fixed

    Revision 1.3  2023/08/06 19:50:52  bomr
    Changed some floating point output formatting

    Revision 1.2  2023/08/06 19:33:46  bomr
    Added algorithm: Internal (but not sure how useful it is)

    Revision 1.1  2023/08/06 19:21:26  bomr
    Intial entry to CVS. Seems to be working.


*/

#include    <stdio.h>
#include    <stdlib.h>
#include    <assert.h>
#include    <math.h>
// #include    "scale.h"


char * scale_c = "$Revision: 1.4 $";

/* The set of potential multipliers for well-tempered numbers.      */
/* (10.0 is included only as a convenience for computing geometric  */
/* means for Lewart's algorithm.)                                   */

// static double pdSet[] = { 1.0, 2.0, 5.0 };
static double pdSet[] = { 1.0, 2.0, 5.0, 10.0 };
// static double pdSet[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 10.0 };

#define SET_LEN (sizeof(pdSet)/sizeof(double)-1)

/* Function prototypes 
*/
static double   scFirstNiceNum( double, int *, double * );
static double   scNextNiceNum( double *, int, int *, double * );
static void     scCalcExtLabel( double, double, double, int *, int * );
static void     scCalcIntLabel( double, double, double, int *, int * );
static double   scPower( double, int );

/* Enhanced Dixon-Kronmal algorithm 
*/
void  scDixonKronmal( double   dDataMin,        /* Data Minimum         */   
                      double   dDataMax,        /* Data Maximum         */
                      int      nExactIntervals, /* Exact Number of intervals to use */
                      double * pdScaleMin,      /* Scale Minimum (result returned)  */
                      double * pdScaleMax,      /* Scale Maximum (result returned)  */
                      int *    pnActual         /* Actual number of intervals used  */ 
                      ){
                        
double      dIntervalSize;
int         iIndex;
double      dPowerOfTen, 
            dNiceNum;
int         nLoMult,
            nAdjLoMult,
            nHiMult,
            nAdjHiMult,
            nActualIntervals,
            nDiffIntervals,
            nAdjIntervals;
            
    assert( dDataMin < dDataMax );
    assert( nExactIntervals >= 2 );
    
    /*  Calculate the smallest possible interval size
    */
    dIntervalSize = (dDataMax - dDataMin ) / nExactIntervals;

    /*  Calculate the smallest nice number not smaller the dInterValSize
    */
    for( dNiceNum = scFirstNiceNum( dIntervalSize, &iIndex, &dPowerOfTen );
            dNiceNum < dIntervalSize;
                dNiceNum = scNextNiceNum( pdSet, SET_LEN, &iIndex, &dPowerOfTen ) ){
                ;
    }

    /*  Produce the scale using the nice number
    */
    scCalcExtLabel( dDataMin, dDataMax, dNiceNum, &nLoMult, &nHiMult );

    
    /*  Continue to rescale the data with new nice numbers until the
    **  requested number of intervals is not exceeded
    */
    while( nHiMult - nLoMult > nExactIntervals ){
        dNiceNum = scNextNiceNum( pdSet, SET_LEN, &iIndex, &dPowerOfTen );
        scCalcExtLabel( dDataMin, dDataMax, dNiceNum, &nLoMult, &nHiMult );
    }
                                    
    /*  Calculate the actual number of intervals spanned by the data
    */
    nActualIntervals = nHiMult - nLoMult;
    
    /*  Adjust lo and hi multiples to account for the additional
    **  intervals required.  Adjust in favor of centering.
    */
    nDiffIntervals = nExactIntervals - nActualIntervals;
    nAdjIntervals = nDiffIntervals / 2;
    if( nDiffIntervals & 1 ){
        
        /* nDiffIntervals is ODD.  Decide where the extra interval should go
        */
        if( (dDataMin - nLoMult * dNiceNum) < (nHiMult * dNiceNum -dDataMax) ){
            nAdjIntervals++;
        }
    }
    nAdjLoMult = nLoMult - nAdjIntervals;
    nAdjHiMult = nAdjLoMult + nExactIntervals;
    
    /*  Avoid adjustments that cause negative scales for non-negative data
    */
    if( nAdjLoMult < 0 && nLoMult >= 0 ){
        nAdjLoMult = 0;
        nAdjHiMult = nExactIntervals;
    }
    
    /*  Avoid adjustments that cause positive scales for non- positive data
    */
    if( nAdjHiMult > 0 && nHiMult <= 0 ){
        nAdjHiMult = 0;
        nAdjLoMult = -nExactIntervals;
    }
    
    /*  Calculate scale limits
    */
    *pdScaleMin = nAdjLoMult * dNiceNum;
    *pdScaleMax = nAdjHiMult * dNiceNum;
    *pnActual = nExactIntervals;
    
}


/*  Lewart's algorithm
*/
void    scLewart( double    dDataMin,   /* Input: Data Minimum          */ 
                  double    dDataMax,   /* Input: Data Maximum          */
                  int       nApproxIntervals,   /* Input: Approximate number of intervasl to use */
                  double *  pdScaleMin, /* Ouptut: scale minimum        */
                  double *  pdScaleMax, /* Output: scale maximum        */
                  int *     pnActual
                  ){
                    
double  dIntervalSize;
int     iIndex;
double  dPowerOfTen,
        dNiceNum;
int     nLoMult,
        nHiMult;
        
    assert( dDataMin < dDataMax );
    assert( nApproxIntervals >= 2 );
                        
    /*  Calculate the smallest possible interval size   
    */
    dIntervalSize = ( dDataMax - dDataMin ) / nApproxIntervals;
    
    /*  Find the nice number closest to the smallest potential
    **  interval size.  Use the geometric means of adjacent multiplier
    **  values as break points.
    */
    for( dNiceNum = scFirstNiceNum( dIntervalSize, &iIndex, &dPowerOfTen );
            sqrt( pdSet[iIndex] * pdSet[iIndex-1] ) * dPowerOfTen < dIntervalSize;
                dNiceNum = scNextNiceNum( pdSet, SET_LEN, &iIndex, &dPowerOfTen ) ){
                    
        ;
    }
                      
    /*  Produce the scale using the specified nice number
    */
    scCalcExtLabel( dDataMin, dDataMax, dNiceNum, &nLoMult, &nHiMult );
    
    /* Calculate scale limits
    */
    *pdScaleMin = nLoMult * dNiceNum;
    *pdScaleMax = nHiMult * dNiceNum;
    *pnActual = nHiMult - nLoMult;
    
}

                  
/*  Algorithm for scaling with a maximum number of intervals
*/
void     scMaxInterval( double  dDataMin,   /* Input: Data Minimum          */
                        double  dDataMax,   /* Input: Data Maximum          */
                        int     nMaxIntervals,  /* Input: Maximum number of intervals to use    */
                        double * pdScaleMin,    /* Output: Scale Minimum    */
                        double * pdScaleMax,    /* Output: Scale Maximum    */
                        int    * pnActual       /* Output: Actual number of intervala used  */
                    ){
double  dIntervalSize;
int     iIndex;
double  dPowerOfTen,
        dNiceNum;
int     nLoMult,
        nHiMult;

    assert( dDataMin < dDataMax );
    assert( nMaxIntervals >= 2 );
    
    /*  Calculate the smallest potential interval size
    */
    dIntervalSize = ( dDataMax - dDataMin ) / nMaxIntervals;
    
    /*  Calculate the smallest nicve number not greater than dIntervalSize
    */
    for( dNiceNum = scFirstNiceNum( dIntervalSize, &iIndex, &dPowerOfTen); 
            dNiceNum < dIntervalSize;
                dNiceNum = scNextNiceNum( pdSet, SET_LEN, &iIndex, &dPowerOfTen ) ){
        ;
    }
    
    
    /*  Produce the scale using the specified nice number
    */
    scCalcExtLabel( dDataMin, dDataMax, dNiceNum, &nLoMult, &nHiMult);
    
    /*  Continue to rescale the data with the new nice numbers until the
    **  requested number of intervals is not exceeded
    */
    while( nHiMult - nLoMult > nMaxIntervals ){
        dNiceNum = scNextNiceNum( pdSet, SET_LEN, &iIndex, &dPowerOfTen );
        /* scCalcExtScale( dDataMin, dDataMax, dNiceNum, &nLoMult, &nHiMult ); */
        scCalcExtLabel( dDataMin, dDataMax, dNiceNum, &nLoMult, &nHiMult );
    }
    
    /*  Calculate scale limits
    */
    *pdScaleMin = nLoMult * dNiceNum;
    *pdScaleMax = nHiMult * dNiceNum;
    *pnActual = nHiMult - nLoMult;
    
}


/*  Algorithm for internal labelling
*/
void    scInternal( double  dDataMin,       /* Input: Data Minimum              */
                    double  dDataMax,       /* Input: Data Maximum              */
                    int     nMaxIntervals,  /* Input: Maximum number of intervals to use    */
                    double  * pdRefMin,     /* Output: Minimum reference value      */
                    double  * pdRefMax,     /* Output: Maximum reference value      */            
                    int     * pnActual      /* Output: Actual number of reference value intervals   */
                    ){

double      dIntervalSize;
int         iIndex;
double      dPowerOfTen,
            dNiceNum;
int         nLoMult,
            nHiMult;
            
            
    assert( dDataMin < dDataMax );
    assert( nMaxIntervals >= 5 );
    
    /*  Calculate the smallest potential interval size
    */
    dIntervalSize = ( dDataMax -dDataMin) / nMaxIntervals;
    
    /*  Calculate the smallest nice number not smaller than interval size
    */
    for( dNiceNum = scFirstNiceNum( dIntervalSize, &iIndex, &dPowerOfTen );
            dNiceNum < dIntervalSize;
                dNiceNum = scNextNiceNum( pdSet, SET_LEN, &iIndex, &dPowerOfTen ) ){
        ;
    }
    
    /*  Produce the internal scaling using the specified nice number
    */
    scCalcIntLabel( dDataMin, dDataMax, dNiceNum, &nLoMult, &nHiMult );
    
    /*  Calculate minimum & maximum reference values
    */
    *pdRefMin = nLoMult * dNiceNum;
    *pdRefMax = nHiMult * dNiceNum;
    *pnActual = nHiMult - nLoMult;
    
}


/*  Calculate an initial value for the nice number
*/
double  scFirstNiceNum( double  dIntervalSize,  /* Input: Interval Size     */
                        int   * piIndex,        /* Output: Index into multiplier array for nice number  */
                        double * pdPowerOfTen   /* Output: Power of ten for nice number     */
                        ){
int     iExponent;

    /*  Calculate an intial power of 10
    */
    iExponent = (int)floor( log10( dIntervalSize ) );
    
    /*  Perform some extra checking 
    */
    *pdPowerOfTen = scPower( 10.0, iExponent );
    if( *pdPowerOfTen * 10.0 <= dIntervalSize ){
        *pdPowerOfTen *= 10.0;
    }
    
    /*  Initial index is always zero
    */
    *piIndex = 0;
    
    return( *pdPowerOfTen );
    
}


/*  Calculate the next nice number
*/
double  scNextNiceNum(  double  * pdSet,    /* Input: Set of multipliers    */
                        int       nSet,     /* Input: Number of elements in set  */
                        int     * piIndex,  /* In/Out: Index into multiplier array for nice number  */
                        double  * pdPowerOfTen  /* In/Out: power of ten for nice number */
                        ){

    /*  Increment the index
    */
    (*piIndex)++;
    
    /*  If the maximum index has been exceeded, reset the index to
    **  zero, and increase the power of 10.
    */
    if( *piIndex >= nSet ){
        *piIndex = 0;
        *pdPowerOfTen *= 10.0;
    }
                             
    return( pdSet[*piIndex] * *pdPowerOfTen );
}

/*  Calculate an externally labeled scale
*/
void    scCalcExtLabel( double      dDataMin,   /* Input: Data Minimum          */
                        double      dDataMax,   /* Input: Data Maximum          */
                        double      dNiceNum,   /* Input: Nice Number           */
                        int       * pnLoMult,   /* Output: Multiplier for scale minimum */
                        int       * pnHiMult    /* Output: Multiplier for scale maximum */
                        ){

    /*  Calculate the low multiple
    */
    *pnLoMult = (int) floor( dDataMin / dNiceNum );
    
    /*  Perform some extra checking
    */
    if( (double) (*pnLoMult + 1) * dNiceNum <= dDataMin ){
        (*pnLoMult)++;
    }
    
    /*  Calculate the high multiple
    */
    *pnHiMult = (int) floor( dDataMax / dNiceNum );
    
    /*  Perform some extra checking
    */
    if( (double) (*pnHiMult - 1) * dNiceNum >= dDataMax ){
        (*pnHiMult)--;
    }
}
    
/*  Calculate an internally labelled scale
*/
void    scCalcIntLabel( double  dDataMin,       /* Input: Data minimum      */
                        double  dDataMax,       /* Input: Data maximum      */
                        double  dNiceNum,       /* Input: Nice number       */
                        int   * pnLoMult,       /* Output: Multiplier for minimum reference value   */
                        int   * pnHiMult        /* Output: Multiplier for maximum reference value   */
                        ){

    /*  Calculate the low multiple
    */
    *pnLoMult = (int) ceil( dDataMin / dNiceNum );

    /*  Perform some extra checking
    */
    if( (double)( *pnLoMult - 1) * dNiceNum >= dDataMin ){
        (*pnLoMult)--;
    }

    /*  Calculate the high multiple
    */
    *pnHiMult = (int) ceil( dDataMax / dNiceNum );

    /*  Perform some extra checking
    */
    if( (double)( *pnHiMult + 1) * dNiceNum <= dDataMax ){
        (*pnHiMult)++;
    }
}

/*  Raise a double to an integer power
**  (adapted from an algorithm described in 'Algorithms' by Robert Sedgewick,
**  First Ed., pp 46-47
*/
static double   scPower( double     dRoot,      /* Input: Root to be raised to a power  */
                         int        iExponent   /* Input: Power to which the root should be raised  */
                         ){

double  dResult;

    /*  For negative exponents, invert root and use a positive exponent
    */
    if( iExponent < 0 ){
        dRoot = 1.0 / dRoot;
        iExponent = -iExponent;
    }
    
    /*  Perform multiple multiplications
    */
    dResult = 1.0;
    while( iExponent ){
        if( iExponent & 1 ){
            dResult *= dRoot;
        }
        
        iExponent >>= 1;
        
        if( iExponent ){
            dRoot *= dRoot;
        }
    }
    return( dResult );
}
