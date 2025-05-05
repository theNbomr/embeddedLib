#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

            
typedef struct linConv{
    // Instance parameters
    u_int16_t   rawLow;
    double   engLow;
    u_int16_t   rawHigh;
    double   engHigh;
    
    // Computed values
    // double   rawRange;
    //double   engRange;
    double   rawOffset;
    double   engSlope;
    double   engOffset;
} linConv_t;
            
static int calcScale( linConv_t * linConv );
static int rawToEng( linConv_t * linConv, u_int16_t rawVal, double * engVal );
            
int main( int argc, char * argv[] ){
int     ret;
linConv_t   linConv;
    
    sscanf( argv[1], "%i",  &linConv.rawLow );
    sscanf( argv[2], "%i",  &linConv.rawHigh );
    sscanf( argv[3], "%lf", &linConv.engLow );
    sscanf( argv[4], "%lf", &linConv.engHigh );

    printf( "Raw Low: %d,  RawHigh: %d,  Eng Low: %g, Eng High: %g\n", linConv.rawLow, linConv.rawHigh, linConv.engLow, linConv.engHigh );
    
    ret = calcScale( &linConv );
    
    if( 0 == ret ){
        printf( "(%d) Eng Slope: %g, Eng Offset: %g rawOffset %g\n\n", ret, linConv.engSlope, linConv.engOffset, linConv.rawOffset );
    }
    else{
        printf( "Error: raw range = 0\n" );
        exit( 1 );
    }
    
    int i, incr;
    double engVal;
    incr = ((u_int16_t)(linConv.rawHigh - linConv.rawLow)>>4);
    
    for( i = linConv.rawLow; i <= linConv.rawHigh; i += incr ){        
        ret = rawToEng( &linConv, i, &engVal );
        printf( "%i\t%g\n", i, engVal);
    }
    i = linConv.rawHigh;
    ret = rawToEng( &linConv, i, &engVal );
    printf( "\n%i\t%g\n", i, engVal);
    
    exit( 0 );
}


/*
 * Compute the linear scaling parameters need to convert
 * raw ADC (unsigned 16-bit integer) to engineering units,
 * using the ranges of values that the ADC maps to the 
 * respective Engineering Units scale.
 * 
 * Used to convert, for example, 10-bit ADC to 0-10VDC real world units,
 * or 16-bit ADC to -50C to +50C temperatures.
 */

static int calcScale( linConv_t * linConv ){

u_int16_t rawRange = linConv->rawHigh - linConv->rawLow;
double    engRange = linConv->engHigh - linConv->engLow;

    if( linConv->rawLow != linConv->rawHigh ){
        linConv->rawOffset = linConv->rawLow;
        linConv->engSlope = engRange/rawRange;
        linConv->engOffset = linConv->engLow;
        return( 0 );
    }
    else{
        return( -1 );   /* Error: div by zero */
    }    
}


/* 
 *  Convert the specified raw value to Engineering Units, 
 *  per the provided scaling paramters.
 * 
 * 
 */
static int rawToEng( linConv_t * linConv, u_int16_t rawVal, double * engVal ){

    *engVal = ( ( ( (double)rawVal - linConv->rawOffset ) * linConv->engSlope ) + linConv->engOffset );
    return( 0 );
    
}

