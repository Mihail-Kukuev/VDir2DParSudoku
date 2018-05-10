#include "task_def.h"


//--------------------//
//    Test solution   //
//--------------------//

double u ( double x1, double x2, double t )
{
    return exp ( x1 + x2 + t );
}




//-------------------//
//     Conditions    //
//-------------------//

double u0 ( double x1, double x2 ) // t = 0
{
    return u ( x1, x2, 0 );
}


double m1 ( double x2, double t )  // x1 = 0
{
    return u ( 0, x2, t );
}


double m2 ( double x2, double t )  // x1 = L1
{
    return u ( ( double ) L1, x2, t );
}


double m3 ( double x1, double t )  // x2 = 0
{
    return u ( x1, 0, t );
}


double m4 ( double x1, double t )  // x2 = L2
{
    return u ( x1, ( double ) L2, t );
}



//---------------------//
//     Coefficients    //
//---------------------//

double k1 ( double x1, double x2, double t )
{
    return exp ( -( x1 + x2 + t ) );
}


double k2 ( double x1, double x2, double t )
{
    return exp ( -( x1 + x2 + t ) );
}


double f ( double x1, double x2, double t )
{
    return exp ( x1 + x2 + t );
}
