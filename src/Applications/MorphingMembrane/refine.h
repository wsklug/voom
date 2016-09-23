#ifndef _REFINE_
#define _REFINE_

/* External Variables */
extern double *X;
extern int    NUMPTS,*V,*E,NUMTRI;
extern int    *RefineList;
/* Funtion Prototypes */
double hfun(double *x);
void   barycenter(int t, double *xcg);
int    longestedge(int t);
double edgelen(int t, int e);
void   LEPPdivide(int t, int *NumDiv);


#endif _REFINE_ 
