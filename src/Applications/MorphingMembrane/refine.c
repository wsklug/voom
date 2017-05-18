#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "refine.h"

double hfun(double *x)
{
  double h;
  double r;

  r = sqrt( pow(0.0-x[1],2) + pow(1.0-x[0],2) );
  h = 0.05*sqrt( sqrt( pow(r,3) ) );
  
  return h;
}

/*--------------------------------------------------------*/

void barycenter(int t, double *xcg)
{
  int i, j;

  for (i=0; i<2; i++)
    {
      xcg[i]=0;
      for (j=0; j<3; j++)
	{
	  xcg[i] = xcg[i] + X[2*(V[3*(t-1)+j]-1) + i ]/3.0;
	}
    }
  /*
  xcg[0] = ( X[2*V[3*(t-1)+0]+0] + X[2*V[3*(t-1)+1]+0] + X[2*V[3*(t-1)+2]+0] )/3.0;
  xcg[1] = ( X[2*V[3*(t-1)+0]+1] + X[2*V[3*(t-1)+1]+1] + X[2*V[3*(t-1)+2]+1] )/3.0;
  */
}

/*--------------------------------------------------------*/

int longestedge(int t)
{
  int LE;
  int e;
  double len;
  double maxlen;
  
  if (t==0)
    {
      return 0;
    }
  else
    {
      maxlen=0.0;
      
      for (e=0; e<3; e++)
	{
	  len = edgelen(t,e);
	  if (maxlen < len)
	    {
	      maxlen = len;
	      LE = e;
	    }
	}
      return LE;
    }	    

}

/*--------------------------------------------------------*/

double edgelen(int t, int e)
{
  int v1, v2;
  double len;

  switch(e)
    {
    case 0:
      v1 = V[3*(t-1)+0];
      v2 = V[3*(t-1)+1]; 
      break;
      
    case 1:
      v1 = V[3*(t-1)+1];
      v2 = V[3*(t-1)+2]; 
      break;
      
    case 2:
      v1 = V[3*(t-1)+2];
      v2 = V[3*(t-1)+0]; 
      break;
    }

  len = sqrt(  pow( X[2*(v1-1)+0] - X[2*(v2-1)+0], 2 ) 
	     + pow( X[2*(v1-1)+1] - X[2*(v2-1)+1], 2 ) );

  return len;

}

/*--------------------------------------------------------*/

void LEPPdivide(int t, int *NumDiv)
{
  int t1;
  int t2;
  int t3;
  int t4;
  int e1[3];
  int e2[3];
  int P;
  int i;

  t1 = t;

  e1[0] = longestedge(t1);
  
  t2 = E[ 3*(t1-1)+e1[0] ];

  e2[0] = longestedge(t2);

/*  while (t2!=0 && E[ 3*(t2-1)+e2[0] ]!=t1) */
  while (t2!=0 && edgelen(t1, e1[0])<edgelen(t2,e2[0]))
    {
      t1    = t2;
      e1[0] = e2[0];
      t2    = E[ 3*(t1-1)+e1[0] ];
      e2[0] = longestedge(t2);
    }

  if (t2!=0 && edgelen(t1, e1[0])==edgelen(t2,e2[0]))
    {
      for (i=0; i<3; i++)
	{
	  if (E[3*(t2-1)+i]==t1)
	    {
	      e2[0] = i;
	      break;
	    }
	}
    }

  /* ************************************************************************************** */
  /* At this point, the following should be true:                                           */
  /*                                                                                        */
  /*   (1)  t2 = E[ 3*(t1-1)+e1[0] ]; i.e., t2 is the id of the triangle adjacent to        */
  /*        the longest edge (e1[0] = 0, 1, or 2, whichever is longest) of the triangle     */
  /*   	    with id t1;                                                                     */
  /*	                                                                                    */
  /*  and                                                                                   */
  /*                                                                                        */
  /*   (2)  either:                                                                         */
  /*                                                                                        */
  /*       (a) t1 = E[ 3*(t2-1)+e2[0] ]; i.e, t1 is the triangle adjacent to the            */
  /*	       longest edge of t2 (t1 and t2 share a common longest edge)                   */
  /*                                                                                        */
  /*	    or                                                                              */
  /*                                                                                        */
  /*	   (b) t2=0, i.e., the longest edge of t1 is along a boundary                       */
  /* ************************************************************************************** */
  
  if (t2==0)
    {
      printf("\n\n  Dividing boundary triangle %d.\n", t1);
      RefineList[t1-1]=-1;
    }
  else
    {
	  printf("\n\n  Dividing adjacent triangles %d and %d.\n", t1, t2);
	  RefineList[t1-1]=-1;
	  RefineList[t2-1]=-1;
    }
  
  for (i=1; i<3; i++)
    {
      if (e1[0]+i>2)
	{
	  e1[i]=e1[0]+i-3;
	}
      else
	{
	  e1[i]=e1[0]+i;
	}
      if (e2[0]+i>2)
	{
	  e2[i]=e2[0]+i-3;
	}
      else
	{
	  e2[i]=e2[0]+i;
	}
    }

  /* insert P at middle of e1[0] */

  X = realloc(X, 2*(NUMPTS+1)*sizeof(double));
  NUMPTS++;
  P = NUMPTS;

  X[2*(P-1)+0] = (0.5)*( X[ 2*(V[ 3*(t1-1)+e1[0] ] - 1)     ]
			+X[ 2*(V[ 3*(t1-1)+e1[1] ] - 1)     ] );
  X[2*(P-1)+1] = (0.5)*( X[ 2*(V[ 3*(t1-1)+e1[0] ] - 1) + 1 ]
			+X[ 2*(V[ 3*(t1-1)+e1[1] ] - 1) + 1 ] );

  /* Add triangle t3 */

  V = realloc(V, 3*(NUMTRI+1)*sizeof(int));
  E = realloc(E, 3*(NUMTRI+1)*sizeof(int));
  RefineList = realloc(RefineList, (NUMTRI+1)*sizeof(int));

  
  NUMTRI++;
  t3 = NUMTRI;

  RefineList[t3-1]=-1;

  printf("\n\n  Creating triangle %d.\n", t3);

  /* Connect vertices of t1 and t3 */

  V[3*(t3-1)+0] = P;
  V[3*(t3-1)+1] = V[ 3*(t1-1)+e1[1] ];
  V[3*(t3-1)+2] = V[ 3*(t1-1)+e1[2] ];

  V[3*(t1-1)+e1[1]] = P;

  if (t2 == 0)
    {
      t4=0;
    }
  else
    {
      /* Add triangle t4 */

      V = realloc(V, 3*(NUMTRI+1)*sizeof(int));
      E = realloc(E, 3*(NUMTRI+1)*sizeof(int));
      RefineList = realloc(RefineList, (NUMTRI+1)*sizeof(int));
 
      NUMTRI++;
      t4 = NUMTRI;
      
      RefineList[t4-1]=-1;

      printf("\n\n  Creating triangle %d.\n", t4);

      /* Connect vertices of t2 and t4 */
      
      V[3*(t4-1)+0] = P;
      V[3*(t4-1)+1] = V[ 3*(t2-1)+e2[1] ];
      V[3*(t4-1)+2] = V[ 3*(t2-1)+e2[2] ];
      
      V[3*(t2-1)+e2[1]] = P;
    }

  /* Fix adjacency for t1 and t3 */

  E[3*(t3-1)+0] = t2;
  E[3*(t3-1)+1] = E[3*(t1-1)+e1[1]];
  E[3*(t3-1)+2] = t1;

  E[3*(t1-1)+e1[0]] = t4;
  E[3*(t1-1)+e1[1]] = t3;

  if (t2!=0)
    {
      E[3*(t4-1)+0] = t1;
      E[3*(t4-1)+1] = E[3*(t2-1)+e2[1]];
      E[3*(t4-1)+2] = t2;
      
      E[3*(t2-1)+e2[0]] = t3;
      E[3*(t2-1)+e2[1]] = t4;
    }

  /* Fix adjacency for element that was adjacent to t1 and is now adjacent to t3
               and for element that was adjacent to t2 and is now adjacent to t4  */

  for (i=0; i<3; i++)
    {
      if ( E[3*( E[3*(t3-1)+1] - 1 ) + i] == t1 )
	{
	  E[3*( E[3*(t3-1)+1] - 1 ) + i] = t3;
	}

      if (t2!=0)
	{
	  if ( E[3*( E[3*(t4-1)+1] - 1 ) + i] == t2 )
	    {
	      E[3*( E[3*(t4-1)+1] - 1 ) + i] = t4;
	    }
	}
    }

  if (t2==0)
    {
      *NumDiv++;
    }
  else
    {
      *NumDiv = *NumDiv + 2;
    }

}

/*--------------------------------------------------------*/
