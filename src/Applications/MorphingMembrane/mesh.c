#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "brep.h"   
#include "mesh.h"

int NumNewElmts(GEdge *edge, double h)
{
  static double N, f, L;
  L  = elLength(edge->HeadElem);
  f  = modf(L/h, &N);
  if (L/N > h)
    {
      return N;
    }
  else
    {
      return N-1;  
    }
}
void MeshEdges(void)
{
  static int ge;      /* Global Edge Index 0<=ge<NGEdges           */
  static int n;       /* # of nodes inserted in Edge[ge] so far and
		         # of elements inserted +1                 */
  static int N;       /* total # of nodes to be inerted in Edge[ge]*/
  EdgeElem *head;     /* pointer to head element of Edge[ge]       */
  EdgeElem *element;  /* pointer to element most recently inserted 
			 in Edge[ge].  If none inserted yet, points
			 to head element                           */
  double h;           /* The actual hsize for elements in an edge.
			 Determined by calculating the number of
			 equilength elements (close to global hsize) 
			 that can fit in edge. Differs from global 
			 hsize.                                    */
  for(ge=0; ge<NGEdges; ge++)
    {
      n=0;
      head = GEdges[ge].HeadElem;

      if (elLength(head)<=hsize) break;

      N = NumNewElmts(GEdges+ge, hsize);
      h = elLength(head)/(N+1);
      element = head;

      /* Insert 1st Node */
      coord = realloc(coord, 2*(NNodes+1)*sizeof(double));
      NNodes++;
      coords_on_edge(&(coord[2*NNodes-2]), &(coord[2*NNodes-1]), head, h);

      /* Initialize Next Element */
      element->next = malloc(sizeof(EdgeElem));
      element->next->NodeA = NNodes;
      normal_to_edge(&(element->next->nxA), &(element->next->nyA), head, h);
      n++;
      element = element->next;

      /* Loop: insert rest of nodes and initialize rest of elements */
      while (n < N)
	{
	  /* Insert Next Node */
	  coord = realloc(coord, 2*(NNodes+1)*sizeof(double));
	  NNodes++; 
	  coords_on_edge(coord+2*NNodes-2, coord+2*NNodes-1, head, (n+1)*h);

	  /* Connect EndB of element to inserted node; 
         Calculate normal to edge at new node.     */
	  element->NodeB = NNodes;
      normal_to_edge(&(element->nxB), &(element->nyB), head, (n+1)*h);

	  /* Initialize Next Element */
	  element->next = malloc(sizeof(EdgeElem));
      element->next->NodeA = NNodes;
      element->next->nxA = element->nxB;
      element->next->nyA = element->nyB;
	  
	  n++;
	  element=element->next;
	}
      
      /* Connect tail of last element to tail of head */
      element->NodeB = head->NodeB;
      element->nxB   = head-> nxB ;
      element->nyB   = head-> nyB ;
      element->next  = NULL;
      /* Connect tail of head to nose of head->next */
      head->NodeB = head->next->NodeA;
      head-> nxB  = head->next-> nxA ;
      head-> nyB  = head->next-> nyA ;
    }

}

void TGridInsert(double h, double base, double x0, double y0)
{
  double  m;   /* Number of points per side of triangular grid          */
  int     i;   /* Index corresp. to row in triangular grid;  0<=i < m   */
  int     j;   /* Index corresp. to diag in triangular grid; 0<=j < m-i */
  int    id;   /* Index of coord array GridX[]; 1<=id<=(1/2)m*(m+1)     */

  modf(base/h, &m);
  GridN = m*(m+1)/2;
  GridX = malloc( 2*GridN*sizeof(double) );
  id = 0;
  for(i=0; i<m; i++)
  {
    for(j=0; j<m-i; j++)
    {
      GridX[2*id  ] = i*h/2 + j*h + x0;
      GridX[2*id+1] = i*h*sqrt(3.0)/2 + y0;
      id++;
    }
  }
}

void PGridInsert(double h, double base, double height, double x0, double y0)
{
  double  m;   /* Number of points along base of parallelagram grid     */
  double  n;   /* Number of points along sides of  parallelagram  grid  */
  int     i;   /* Index corresp. to row in triangular grid;  0<=i < m   */
  int     j;   /* Index corresp. to diag in triangular grid; 0<=j < n   */
  int    id;   /* Index of coord array GridX[]; 1<=id<=(n-1)*(m-1)      */

  modf(base/h, &m);
  modf(2*height/(sqrt(3.0)*h), &n);
  GridN = n*m;
  GridX = malloc( 2*GridN*sizeof(double) );
  id = 0;
  for(i=0; i<n; i++)
  {
    for(j=0; j<m; j++)
    {
      GridX[2*id  ] = i*h/2 + j*h + x0;
      GridX[2*id+1] = i*h*sqrt(3.0)/2 + y0;
      id++;
    }
  }
}

void PrintGrid(void)
{
  int i;
  FILE *fl;
  fl = fopen("points.tec","w");
  fprintf(fl,"ZONE T=AllPoints I=%d\n",GridN);
  for (i=0; i<GridN; i++)
    fprintf(fl,"%30.30f %30.30f\n",GridX[2*i],GridX[2*i+1]); 
  fclose(fl);
}

int NearBndry(double x, double y)
{
  int i;
  for (i=0; i<NNodes; i++)
  {
    if ( fabs(x-coord[2*i])<hsize/2 && fabs(y-coord[2*i+1])<hsize/2)
    {
      return 1;
    }
  }
  return 0;
}

void InsertNodes(void)
{
  int sb;    /* sub body counter              */
  int nd;    /* # nodes inserted so far       */
  int lc;    /* list counter                  */
  int gn;    /* grid node counter             */
  int lp;    /* loop counter                  */
  int ge;    /* global edge counter           */
  int *nflg; /* nflg[i]=1 if coord[i] inserted
		        0 otherwise           */
  int *nloc; /* nloc[i]=location of coord[i]
		in X array if already inserted*/

  EdgeElem *ee;   /* edge element pointer     */

  nflg = calloc(NNodes,sizeof(int));
  nloc = malloc(NNodes*sizeof(int));
  X    = malloc(2*(NNodes+GridN+3)*sizeof(double));
  LIST = malloc(2*(NNodes+GridN+3)*sizeof(int   ));
  nd=0;
  lc=0;

  for (sb=0; sb<NSubBodies; sb++)
    {
      N[sb]=0;
      /* get nodes on boundary of sb; insert in X if not already inserted in another sb*/
      for (lp=0; lp<SubBodies[sb].NL; lp++)
	{
	  for (ge=0; ge<SubBodies[sb].Lp[lp]->NE; ge++)
	    {
	      for (ee=SubBodies[sb].Lp[lp]->Edge[ge]->HeadElem; ee!=NULL; ee=ee->next)
		{
		  if (SubBodies[sb].Lp[lp]->EdgeDir[ge]==1)
		    {
		      if (nflg[ee->NodeA-1]==0)
			{
			  X[2*nd  ] = coord[2*(ee->NodeA-1)  ];
			  X[2*nd+1] = coord[2*(ee->NodeA-1)+1];
			  nflg[ee->NodeA-1]=1;
			  nloc[ee->NodeA-1]=nd;
			  LIST[lc] = nd+1;
			  nd++;
			}
		      else
			{
			  LIST[lc] = nloc[ee->NodeA-1]+1;
			}
		    }
		  else
		    {
		      if (nflg[ee->NodeB-1]==0)
			{
			  X[2*nd  ] = coord[2*(ee->NodeB-1)  ];
			  X[2*nd+1] = coord[2*(ee->NodeB-1)+1];
			  nflg[ee->NodeB-1]=1;
			  nloc[ee->NodeB-1]=nd;
			  LIST[lc] = nd+1;
			  nd++;
			}
		      else
			{
			  LIST[lc] = nloc[ee->NodeB-1]+1;
			}
		    }
		  N[sb]++;
		  lc++;
		}
	    }
	}
      /* get nodes inside boundary of sb; insert in X */
      for (gn=0; gn<GridN; gn++)
	{
	  if (InsideSb(sb, GridX[2*gn], GridX[2*gn+1]))
	    {
	      if (NearBndry(GridX[2*gn],GridX[2*gn+1])==0)
		{
		  X[2*nd  ] = GridX[2*gn  ];
		  X[2*nd+1] = GridX[2*gn+1];
		  LIST[lc] = nd+1;
		  N[sb]++;
		  nd++;
		  lc++;
		}
	    }
	}
    }
  NUMPTS = nd;
  X      = realloc(X,2*(NUMPTS+3)*sizeof(double));
  LIST   = realloc(LIST,lc*sizeof(int));
}

void rmExtElmts(void)
{
  int i;
  double x,y,x1,x2,x3,y1,y2,y3;
  for (i=0; i<NUMTRI; i++)
    {
      x1 = X[2*(V[3*i  ]-1)  ];
      y1 = X[2*(V[3*i  ]-1)+1];
      x2 = X[2*(V[3*i+1]-1)  ];
      y2 = X[2*(V[3*i+1]-1)+1];
      x3 = X[2*(V[3*i+2]-1)  ];
      y3 = X[2*(V[3*i+2]-1)+1];

      x = (x1+x2+x3)/3;
      y = (y1+y2+y3)/3;

      if (inside(x,y)!=1)
	{
	  rmTriang(i);
	}
    }
  return;
}
      
      

void rmTriang(int T)
{

  int i, j, k, adj;

  for (i=0; i<3; i++)
    {
      /* fix triangles adjacent to T */
      adj = E[3*T+i];
      if (adj != 0)
	{
	  for (j=0; j<3; j++)
	    {
	      if (E[3*(adj-1)+j]==T+1)
		{
		  E[3*(adj-1)+j]=0;
		  break;
		}
	    }
	}
      /* move triangle NUMTRI-1 to location T */
      V[3*T+i] = V[3*(NUMTRI-1)+i];
      E[3*T+i] = E[3*(NUMTRI-1)+i];

      /* inform moved triangles's neighbors of the move */
      adj = E[3*(NUMTRI-1)+i];
      if (adj != 0)
	{
	  for (j=0; j<3; j++)
	    {
	      if (E[3*(adj-1)+j]==NUMTRI)
		{
		  E[3*(adj-1)+j]=T+1;
		  break;
		}
	    }
	
	}
    }
  NUMTRI--;
  V = realloc(V, 3*NUMTRI*sizeof(int));
  E = realloc(E, 3*NUMTRI*sizeof(int));

} 

void connectSBmeshes(void)
{
  int i, j, ei, ai, bi, ej, aj, bj;

  for (i=0; i<NUMTRI; i++)
    {
      if (E[3*i]==0 || E[3*i+1]==0 || E[3*i+2]==0)
	{
	  for (j=i+1; j<NUMTRI; j++)
	    {
	      if (E[3*j]==0 || E[3*j+1]==0 || E[3*j+2]==0)
		{
		  if (E[3*i]==0)
		    {
		      ei=0;
		      ai=0;
		      bi=1;
		    }
		  if (E[3*i+1]==0)
		    {
		      ei=1;
		      ai=1;
		      bi=2;
		    }
      		  if (E[3*i+2]==0)
		    {
		      ei=2;
		      ai=2;
		      bi=0;
		    }
		  
		  if (E[3*j]==0)
		    {
		      ej=0;
		      aj=0;
		      bj=1;
		    }
		  if (E[3*j+1]==0)
		    {
		      ej=1;
		      aj=1;
		      bj=2;
		    }
		  if (E[3*j+2]==0)
		    {
		      ej=2;
		      aj=2;
		      bj=0;
		    }
		  
		  if ( V[3*i+ai]==V[3*j+bj] && V[3*i+bi]==V[3*j+aj] )
		    {
		      E[3*i+ei]=j+1;
		      E[3*j+ej]=i+1;
		      break;
		    }
		}
	    }
	}
    }
}

void PrintNodes(void)
{
  int i;
  FILE *fl;
  fl = fopen("nodes.plt","w");
  fprintf(fl,"ZONE T=Nodes I=%d\n",NUMPTS);
  for (i=0; i<NUMPTS; i++)
  {
    fprintf(fl,"%30.30f %30.30f\n",X[2*i],X[2*i+1]); 
  }
  fclose(fl);
}

void PrintElements(void)
{
  int i, v1, v2, v3;
  FILE *fl;
  fl = fopen("elmts.plt","w");
  for (i=0; i<NUMTRI; i++)
  {
    fprintf(fl,"ZONE T=Elmt%d I=4\n",i+1);
    v1 = V[3*i  ]-1;
    v2 = V[3*i+1]-1;
    v3 = V[3*i+2]-1;
    fprintf(fl,"%30.30f %30.30f\n",X[2*v1],X[2*v1+1]); 
    fprintf(fl,"%30.30f %30.30f\n",X[2*v2],X[2*v2+1]); 
    fprintf(fl,"%30.30f %30.30f\n",X[2*v3],X[2*v3+1]); 
    fprintf(fl,"%30.30f %30.30f\n",X[2*v1],X[2*v1+1]); 
  }
  fclose(fl);
}

void PrintElements2(void)
{
  int i;
  FILE *fl;
  static char t1[]="TITLE = \"Example: 2D DELAUNAY MESH\" ";

  
  if((fl=fopen("elmts.plt","w")) != NULL)
    {
      fprintf(fl,t1);
      fprintf(fl,"\nZONE T = \"TRIANGLES\", N = %d"
	      " , E=%d, F= FEPOINT , ET=TRIANGLE\n",NUMPTS,NUMTRI);
      /* save coordinates */
      for(i=0;i<NUMPTS;i++)
	fprintf(fl,"%10.8f %10.8f\n",X[2*i],X[2*i+1]);
      /* save connectivity */
      for(i=0;i<NUMTRI;i++)
	fprintf(fl,"%d %d %d\n",V[3*i],V[3*i+1],V[3*i+2]);
      fclose(fl);
    }
  else
    printf("Error opening file elmts.plt");
}
