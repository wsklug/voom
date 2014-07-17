#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "brep.h"
#include "mesh.h"
#include "refine.h"

double hsize;
double *GridX; /*  GridX[] and GridY[] are the global x- and y-coord.    */
int     GridN; /*  GridN = number of pts inserted inside the boundary.   */
double *X,*Y;
int     NUMPTS,*N,*LIST,*BIN,*V,*E,NUMTRI,*NT;
int    *RefineList;
double    *aspect;
int NNodes, NGEdges, NLoops, NSubBodies, NBodies;

double	*coord;
GEdge	*GEdges;
Loop	*Loops;
SubBody *SubBodies;
Body	*Bodies;

void mesh(void);
int  refinemesh(void);
int  GetRefineList(void);
double AspectRatio(int t);
void PrintAspectRatios(char opt[]);
void PrintVTK(char * fname);
void PrintObj(char * fname);
void FreeMesh(void);

int main(int argc, char* argv[])
{
 
  int i, NumDivisions;
  char brepName[50], vtkName[50], objName[50];
  hsize=atof(argv[2]);/*0.1;*/

  printf("Input hsize = %g\n",hsize);

  sprintf(brepName,"%s.brep",argv[1]);
  sprintf(vtkName, "%s.vtk", argv[1]);
  sprintf(objName, "%s.obj", argv[1]);

  ReadBrep(brepName);
  MeshEdges();
  PrintEdges();
  mesh();
  printf("Meshed brep with %d nodes and %d elements.\n",NUMPTS,NUMTRI);
/*   PrintAspectRatios("w"); */
  
/*   RefineList = malloc(NUMTRI*sizeof(int)); */
/*   while(GetRefineList()>0) */
/*     { */
/*       NumDivisions=refinemesh(); */
/*     } */
/*   free(RefineList); */
  
/*   PrintElements2(); */

/*   PrintAspectRatios("a"); */
  PrintVTK(vtkName);
  PrintObj(objName);
  FreeMesh();
  return 0;
}

/*--------------------------------------------------------*/

void mesh(void)
{
  int sb, i, j;
  int Tstart, Lstart;

  double xmin,xmax,ymin,ymax;
  double height,base,x0,y0;
  xmin=xmax=coord[0];
  ymin=ymax=coord[1];
  for (i=0; i<NNodes; i++) {
    if( coord[2*i  ] < xmin ) xmin = coord[2*i  ];
    if( coord[2*i  ] > xmax ) xmax = coord[2*i  ];
    if( coord[2*i+1] < ymin ) ymin = coord[2*i+1];
    if( coord[2*i+1] > ymax ) ymax = coord[2*i+1];
  }
  height = ymax-ymin;
  y0 = ymin;
  x0 = xmin - height/sqrt(3.0);
  base = xmax-xmin + 2.0*height/sqrt(3.0);
  printf("Creating grid.\n");
  printf("x min/max (%g,%g)\n",xmin,xmax);
  printf("y min/max (%g,%g)\n",ymin,ymax);
  printf("hsize = %g\n",hsize);
  printf("height = %g\n",height);
  printf("base = %g\n",base);
  printf("x0 = %g\n",x0);
  printf("y0 = %g\n",y0);
/*   PGridInsert(hsize, 5.0, 2.0,-1.1547 ,-1.0); */
  PGridInsert(hsize, base, height, x0 , y0);
  N = malloc(NSubBodies*sizeof(int));
  InsertNodes();
  PrintNodes();
  BIN = malloc(NUMPTS*sizeof(int));
  V   = malloc(3*(2*NUMPTS+1)*sizeof(int));
  E   = malloc(3*(2*NUMPTS+1)*sizeof(int));

  Lstart=0;
  Tstart=0;
  NUMTRI=0;
  NT = malloc(NSubBodies*sizeof(int));
  for (sb=0; sb<NSubBodies; sb++)
    {
      deltri_(&NUMPTS,N+sb,X,LIST+Lstart,BIN,V+Tstart,E+Tstart,NT+sb);
      for (i=0; i<NT[sb]; i++)
	{
	  for (j=0; j<3; j++)
	    {
	      if (E[3*(NUMTRI+i)+j]>0)
		{
		  E[3*(NUMTRI+i)+j]=E[3*(NUMTRI+i)+j]+NUMTRI;
		}
	    }
	}
      Lstart=Lstart+N[sb];
      Tstart=Tstart+3*NT[sb];
      NUMTRI=NUMTRI+NT[sb];
    }
  
  rmExtElmts();
  connectSBmeshes();
  
  free(coord);
  free(GEdges);
  free(Loops);
  free(SubBodies);
  free(Bodies);
  free(GridX);
  free(BIN);
  free(LIST);
}  

/*--------------------------------------------------------*/

int refinemesh(void)
{
  int t;
  int NumDiv;
  double xcg[2];
  double len;
  int m = NUMTRI;
  NumDiv=0;
  for (t=1; t<=m; t++)
    {
      if (RefineList[t-1]==1)
	{
	  LEPPdivide(t,&NumDiv);
	  printf("\n\n Divided %d triangles.\n\n", NumDiv);
	}
    }
  return NumDiv;
}

int GetRefineList(void)
{
  int t;
  int N2refine;
  double xcg[2];
  double len;
  
  N2refine=0;

  RefineList = realloc(RefineList, NUMTRI*sizeof(int));

  for (t=1; t<=NUMTRI; t++)
    {
      barycenter(t, xcg);
      len=edgelen(t,longestedge(t));
      if (len > hfun(xcg) && len > 0.001)
	{
	  RefineList[t-1]=1;
	  N2refine++;
	}
      else
	{
      	  RefineList[t-1]=0;
	}
    }
  return N2refine;
}

/* ************************************************ */
  
double AspectRatio(int t)
{
  double a, b, c;
  double s;
  double h;
  double rho;
  
  a = edgelen(t, 0);
  b = edgelen(t, 1);
  c = edgelen(t, 2);
  
  s = 0.5*(a + b + c);
  
  rho = sqrt(s*(s-a)*(s-b)*(s-c))/s;
  
  h = a*b*c/(4*sqrt(s*(s-a)*(s-b)*(s-c)));
  
  return h/rho;
  
}

void PrintAspectRatios(char opt[])
{
  int i, j;
  double max, min;
  int *bin;
  int Nbins;
  FILE *fl;
 
  
  if ((fl=fopen("hist.dat",opt))==NULL)
    {
      fprintf(fl,"\nCan't open hist.dat!  Terminating.\n\n");
      exit(0);
    }
  
  Nbins = 10;
  bin = calloc(Nbins,sizeof(int));

  aspect=malloc(NUMTRI*sizeof(double));
  for (i=0; i<NUMTRI; i++)
    {
      aspect[i]=AspectRatio(i+1);
    }

  /* get max of data */
  
  max=aspect[0];
  min=aspect[0];
  for (i=0; i<NUMTRI; i++)
    {
      if (aspect[i]>10.0)
	{
	  printf("\nWarning! Element %d has a high aspect ratio of %8.8g\n",
		 i,aspect[i]);
	  break;
        }
      if (aspect[i]>max)
        {
          max=aspect[i];
        }
      if  (aspect[i]<min)
        {
          min=aspect[i];
        }
    }
  fprintf(fl,"max = %8.8g\n",max);
  fprintf(fl,"min = %8.8g\n",min);
  for (i=0; i<NUMTRI; i++)
    {
      for (j=1; j<=Nbins; j++)
        {
	  if (aspect[i]<=min+j*(max-min)/Nbins)
	    {
              bin[j-1]++;
	      break;
	    }
        }
    }

  for (i=0; i<Nbins; i++)
    {
      fprintf(fl,"%d %8.4g %d \n", i+1,(0.5)*(min+(i+1)*(max-min)/Nbins),bin[i]);
    }
  free(bin);
  free(aspect);
  fclose(fl);
}


void FreeMesh(void)
{
  free(V);
  free(E);
  free(X);
}


void PrintVTK(char * fname)
{
  int i;
  FILE *fl;
  static char t1[]="TITLE = \"Example: 2D DELAUNAY MESH\" ";

  
  if((fl=fopen(fname,"w")) != NULL)
    {
      fprintf(fl,"# vtk DataFile Version 2.0\n");
      fprintf(fl,"Test example\n");
      fprintf(fl,"ASCII\n");
      fprintf(fl,"DATASET POLYDATA\n");
      fprintf(fl,"POINTS  %d  double\n",NUMPTS);
    
      /* save coordinates */
      for(i=0;i<NUMPTS;i++)
	fprintf(fl,"%10.8f %10.8f %10.8f\n",X[2*i],X[2*i+1],0.0);

      fprintf(fl,"POLYGONS %d  %d \n",NUMTRI,4*NUMTRI);
      /* save connectivity */
      for(i=0;i<NUMTRI;i++)
	fprintf(fl,"3 %d %d %d\n",V[3*i]-1,V[3*i+1]-1,V[3*i+2]-1);
      fclose(fl);
    }
  else
    printf("Error opening file %s",fname);
}

void PrintObj(char * fname)
{
  int i;
  FILE *fl;
  static char t1[]="TITLE = \"Example: 2D DELAUNAY MESH\" ";

  
  if((fl=fopen(fname,"w")) != NULL)
    {    
      /* save coordinates */
      for(i=0;i<NUMPTS;i++)
	fprintf(fl,"v %10.8f %10.8f %10.8f\n",X[2*i],X[2*i+1],0.0);

      /* save connectivity */
      for(i=0;i<NUMTRI;i++)
	fprintf(fl,"f %d %d %d\n",V[3*i],V[3*i+1],V[3*i+2]);
      fclose(fl);
    }
  else
    printf("Error opening file %s",fname);
}
