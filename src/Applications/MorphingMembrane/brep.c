#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "brep.h"   


int i, j, k, l, i_tmp, ID;
static char line[80], tmp[20];

void ReadBrep(char * fname)
{
	FILE *inp;
	inp = fopen(fname,"r"); 
    if (inp == 0)
    {
      printf("\nCan't open %s!  Terminating.\n\n",fname);
      exit(EXIT_FAILURE);
    }

  /* Read in Nodal Coords */

	fscanf(inp, "%d", &NNodes);              
	
    coord = malloc( 2*NNodes*sizeof(double) );
	for (i=0; i<NNodes; i++)
		fscanf(inp, "%d%lg%lg", &ID, &(coord[2*i]), &(coord[2*i+1]) );


  /* Read in Edge data */

	fscanf(inp, "%d", &NGEdges);

	GEdges = malloc( NGEdges*sizeof(GEdge) );
	for (i=0; i<NGEdges; i++)
	{
	  GEdges[i].HeadElem=malloc(sizeof(EdgeElem));

	  fscanf(inp, "%d%d%d%lg%lg%lg%lg", &ID, 
		 &(GEdges[i].HeadElem->NodeA), &(GEdges[i].HeadElem->NodeB),
		 &(GEdges[i].HeadElem-> nxA ), &(GEdges[i].HeadElem-> nyA ), 
		 &(GEdges[i].HeadElem-> nxB ), &(GEdges[i].HeadElem-> nyB )  );

	  GEdges[i].HeadElem->next=NULL;
	}

  /* Read in Loop data */


	fscanf(inp, "%d", &NLoops);

	Loops = malloc( NLoops*sizeof(Loop) );

	for (i = 0; i < NLoops; i++) 
	{
	  fscanf(inp, "%d%d\n", &ID, &(Loops[i].NE));

	  Loops[i].Edge=malloc( (Loops[i].NE)*sizeof(GEdge) );
	  Loops[i].EdgeDir = malloc( (Loops[i].NE)*sizeof(int) );
	  
	  fgets(line, 80, inp);

	  j=0;
	  for (l=0; l < Loops[i].NE; l++)
	  {
		  k=0;
		  while ((line[j]!=' ') && (line[j]!='\n') && (line[j]!='\0')) 
		  {
			  tmp[k]=line[j]; 
			  j++; k++;
		  }
		  tmp[k]='\0'; j++;
		  i_tmp=atoi(tmp);
		  Loops[i].Edge[l]=&GEdges[abs(i_tmp) - 1];
		  Loops[i].EdgeDir[l]=(i_tmp/(abs(i_tmp)));
	  }
    }

  /* Read in SubBody Data */

	fscanf(inp, "%d", &NSubBodies);

	SubBodies = malloc( NSubBodies*sizeof(SubBody) );

	for (i = 0; i < NSubBodies; i++) 
	{
	  fscanf(inp, "%d%d\n", &ID, &(SubBodies[i].NL));
	  
	  SubBodies[i].Lp=malloc( (SubBodies[i].NL)*sizeof(Loop) );
	  
	  fgets(line, 80, inp);

	  j=0;
	  for (l=0; l < SubBodies[i].NL; l++)
	  {
		  k=0;
		  while ((line[j]!=' ') && (line[j]!='\n') && (line[j]!='\0')) 
		  {
			  tmp[k]=line[j]; 
			  j++; k++;
		  }
		  tmp[k]='\0'; j++;
		  i_tmp=atoi(tmp);
		  SubBodies[i].Lp[l]=&Loops[abs(i_tmp) - 1];
	  }
    }

  /* Read in Body Data */

	fscanf(inp, "%d", &NBodies);

	Bodies = malloc( NBodies*sizeof(Body) );

	for (i = 0; i < NBodies; i++) 
	{
	  fscanf(inp, "%d%d\n", &ID, &(Bodies[i].NSB));

	  Bodies[i].SBody=malloc( (Bodies[i].NSB)*sizeof(SubBody) );
	  
	  fgets(line, 80, inp);

	  j=0;
	  for (l=0; l < Bodies[i].NSB; l++)
	  {
		  k=0;
		  while ((line[j]!=' ') && (line[j]!='\n') && (line[j]!='\0')) 
		  {
			  tmp[k]=line[j]; 
			  j++; k++;
		  }
		  tmp[k]='\0'; j++;
		  i_tmp=atoi(tmp);
		  Bodies[i].SBody[l]=&SubBodies[abs(i_tmp) - 1];
	  }
    }

   fclose(inp);

}


void coords_on_edge(double *x, double *y, EdgeElem *Edge, double s)
{
	double xa = coord[2*(Edge->NodeA)-2];
	double ya = coord[2*(Edge->NodeA)-1];
	double xb = coord[2*(Edge->NodeB)-2];
	double yb = coord[2*(Edge->NodeB)-1];

	double dx = xb - xa;
	double dy = yb - ya;

/* Rotate Normals from in x,y coords. to xi,eta coords. 
   Note that Normal Vectors aren't normalized, but that's
   ok since we only need ration of n_xi/n_eta for tangents */

	double n_xi_A  =  dx*(Edge->nxA) + dy*(Edge->nyA);
	double n_eta_A = -dy*(Edge->nxA) + dx*(Edge->nyA);
	double n_xi_B  =  dx*(Edge->nxB) + dy*(Edge->nyB);
	double n_eta_B = -dy*(Edge->nxB) + dx*(Edge->nyB);
	
	double tanA = -(n_xi_A)/(n_eta_A);
	double tanB = -(n_xi_B)/(n_eta_B);

	double xi  = s/elLength(Edge);
	double eta = tanA*( pow(xi,3) - 2*pow(xi,2) + xi )
		   + tanB*( pow(xi,3) -   pow(xi,2)	 );
	*x = dx*xi - dy*eta + xa ;
	*y = dy*xi + dx*eta + ya ;

}

void normal_to_edge(double *nx, double *ny, EdgeElem *Edge, double s)
{
	double xa = coord[2*(Edge->NodeA)-2];
	double ya = coord[2*(Edge->NodeA)-1];
	double xb = coord[2*(Edge->NodeB)-2];
	double yb = coord[2*(Edge->NodeB)-1];

	double dx = xb - xa;
	double dy = yb - ya;

/* Rotate Normals from in x,y coords. to xi,eta coords. 
   Note that Normal Vectors aren't normalized, but that's
   ok since we only need ration of n_xi/n_eta for tangents */

	double n_xi_A  =  dx*(Edge->nxA) + dy*(Edge->nyA);
	double n_eta_A = -dy*(Edge->nxA) + dx*(Edge->nyA);
	double n_xi_B  =  dx*(Edge->nxB) + dy*(Edge->nyB);
	double n_eta_B = -dy*(Edge->nxB) + dx*(Edge->nyB);
	
	double tanA = -(n_xi_A)/(n_eta_A);
	double tanB = -(n_xi_B)/(n_eta_B);

	double xi  = s/elLength(Edge);
	double theta = tanA*( 3*pow(xi,2) - 4*xi + 1 )
		     + tanB*( 3*pow(xi,2) - 2*xi     );

	/*  n_xi  = -1*sin(theta) */
	/*  n_eta =    cos(theta) */
	*nx = -dx*sin(theta) - dy*cos(theta);
	*ny = -dy*sin(theta) + dx*cos(theta);
  
}

double intersect(double yq, EdgeElem *Edge)
{
	double ya = coord[2*(Edge->NodeA)-1];
	double yb = coord[2*(Edge->NodeB)-1];
	double s  = elLength(Edge)*(yq-ya)/(yb-ya);
	double x;
    double y;
	coords_on_edge(&x, &y, Edge, s);
	return x;
}

int inside(double xq, double yq)
{
  int i, j, k;
  int counter = 0;
  GEdge *GE;
  EdgeElem *EE;
  double xa, ya, xb, yb, xx;
  for(i=0; i<NLoops; i++)
    {
       for(j=0; j<Loops[i].NE; j++)
         {
          GE = Loops[i].Edge[j];
          k=0;
          for (EE=GE->HeadElem;EE!=NULL;EE=EE->next)
            {
         	  xa = coord[2*(EE->NodeA)-2];
         	  ya = coord[2*(EE->NodeA)-1];
         	  xb = coord[2*(EE->NodeB)-2];
         	  yb = coord[2*(EE->NodeB)-1];
            
         	  if( (ya<=yq && yq<yb) || (yb<=yq && yq<ya) )
                {
         		  xx=intersect(yq, EE);
                  if( xq < xx) 
                    {
                      counter++;
                    }
              k++;
             }
           }
         }
    }
  if( (counter % 2)==1 ) 
  {
    return 1;
  }
  else 
  {
    return 0;
  }
}

int InsideSb(int sb, double xq, double yq)
{
  int i, j, k;
  int counter = 0;
  GEdge *GE;
  EdgeElem *EE;
  double xa, ya, xb, yb, xx;
  for(i=0; i<SubBodies[sb].NL; i++)
    {
      for(j=0; j<SubBodies[sb].Lp[i]->NE; j++)
	{
          GE = SubBodies[sb].Lp[i]->Edge[j];
          k=0;
          for (EE=GE->HeadElem;EE!=NULL;EE=EE->next)
            {
	      xa = coord[2*(EE->NodeA)-2];
	      ya = coord[2*(EE->NodeA)-1];
	      xb = coord[2*(EE->NodeB)-2];
	      yb = coord[2*(EE->NodeB)-1];
	      
	      if( (ya<=yq && yq<yb) || (yb<=yq && yq<ya) )
                {
		  xx=intersect(yq, EE);
                  if( xq < xx) 
                    {
                      counter++;
                    }
		  k++;
		}
	    }
	}
    }
  if( (counter % 2)==1 ) 
    {
      return 1;
    }
  else 
    {
      return 0;
    }
}


void box(double NumPoints, double xmin, double xmax, double ymin, double ymax)
{
    FILE *fl_in, *fl_ex;
    int i, j=0, k=0;
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    double x, y;

    srand(time(NULL));

    fl_in = fopen("in_pts.tec","w"); 
    fl_ex = fopen("ex_pts.tec","w"); 

    fprintf(fl_in,"                                  \n");
    fprintf(fl_ex,"                                  \n");

    for (i=0; i < NumPoints; i++)
    {
        x = rand()*dx/RAND_MAX + xmin;
        y = rand()*dy/RAND_MAX + ymin;

        if( inside(x,y)==1 )
        {
            fprintf(fl_in,"%30.30f %30.30f\n", x, y);
            j++;
        }
        else
        {
            fprintf(fl_ex,"%30.30f %30.30f\n", x, y);
            k++;
        }
    }

    printf("\nTotal No. of Points Checked: %d\n", i);

    fseek(fl_in, 0, SEEK_SET);
    fseek(fl_ex, 0, SEEK_SET);
    fprintf(fl_in,"ZONE T=INTERIOR I=%d",j);
    fprintf(fl_ex,"ZONE T=EXTERIOR I=%d",k);
    fseek(fl_in, 0, SEEK_END);
    fseek(fl_ex, 0, SEEK_END);
    fclose(fl_in);
    fclose(fl_ex);
}


void PrintEdges()
{
  int l=0;
  int ge;
  double x,y;
  GEdge *GE;
  EdgeElem *EE;
  FILE *fl;
  fl = fopen("edges.plt","w");

  for (ge=0;ge<NGEdges;ge++)
    {
      GE = &GEdges[ge];
      for (EE=GE->HeadElem;EE!=NULL;EE=EE->next)
	{
	  l++;
 	  fprintf(fl,"ZONE T=Elem%d I=2\n",l); 
          x = coord[2*(EE->NodeA)-2];
          y = coord[2*(EE->NodeA)-1];
 	  fprintf(fl,"%30.30f %30.30f\n",x,y); 
          x = coord[2*(EE->NodeB)-2];
          y = coord[2*(EE->NodeB)-1];
	  fprintf(fl,"%30.30f %30.30f\n",x,y); 
	   
	}
    }
  fclose(fl);
}

double elLength(EdgeElem *elem)
{
  double xa = coord[2*(elem->NodeA)-2];
  double ya = coord[2*(elem->NodeA)-1];
  double xb = coord[2*(elem->NodeB)-2];
  double yb = coord[2*(elem->NodeB)-1];
  double dx = xb - xa;
  double dy = yb - ya;
  return sqrt( pow(dx,2) + pow(dy,2) );
}
