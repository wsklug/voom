#ifndef _MESH_
#define _MESH_

/* External Variables */
extern int       NNodes,NGEdges,NLoops,NSubBodies,NBodies;
extern int       NUMPTS,*N,*LIST,*BIN,*V,*E,NUMTRI,*NT;
extern double    hsize;
extern double   *X;
extern double	*coord;
extern GEdge	*GEdges;
extern Loop	*Loops;
extern SubBody	*SubBodies;
extern Body	*Bodies;
extern double   *GridX; /*  GridX[] is the global coord. array for */
		        /*  nodes inserted inside the boundary.    */
extern int       GridN; /*  GridN = number of pts inserted inside  */
                        /*  the boundary.                          */


/* Funtion Prototypes */

void   MeshEdges(void);
int     NumNewElmts(GEdge *edge, double h);
void   TGridInsert(double h, double base, double x0, double y0);
void   PGridInsert(double h, double base, double height, double x0, double y0);
void   PrintGrid(void);
int    NearBndry(double x, double y);
void   InsertNodes(void);
void   PrintNodes(void);
void   PrintElements(void);
void   PrintElements2(void);
void   rmTriang(int T);
void   rmExtElmts(void);
void   connectSBmeshes(void);
extern void deltri_(int *, int *, double *,
		    int *, int *, int *, int *, int *);

#endif _MESH_ 
