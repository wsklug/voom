#ifndef _BREP_
#define _BREP_

/* Data Structure Declarations */
/* --------------------------- */

typedef struct EdgeElem	EdgeElem;
typedef struct GEdge	GEdge;
typedef struct Loop	Loop;
typedef struct SubBody  SubBody;
typedef struct Body	Body;


struct EdgeElem{
	int NodeA, NodeB;
	double nxA, nyA, nxB, nyB;
	struct EdgeElem *next;
};

struct GEdge{
	struct EdgeElem *HeadElem;
};

struct Loop{
	int NE;
	int *EdgeDir;
	struct GEdge **Edge;
};

struct SubBody{
	int NL;
	struct Loop **Lp;
};

struct Body{
	int NSB;
	struct SubBody **SBody;
};



/* Function Prototypes */
/* ------------------- */
void   ReadBrep(char * fname);
void   coords_on_edge(double *x, double *y, EdgeElem *Edge, double s);
void   normal_to_edge(double *nx, double *ny, EdgeElem *Edge, double s);
double intersect(double yq, EdgeElem *Edge);
int    inside(double xq, double yq);
int    InsideSb(int sb, double xq, double yq);
void   fill_bb(double xmin, double xmax, double ymin, double ymax);
void   separate_pts(double Points[]);
void   PrintEdges();
void   box(double NumPoints, double xmin, double xmax, double ymin, double ymax);
double elLength(EdgeElem *elem);


/* External Variables */
/* ------------------ */

extern double	*coord;
extern GEdge	*GEdges;
extern Loop	*Loops;
extern SubBody *SubBodies;
extern Body	*Bodies;
extern int NNodes, NGEdges, NLoops, NSubBodies, NBodies;


#endif  /* _BREP_ */
