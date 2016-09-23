unsigned NumberOfNodes=5;
double CoordinateArray[][3] = { {0.0,0.0,0.0},
			       {0.0,1.0,0.0},
			       {0.5,0.5,0.0},
			       {1.0,0.0,0.0},
			       {1.0,1.0,0.0}};

double NormalArray[][3]={ {0.0,0.0,1.0}, 
			 {0.0,0.0,1.0}, 
			 {0.0,0.0,1.0}, 
			 {0.0,0.0,1.0}, 
			 {0.0,0.0,1.0}}; 

unsigned NumberOfElements=4;
int Connectivity[][3]={ {0, 2, 1},
		       {0, 3, 2},
		       {2, 3, 4},
		       {1, 2, 4} };
