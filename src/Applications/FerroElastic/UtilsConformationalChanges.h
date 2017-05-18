/* FUNCTIONS IMPLEMENTED IN UTILS FILE
void printAllParaview
void printStretchDirectioParaview
ReadInputFile
*/

struct InputData {
  // File names
  string Model;
  string Outline;
  // General analysis settings
  int min;            // Minimize with respect to eta/theta flag (0-> don't do it, 1-> min wrt eta, 2-> min wrt eta and theta)
  int Ns;             // Number of stretch/theta dof
  double in_eta;      // ta initial guess
  double in_theta;    // theta initial guess

  int NPout;          // Nodes needed to build the outline
  double HexonShift;  // Scale factor for plotting shear and stretch directions
  // Bending
  double KC;
  double KG;
  double C0;
  // Stretching
  double mu;
  double kS;
  double AngleShift; // Uniform angle shift, MUST BE IN RADIANT !!!
  int WcType;        // Type of conformational energy (if 0 ->  no conformational energy added) 
  double WcConst0;   // K eta
  double WcConst1;   // K theta
  double WcConst2;   // Assigned \bar{eta}
};



InputData ReadInputFile(const string InputFileName)
{
  string line;
  ifstream InputFile(InputFileName.c_str());

  if (!InputFile) {
    std::cout << "Cannot open input file: "
	      << InputFileName << std::endl;
    exit(0);
  }

  InputData InputValues;
  
  while ( !InputFile.eof() )
  {
    getline(InputFile,line);
    switch (line) {
    case "Model":
      getline(InputFile, InputValues.Model);
      break;
    case "in_eta":
      getline(InputFile, line);
      InputValues.in_eta = atof(line.c_str());
      break;
    case "in_theta":
      getline(InputFile, line);
      InputValues.in_theta = atof(line.c_str());
      break;
    case "KC":
      getline(InputFile, line);
      InputValues.KC = atof(line.c_str());
      break;
    case "KG":
      getline(InputFile, line);
      InputValues.KG = atof(line.c_str());
      break;
    case "C0":
      getline(InputFile, line);
      InputValues.C0 = atof(line.c_str());
      break;
    case "mu":
      getline(InputFile, line);
      InputValues.mu = atof(line.c_str());
      break;
    case "kS":
      getline(InputFile, line);
      InputValues.kS = atof(line.c_str());
      break;
    case "min":
      getline(InputFile, line);
      InputValues.min = atoi(line.c_str());
      break;
    case "Ns":
      getline(InputFile, line);
      InputValues.Ns = atoi(line.c_str());
      break;
    case "AngleShift":
      getline(InputFile, line);
      InputValues.AngleShift = atof(line.c_str());
      break;
    case "OutlineName":
      getline(InputFile, InputValues.OutlineName);
      break;
    case "NPout":
      getline(InputFile, line);
      InputValues.NPout = atoi(line.c_str());
      break;
    case "HexonShift":
      getline(InputFile, line);
      InputValues.HexonShift = atof(line.c_str());
      break;
    case "WcType":
      getline(InputFile, line);
      InputValues.WcType = atoi(line.c_str());
      break;
    case "WcConst0":
      getline(InputFile, line);
      InputValues.WcConst0 = atof(line.c_str());
      break;
    case "WcConst1":
      getline(InputFile, line);
      InputValues.WcConst1 = atof(line.c_str());
      break;
    case "WcConst2":
      getline(InputFile, line);
      InputValues.WcConst2 = atof(line.c_str());
      break;
   
    default:
      cout << "Input " << line << " unknown";
  }
  
  InputFile.close();
  
  // List input parameters
  cout << " in_eta     : " << InputValues.in_eta      << endl
       << " in_theta   : " << InputValues.in_theta    << endl
       << " KC         : " << InputValues.KC          << endl
       << " KG         : " << InputValues.KG          << endl
       << " C0         : " << InputValues.C0          << endl
       << " mu         : " << InputValues.mu          << endl
       << " kS         : " << InputValues.kS          << endl
       << " min        : " << InputValues.min         << endl
       << " NstretchDOF: " << InputValues.Ns          << endl
       << " AngleShift : " << InputValues.AngleShift  << endl
       << " OutlineName: " << InputValues.OutlineName << endl
       << " NPout      : " << InputValues.NPout       << endl
       << " HexonShift : " << InputValues.HexonShift  << endl
       << " WcType     : " << InputValues.WcType      << endl
       << " keta       : " << InputValues.WcConst0    << endl
       << " ktheta     : " << InputValues.WcConst1    << endl
       << " eta bar    : " << InputValues.WcConst2    << endl;

  return InputValues; 
}



void printAllParaview(const std::string fileNameA,
		      const std::string fileNameB,
		      const std::string OutlineName,
		      const unsigned int NPout,
		      const Model::BodyContainer & bdc, 
		      const vector<ScalarFieldNode<3>* > & stretchNodes, 
		      const vector<ScalarFieldNode<3>* > & directionNodes,  
		      const vector<unsigned int> & CapsomersNum)
{
  // OutPut bending and stretching energy in current configuration
  ofstream ofs(fileNameA.c_str());
  if (!ofs) {
    std::cout << "can not open output ("
	      << fileNameA
	      << ") file." << std::endl;
    exit(0);
  }

  std::string OutlineFile = OutlineName+".vtk";
  ofstream ofs2(fileNameB.c_str());
  if (!ofs2) {
    std::cout << "can not open output ("
	      << fileNameB
	      << ") file." << std::endl;
    exit(0);
  }
  unsigned int ind = 0;

  // Start from bending body which has all elements
  Body::ElementContainer BendingElements = bdc[0]->elements();
  Body::NodeContainer BendingNodes = bdc[0]->nodes();

  // Node Section
  ofs << "# vtk DataFile Version 3.0" << endl
	<< "Test example" << endl
	<< "ASCII" << endl
	<< "DATASET POLYDATA" << endl
        << "POINTS  " << BendingNodes.size() << "  double" << endl;

  ofs2 << "# vtk DataFile Version 3.0" << endl
	 << "vtk output" << endl
	 << "ASCII" << endl
 	 << "DATASET POLYDATA" << endl
         << "POINTS " << NPout << " float" << endl;

  // Output nodal postions
  Body::NodeIterator pn; 
  for (pn = BendingNodes.begin(); pn!= BendingNodes.end(); pn ++)
  {
    DeformationNode<3> * node = dynamic_cast<DeformationNode<3>* >(*pn);

    if (node == NULL) exit(0);

    const Vector3D & nodalPos =  node->point();
    ofs << std::setprecision(16) 
	<< nodalPos(0) << "  "
	<< nodalPos(1) << "  "
	<< nodalPos(2) << std::endl;

    // Update nodal position for the outline (up to 212 nodes)
    // and considering that the nodes are always in the same order (all the mesh need to be consistent)
    if (ind < NPout)
    {
      ofs2 << std::setprecision(16) 
	   << nodalPos(0) << "  "
	   << nodalPos(1) << "  "
	   << nodalPos(2) << std::endl;
      ind++;
    } 
  }

  // Element Section
  int nel = BendingElements.size();
  ofs << "POLYGONS  " << nel << "  "
      << 4*nel << endl;
  for(int e = 0; e < nel; e++)
  {
    if(!(bdc[0]->active(e)) ) exit(0); // All elements should be active in this application
    // ::FElement_t * pe = dynamic_cast<LoopShellBody::FElement_t *>(BendingElements[e]);
    ofs << 3 << "  "
	<< setw(10) << BendingElements[e]->baseNodes()[0] -> id()
	<< setw(10) << BendingElements[e]->baseNodes()[1] -> id()
	<< setw(10) << BendingElements[e]->baseNodes()[2] -> id()
	<< endl;
  }





  // Outline for pentamers and examers
  ifstream ifs(OutlineFile.c_str());
  if (!ifs) {
    cout << "Cannot open input file: " << OutlineFile << endl;
    exit(0);
  }

  string token;
  ifs >> token; 
  while( token != "LINES") ifs >> token; 
  ofs2 << token << " ";
  ifs >> token;
  ofs2 << token << " ";
  ifs >> token;
  ofs2 << token << endl;
  ifs >> token;
  while(!ifs.eof())
  {
    ofs2 << token << " ";
    ifs >> token;
  }
  ifs.close();





  // Cell section
  ofs << "CELL_DATA    " << nel << endl;
  // Bending energy
  ofs << "SCALARS    BendingEnergy    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  for(int e = 0; e < nel; e++) {
      ofs << BendingElements[e]->energy() << endl;
  }
  ofs << std::endl;

  // Stretching energy
  Body::ElementContainer PentElements = bdc[1]->elements();
  Body::ElementContainer HexElements = bdc[2]->elements();
  ofs << "SCALARS    StretchingEnergy    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  
  for(int e = 0; e < HexElements.size(); e++) {
    ofs << HexElements[e]->energy() << endl;
  }

  for(int e = 0; e < PentElements.size(); e++) {
      ofs << PentElements[e]->energy() << endl;
  }
  ofs << endl;

 // Conformational energy
  ofs << "SCALARS    ConformationalEnergy    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  
  for(int e = 0; e < HexElements.size(); e++) {
      C0MembraneStretch * TempE = dynamic_cast<C0MembraneStretch *>(HexElements[e]);
      ofs << TempE->conformationalEnergy() << endl;
  }

  for(int e = 0; e < PentElements.size(); e++) {
      ofs << 0.0 << endl;
  }
  ofs << endl;



  // Stretch magnitude
  ofs << "SCALARS    MaxPrStr    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  double etaPrint = 0.0; 
  
  if (stretchNodes.size() == 1) {
    double etaPrint = stretchNodes[0]->point(); 
    for(int e = 0; e < HexElements.size(); e++) {
      if (etaPrint > 1.0) {
	ofs << etaPrint << endl;
      }
      else {
	ofs << 1.0/etaPrint << endl;
      }
    }
  }
  else {
    for(int e = 0; e < HexElements.size(); e++) {
      double etaPrint = stretchNodes[CapsomersNum[e]]->point(); 
      if (etaPrint > 1.0) {
	ofs << etaPrint << endl;
      }
      else {
	ofs << 1.0/etaPrint << endl;
      }
    }
  }

  for(int e = 0; e < PentElements.size(); e++) {
      ofs << 0.0 << endl;
  }
  ofs << endl;

  ofs << "SCALARS    MinPrStr    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  
 if (stretchNodes.size() == 1) {
    double etaPrint = stretchNodes[0]->point(); 
    for(int e = 0; e < HexElements.size(); e++) {
      if (etaPrint < 1.0) {
	ofs << etaPrint << endl;
      }
      else {
	ofs << 1.0/etaPrint << endl;
      }
    }
  }
  else {
    for(int e = 0; e < HexElements.size(); e++) {
      double etaPrint = stretchNodes[CapsomersNum[e]]->point(); 
      if (etaPrint < 1.0) {
	ofs << etaPrint << endl;
      }
      else {
	ofs << 1.0/etaPrint << endl;
      }
    }
  }

  for(int e = 0; e < PentElements.size(); e++) {
      ofs << 0.0 << endl;
  }
  ofs << endl;



  // Stretch direction
  ofs << "SCALARS    StretchDirection    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  
  if (directionNodes.size() == 1) {
    for(int e = 0; e < HexElements.size(); e++) {
      ofs << directionNodes[0]->point() << endl;
    }
  }
  else {
    for(int e = 0; e < HexElements.size(); e++) {
      ofs << directionNodes[CapsomersNum[e]]->point() << endl;
    }
  }

  for(int e = 0; e < PentElements.size(); e++) {
      ofs << 0.0 << endl;
  }
  ofs << endl;

  ofs.close();
  ofs2.close();

}



void printStretchDirectioParaview(const std::string fileName,
				const std::string meshFileName,
				const Model::BodyContainer & bdc, 
				const double AngleShift,
				const double HexonShift,
				const vector<ScalarFieldNode<3>* > & stretchNodes,
				const vector<ScalarFieldNode<3>* > & directionNodes,
				const vector<unsigned int> & CapsomersNum,
				bool final)
{
  const double PI = 3.14159265;
  double alpha = AngleShift*PI/180;
  // OutPut bending and stretching energy in current configuration
  ofstream ofs(fileName.c_str());
  if (!ofs) {
    std::cout << "can not open output ("
	      << fileName
	      << ") file." << std::endl;
    exit(0);
  }

  // Start from bending body which has all elements
  Body::ElementContainer BendingElements = bdc[0]->elements();
  Body::NodeContainer BendingNodes = bdc[0]->nodes();

  // As many points as are the hexons 
  int NumHexon = CapsomersNum.back()-11;
  ofs << "# vtk DataFile Version 2.0\n"
      << "Test example" << endl
      << "ASCII" << endl
      << "DATASET POLYDATA" << endl
      << "POINTS  " << 4*NumHexon << "  double" << endl;

  // Compute geometric center of each hexon
  int CurrentCapsomer =  CapsomersNum.front();
  int i = -1, j = 0;
  double ind = 0.0;
  vector<Vector3D > HexonCenters;
  Vector3D Hexon(0.0);
  while (CurrentCapsomer < NumHexon)
  {
    i++;
    DeformationNode<3> * nodeA = dynamic_cast<DeformationNode<3>* >(BendingElements[i]->baseNodes()[0]);
    DeformationNode<3> * nodeB = dynamic_cast<DeformationNode<3>* >(BendingElements[i]->baseNodes()[1]);
    DeformationNode<3> * nodeC = dynamic_cast<DeformationNode<3>* >(BendingElements[i]->baseNodes()[2]);

    const Vector3D & nodalPosA =  nodeA->point();
    const Vector3D & nodalPosB =  nodeB->point();
    const Vector3D & nodalPosC =  nodeC->point();
    
    if (CapsomersNum[i] == CurrentCapsomer)
    {
      Hexon(0) += nodalPosA(0) + nodalPosB(0) + nodalPosC(0);
      Hexon(1) += nodalPosA(1) + nodalPosB(1) + nodalPosC(1);
      Hexon(2) += nodalPosA(2) + nodalPosB(2) + nodalPosC(2);
      ind += 3.0;
      
    }
    else
    {
      Hexon(0) /= ind;
      Hexon(1) /= ind;
      Hexon(2) /= ind;
      HexonCenters.push_back(Hexon);
      
      Hexon(0) = nodalPosA(0) + nodalPosB(0) + nodalPosC(0);
      Hexon(1) = nodalPosA(1) + nodalPosB(1) + nodalPosC(1);
      Hexon(2) = nodalPosA(2) + nodalPosB(2) + nodalPosC(2);
      ind = 3.0;
      CurrentCapsomer = CapsomersNum[i];
    }
  }

  



  // Read in the initial stretch directions\
  // Create input stream
  ifstream ifs;
  ifs.open( meshFileName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << meshFileName << endl;
    exit(0);
  }
  // Create vector of initial stretch directions
  vector<tvmet::Vector<double,3> > InitialStretchDirections;
  vector<tvmet::Vector<double,3> > RotationAxis;
  
  string token;
  ifs >> token; 
  while( token != "shear_direction") ifs >> token; 
  ifs >> token;

  int SavedCapsomer = -1;
  CurrentCapsomer = CapsomersNum.front();
  i = 0;
  while (CurrentCapsomer < NumHexon)
  {
    if (CurrentCapsomer != SavedCapsomer)
    {
      tvmet::Vector<double,3> StretchDir;
      ifs >>  StretchDir(0) >> StretchDir(1) >> StretchDir(2);
      StretchDir = StretchDir/tvmet::norm2(StretchDir);
      InitialStretchDirections.push_back(StretchDir);
      // cout << StretchDir << CurrentCapsomer << endl;
      SavedCapsomer = CurrentCapsomer;

      // Computing axis of rotation
      DeformationNode<3> * nodeA = dynamic_cast<DeformationNode<3>* >(BendingElements[i]->baseNodes()[0]);
      DeformationNode<3> * nodeB = dynamic_cast<DeformationNode<3>* >(BendingElements[i]->baseNodes()[1]);
      DeformationNode<3> * nodeC = dynamic_cast<DeformationNode<3>* >(BendingElements[i]->baseNodes()[2]);

      tvmet::Vector<double,3> e31(nodeA->position()-nodeC->position()), 
 	                      e32(nodeB->position()-nodeC->position());
      
      tvmet::Vector<double,3> n = tvmet::cross(e31,e32);
      n = n/tvmet::norm2(n);

      RotationAxis.push_back(n);
    }
    else
    {
      ifs >> token >> token >> token;
    }
    i++;
    CurrentCapsomer = CapsomersNum[i];
    // cout << i << " " << CurrentCapsomer << " " << SavedCapsomer << endl;
  }





  // Compute stretch directions and positions
  vector<Vector3D > HexonPoint(4*NumHexon, Vector3D(0.0)),  Stretch(4*NumHexon, Vector3D(0.0));
  double beta = 0.0, eta = 0.0;
  for(i = 0; i < NumHexon; i++) {
    // Imposed initial angle shift (IN RADIANTS!!)
    if (final) {
      if (directionNodes.size() == 1) {
	beta = alpha + directionNodes[0]->point();
      }
      else {
	beta = alpha + directionNodes[i]->point();
      }
    }
    else {
      beta = alpha;
    }

    tvmet::Vector<double, 3> StretchLoc, StretchLocP;
    StretchLoc = InitialStretchDirections[i]*cos(beta) + tvmet::cross(RotationAxis[i],InitialStretchDirections[i])*sin(beta);
    StretchLoc /= tvmet::norm2(StretchLoc);

    // Move point of application of InitStretch in the Yloc direction
    StretchLocP = tvmet::cross(RotationAxis[i],StretchLoc);
    if (stretchNodes.size() == 1) {
      eta = stretchNodes[0]->point();
    }
    else {
      eta = stretchNodes[i]->point();
    }
    
    StretchLoc *= eta;
    StretchLocP /= eta;

    if (eta > 1.0-1.0e-12) {
      HexonPoint[i] = HexonCenters[i];
      HexonPoint[i+NumHexon] = HexonCenters[i];
      HexonPoint[i+2*NumHexon] = HexonCenters[i] + HexonShift*StretchLocP;
      HexonPoint[i+3*NumHexon] = HexonCenters[i] - HexonShift*StretchLocP;

      Stretch[i] = StretchLoc;
      Stretch[i+NumHexon] = -StretchLoc;
      Stretch[i+2*NumHexon] = -StretchLocP;
      Stretch[i+3*NumHexon] =  StretchLocP;
    }
    else {
      HexonPoint[i] = HexonCenters[i];
      HexonPoint[i+NumHexon] = HexonCenters[i];
      HexonPoint[i+2*NumHexon] = HexonCenters[i] + HexonShift*StretchLoc;
      HexonPoint[i+3*NumHexon] = HexonCenters[i] - HexonShift*StretchLoc;

      Stretch[i] = StretchLocP;
      Stretch[i+NumHexon] = -StretchLocP;
      Stretch[i+2*NumHexon] = -StretchLoc;
      Stretch[i+3*NumHexon] =  StretchLoc;

    }
    
  }



  // Print hexon position and stretch directions
  for(i = 0; i < 4*NumHexon; i++) {
    ofs << HexonPoint[i](0) << " " 
	<< HexonPoint[i](1) << " "
	<< HexonPoint[i](2) << endl;
  }
  ofs << std::endl;
 
  // Cell section
  ofs << "POINT_DATA    " << 4*NumHexon << endl;
  // Initial stretch directions
  ofs << "VECTORS    StretchDirections    double" << endl;

  for(i = 0; i < 4*NumHexon; i++) {
    ofs << Stretch[i](0) << " " 
	<< Stretch[i](1) << " "
	<< Stretch[i](2) << endl;
  }
  ofs << endl;
  
  ofs.close();

}





    
    
   
