#include "OPSBody.h"

namespace voom {

OPSBody::OPSBody(const vector<OPSNode*> & nodes, OPSParams &p, double r) :
    _opsNodes(nodes), _prop(p), _searchR(r) {
    _dof = 0;
    // Initialize Body.h containers
    _nodes.insert(_nodes.begin(), nodes.begin(), nodes.end());
    for (ConstNodeIterator n = _nodes.begin(); n != _nodes.end(); n++) {
        _dof += (*n)->dof();
    }

    _numNodes = nodes.size();
    _radius = 0.0;
    _volume = 0.0;
    _volConstraintOn = false;
    _volConstraint = 0.0;
    _PVen = 0.0;

    std::vector<Vector3D> volDiff(_numNodes, Vector3D(0.0));
    _volDiff = volDiff;

    // Construct vtkPolyData from nodes
    _polyData = vtkSmartPointer<vtkPolyData>::New();
    //Construct the kd-tree
    _kdTree = vtkSmartPointer<vtkKdTree>::New();
    updatePolyDataAndKdTree();

    //Initialize the _neighbors vector
    for(int i=0; i < _numNodes; i++){
        _neighbors.push_back( vtkSmartPointer<vtkIdList>::New() );
    }
    updateNeighbors();

    //Set the initial nearest neighbor vector
    for (int i = 0; i < _numNodes; i++) {
        Vector3D currPos(0.0);
        _polyData->GetPoint(i, &(currPos[0]));
        vtkSmartPointer<vtkIdList> neighbors =
                vtkSmartPointer<vtkIdList>::New();
        _kdTree->FindClosestNPoints(2, &(currPos[0]), neighbors);
        neighbors->DeleteId(i);
        _initialNearestNeighbor.push_back( neighbors->GetId(0) );
    }
}

/*
*Updates _polyData with latest deformed node positions and normals
 */
void OPSBody::updatePolyDataAndKdTree() {
    vtkSmartPointer<vtkPoints> pts =
            vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> pts2 =
            vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> normals =
            vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkPolyData> unitSphere =
            vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> final;
    vtkSmartPointer<vtkDataSetSurfaceFilter> dssf =
            vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    vtkSmartPointer<vtkIdFilter> idf =
            vtkSmartPointer<vtkIdFilter>::New();
    vtkSmartPointer<vtkDelaunay3D> d3D =
            vtkSmartPointer<vtkDelaunay3D>::New();
    vtkSmartPointer<vtkCellArray> finalCells =
            vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkCellArray> interim;
    vtkSmartPointer<vtkIdTypeArray> origIds;
    vtkSmartPointer<vtkIdList> pointIds =
            vtkSmartPointer<vtkIdList>::New();

    normals->SetNumberOfComponents(3);
    normals->SetName("PointNormals");

    for (opsNodeIterator n = _opsNodes.begin(); n != _opsNodes.end(); ++n) {
        Vector3D coords, currRotVec, currNormal;
        coords = (*n)->deformedPosition();
        currRotVec = (*n)->deformedRotationVector();
        currNormal = OPSNode::convertRotVecToNormal( currRotVec );
        pts->InsertNextPoint(&(coords[0]));
        normals->InsertNextTuple(&(currNormal[0]));
    }
    _polyData->SetPoints(pts);
    _polyData->GetPointData()->SetNormals(normals);
    _kdTree->BuildLocatorFromPoints( _polyData );

    for(int i=0; i < _polyData->GetNumberOfPoints(); i++){
        Vector3D x(0.0), temp(0.0);
        _polyData->GetPoint(i,&(temp[0]));
        x = temp / tvmet::norm2(temp);
        //We are rounding off to 2 decimals as a hack to aid vtkDelaunay3D
        //Otherwise sometimes we don't get a convex hull
        x[0] = std::round( x[0]*100 )/100;
        x[1] = std::round( x[1]*100 )/100;
        x[2] = std::round( x[2]*100 )/100;
        pts2->InsertNextPoint(&(x[0]));
    }
    unitSphere->SetPoints(pts2);

    idf->SetIdsArrayName("PointIds");
    idf->PointIdsOn();
    idf->SetInputData(unitSphere);

    //Calculate ideal number of triangles.
    int idealTriCount, pentCount = 12, hexCount;
    hexCount = _numNodes - pentCount;
    idealTriCount = (6*hexCount + 5*pentCount)/3;

    d3D->SetInputConnection(idf->GetOutputPort());
    dssf->SetInputConnection(d3D->GetOutputPort());
    dssf->Update();
    final = dssf->GetOutput();
    if( final->GetNumberOfPolys() != idealTriCount){
        cout<<"The mesh has " << final->GetNumberOfPolys() << " triangles." << std::endl;
        cout<< "Bad Delaunay triangulation detected!" <<std::endl;
        vtkSmartPointer<vtkPolyDataWriter> writer =
                vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetInputData(final);
        writer->SetFileName("BadMesh.vtk");
        writer->Write();
        exit(EXIT_FAILURE);
    }
    interim = final->GetPolys();
    interim->InitTraversal();
    origIds = vtkIdTypeArray::SafeDownCast(
                final->GetPointData()->GetArray("PointIds"));
    while(interim->GetNextCell(pointIds)){
        int numIds = pointIds->GetNumberOfIds();
        finalCells->InsertNextCell(numIds);
        for(int j=0; j < numIds; j++ ){
            int id = (int)origIds->GetTuple1( pointIds->GetId(j) );
            finalCells->InsertCellPoint(id);
        }
    }
    _polyData->SetPolys(finalCells);
/*
    vtkSmartPointer<vtkMassProperties> props =
            vtkSmartPointer<vtkMassProperties>::New();
    props->SetInputData(_polyData);
    props->Update();
    _volume = props->GetVolume();
*/
}

/*
*Store the neighbors information for each node
 */
void OPSBody::updateNeighbors(){
    for(int i=0; i < _numNodes; i++){
        Vector3D pos(0.0);
        _polyData->GetPoint( i, &(pos[0]) );
        _neighbors[i]->Reset();
        _kdTree->FindPointsWithinRadius( _searchR, &(pos[0]), _neighbors[i] );
        _neighbors[i]->DeleteId(i);
    }
}

/*
*Compute the OPSBody energy and energy contribution due to some other elements
*inherited from the parent Body class _elements container
 */
void OPSBody::compute(bool f0, bool f1, bool f2) {

    double E, re, a, b;
    double A, B, G;

    E = _prop.D_e; re = _prop.r_e; a = _prop.a; b = _prop.b;
    A = _prop.alpha; B = _prop.beta; G = _prop.gamma;

    // Initialize energy to be zero
    if (f0){
        _energy = 0.0;
        _morseEn = 0.0;        
        _normalEn = 0.0;
        _circEn = 0.0;
        _PVen = 0.0;
    }
    if(f0 && !f1){
        if(_volConstraintOn){
            double volume = calcAvgVolume();
            double dV = volume - _volConstraint; //Change in volume
            double P = _PV->point(); //Pressure
            _PVen = -P*dV;
            _energy += _PVen ;
        }
        for(int i=0; i < _numNodes; i++){
            Vector3D vi(0.0), xi(0.0), p(0.0);
            Matrix3X3 M(0.0);

            xi = _opsNodes[i]->deformedPosition();
            vi = _opsNodes[i]->deformedRotationVector();
            p = OPSNode::convertRotVecToNormal(vi);
            M = diffNormalRotVec(vi);

            for(int j=0; j < _neighbors[i]->GetNumberOfIds(); j++){
                double r, n_dot_rij, exp1, exp2;
                double morseEn, Ker, Phi_n, Phi_c;
                Vector3D vj(0.0), xj(0.0), q(0.0), m(0.0), n(0.0), rij(0.0);
                Matrix3X3 N(0.0);
                vtkIdType currId = _neighbors[i]->GetId(j);

                xj = _opsNodes[currId]->deformedPosition();
                vj = _opsNodes[currId]->deformedRotationVector();
                q = OPSNode::convertRotVecToNormal(vj);
                N = diffNormalRotVec(vj);
                rij = xj - xi;
                m = p - q;
                n = p + q;
                r = tvmet::norm2(rij);
                n_dot_rij = tvmet::dot(n,rij);

                exp1 = exp( -2*a*(r - re) );
                exp2 = exp( -a*(r - re) );

                morseEn = E*( exp1 - 2*exp2 );
                Ker = (E/G)*exp( -r*r/(2*b*b) );
                Phi_n = tvmet::norm2(m);
                Phi_n *= Phi_n;
                Phi_c = n_dot_rij/r;
                Phi_c *= Phi_c;

                _morseEn += morseEn;                
                _normalEn += Ker*Phi_n;
                _circEn += Ker*Phi_c;
                _energy += morseEn + Ker*(Phi_n + Phi_c);
            }
        }
    }
    else{
        std::vector<Vector3D> volDiff(_numNodes, Vector3D(0.0));
        if(_volConstraintOn){
            //Take care of the volume constraint first
            volDiff = calcAvgVolDerivative();
            double volume = calcAvgVolume();
            double dV = volume - _volConstraint; //Change in volume
            double P = _PV->point(); //Pressure
            _PVen = -P*dV;
            _energy += _PVen ;
            _PV->setForce(0, -dV );
        }
        for(int i=0; i < _numNodes; i++){
            Vector3D vi(0.0), xi(0.0), p(0.0);
            Matrix3X3 M(0.0);

            if(_volConstraintOn){
                //Add the Volume constraint force to the node
                Vector6D pvF(volDiff[i][0], volDiff[i][1], volDiff[i][2],0,0,0);
                pvF = _PV->point()*pvF;
                _opsNodes[i]->updateForce(pvF);
            }

            xi = _opsNodes[i]->deformedPosition();
            vi = _opsNodes[i]->deformedRotationVector();
            p = OPSNode::convertRotVecToNormal(vi);
            M = diffNormalRotVec(vi);

            for(int j=0; j < _neighbors[i]->GetNumberOfIds(); j++){
                double r, n_dot_rij, exp1, exp2;
                double morseEn, Ker, Phi_n, Phi_c;
                Vector3D vj(0.0), xj(0.0), q(0.0), m(0.0), n(0.0), rij(0.0);
                Vector3D dMorseXi, dMorseXj, dKerXi, dKerXj;
                Vector3D dPhi_nVi(0.0), dPhi_nVj(0.0);
                Vector3D dPhi_cVi(0.0), dPhi_cVj(0.0);
                Vector3D dPhi_cXi(0.0), dPhi_cXj(0.0);
                Matrix3X3 N(0.0);
                vtkIdType currId = _neighbors[i]->GetId(j);

                xj = _opsNodes[currId]->deformedPosition();
                vj = _opsNodes[currId]->deformedRotationVector();
                q = OPSNode::convertRotVecToNormal(vj);
                N = diffNormalRotVec(vj);
                rij = xj - xi;
                m = p - q;
                n = p + q;
                r = tvmet::norm2(rij);
                n_dot_rij = tvmet::dot(n,rij);                

                exp1 = exp( -2*a*(r - re) );
                exp2 = exp( -a*(r - re) );

                morseEn = E*( exp1 - 2*exp2 );
                Ker = (E/G)*exp( -r*r/(2*b*b) );
                Phi_n = tvmet::norm2(m);
                Phi_n *= Phi_n;
                Phi_c = n_dot_rij/r;
                Phi_c *= Phi_c;

                // Evaluate morse derivatives
                dMorseXi = (2*E*a/r)*( exp1 - exp2 )*rij;
                dMorseXj = -dMorseXi;

                // Evaluate kernel derivatives
                dKerXi = (Ker/(b*b))*rij;
                dKerXj = -dKerXi;

                //Evaluate co-normality derivatives
                dPhi_nVi = 2*M*m;
                dPhi_nVj = -2*N*m;

                //Evaluate co-circularity derivatives
                dPhi_cXi = (2*n_dot_rij/(r*r*r*r))*( n_dot_rij*rij -r*r*n );
                dPhi_cXj = -dPhi_cXi;
                dPhi_cVi = (2*n_dot_rij/(r*r))*M*rij;
                dPhi_cVj = (2*n_dot_rij/(r*r))*N*rij;

                Vector3D centerDx(0.0), centerDv(0.0), neighborDx(0.0), neighborDv(0.0);

                centerDx = dMorseXi + Ker*dPhi_cXi
                        + dKerXi*(Phi_n + Phi_c );
                centerDv = Ker*(dPhi_nVi + dPhi_cVi );

                neighborDx = dMorseXj + Ker*(dPhi_cXj)
                        + dKerXj*(Phi_n + Phi_c );
                neighborDv = Ker*(dPhi_nVj + dPhi_cVj);

                Vector6D centerForce(centerDx[0], centerDx[1], centerDx[2],
                        centerDv[0], centerDv[1], centerDv[2]);
                Vector6D neighborForce(neighborDx[0], neighborDx[1],
                        neighborDx[2], neighborDv[0], neighborDv[1], neighborDv[2]);

                for(int k=0; k < 6; k++){
                    _opsNodes[i]->addForce( k, centerForce(k) );
                    _opsNodes[currId]->addForce( k, neighborForce(k) );
                }

                _morseEn += morseEn;
                _normalEn += Ker*Phi_n;
                _circEn += Ker*Phi_c;
                _energy += morseEn + Ker*(Phi_n + Phi_c);
            }
        }        
    }

    //If the inherited _elements container is non-empty
    for (int i = 0; i < _elements.size(); i++) {
        _elements[i]->compute(f0, f1, f2);
        _energy += _elements[i]->energy();
    }

    return;

}

/*
 * Print a Paraview file by generating a triangulation
 */
void OPSBody::printParaview( const string fileName ) const{
    vtkSmartPointer<vtkPolyDataWriter> writer =
            vtkSmartPointer<vtkPolyDataWriter>::New();
    vtkSmartPointer<vtkIdFilter> idf =
            vtkSmartPointer<vtkIdFilter>::New();
    idf->SetInputData(_polyData);
    idf->PointIdsOn();
    idf->SetIdsArrayName("PointIds");
    idf->Update();

    writer->SetFileName( fileName.c_str() );
    writer->SetInputConnection(idf->GetOutputPort());
    writer->Write();
}

/*
 * Calculate average edge length as if the particles were
 * triangulated
 */
double OPSBody::getAverageEdgeLength(){
    double avg = 0;
    int numEdges = 0;
    for(int i=0; i < _numNodes; i++){
        Vector3D center;
        center = _opsNodes[i]->deformedPosition();
        for(int j=0; j < _neighbors[i]->GetNumberOfIds(); j++){
            vtkIdType currId = _neighbors[i]->GetId(j);
            Vector3D neighbor;
            neighbor = _opsNodes[currId]->deformedPosition();
            avg += tvmet::norm2( center - neighbor );
            numEdges++;
        }
    }
    avg /= numEdges;
    return avg;
}

/*
 * Get average radius
 */
double OPSBody::getAverageRadius(){
    _radius = 0.0;
    for(int i=0; i < _numNodes; i++){
        Vector3D x;
        x = _opsNodes[i]->deformedPosition();
        _radius += tvmet::norm2(x);
    }
    _radius /= _numNodes;
    return _radius;
}

/*
 * Calculate volume of the polydata and derivative
 */
void OPSBody::calcVolumeAndDerivative(){
    _volume = 0.0;
    vtkSmartPointer<vtkCellArray> cells = _polyData->GetPolys();
    vtkSmartPointer<vtkIdList> verts =
            vtkSmartPointer<vtkIdList>::New();
    cells->InitTraversal();
    while( cells->GetNextCell(verts) ){
        assert(verts->GetNumberOfIds() == 3);
        int ida, idb, idc;
        ida = verts->GetId(0);
        idb = verts->GetId(1);
        idc = verts->GetId(2);
        Vector3D a(0.0), b(0.0), c(0.0), crossTemp(0.0);
        _polyData->GetPoint( ida, &a[0] );
        _polyData->GetPoint( idb, &b[0] );
        _polyData->GetPoint( idc, &c[0] );

        //Calculate volume
        crossTemp = tvmet::cross(b,c);
        _volume += (tvmet::dot( a, crossTemp ))/6.0;

        //Calculate derivative wrt volume
        Vector3D dVa( b[1]*c[2]-b[2]*c[1],b[2]*c[0]-b[0]*c[2],b[0]*c[1]-b[1]*c[0] );
        Vector3D dVb( a[2]*c[1]-a[1]*c[2],a[0]*c[2]-a[2]*c[0],a[1]*c[0]-a[0]*c[1] );
        Vector3D dVc( a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0] );
        Vector3D temp;

        temp = _volDiff[ida];
        _volDiff[ida] = temp + (1/6)*dVa;
        temp = _volDiff[idb];
        _volDiff[idb] = temp + (1/6)*dVb;
        temp = _volDiff[idc];
        _volDiff[idc] = temp + (1/6)*dVc;
    }
}

/*
 * Calculate volume of the polydata
 */
double OPSBody::calcVolume(){
    double volume = 0.0;
    vtkSmartPointer<vtkCellArray> cells = _polyData->GetPolys();
    vtkSmartPointer<vtkIdList> verts =
            vtkSmartPointer<vtkIdList>::New();
    cells->InitTraversal();
    while( cells->GetNextCell(verts) ){
        assert(verts->GetNumberOfIds() == 3);
        int ida, idb, idc;
        ida = verts->GetId(0);
        idb = verts->GetId(1);
        idc = verts->GetId(2);
        Vector3D a(0.0), b(0.0), c(0.0), crossTemp(0.0);
        _polyData->GetPoint( ida, &a[0] );
        _polyData->GetPoint( idb, &b[0] );
        _polyData->GetPoint( idc, &c[0] );

        //Calculate volume
        crossTemp = tvmet::cross(b,c);
        volume += (tvmet::dot( a, crossTemp ))/6.0;
    }
    return volume;
}

/*
 * Calculate volume derivative
 */
std::vector<Vector3D> OPSBody::calcVolumeDerivative(){
    vtkSmartPointer<vtkCellArray> cells = _polyData->GetPolys();
    vtkSmartPointer<vtkIdList> verts =
            vtkSmartPointer<vtkIdList>::New();
    std::vector<Vector3D> dV( _polyData->GetNumberOfPoints(), Vector3D(0.0) );
    cells->InitTraversal();
    while( cells->GetNextCell(verts) ){
        assert(verts->GetNumberOfIds() == 3);
        int ida, idb, idc;
        ida = verts->GetId(0);
        idb = verts->GetId(1);
        idc = verts->GetId(2);
        Vector3D a(0.0), b(0.0), c(0.0);
        _polyData->GetPoint( ida, &a[0] );
        _polyData->GetPoint( idb, &b[0] );
        _polyData->GetPoint( idc, &c[0] );

        //Calculate derivative wrt volume
        Vector3D dVa( b[1]*c[2]-b[2]*c[1],b[2]*c[0]-b[0]*c[2],b[0]*c[1]-b[1]*c[0] );
        Vector3D dVb( a[2]*c[1]-a[1]*c[2],a[0]*c[2]-a[2]*c[0],a[1]*c[0]-a[0]*c[1] );
        Vector3D dVc( a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0] );        

        dV[ida] += (1/6)*dVa;
        dV[idb] += (1/6)*dVb;
        dV[idc] += (1/6)*dVc;
    }
    return dV;
}
/*
 * Calculate derivative of the average volume that does not need triangulation
 */
std::vector<Vector3D> OPSBody::calcAvgVolDerivative(){
    std::vector<Vector3D> dV( _numNodes, Vector3D(0.0) );
    double R = getAverageRadius();
    double factor = 12.5663706144*R*R/_numNodes; //4*pi*R^2/N

    for(int i=0; i < _numNodes; i++){
        double x_norm;
        Vector3D x;
        x = _opsNodes[i]->deformedPosition();
        x_norm = tvmet::norm2(x);
        dV[i] = (factor/x_norm)*x;
    }
    return dV;
}

/*
 *  Calculate asphericity
 */
double OPSBody::getAsphericity(){
    double asphericity = 0.0;
    double R0 = _radius;
    for(int i=0; i < _numNodes; i++){
        Vector3D x;
        x = _opsNodes[i]->deformedPosition();
        double R = tvmet::norm2(x);
        asphericity += (R - R0)*(R - R0);
    }
    asphericity /= (_numNodes*R0*R0);
    return asphericity;
}
/*
 * Get Loop subdivision surface asphericity
 */
double OPSBody::getLoopAsphericity(){
    double asphericity = 0;
    vtkSmartPointer<vtkLoopSubdivisionFilter> lsf =
            vtkSmartPointer<vtkLoopSubdivisionFilter>::New();
    lsf->SetInputData(_polyData);
    lsf->SetNumberOfSubdivisions(1);
    lsf->Update();
    vtkSmartPointer<vtkPolyData> sub = lsf->GetOutput();
    //Quadrature order 2 for a triangle has 3 edge mid-points as Gauss points
    vtkSmartPointer<vtkCellArray> polys = sub->GetPolys();
    double R0 = 0.0;
    int numQuadPoints = 0;
    polys->InitTraversal();
    vtkSmartPointer<vtkIdList> verts =
            vtkSmartPointer<vtkIdList>::New();
    int count = 0;
    while(polys->GetNextCell(verts)){
        int numVerts = verts->GetNumberOfIds();
        if(numVerts == 3){
            Vector3D A(0.0), B(0.0), C(0.0), D(0.0), E(0.0), F(0.0);
            sub->GetPoint( verts->GetId(0), &(A[0]) );
            sub->GetPoint( verts->GetId(1), &(B[0]) );
            sub->GetPoint( verts->GetId(2), &(C[0]) );
            D = 0.5*(A+B);
            E = 0.5*(B+C);
            F = 0.5*(C+A);
            R0 += (tvmet::norm2(D) + tvmet::norm2(E) + tvmet::norm2(F));
            numQuadPoints += 3;
        }
    }
    R0 /= numQuadPoints;
    verts->Reset();
    polys->InitTraversal();
    while(polys->GetNextCell(verts)){
        int numVerts = verts->GetNumberOfIds();
        double term1, term2, term3;
        if(numVerts == 3){
            Vector3D A(0.0), B(0.0), C(0.0), D(0.0), E(0.0), F(0.0);
            sub->GetPoint( verts->GetId(0), &(A[0]) );
            sub->GetPoint( verts->GetId(1), &(B[0]) );
            sub->GetPoint( verts->GetId(2), &(C[0]) );
            D = 0.5*(A+B);
            E = 0.5*(B+C);
            F = 0.5*(C+A);

            term1 = (tvmet::norm2(D) - R0);
            term2 = (tvmet::norm2(E) - R0);
            term3 = (tvmet::norm2(F) - R0);

            asphericity += (term1*term1 + term2*term2 + term3*term3);
        }
    }
    asphericity /= (numQuadPoints/R0*R0);
    return asphericity;
}

/*
 * The mean squared displacement
 */
double OPSBody::msd() {
    double msd = 0;
    int nn = -1;
    // We will subtract off the radial displacement.
    for (int i = 0; i < _numNodes; i++) {
        Vector3D xi, xj, diff, xi_diff, xj_diff;
        Vector3D xi0, xj0, xi1, xj1, ni0, nj0;

        nn = _initialNearestNeighbor[i];

        xi0 = _opsNodes[i]->referencePosition();
        xi1 = _opsNodes[i]->deformedPosition();
        ni0 = xi0 / norm2( xi0 );

        xj0 = _opsNodes[nn]->referencePosition();
        xj1 = _opsNodes[nn]->deformedPosition();
        nj0 = xj0 / norm2( xj0 );

        xi_diff = (xi1 - xi0);
        xj_diff = (xj1 - xj0);

        xi = xi_diff - (tvmet::dot(xi_diff, ni0)*ni0);
        xj = xj_diff - (tvmet::dot(xj_diff, nj0)*nj0);

        diff = xi - xj;
        msd += tvmet::dot(diff, diff);
    }
    msd = msd / (2*_numNodes);
    return msd;
}

/*
 * Update property of the OPSBody
 */
void OPSBody::updateProperty(Property p, double val) {
    switch (p) {
    case E:
        _prop.D_e = val;
        break;
    case r:
        _prop.r_e = val;
        break;
    case av:
        _prop.a = val;
        break;    
    case bv:
        _prop.b = val;
        break;
    case A:
        _prop.alpha = val;
        break;
    case B:
        _prop.beta = val;
        break;
    case G:
        _prop.gamma = val;
        break;
    }
}

/*
 * Convert rotation vector to unit normal. A similar function
 * exists in OPSNode. We need to verify if they are equal.
 */
Vector3D OPSBody::rotVecToNormal(Vector3D U){
    // Assume z-axis of Global Coord Sys is the reference for ptNormal rotation
    Vector3D zaxis(0.0);
    zaxis[2] = 1.0;
    Vector3D p(0.0);
    double u = tvmet::norm2(U); // Angle of rotation
    // Check 0-angle case
    if (u < 1e-10) {
        p[2] = 1.0;
    }
    else if (std::abs(u - M_PI) < 1e-10) {
        p[2] = -1.0;
    }
    else {
        double sinA = sin( 0.5*u );
        double cosA = cos( 0.5*u );
        p[0] = (2*sinA/(u*u))*( U[1]*u*cosA + U[0]*U[2]*sinA );
        p[1] = (2*sinA/(u*u))*( U[1]*U[2]*sinA - U[0]*u*cosA );
        p[2] = cosA*cosA + (sinA*sinA/(u*u))*( U[2]*U[2] - U[1]*U[1] - U[0]*U[0] );
    }
    return p;
}

/*
 * Differentiate normal vector with respect to components of
 * rotation vector
 */
OPSBody::Matrix3X3 OPSBody::diffNormalRotVec(Vector3D vi){
    double u0, u1, u2, u0_2, u1_2, u2_2;
    double u, u_2, u_3, u_4;
    double sin_alpha, cos_alpha, sin_alpha_2;
    double dp0du0, dp1du0, dp2du0;
    double dp0du1, dp1du1, dp2du1;
    double dp0du2, dp1du2, dp2du2;
    Matrix3X3 finalMat(0.0);

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    u0_2 = u0*u0; u1_2 = u1*u1; u2_2 = u2*u2;
    u_2 = u0_2 + u1_2 + u2_2;
    u = sqrt(u_2);
    u_3 = u_2*u;
    u_4 = u_3*u;
    sin_alpha = sin( 0.5*u );
    sin_alpha_2 = sin_alpha*sin_alpha;
    cos_alpha = cos( 0.5*u );

    dp0du0 = 2*(-u0*u1*sin_alpha/2 + u2*sin_alpha + u0_2*u2*cos_alpha/(2*u) +
                u0*u1*cos_alpha/u)*sin_alpha/u_2 +
            u0*(u*u1*cos_alpha + u0*u2*sin_alpha)*cos_alpha/u_3 -
            4*u0*(u*u1*cos_alpha + u0*u2*sin_alpha)*sin_alpha/u_4;

    dp1du0 = 2*(-u*cos_alpha + u0_2*sin_alpha/2 - u0_2*cos_alpha/u +
                u0*u1*u2*cos_alpha/(2*u))*sin_alpha/u_2 +
            u0*(-u*u0*cos_alpha + u1*u2*sin_alpha)*cos_alpha/u_3 -
            4*u0*(-u*u0*cos_alpha + u1*u2*sin_alpha)*sin_alpha/u_4;

    dp2du0 = -u0*sin_alpha*cos_alpha/u - 2*u0*sin_alpha_2/u_2 +
            u0*(-u0_2 - u1_2 + u2_2)*sin_alpha*cos_alpha/u_3 -
            2*u0*(-u0_2 - u1_2 + u2_2)*sin_alpha_2/u_4;

    dp0du1 = 2*(u*cos_alpha - u1_2*sin_alpha/2 + u0*u1*u2*cos_alpha/(2*u) +
                u1_2*cos_alpha/u)*sin_alpha/u_2 +
            u1*(u*u1*cos_alpha + u0*u2*sin_alpha)*cos_alpha/u_3 -
            4*u1*(u*u1*cos_alpha + u0*u2*sin_alpha)*sin_alpha/u_4;

    dp1du1 = 2*(u0*u1*sin_alpha/2 + u2*sin_alpha - u0*u1*cos_alpha/u +
                u1_2*u2*cos_alpha/(2*u))*sin_alpha/u_2 +
            u1*(-u*u0*cos_alpha + u1*u2*sin_alpha)*cos_alpha/u_3 -
            4*u1*(-u*u0*cos_alpha + u1*u2*sin_alpha)*sin_alpha/u_4;

    dp2du1 = -u1*sin_alpha*cos_alpha/u - 2*u1*sin_alpha_2/u_2 +
            u1*(-u0_2 - u1_2 + u2_2)*sin_alpha*cos_alpha/u_3 -
            2*u1*(-u0_2 - u1_2 + u2_2)*sin_alpha_2/u_4;

    dp0du2 = 2*(u0*sin_alpha - u1*u2*sin_alpha/2 + u0*u2_2*cos_alpha/(2*u) +
                u1*u2*cos_alpha/u)*sin_alpha/u_2 +
            u2*(u*u1*cos_alpha + u0*u2*sin_alpha)*cos_alpha/u_3 -
            4*u2*(u*u1*cos_alpha + u0*u2*sin_alpha)*sin_alpha/u_4;

    dp1du2 = 2*(u0*u2*sin_alpha/2 + u1*sin_alpha - u0*u2*cos_alpha/u +
                u1*u2_2*cos_alpha/(2*u))*sin_alpha/u_2 +
            u2*(-u*u0*cos_alpha + u1*u2*sin_alpha)*cos_alpha/u_3 -
            4*u2*(-u*u0*cos_alpha + u1*u2*sin_alpha)*sin_alpha/u_4;

    dp2du2 = -u2*sin_alpha*cos_alpha/u + 2*u2*sin_alpha_2/u_2 +
            u2*(-u0_2 - u1_2 + u2_2)*sin_alpha*cos_alpha/u_3 -
            2*u2*(-u0_2 - u1_2 + u2_2)*sin_alpha_2/u_4;

    finalMat = dp0du0, dp1du0, dp2du0,
            dp0du1, dp1du1, dp2du1,
            dp0du2, dp1du2, dp2du2;
    return finalMat;
}

/*
*Morse function
 */
double OPSBody::morse(Vector3D xi, Vector3D xj){
    double E, re, a, r, potential;
    E = _prop.D_e; re = _prop.r_e; a = _prop.a;
    r = tvmet::norm2(xi - xj);
    potential = E*( -2*exp( -a*(r-re) ) + exp( -2*a*(r-re) ) );
    return potential;
}

/*
*The kernel function
 */
double OPSBody::psi(Vector3D xi, Vector3D xj){
    double K, b, r, psi0;
    Vector3D rij(0.0);

    rij = xj - xi;
    r = tvmet::norm2(rij);
    K = (_prop.D_e/_prop.gamma); b = _prop.b;
    psi0 = K*exp(- (r*r) / (2*b*b) );

    return psi0;
}

/*
 * The Co-Circularity function

double OPSBody::phi_p(Vector3D vi, Vector3D xi, Vector3D xj){
    double r, phi_p0;
    Vector3D p(0.0), rij(0.0);

    p = OPSNode::convertRotVecToNormal(vi);
    rij = xj - xi;
    r = tvmet::norm2(rij);
    phi_p0 = tvmet::dot(rij, p)/r;
    phi_p0 *= phi_p0;

    return phi_p0;
}
*/

/*
 * The Co-Normality function
 */
double OPSBody::phi_n(Vector3D vi, Vector3D vj){
    double phi_n0;
    Vector3D p(0.0), q(0.0);
    p = OPSNode::convertRotVecToNormal(vi);
    q = OPSNode::convertRotVecToNormal(vj);
    phi_n0 = tvmet::norm2( p - q );
    phi_n0 *= phi_n0;
    return phi_n0;
}

/*
 * The Co-Circularity function
 */
double OPSBody::phi_c(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj){
    double r, phi_c0;
    Vector3D p(0.0), q(0.0), n(0.0), rij(0.0);

    p = OPSNode::convertRotVecToNormal(vi);
    q = OPSNode::convertRotVecToNormal(vj);
    n = p + q;
    rij = xj - xi;
    r = tvmet::norm2(rij);
    phi_c0 = tvmet::dot(rij, n)/r;
    phi_c0 *= phi_c0;

    return phi_c0;
}

/*
 * Derivative of Morse with respect to xi
 */
Vector3D OPSBody::DmorseDxi(Vector3D xi, Vector3D xj){
    double epsilon, l0, a, r;
    Vector3D rij(0.0), DmorseDxi( 0.0 );;
    epsilon = _prop.D_e; l0 = _prop.r_e; a = _prop.a;

    rij = xj-xi;
    r = tvmet::norm2( rij );
    DmorseDxi = -(2/r)*epsilon*a*( exp(-a*(-l0+r)) - exp(-2*a*(-l0+r)) )*rij;

    return DmorseDxi;
}

/*
 * Derivative psi wrt xi
 */
Vector3D OPSBody::DpsiDxi(Vector3D xi, Vector3D xj){
    double K, b, r;
    Vector3D DpsiDxi0(0.0), rij(0.0);

    K = (_prop.D_e/_prop.gamma); b = _prop.b;
    rij = xj - xi;
    r = tvmet::norm2(rij);
    DpsiDxi0 = (K/(b*b))*exp(-r*r/(2*b*b))*rij;
    return DpsiDxi0;
}

/*
 * Derivative of phi_p wrt xi

Vector3D OPSBody::Dphi_pDxi(Vector3D vi, Vector3D xi, Vector3D xj){
    double r, p_dot_rij;
    Vector3D p(0.0), rij(0.0), ans(0.0);

    p = OPSNode::convertRotVecToNormal(vi);
    rij = xj - xi;
    r = tvmet::norm2(rij);
    p_dot_rij = tvmet::dot(p,rij);
    ans = (2*p_dot_rij)/(r*r*r*r)*( p_dot_rij*rij - r*r*p );

    return ans;
}
 */
/*
 * Derivative of phi_p wrt xj

Vector3D OPSBody::Dphi_pDxj(Vector3D vi, Vector3D xi, Vector3D xj){
    Vector3D ans(0.0);
    ans = -Dphi_pDxi(vi,xi,xj);
    return ans;
}
 */
/*
 * Derivative of phi_p wrt vi
Vector3D OPSBody::Dphi_pDvi(Vector3D vi, Vector3D xi, Vector3D xj){
    double r, p_dot_rij;
    Vector3D p(0.0), rij(0.0), ans(0.0);
    Matrix3X3 M(0.0);

    p = OPSNode::convertRotVecToNormal(vi);
    rij = xj - xi;
    r = tvmet::norm2(rij);
    p_dot_rij = tvmet::dot(p,rij);
    M = diffNormalRotVec(vi);
    ans = (2*p_dot_rij)/(r*r)*( M*rij );

    return ans;
}
 */
/*
 * Derivative of phi_n wrt vi
 */
Vector3D OPSBody::Dphi_nDvi(Vector3D vi, Vector3D vj){
    Vector3D p(0.0), q(0.0), m(0.0), ans(0.0);
    Matrix3X3 M(0.0);

    p = OPSNode::convertRotVecToNormal(vi);
    q = OPSNode::convertRotVecToNormal(vj);
    m = p - q;
    M = diffNormalRotVec(vi);
    ans = 2*M*m;

    return ans;
}

/*
 * Derivative of phi_n wrt vj
 */
Vector3D OPSBody::Dphi_nDvj(Vector3D vi, Vector3D vj){
    Vector3D p(0.0), q(0.0), m(0.0), ans(0.0);
    Matrix3X3 N(0.0);

    p = OPSNode::convertRotVecToNormal(vi);
    q = OPSNode::convertRotVecToNormal(vj);
    m = p - q;
    N = diffNormalRotVec(vj);
    ans = -2*N*m;

    return ans;
}

/*
 * Derivative of phi_c wrt xi
 */
Vector3D OPSBody::Dphi_cDxi(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj){
    double r, n_dot_rij;
    Vector3D p(0.0), q(0.0), n(0.0), rij(0.0), ans(0.0);

    p = OPSNode::convertRotVecToNormal(vi);
    q = OPSNode::convertRotVecToNormal(vj);
    n = p + q;
    rij = xj - xi;
    r = tvmet::norm2(rij);
    n_dot_rij = tvmet::dot(n,rij);
    ans = (2*n_dot_rij)/(r*r*r*r)*( n_dot_rij*rij - r*r*n );

    return ans;
}

/*
 * Derivative of phi_c wrt xj
 */
Vector3D OPSBody::Dphi_cDxj(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj){
    Vector3D ans(0.0);
    ans = -Dphi_cDxi(vi,vj,xi,xj);
    return ans;
}

/*
 * Derivative of phi_c wrt vi
 */
Vector3D OPSBody::Dphi_cDvi(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj){
    double r, n_dot_rij;
    Vector3D p(0.0), q(0.0), n(0.0), rij(0.0), ans(0.0);
    Matrix3X3 M(0.0);

    p = OPSNode::convertRotVecToNormal(vi);
    q = OPSNode::convertRotVecToNormal(vj);
    n = p + q;
    rij = xj - xi;
    r = tvmet::norm2(rij);
    n_dot_rij = tvmet::dot(n,rij);
    M = diffNormalRotVec(vi);
    ans = (2*n_dot_rij)/(r*r)*( M*rij );

    return ans;
}

/*
 * Derivative of phi_c wrt vj
 */
Vector3D OPSBody::Dphi_cDvj(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj){
    double r, n_dot_rij;
    Vector3D p(0.0), q(0.0), n(0.0), rij(0.0), ans(0.0);
    Matrix3X3 N(0.0);

    p = OPSNode::convertRotVecToNormal(vi);
    q = OPSNode::convertRotVecToNormal(vj);
    n = p + q;
    rij = xj - xi;
    r = tvmet::norm2(rij);
    n_dot_rij = tvmet::dot(n,rij);
    N = diffNormalRotVec(vj);
    ans = (2*n_dot_rij)/(r*r)*( N*rij );

    return ans;
}

}
