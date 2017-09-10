/*
*OPSBody.cc
** Created on: Aug 5, 2017
*     Author: amit
 */
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
    _numBadTri = 0;

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

    int ptIdx = 0;
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
        pts2->InsertNextPoint(&(x[0]));
    }
    unitSphere->SetPoints(pts2);

    idf->SetIdsArrayName("PointIds");
    idf->PointIdsOn();
    idf->SetInputData(unitSphere);

    //Calculate ideal number of triangles. If after delaunay triangulation we get a
    //different number of triangles we will print out some intermediate files and try to re-run.
    int idealTriCount, pentCount = 12, hexCount;
    hexCount = _numNodes - pentCount;
    idealTriCount = (6*hexCount + 5*pentCount)/3;

    d3D->SetInputConnection(idf->GetOutputPort());
    dssf->SetInputConnection(d3D->GetOutputPort());
    dssf->Update();
    final = dssf->GetOutput();    
    if( final->GetNumberOfPolys() != idealTriCount){
        cout<< "Bad Delaunay triangulation detected!" <<std::endl;
        cout<< "Printing the unit sphere generated." <<std::endl;        
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
*Compute the OPSbody energy and energy contribution due to some other elements
*inherited from the parent Body class _elements container
 */
void OPSBody::compute(bool f0, bool f1, bool f2) {

    double aM, aP, aN, aC;
    double E, re, s;
    double K, a, b;

    aM = _prop.alphaM; aP = _prop.alphaP; aN = _prop.alphaN; aC = _prop.alphaC;
    E = _prop.epsilon; re = _prop.r_e; s = _prop.s;
    K = _prop.K; a = _prop.a; b = _prop.b;

    // Initialize energy to be zero
    if (f0){
        _energy = 0.0;
        _morseEn = 0.0;
        _planarEn = 0.0;
        _normalEn = 0.0;
        _circularEn = 0.0;
    }
    if(f0 && !f1){
        for(int i=0; i < _numNodes; i++){
            Vector3D vi(0.0), xi(0.0);
            xi = _opsNodes[i]->deformedPosition();
            vi = _opsNodes[i]->deformedRotationVector();

            for(int j=0; j < _neighbors[i]->GetNumberOfIds(); j++){
                Vector3D vj(0.0), xj(0.0);
                double morseEn, KerXi, Phi_pXi, Phi_nXi, Phi_cXi;
                vtkIdType currId = _neighbors[i]->GetId(j);

                xj = _opsNodes[currId]->deformedPosition();
                vj = _opsNodes[currId]->deformedRotationVector();

                morseEn = morse(xi,xj);
                KerXi = psi( vi, xi, xj);
                Phi_pXi = phi_p( vi, xi, xj );
                Phi_nXi = phi_n( vi, vj );
                Phi_cXi = phi_c( vi, vj, xi, xj );
                _morseEn += aM*morseEn;
                _planarEn += aP*KerXi*Phi_pXi;
                _normalEn += aN*KerXi*Phi_nXi;
                _circularEn += aC*KerXi*Phi_cXi;
                _energy += aM*morseEn + ( aP*Phi_pXi + aN*Phi_nXi + aC*Phi_cXi )*KerXi;
            }
        }
    }
    else{
        for(int i=0; i < _numNodes; i++){
            Vector3D vi(0.0), xi(0.0);
            xi = _opsNodes[i]->deformedPosition();
            vi = _opsNodes[i]->deformedRotationVector();

            for(int j=0; j < _neighbors[i]->GetNumberOfIds(); j++){
                Vector3D vj(0.0), xj(0.0);
                double morseEn, KerXi, Phi_pXi, Phi_nXi, Phi_cXi;
                Vector3D dMorseXi, dMorseXj, dKerXi, dKerXj, dKerVi;
                Vector3D dPhi_pXi, dPhi_pXj, dPhi_pVi, dPhi_nVi, dPhi_nVj;
                Vector3D dPhi_cXi, dPhi_cXj, dPhi_cVi, dPhi_cVj;
                vtkIdType currId = _neighbors[i]->GetId(j);

                xj = _opsNodes[currId]->deformedPosition();
                vj = _opsNodes[currId]->deformedRotationVector();

                // Evaluate morse potential derivatives
                morseEn = morse(xi,xj);
                dMorseXi = DmorseDxi( xi, xj );
                dMorseXj = DmorseDxj( xi, xj );

                // Evaluate kernel and its derivatives
                KerXi = psi( vi, xi, xj );
                dKerXi = DpsiDxi( vi, xi, xj );
                dKerXj = DpsiDxj( vi, xi, xj );
                dKerVi = DpsiDvi( vi, xi, xj );

                // Evaluate co-planarity potential and its derivatives
                Phi_pXi = phi_p( vi, xi, xj );
                dPhi_pXi = Dphi_pDxi( vi, xi, xj );
                dPhi_pXj = Dphi_pDxj( vi, xi, xj );
                dPhi_pVi = Dphi_pDvi( vi, xi, xj );

                // Evaluate co-normality potential and its derivatives
                Phi_nXi = phi_n( vi, vj );
                dPhi_nVi = Dphi_nDvi( vi, vj );
                dPhi_nVj = Dphi_nDvj( vi, vj );

                // Evaluate co-circularity potential and its derivatives
                Phi_cXi = phi_c( vi, vj, xi, xj );
                dPhi_cXi = Dphi_cDxi( vi, vj, xi, xj );
                dPhi_cXj = Dphi_cDxj( vi, vj, xi, xj );
                dPhi_cVi = Dphi_cDvi( vi, vj, xi, xj );
                dPhi_cVj = Dphi_cDvj( vi, vj, xi, xj );

                Vector3D centerDx(0.0), centerDv(0.0), neighborDx(0.0), neighborDv(0.0);

                centerDx = aM*dMorseXi + aP*( dPhi_pXi*KerXi + Phi_pXi*dKerXi )
                        + aN*( Phi_nXi*dKerXi ) + aC*( dPhi_cXi*KerXi
                                                       + Phi_cXi*dKerXi );
                centerDv = aP*( dPhi_pVi*KerXi + Phi_pXi*dKerVi ) +
                        aN*( dPhi_nVi*KerXi + Phi_nXi*dKerVi ) +
                        aC*( dPhi_cVi*KerXi + Phi_cXi*dKerVi );

                neighborDx = aM*dMorseXj + aP*( dPhi_pXj*KerXi +
                                                Phi_pXi*dKerXj ) + aN*( Phi_nXi*dKerXj )
                        + aC*( dPhi_cXj*KerXi + Phi_cXi*dKerXj );
                neighborDv = KerXi*( aN*dPhi_nVj + aC*dPhi_cVj );

                Vector6D centerForce(centerDx[0], centerDx[1], centerDx[2],
                        centerDv[0], centerDv[1], centerDv[2]);
                Vector6D neighborForce(neighborDx[0], neighborDx[1],
                        neighborDx[2], neighborDv[0], neighborDv[1], neighborDv[2]);

                for(int k=0; k < 6; k++){
                    _opsNodes[i]->addForce( k, centerForce(k) );
                    _opsNodes[currId]->addForce( k, neighborForce(k) );
                }

                _morseEn += aM*morseEn;
                _planarEn += aP*KerXi*Phi_pXi;
                _normalEn += aN*KerXi*Phi_nXi;
                _circularEn += aC*KerXi*Phi_cXi;
                _energy += aM*morseEn + ( aP*Phi_pXi + aN*Phi_nXi + aC*Phi_cXi )*KerXi;
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

    writer->SetFileName( fileName.c_str() );
    writer->SetInputData( _polyData );
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
    double radius = 0.0;
    for(int i=0; i < _numNodes; i++){
        Vector3D x;
        x = _opsNodes[i]->deformedPosition();
        radius += tvmet::norm2(x);
    }
    radius /= _numNodes;
    return radius;
}

/*
 *  Calculate asphericity
 */
double OPSBody::getAsphericity(){
    double asphericity = 0.0;
    double R0 = getAverageRadius();
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
    case aM:
        _prop.alphaM = val;
        break;
    case aP:
        _prop.alphaP = val;
        break;
    case aN:
        _prop.alphaN = val;
        break;
    case aC:
        _prop.alphaC = val;
        break;
    case E:
        _prop.epsilon = val;
        break;
    case r:
        _prop.r_e = val;
        break;
    case sv:
        _prop.s = val;
        break;
    case Kv:
        _prop.K = val;
        break;
    case av:
        _prop.a = val;
        break;
    case bv:
        _prop.b = val;
        break;
    }
}

/*
*Morse function
 */
double OPSBody::morse(Vector3D xi, Vector3D xj){
    double E, re, s, r, potential;
    E = _prop.epsilon; re = _prop.r_e; s = _prop.s;
    r = tvmet::norm2(xi - xj);
    potential = E*( -2*exp( -s*(r-re) ) + exp( -2*s*(r-re) ) );
    return potential;
}

/*
*The kernel function
 */
double OPSBody::psi(Vector3D vi, Vector3D xi, Vector3D xj){
    double K, a, b, psi0;
    double u0, u1, u2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, sin_alpha_i, cos_alpha_i;
    double sin_alpha_i_2, cos_alpha_i_2;
    double term1, term2, term3;
    
    K = _prop.K; a = _prop.a; b = _prop.b;
    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vi_mag = sqrt(vi_mag_sqr);
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;

    term1 = cos_alpha_i_2*(-x2 + y2) + cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1) / vi_mag
            + cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0) / vi_mag	- sin_alpha_i_2*u0*u0*(-x2 + y2)/ vi_mag_sqr
            + sin_alpha_i_2*u0*u2*(-2*x0 + 2*y0) / vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/ vi_mag_sqr
            + sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1) / vi_mag_sqr + sin_alpha_i_2*u2*u2*(-x2 + y2)/ vi_mag_sqr;
    
    term2 = cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 + 2*y2)/ vi_mag
            - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)	/ vi_mag
            + sin_alpha_i_2*u0*u0 *(-x0 + y0) / vi_mag_sqr + sin_alpha_i_2*u0*u1 *(-2*x1 + 2*y1)/ vi_mag_sqr
            + sin_alpha_i_2*u0*u2 *(-2*x2 + 2*y2)/ vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0) / vi_mag_sqr
            - sin_alpha_i_2*u2*u2*(-x0 + y0) / vi_mag_sqr;

    term3 = cos_alpha_i_2*(-x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/ vi_mag
            + cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/ vi_mag
            - sin_alpha_i_2*u0*u0*(-x1 + y1)/ vi_mag_sqr + sin_alpha_i_2*u0*u1*(-2*x0 + 2*y0)/ vi_mag_sqr
            + sin_alpha_i_2*u1*u1*(-x1 + y1) / vi_mag_sqr + sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/ vi_mag_sqr
            - sin_alpha_i_2*u2*u2*(-x1 + y1)/ vi_mag_sqr;

    psi0 = K*exp(-(term1*term1) / (2*b*b) + (-(term2*term2) - (term3*term3))/ (2*a*a));
    return psi0;

}
/*
*The co-planarity function
 */
double OPSBody::phi_p(Vector3D vi, Vector3D xi, Vector3D xj){
    double u0, u1, u2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, sin_alpha_i, cos_alpha_i;
    double sin_alpha_i_2, cos_alpha_i_2, phi_p0;
    double term1;

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vi_mag = sqrt(vi_mag_sqr);
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;

    term1 = (-x0 + y0)*(2*cos_alpha_i*sin_alpha_i*u1 / vi_mag + 2*sin_alpha_i_2*u0*u2 / vi_mag_sqr)
            + (-x1 + y1)*(-2*cos_alpha_i*sin_alpha_i*u0 / vi_mag + 2*sin_alpha_i_2*u1*u2 / vi_mag_sqr)
            + (-x2 + y2)*(cos_alpha_i_2	- sin_alpha_i_2*u0*u0 / vi_mag_sqr
                          - sin_alpha_i_2*u1*u1 / vi_mag_sqr	+ sin_alpha_i_2*u2*u2 / vi_mag_sqr);
    
    phi_p0 = term1*term1;
    return phi_p0;
}
/*
*The co-normality function
 */
double OPSBody::phi_n(Vector3D vi, Vector3D vj){
    double u0, u1, u2, v0, v1, v2;
    double vi_mag_sqr, vi_mag, vj_mag_sqr, vj_mag;
    double sin_alpha_i, cos_alpha_i, sin_alpha_j, cos_alpha_j;
    double sin_alpha_i_2, cos_alpha_i_2, sin_alpha_j_2, cos_alpha_j_2, phi_n0;
    double term1, term2, term3;
    
    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    v0 = vj[0]; v1 = vj[1]; v2 = vj[2];

    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vj_mag_sqr = v0*v0 + v1*v1 + v2*v2;
    vi_mag = sqrt(vi_mag_sqr);
    vj_mag = sqrt(vj_mag_sqr);
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;
    sin_alpha_j = sin( 0.5*vj_mag );
    cos_alpha_j = cos( 0.5*vj_mag );
    sin_alpha_j_2 = sin_alpha_j*sin_alpha_j;
    cos_alpha_j_2 = cos_alpha_j*cos_alpha_j;

    term1 = -2*cos_alpha_i*sin_alpha_i*u0 / vi_mag
            + 2*cos_alpha_j*sin_alpha_j*v0 / vj_mag
            + 2*sin_alpha_i_2*u1*u2 / vi_mag_sqr
            - 2*sin_alpha_j_2*v1*v2 / vj_mag_sqr;
    
    term2 = 2*cos_alpha_i*sin_alpha_i*u1 / vi_mag
            - 2*cos_alpha_j*sin_alpha_j*v1 / vj_mag
            + 2*sin_alpha_i_2*u0*u2 / vi_mag_sqr
            - 2*sin_alpha_j_2*v0*v2 / vj_mag_sqr;
    
    term3 = cos_alpha_i_2 - cos_alpha_j_2
            - sin_alpha_i_2*u0*u0 / vi_mag_sqr
            - sin_alpha_i_2*u1*u1 / vi_mag_sqr
            + sin_alpha_i_2*u2*u2 / vi_mag_sqr
            + sin_alpha_j_2*v0*v0 / vj_mag_sqr
            + sin_alpha_j_2*v1*v1 / vj_mag_sqr
            - sin_alpha_j_2*v2*v2 / vj_mag_sqr;
    
    phi_n0 = term1*term1 + term2*term2 + term3*term3;
    return phi_n0;
}

/*
*The co-circularity function
 */
double OPSBody::phi_c(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj){
    double u0, u1, u2, v0, v1, v2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, vj_mag_sqr, vj_mag;
    double sin_alpha_i, cos_alpha_i, sin_alpha_j, cos_alpha_j;
    double sin_alpha_i_2, cos_alpha_i_2, sin_alpha_j_2, cos_alpha_j_2, phi_c0;
    double term1;

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    v0 = vj[0]; v1 = vj[1]; v2 = vj[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vj_mag_sqr = v0*v0 + v1*v1 + v2*v2;
    vi_mag = sqrt(vi_mag_sqr);
    vj_mag = sqrt(vj_mag_sqr);
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;
    sin_alpha_j = sin( 0.5*vj_mag );
    cos_alpha_j = cos( 0.5*vj_mag );
    sin_alpha_j_2 = sin_alpha_j*sin_alpha_j;
    cos_alpha_j_2 = cos_alpha_j*cos_alpha_j;
    
    term1 = (-x0 + y0)*(2*cos_alpha_i*sin_alpha_i*u1 / vi_mag
                        + 2*cos_alpha_j*sin_alpha_j*v1 / vj_mag
                        + 2*sin_alpha_i_2*u0*u2 / vi_mag_sqr
                        + 2*sin_alpha_j_2*v0*v2 / vj_mag_sqr)
            + (-x1 + y1)*(-2*cos_alpha_i*sin_alpha_i*u0 / vi_mag
                          - 2*cos_alpha_j*sin_alpha_j*v0/ vj_mag
                          + 2*sin_alpha_i_2*u1*u2 / vi_mag_sqr
                          + 2*sin_alpha_j_2*v1*v2 / vj_mag_sqr)
            + (-x2 + y2)*(cos_alpha_i_2 + cos_alpha_j_2
                          - sin_alpha_i_2*u0*u0 / vi_mag_sqr
                          - sin_alpha_i_2*u1*u1 / vi_mag_sqr
                          + sin_alpha_i_2*u2*u2 / vi_mag_sqr
                          - sin_alpha_j_2*v0*v0 / vj_mag_sqr
                          - sin_alpha_j_2*v1*v1 / vj_mag_sqr
                          + sin_alpha_j_2*v2*v2 / vj_mag_sqr);

    phi_c0 = term1*term1;

    return phi_c0;
}

/*
*Derivative of Morse with respect to xi
 */
Vector3D OPSBody::DmorseDxi(Vector3D xi, Vector3D xj){
    double epsilon, l0, a, rij, DmorseDxi0, DmorseDxi1, DmorseDxi2;
    double x0, x1, x2, y0, y1, y2;

    epsilon = _prop.epsilon; l0 = _prop.r_e; a = _prop.s;
    rij = tvmet::norm2( xi-xj );

    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    DmorseDxi0 = epsilon
            *(2*a*(x0 - y0)*exp(-a*(-l0 + rij)) / rij
              - 2*a*(x0 - y0)*exp(-2*a*(-l0 + rij)) / rij);

    DmorseDxi1 = epsilon
            *(2*a*(x1 - y1)*exp(-a*(-l0 + rij)) / rij
              - 2*a*(x1 - y1)*exp(-2*a*(-l0 + rij)) / rij);

    DmorseDxi2 = epsilon
            *(2*a*(x2 - y2)*exp(-a*(-l0 + rij)) / rij
              - 2*a*(x2 - y2)*exp(-2*a*(-l0 + rij)) / rij);

    Vector3D DmorseDxi( DmorseDxi0, DmorseDxi1, DmorseDxi2 );

    return DmorseDxi;
}

/*
*Derivative of Morse with respect to xj
 */
Vector3D OPSBody::DmorseDxj(Vector3D xi, Vector3D xj){
    double epsilon, l0, a, rij, DmorseDxj0, DmorseDxj1, DmorseDxj2;
    double x0, x1, x2, y0, y1, y2;

    epsilon = _prop.epsilon; l0 = _prop.r_e; a = _prop.s;
    rij = tvmet::norm2( xi-xj );

    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    DmorseDxj0 = epsilon
            *(2*a*(-x0 + y0)*exp(-a*(-l0 + rij)) / rij
              - 2*a*(-x0 + y0)*exp(-2*a*(-l0 + rij)) / rij);

    DmorseDxj1 = epsilon
            *(2*a*(-x1 + y1)*exp(-a*(-l0 + rij)) / rij
              - 2*a*(-x1 + y1)*exp(-2*a*(-l0 + rij)) / rij);

    DmorseDxj2 = epsilon
            *(2*a*(-x2 + y2)*exp(-a*(-l0 + rij)) / rij
              - 2*a*(-x2 + y2)*exp(-2*a*(-l0 + rij)) / rij);

    Vector3D DmorseDxj( DmorseDxj0, DmorseDxj1, DmorseDxj2 );

    return DmorseDxj;
}

/*
*Derivative of kernel psi wrt vi
 */
Vector3D OPSBody::DpsiDvi(Vector3D vi, Vector3D xi, Vector3D xj){
    double K, a, b, psi0;
    double u0, u1, u2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, sin_alpha_i, cos_alpha_i;
    double sin_alpha_i_2, cos_alpha_i_2;
    double DpsiDu0, DpsiDu1, DpsiDu2;
    double vi_mag_3, vi_mag_4;
    double u0_3, u1_3, u2_3;
    double t1,t2,t3,t4,t5,t6,t7,t8,t9;

    K = _prop.K; a = _prop.a; b = _prop.b;
    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];
    
    u0_3 = u0*u0*u0;
    u1_3 = u1*u1*u1;
    u2_3 = u2*u2*u2;
    
    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vi_mag = sqrt(vi_mag_sqr);
    vi_mag_3 = vi_mag*vi_mag*vi_mag;
    vi_mag_4 = vi_mag*vi_mag*vi_mag*vi_mag;
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;
    
    t1 = cos_alpha_i_2*(-x2 + y2) +
            cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
            cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
            sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr;
    t2 = cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 +
                                                                2*y2)/vi_mag - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
            sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                            2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr;
    t3 = cos_alpha_i_2*(-x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
            cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x1 +
                                 y1)/vi_mag_sqr;

    DpsiDu0 = K*(-(cos_alpha_i_2*(-x2 + y2) +
                   cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
                   cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
                   sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                                    2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
                   sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                   sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr)*(cos_alpha_i_2*u0*u0*(-
                                                                                    2*x1 + 2*y1)/vi_mag_sqr + cos_alpha_i_2*u0*u1*(2*x0 -
                                                                                                                                   2*y0)/vi_mag_sqr - 2*cos_alpha_i*sin_alpha_i*u0_3*(-x2 +
                                                                                                                                                                                      y2)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u0*u2*(-2*x0 +
                                                                                                                                                                                                                                           2*y0)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*u0*(-2*x1 +
                                                                                                                                                                                                                                                                                               2*y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*u1*u1*(-x2 +
                                                                                                                                                                                                                                                                                                                                                      y2)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u1*u2*(-2*x1 +
                                                                                                                                                                                                                                                                                                                                                                                                           2*y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*u1*(2*x0 -
                                                                                                                                                                                                                                                                                                                                                                                                                                                               2*y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u2*u2*(-x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*(-x2 + y2)/vi_mag +
                                                               2*cos_alpha_i*sin_alpha_i*(-2*x1 + 2*y1)/vi_mag +
                                                               4*sin_alpha_i_2*u0_3*(-x2 + y2)/(vi_mag_4) -
                                                               4*sin_alpha_i_2*u0*u0*u2*(-2*x0 + 2*y0)/(vi_mag_4) -
                                                               sin_alpha_i_2*u0*u0*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                               4*sin_alpha_i_2*u0*u1*u1*(-x2 + y2)/(vi_mag_4) -
                                                               4*sin_alpha_i_2*u0*u1*u2*(-2*x1 + 2*y1)/(vi_mag_4) -
                                                               sin_alpha_i_2*u0*u1*(2*x0 - 2*y0)/vi_mag_sqr -
                                                               4*sin_alpha_i_2*u0*u2*u2*(-x2 + y2)/(vi_mag_4) -
                                                               4*sin_alpha_i_2*u0*(-x2 + y2)/vi_mag_sqr + 2*sin_alpha_i_2*u2*(-
                                                                                                                              2*x0 + 2*y0)/vi_mag_sqr)/(2*b*b) + (-(cos_alpha_i_2*(-x0 + y0) +
                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u1*(-2*x2 + 2*y2)/vi_mag -
                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
                                                                                                                                                                    sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                     2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                    2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
                                                                                                                                                                    sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr)*(cos_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                     2*x2 + 2*y2)/vi_mag_sqr - cos_alpha_i_2*u0*u2*(-2*x1 +
                                                                                                                                                                                                                                                                                    2*y1)/vi_mag_sqr + 2*cos_alpha_i*sin_alpha_i*u0_3*(-x0 +
                                                                                                                                                                                                                                                                                                                                       y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u0*u1*(-2*x1 +
                                                                                                                                                                                                                                                                                                                                                                                            2*y1)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                   2*y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*u1*u1*(-x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          y0)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*u1*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            2*y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*u2*u2*(-x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u2*(-2*x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     2*y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*(-x0 + y0)/vi_mag -
                                                                                                                                                                                                                4*sin_alpha_i_2*u0_3*(-x0 + y0)/(vi_mag_4) -
                                                                                                                                                                                                                4*sin_alpha_i_2*u0*u0*u1*(-2*x1 + 2*y1)/(vi_mag_4) -
                                                                                                                                                                                                                4*sin_alpha_i_2*u0*u0*u2*(-2*x2 + 2*y2)/(vi_mag_4) +
                                                                                                                                                                                                                4*sin_alpha_i_2*u0*u1*u1*(-x0 + y0)/(vi_mag_4) -
                                                                                                                                                                                                                sin_alpha_i_2*u0*u1*(-2*x2 + 2*y2)/vi_mag_sqr +
                                                                                                                                                                                                                4*sin_alpha_i_2*u0*u2*u2*(-x0 + y0)/(vi_mag_4) +
                                                                                                                                                                                                                sin_alpha_i_2*u0*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                                                                                                                                                                                4*sin_alpha_i_2*u0*(-x0 + y0)/vi_mag_sqr + 2*sin_alpha_i_2*u1*(-
                                                                                                                                                                                                                                                                               2*x1 + 2*y1)/vi_mag_sqr + 2*sin_alpha_i_2*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                             2*y2)/vi_mag_sqr) - (cos_alpha_i_2*(-x1 + y1) -
                                                                                                                                                                                                                                                                                                                                                  cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
                                                                                                                                                                                                                                                                                                                                                  cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
                                                                                                                                                                                                                                                                                                                                                  sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                                                                                                                                                                                                   2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
                                                                                                                                                                                                                                                                                                                                                  sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                                                                                                                                                                                                                                                  sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr)*(-cos_alpha_i_2*u0*u0*(-
                                                                                                                                                                                                                                                                                                                                                                                                                    2*x2 + 2*y2)/vi_mag_sqr + cos_alpha_i_2*u0*u2*(-2*x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                   2*y0)/vi_mag_sqr - 2*cos_alpha_i*sin_alpha_i*u0_3*(-x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      y1)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u0*u1*(-2*x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           2*y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u0*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               2*y2)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u1*u1*(-x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      y1)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u1*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           2*y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*u2*u2*(-x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*u2*(-2*x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    2*y0)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*(-x1 + y1)/vi_mag -
                                                                                                                                                                                                                                                                                                                                                                                              2*cos_alpha_i*sin_alpha_i*(-2*x2 + 2*y2)/vi_mag +
                                                                                                                                                                                                                                                                                                                                                                                              4*sin_alpha_i_2*u0_3*(-x1 + y1)/(vi_mag_4) -
                                                                                                                                                                                                                                                                                                                                                                                              4*sin_alpha_i_2*u0*u0*u1*(-2*x0 + 2*y0)/(vi_mag_4) +
                                                                                                                                                                                                                                                                                                                                                                                              sin_alpha_i_2*u0*u0*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                                                                                                                                                                                                                                                                                              4*sin_alpha_i_2*u0*u1*u1*(-x1 + y1)/(vi_mag_4) -
                                                                                                                                                                                                                                                                                                                                                                                              4*sin_alpha_i_2*u0*u1*u2*(-2*x2 + 2*y2)/(vi_mag_4) +
                                                                                                                                                                                                                                                                                                                                                                                              4*sin_alpha_i_2*u0*u2*u2*(-x1 + y1)/(vi_mag_4) -
                                                                                                                                                                                                                                                                                                                                                                                              sin_alpha_i_2*u0*u2*(-2*x0 + 2*y0)/vi_mag_sqr -
                                                                                                                                                                                                                                                                                                                                                                                              4*sin_alpha_i_2*u0*(-x1 + y1)/vi_mag_sqr + 2*sin_alpha_i_2*u1*(-
                                                                                                                                                                                                                                                                                                                                                                                                                                                             2*x0 + 2*y0)/vi_mag_sqr))/(2*a*a))*exp(-(t1*t1)/(2*b*b) + (-(t2*t2)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        - (t3*t3))/(2*a*a));

    t4 = cos_alpha_i_2*(-x2 + y2) +
            cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
            cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
            sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr;
    t5 = cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 +
                                                                2*y2)/vi_mag - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
            sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                            2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr;
    t6 = cos_alpha_i_2*(-x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
            cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x1 +	y1)/vi_mag_sqr;

    DpsiDu1 = K*(-(cos_alpha_i_2*(-
                                  x2 + y2) + cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
                   cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
                   sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                                    2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
                   sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                   sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr)*(cos_alpha_i_2*u0*u1*(-
                                                                                    2*x1 + 2*y1)/vi_mag_sqr + cos_alpha_i_2*u1*u1*(2*x0 -
                                                                                                                                   2*y0)/vi_mag_sqr - 2*cos_alpha_i*sin_alpha_i*u0*u0*u1*(-x2 +
                                                                                                                                                                                          y2)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u1*u2*(-2*x0 +
                                                                                                                                                                                                                                               2*y0)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*u1*(-2*x1 +
                                                                                                                                                                                                                                                                                                   2*y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1_3*(-x2 +
                                                                                                                                                                                                                                                                                                                                                      y2)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u1*u1*u2*(-2*x1 +
                                                                                                                                                                                                                                                                                                                                                                                                           2*y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*u1*(2*x0 -
                                                                                                                                                                                                                                                                                                                                                                                                                                                               2*y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u1*u2*u2*(-x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*(-x2 + y2)/vi_mag +
                                                               2*cos_alpha_i*sin_alpha_i*(2*x0 - 2*y0)/vi_mag +
                                                               4*sin_alpha_i_2*u0*u0*u1*(-x2 + y2)/(vi_mag_4) -
                                                               4*sin_alpha_i_2*u0*u1*u2*(-2*x0 + 2*y0)/(vi_mag_4) -
                                                               sin_alpha_i_2*u0*u1*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                               4*sin_alpha_i_2*u1_3*(-x2 + y2)/(vi_mag_4) -
                                                               4*sin_alpha_i_2*u1*u1*u2*(-2*x1 + 2*y1)/(vi_mag_4) -
                                                               sin_alpha_i_2*u1*u1*(2*x0 - 2*y0)/vi_mag_sqr -
                                                               4*sin_alpha_i_2*u1*u2*u2*(-x2 + y2)/(vi_mag_4) -
                                                               4*sin_alpha_i_2*u1*(-x2 + y2)/vi_mag_sqr + 2*sin_alpha_i_2*u2*(-
                                                                                                                              2*x1 + 2*y1)/vi_mag_sqr)/(2*b*b) + (-(cos_alpha_i_2*(-x0 + y0) +
                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u1*(-2*x2 + 2*y2)/vi_mag -
                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
                                                                                                                                                                    sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                     2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                    2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
                                                                                                                                                                    sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr)*(cos_alpha_i_2*u1*u1*(-
                                                                                                                                                                                                                                     2*x2 + 2*y2)/vi_mag_sqr - cos_alpha_i_2*u1*u2*(-2*x1 +
                                                                                                                                                                                                                                                                                    2*y1)/vi_mag_sqr + 2*cos_alpha_i*sin_alpha_i*u0*u0*u1*(-x0 +
                                                                                                                                                                                                                                                                                                                                           y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u1*u1*(-2*x1 +
                                                                                                                                                                                                                                                                                                                                                                                                2*y1)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u1*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                       2*y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1_3*(-x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          y0)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*u1*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            2*y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*u2*u2*(-x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u1*u2*(-2*x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     2*y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*(-x0 + y0)/vi_mag +
                                                                                                                                                                                                                2*cos_alpha_i*sin_alpha_i*(-2*x2 + 2*y2)/vi_mag -
                                                                                                                                                                                                                4*sin_alpha_i_2*u0*u0*u1*(-x0 + y0)/(vi_mag_4) -
                                                                                                                                                                                                                4*sin_alpha_i_2*u0*u1*u1*(-2*x1 + 2*y1)/(vi_mag_4) -
                                                                                                                                                                                                                4*sin_alpha_i_2*u0*u1*u2*(-2*x2 + 2*y2)/(vi_mag_4) +
                                                                                                                                                                                                                2*sin_alpha_i_2*u0*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                                                                                                                                                                                4*sin_alpha_i_2*u1_3*(-x0 + y0)/(vi_mag_4) - sin_alpha_i_2*u1*u1*(-
                                                                                                                                                                                                                                                                                  2*x2 + 2*y2)/vi_mag_sqr + 4*sin_alpha_i_2*u1*u2*u2*(-x0 +
                                                                                                                                                                                                                                                                                                                                      y0)/(vi_mag_4) + sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr -
                                                                                                                                                                                                                4*sin_alpha_i_2*u1*(-x0 + y0)/vi_mag_sqr) - (cos_alpha_i_2*(-x1 +
                                                                                                                                                                                                                                                                            y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
                                                                                                                                                                                                                                                             cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
                                                                                                                                                                                                                                                             sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                                                                                                              2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
                                                                                                                                                                                                                                                             sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                                                                                                                                                             sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr)*(-cos_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                                                                                                               2*x2 + 2*y2)/vi_mag_sqr + cos_alpha_i_2*u1*u2*(-2*x0 +
                                                                                                                                                                                                                                                                                                                                                                              2*y0)/vi_mag_sqr - 2*cos_alpha_i*sin_alpha_i*u0*u0*u1*(-x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                     y1)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u1*u1*(-2*x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          2*y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u1*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              2*y2)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u1_3*(-x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 y1)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u1*u1*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      2*y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*u2*u2*(-x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*u2*(-2*x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               2*y0)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*(-x1 + y1)/vi_mag +
                                                                                                                                                                                                                                                                                                         4*sin_alpha_i_2*u0*u0*u1*(-x1 + y1)/(vi_mag_4) -
                                                                                                                                                                                                                                                                                                         4*sin_alpha_i_2*u0*u1*u1*(-2*x0 + 2*y0)/(vi_mag_4) +
                                                                                                                                                                                                                                                                                                         sin_alpha_i_2*u0*u1*(-2*x2 + 2*y2)/vi_mag_sqr +
                                                                                                                                                                                                                                                                                                         2*sin_alpha_i_2*u0*(-2*x0 + 2*y0)/vi_mag_sqr -
                                                                                                                                                                                                                                                                                                         4*sin_alpha_i_2*u1_3*(-x1 + y1)/(vi_mag_4) -
                                                                                                                                                                                                                                                                                                         4*sin_alpha_i_2*u1*u1*u2*(-2*x2 + 2*y2)/(vi_mag_4) +
                                                                                                                                                                                                                                                                                                         4*sin_alpha_i_2*u1*u2*u2*(-x1 + y1)/(vi_mag_4) -
                                                                                                                                                                                                                                                                                                         sin_alpha_i_2*u1*u2*(-2*x0 + 2*y0)/vi_mag_sqr +
                                                                                                                                                                                                                                                                                                         4*sin_alpha_i_2*u1*(-x1 + y1)/vi_mag_sqr + 2*sin_alpha_i_2*u2*(-
                                                                                                                                                                                                                                                                                                                                                                        2*x2 + 2*y2)/vi_mag_sqr))/(2*a*a))*exp(-(t4*t4)/(2*b*b) + (-
                                                                                                                                                                                                                                                                                                                                                                                                                                   (t5*t5) - (t6*t6))/(2*a*a));
    
    t7 = cos_alpha_i_2*(-x2 + y2) +
            cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
            cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
            sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr;
    t8 = cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 +
                                                                2*y2)/vi_mag - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
            sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                            2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr;
    t9 = cos_alpha_i_2*(-x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
            cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr;

    DpsiDu2 = K*(-(cos_alpha_i_2*(-
                                  x2 + y2) + cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
                   cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
                   sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                                    2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
                   sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                   sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr)*(cos_alpha_i_2*u0*u2*(-
                                                                                    2*x1 + 2*y1)/vi_mag_sqr + cos_alpha_i_2*u1*u2*(2*x0 -
                                                                                                                                   2*y0)/vi_mag_sqr - 2*cos_alpha_i*sin_alpha_i*u0*u0*u2*(-x2 +
                                                                                                                                                                                          y2)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u2*u2*(-2*x0 +
                                                                                                                                                                                                                                               2*y0)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u0*u2*(-2*x1 +
                                                                                                                                                                                                                                                                                                   2*y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*u1*u2*(-x2 +
                                                                                                                                                                                                                                                                                                                                                          y2)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u1*u2*u2*(-2*x1 +
                                                                                                                                                                                                                                                                                                                                                                                                               2*y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*u2*(2*x0 -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                   2*y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u2_3*(-x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u2*(-x2 + y2)/vi_mag +
                                                               4*sin_alpha_i_2*u0*u0*u2*(-x2 + y2)/(vi_mag_4) -
                                                               4*sin_alpha_i_2*u0*u2*u2*(-2*x0 + 2*y0)/(vi_mag_4) -
                                                               sin_alpha_i_2*u0*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                               2*sin_alpha_i_2*u0*(-2*x0 + 2*y0)/vi_mag_sqr +
                                                               4*sin_alpha_i_2*u1*u1*u2*(-x2 + y2)/(vi_mag_4) -
                                                               4*sin_alpha_i_2*u1*u2*u2*(-2*x1 + 2*y1)/(vi_mag_4) -
                                                               sin_alpha_i_2*u1*u2*(2*x0 - 2*y0)/vi_mag_sqr + 2*sin_alpha_i_2*u1*(-
                                                                                                                                  2*x1 + 2*y1)/vi_mag_sqr - 4*sin_alpha_i_2*u2_3*(-x2 + y2)/(vi_mag_4) +
                                                               4*sin_alpha_i_2*u2*(-x2 + y2)/vi_mag_sqr)/(2*b*b) + (-
                                                                                                                    (cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 +
                                                                                                                                                                            2*y2)/vi_mag - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
                                                                                                                     sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                      2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                     2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
                                                                                                                     sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr)*(cos_alpha_i_2*u1*u2*(-
                                                                                                                                                                                      2*x2 + 2*y2)/vi_mag_sqr - cos_alpha_i_2*u2*u2*(-2*x1 +
                                                                                                                                                                                                                                     2*y1)/vi_mag_sqr + 2*cos_alpha_i*sin_alpha_i*u0*u0*u2*(-x0 +
                                                                                                                                                                                                                                                                                            y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u1*u2*(-2*x1 +
                                                                                                                                                                                                                                                                                                                                                 2*y1)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u2*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                        2*y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*u1*u2*(-x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                               y0)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u1*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 2*y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u2_3*(-x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u2*u2*(-2*x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      2*y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u2*(-x0 + y0)/vi_mag -
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*(-2*x1 + 2*y1)/vi_mag -
                                                                                                                                                                 4*sin_alpha_i_2*u0*u0*u2*(-x0 + y0)/(vi_mag_4) -
                                                                                                                                                                 4*sin_alpha_i_2*u0*u1*u2*(-2*x1 + 2*y1)/(vi_mag_4) -
                                                                                                                                                                 4*sin_alpha_i_2*u0*u2*u2*(-2*x2 + 2*y2)/(vi_mag_4) +
                                                                                                                                                                 2*sin_alpha_i_2*u0*(-2*x2 + 2*y2)/vi_mag_sqr +
                                                                                                                                                                 4*sin_alpha_i_2*u1*u1*u2*(-x0 + y0)/(vi_mag_4) -
                                                                                                                                                                 sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr +
                                                                                                                                                                 4*sin_alpha_i_2*u2_3*(-x0 + y0)/(vi_mag_4) + sin_alpha_i_2*u2*u2*(-
                                                                                                                                                                                                                                   2*x1 + 2*y1)/vi_mag_sqr - 4*sin_alpha_i_2*u2*(-x0 + y0)/vi_mag_sqr) -
                                                                                                                    (cos_alpha_i_2*(-x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 +
                                                                                                                                                                            2*y2)/vi_mag + cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
                                                                                                                     sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                      2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
                                                                                                                     sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                     sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr)*(-cos_alpha_i_2*u0*u2*(-
                                                                                                                                                                                       2*x2 + 2*y2)/vi_mag_sqr + cos_alpha_i_2*u2*u2*(-2*x0 +
                                                                                                                                                                                                                                      2*y0)/vi_mag_sqr - 2*cos_alpha_i*sin_alpha_i*u0*u0*u2*(-x1 +
                                                                                                                                                                                                                                                                                             y1)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u1*u2*(-2*x0 +
                                                                                                                                                                                                                                                                                                                                                  2*y0)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                      2*y2)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u1*u1*u2*(-x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                             y1)/(vi_mag_3) + 2*cos_alpha_i*sin_alpha_i*u1*u2*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  2*y2)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u2_3*(-x1 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     y1)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u2*u2*(-2*x0 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       2*y0)/(vi_mag_3) - 2*cos_alpha_i*sin_alpha_i*u2*(-x1 + y1)/vi_mag +
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*(-2*x0 + 2*y0)/vi_mag +
                                                                                                                                                                 4*sin_alpha_i_2*u0*u0*u2*(-x1 + y1)/(vi_mag_4) -
                                                                                                                                                                 4*sin_alpha_i_2*u0*u1*u2*(-2*x0 + 2*y0)/(vi_mag_4) +
                                                                                                                                                                 sin_alpha_i_2*u0*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                                                                 4*sin_alpha_i_2*u1*u1*u2*(-x1 + y1)/(vi_mag_4) -
                                                                                                                                                                 4*sin_alpha_i_2*u1*u2*u2*(-2*x2 + 2*y2)/(vi_mag_4) +
                                                                                                                                                                 2*sin_alpha_i_2*u1*(-2*x2 + 2*y2)/vi_mag_sqr +
                                                                                                                                                                 4*sin_alpha_i_2*u2_3*(-x1 + y1)/(vi_mag_4) - sin_alpha_i_2*u2*u2*(-
                                                                                                                                                                                                                                   2*x0 + 2*y0)/vi_mag_sqr - 4*sin_alpha_i_2*u2*(-x1 +
                                                                                                                                                                                                                                                                                 y1)/vi_mag_sqr))/(2*a*a))*exp(-(t7*t7)/(2*b*b) + (-
                                                                                                                                                                                                                                                                                                                                   (t8*t8) - (t9*t9))/(2*a*a));

    Vector3D DpsiDvi( DpsiDu0, DpsiDu1, DpsiDu2 );
    return DpsiDvi;
}

/*
*Derivative psi xi
 */
Vector3D OPSBody::DpsiDxi(Vector3D vi, Vector3D xi, Vector3D xj){
    double K, a, b, psi0;
    double u0, u1, u2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, sin_alpha_i, cos_alpha_i;
    double sin_alpha_i_2, cos_alpha_i_2;
    double DpsiDx0, DpsiDx1, DpsiDx2;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9;

    K = _prop.K; a = _prop.a; b = _prop.b;
    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vi_mag = sqrt(vi_mag_sqr);
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;

    t1 = cos_alpha_i_2*(-x2 + y2) + cos_alpha_i*sin_alpha_i*u0*(-2*x1 +
                                                                2*y1)/vi_mag + cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
            sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr;
    t2 = cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 +
                                                                2*y2)/vi_mag - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
            sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                            2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr;
    t3 = cos_alpha_i_2*(-
                        x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
            cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x1 +
                                 y1)/vi_mag_sqr;
    
    DpsiDx0 = K*(-
                 (4*cos_alpha_i*sin_alpha_i*u1/vi_mag -
                  4*sin_alpha_i_2*u0*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x2 + y2) +
                                                     cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
                                                     cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
                                                     sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                                                                      2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
                                                     sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                     sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr)/(2*b*b) + (-(-
                                                                                                             4*cos_alpha_i*sin_alpha_i*u2/vi_mag -
                                                                                                             4*sin_alpha_i_2*u0*u1/vi_mag_sqr)*(cos_alpha_i_2*(-x1 + y1) -
                                                                                                                                                cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
                                                                                                                                                cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
                                                                                                                                                sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                 2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
                                                                                                                                                sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                                                sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr) - (-2*cos_alpha_i_2 -
                                                                                                                                                                                              2*sin_alpha_i_2*u0*u0/vi_mag_sqr + 2*sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                                                                                              2*sin_alpha_i_2*u2*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x0 + y0) +
                                                                                                                                                                                                                                 cos_alpha_i*sin_alpha_i*u1*(-2*x2 + 2*y2)/vi_mag -
                                                                                                                                                                                                                                 cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
                                                                                                                                                                                                                                 sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                                                                                  2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                 2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
                                                                                                                                                                                                                                 sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr))/(2*a*a))*exp(-
                                                                                                                                                                                                                                                                                          (t1*t1)/(2*b*b) + (-(t2*t2) - (t3*t3))/(2*a*a));
    
    t4 = cos_alpha_i_2*(-x2 + y2) + cos_alpha_i*sin_alpha_i*u0*(-2*x1 +
                                                                2*y1)/vi_mag + cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
            sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr;
    t5 = cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 +
                                                                2*y2)/vi_mag - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
            sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                            2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr;
    t6 = cos_alpha_i_2*(-
                        x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
            cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x1 +
                                 y1)/vi_mag_sqr;

    DpsiDx1 = K*(-(-
                   4*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                   4*sin_alpha_i_2*u1*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x2 + y2) +
                                                      cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
                                                      cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
                                                      sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                                                                       2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
                                                      sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                      sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr)/(2*b*b) + (-
                                                                                                            (4*cos_alpha_i*sin_alpha_i*u2/vi_mag -
                                                                                                             4*sin_alpha_i_2*u0*u1/vi_mag_sqr)*(cos_alpha_i_2*(-x0 + y0) +
                                                                                                                                                cos_alpha_i*sin_alpha_i*u1*(-2*x2 + 2*y2)/vi_mag -
                                                                                                                                                cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
                                                                                                                                                sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                 2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                                                2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
                                                                                                                                                sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr) - (-2*cos_alpha_i_2 +
                                                                                                                                                                                              2*sin_alpha_i_2*u0*u0/vi_mag_sqr - 2*sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                                                                                              2*sin_alpha_i_2*u2*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x1 + y1) -
                                                                                                                                                                                                                                 cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
                                                                                                                                                                                                                                 cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
                                                                                                                                                                                                                                 sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                                                                                  2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
                                                                                                                                                                                                                                 sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                                                                                                                                 sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr))/(2*a*a))*exp(-
                                                                                                                                                                                                                                                                                          (t4*t4)/(2*b*b) + (-(t5*t5) - (t6*t6))/(2*a*a));
    
    t7 = cos_alpha_i_2*(-x2 + y2) + cos_alpha_i*sin_alpha_i*u0*(-2*x1 +
                                                                2*y1)/vi_mag + cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
            sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr;
    t8 = cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 +
                                                                2*y2)/vi_mag - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
            sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                            2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr;
    t9 = cos_alpha_i_2*(-x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
            cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr;

    DpsiDx2 = K*(-(-2*cos_alpha_i_2 +
                   2*sin_alpha_i_2*u0*u0/vi_mag_sqr + 2*sin_alpha_i_2*u1*u1/vi_mag_sqr -
                   2*sin_alpha_i_2*u2*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x2 + y2) +
                                                      cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
                                                      cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
                                                      sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                                                                       2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
                                                      sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                      sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr)/(2*b*b) + (-
                                                                                                            (4*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                                                             4*sin_alpha_i_2*u1*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x1 + y1) -
                                                                                                                                                cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
                                                                                                                                                cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
                                                                                                                                                sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                 2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
                                                                                                                                                sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                                                sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr) - (-
                                                                                                                                                                                              4*cos_alpha_i*sin_alpha_i*u1/vi_mag -
                                                                                                                                                                                              4*sin_alpha_i_2*u0*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x0 + y0) +
                                                                                                                                                                                                                                 cos_alpha_i*sin_alpha_i*u1*(-2*x2 + 2*y2)/vi_mag -
                                                                                                                                                                                                                                 cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
                                                                                                                                                                                                                                 sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                                                                                  2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                 2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
                                                                                                                                                                                                                                 sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr))/(2*a*a))*exp(-
                                                                                                                                                                                                                                                                                          (t7*t7)/(2*b*b) + (-(t8*t8) - (t9*t9))/(2*a*a));

    Vector3D DpsiDxi( DpsiDx0, DpsiDx1, DpsiDx2 );
    return DpsiDxi;
}

/*
*Derivative of psi wrt xj
 */
Vector3D OPSBody::DpsiDxj(Vector3D vi, Vector3D xi, Vector3D xj){
    double K, a, b, psi0;
    double u0, u1, u2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, sin_alpha_i, cos_alpha_i;
    double sin_alpha_i_2, cos_alpha_i_2;
    double DpsiDy0, DpsiDy1, DpsiDy2;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9;

    K = _prop.K; a = _prop.a; b = _prop.b;
    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vi_mag = sqrt(vi_mag_sqr);
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;

    t1 = cos_alpha_i_2*(-x2 + y2) + cos_alpha_i*sin_alpha_i*u0*(-2*x1 +
                                                                2*y1)/vi_mag + cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
            sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr;
    t2 = cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 +
                                                                2*y2)/vi_mag - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
            sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                            2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr;
    t3 = cos_alpha_i_2*(-x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
            cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x1 +	y1)/vi_mag_sqr;
    DpsiDy0 = K*(-(-
                   4*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                   4*sin_alpha_i_2*u0*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x2 + y2) +
                                                      cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
                                                      cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
                                                      sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                                                                       2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
                                                      sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                      sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr)/(2*b*b) + (-
                                                                                                            (4*cos_alpha_i*sin_alpha_i*u2/vi_mag +
                                                                                                             4*sin_alpha_i_2*u0*u1/vi_mag_sqr)*(cos_alpha_i_2*(-x1 + y1) -
                                                                                                                                                cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
                                                                                                                                                cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
                                                                                                                                                sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                 2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
                                                                                                                                                sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                                                sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr) - (2*cos_alpha_i_2 +
                                                                                                                                                                                              2*sin_alpha_i_2*u0*u0/vi_mag_sqr - 2*sin_alpha_i_2*u1*u1/vi_mag_sqr -
                                                                                                                                                                                              2*sin_alpha_i_2*u2*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x0 + y0) +
                                                                                                                                                                                                                                 cos_alpha_i*sin_alpha_i*u1*(-2*x2 + 2*y2)/vi_mag -
                                                                                                                                                                                                                                 cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
                                                                                                                                                                                                                                 sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                                                                                  2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                                                                                                                                 2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
                                                                                                                                                                                                                                 sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr))/(2*a*a))*exp(-
                                                                                                                                                                                                                                                                                          (t1*t1)/(2*b*b) + (-(t2*t2) - (t3*t3))/(2*a*a));
    
    t4 = cos_alpha_i_2*(-x2 + y2) + cos_alpha_i*sin_alpha_i*u0*(-2*x1 +
                                                                2*y1)/vi_mag + cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
            sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr;
    t5 = cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 +
                                                                2*y2)/vi_mag - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
            sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                            2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr;
    t6 = cos_alpha_i_2*(-x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
            cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x1 +	y1)/vi_mag_sqr;
    DpsiDy1 = K*(-
                 (4*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                  4*sin_alpha_i_2*u1*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x2 + y2) +
                                                     cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
                                                     cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
                                                     sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                                                                      2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
                                                     sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                     sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr)/(2*b*b) + (-(-
                                                                                                             4*cos_alpha_i*sin_alpha_i*u2/vi_mag +
                                                                                                             4*sin_alpha_i_2*u0*u1/vi_mag_sqr)*(cos_alpha_i_2*(-x0 + y0) +
                                                                                                                                                cos_alpha_i*sin_alpha_i*u1*(-2*x2 + 2*y2)/vi_mag -
                                                                                                                                                cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
                                                                                                                                                sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                 2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                                                2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
                                                                                                                                                sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr) - (2*cos_alpha_i_2 -
                                                                                                                                                                                              2*sin_alpha_i_2*u0*u0/vi_mag_sqr + 2*sin_alpha_i_2*u1*u1/vi_mag_sqr -
                                                                                                                                                                                              2*sin_alpha_i_2*u2*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x1 + y1) -
                                                                                                                                                                                                                                 cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
                                                                                                                                                                                                                                 cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
                                                                                                                                                                                                                                 sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                                                                                                  2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
                                                                                                                                                                                                                                 sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                                                                                                                                 sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr))/(2*a*a))*exp(-
                                                                                                                                                                                                                                                                                          (t4*t4)/(2*b*b) + (-(t5*t5) - (t6*t6))/(2*a*a));
    
    t7 = cos_alpha_i_2*(-x2 + y2) + cos_alpha_i*sin_alpha_i*u0*(-2*x1 +
                                                                2*y1)/vi_mag + cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
            sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr;
    t8 = cos_alpha_i_2*(-x0 + y0) + cos_alpha_i*sin_alpha_i*u1*(-2*x2 +
                                                                2*y2)/vi_mag - cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
            sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                            2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr;
    t9 = cos_alpha_i_2*(-
                        x1 + y1) - cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
            cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
            sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                             2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
            sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
            sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr;
    DpsiDy2 = K*(-(2*cos_alpha_i_2 -
                   2*sin_alpha_i_2*u0*u0/vi_mag_sqr - 2*sin_alpha_i_2*u1*u1/vi_mag_sqr +
                   2*sin_alpha_i_2*u2*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x2 + y2) +
                                                      cos_alpha_i*sin_alpha_i*u0*(-2*x1 + 2*y1)/vi_mag +
                                                      cos_alpha_i*sin_alpha_i*u1*(2*x0 - 2*y0)/vi_mag -
                                                      sin_alpha_i_2*u0*u0*(-x2 + y2)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-
                                                                                                                       2*x0 + 2*y0)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x2 + y2)/vi_mag_sqr +
                                                      sin_alpha_i_2*u1*u2*(-2*x1 + 2*y1)/vi_mag_sqr +
                                                      sin_alpha_i_2*u2*u2*(-x2 + y2)/vi_mag_sqr)/(2*b*b) + (-(-
                                                                                                              4*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                                                              4*sin_alpha_i_2*u1*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x1 + y1) -
                                                                                                                                                 cos_alpha_i*sin_alpha_i*u0*(-2*x2 + 2*y2)/vi_mag +
                                                                                                                                                 cos_alpha_i*sin_alpha_i*u2*(-2*x0 + 2*y0)/vi_mag -
                                                                                                                                                 sin_alpha_i_2*u0*u0*(-x1 + y1)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                  2*x0 + 2*y0)/vi_mag_sqr + sin_alpha_i_2*u1*u1*(-x1 + y1)/vi_mag_sqr +
                                                                                                                                                 sin_alpha_i_2*u1*u2*(-2*x2 + 2*y2)/vi_mag_sqr -
                                                                                                                                                 sin_alpha_i_2*u2*u2*(-x1 + y1)/vi_mag_sqr) -
                                                                                                            (4*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                                                                                                             4*sin_alpha_i_2*u0*u2/vi_mag_sqr)*(cos_alpha_i_2*(-x0 + y0) +
                                                                                                                                                cos_alpha_i*sin_alpha_i*u1*(-2*x2 + 2*y2)/vi_mag -
                                                                                                                                                cos_alpha_i*sin_alpha_i*u2*(-2*x1 + 2*y1)/vi_mag +
                                                                                                                                                sin_alpha_i_2*u0*u0*(-x0 + y0)/vi_mag_sqr + sin_alpha_i_2*u0*u1*(-
                                                                                                                                                                                                                 2*x1 + 2*y1)/vi_mag_sqr + sin_alpha_i_2*u0*u2*(-2*x2 +
                                                                                                                                                                                                                                                                2*y2)/vi_mag_sqr - sin_alpha_i_2*u1*u1*(-x0 + y0)/vi_mag_sqr -
                                                                                                                                                sin_alpha_i_2*u2*u2*(-x0 + y0)/vi_mag_sqr))/(2*a*a))*exp(-
                                                                                                                                                                                                         (t7*t7)/(2*b*b) + (-(t8*t8) - (t9*t9))/(2*a*a));

    Vector3D DpsiDxj( DpsiDy0, DpsiDy1, DpsiDy2 );
    return DpsiDxj;
}

/*
*Derivative of phi_p wrt vi
 */
Vector3D OPSBody::Dphi_pDvi(Vector3D vi, Vector3D xi, Vector3D xj){
    double u0, u1, u2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, sin_alpha_i, cos_alpha_i;
    double sin_alpha_i_2, cos_alpha_i_2, phi_p0;
    double Dphi_pDu0, Dphi_pDu1, Dphi_pDu2;
    double vi_mag_3, vi_mag_4;
    double u0_3, u1_3, u2_3;
    
    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    u0_3 = u0*u0*u0;
    u1_3 = u1*u1*u1;
    u2_3 = u2*u2*u2;
    
    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vi_mag = sqrt(vi_mag_sqr);
    vi_mag_3 = vi_mag*vi_mag*vi_mag;
    vi_mag_4 = vi_mag*vi_mag*vi_mag*vi_mag;
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;

    Dphi_pDu0 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 -
                                                                                                                       sin_alpha_i_2*u0*u0/vi_mag_sqr - sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                       sin_alpha_i_2*u2*u2/vi_mag_sqr))*(2*(-x0 +
                                                                                                                                                            y0)*(cos_alpha_i_2*u0*u1/vi_mag_sqr +
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u0*u2/vi_mag_3 -
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u1/vi_mag_3 -
                                                                                                                                                                 4*sin_alpha_i_2*u0*u0*u2/vi_mag_4 -
                                                                                                                                                                 sin_alpha_i_2*u0*u1/vi_mag_sqr + 2*sin_alpha_i_2*u2/vi_mag_sqr) +
                                                                                                                                                         2*(-x1 + y1)*(-cos_alpha_i_2*u0*u0/vi_mag_sqr +
                                                                                                                                                                       2*cos_alpha_i*sin_alpha_i*u0*u0/vi_mag_3 +
                                                                                                                                                                       2*cos_alpha_i*sin_alpha_i*u0*u1*u2/vi_mag_3 -
                                                                                                                                                                       2*cos_alpha_i*sin_alpha_i/vi_mag + sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                                                                       4*sin_alpha_i_2*u0*u1*u2/vi_mag_4) + 2*(-x2 + y2)*(-
                                                                                                                                                                                                                          cos_alpha_i*sin_alpha_i*u0_3/vi_mag_3 -
                                                                                                                                                                                                                          cos_alpha_i*sin_alpha_i*u0*u1*u1/vi_mag_3 +
                                                                                                                                                                                                                          cos_alpha_i*sin_alpha_i*u0*u2*u2/vi_mag_3 -
                                                                                                                                                                                                                          cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                                                                                                                                                                          2*sin_alpha_i_2*u0_3/vi_mag_4 +
                                                                                                                                                                                                                          2*sin_alpha_i_2*u0*u1*u1/vi_mag_4 -
                                                                                                                                                                                                                          2*sin_alpha_i_2*u0*u2*u2/vi_mag_4 -
                                                                                                                                                                                                                          2*sin_alpha_i_2*u0/vi_mag_sqr));

    Dphi_pDu1 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 -
                                                                                                                       sin_alpha_i_2*u0*u0/vi_mag_sqr - sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                       sin_alpha_i_2*u2*u2/vi_mag_sqr))*(2*(-x0 +
                                                                                                                                                            y0)*(cos_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u1*u2/vi_mag_3 -
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u1*u1/vi_mag_3 +
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i/vi_mag -
                                                                                                                                                                 4*sin_alpha_i_2*u0*u1*u2/vi_mag_4 -
                                                                                                                                                                 sin_alpha_i_2*u1*u1/vi_mag_sqr) + 2*(-x1 + y1)*(-
                                                                                                                                                                                                                 cos_alpha_i_2*u0*u1/vi_mag_sqr +
                                                                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u1/vi_mag_3 +
                                                                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u1*u1*u2/vi_mag_3 +
                                                                                                                                                                                                                 sin_alpha_i_2*u0*u1/vi_mag_sqr -
                                                                                                                                                                                                                 4*sin_alpha_i_2*u1*u1*u2/vi_mag_4 + 2*sin_alpha_i_2*u2/vi_mag_sqr) + 2*(-x2 + y2)*(-
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u0*u0*u1/vi_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u1_3/vi_mag_3 +
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u1*u2*u2/vi_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u1/vi_mag +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u0*u0*u1/vi_mag_4 +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u1_3/vi_mag_4 -
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u1*u2*u2/vi_mag_4 -
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u1/vi_mag_sqr));

    Dphi_pDu2 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 -
                                                                                                                       sin_alpha_i_2*u0*u0/vi_mag_sqr - sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                       sin_alpha_i_2*u2*u2/vi_mag_sqr))*(2*(-x0 +
                                                                                                                                                            y0)*(cos_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u2*u2/vi_mag_3 -
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u1*u2/vi_mag_3 -
                                                                                                                                                                 4*sin_alpha_i_2*u0*u2*u2/vi_mag_4 + 2*sin_alpha_i_2*u0/vi_mag_sqr -
                                                                                                                                                                 sin_alpha_i_2*u1*u2/vi_mag_sqr) + 2*(-x1 + y1)*(-
                                                                                                                                                                                                                 cos_alpha_i_2*u0*u2/vi_mag_sqr +
                                                                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u2/vi_mag_3 +
                                                                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u1*u2*u2/vi_mag_3 +
                                                                                                                                                                                                                 sin_alpha_i_2*u0*u2/vi_mag_sqr -
                                                                                                                                                                                                                 4*sin_alpha_i_2*u1*u2*u2/vi_mag_4 + 2*sin_alpha_i_2*u1/vi_mag_sqr) + 2*(-x2 + y2)*(-
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u0*u0*u2/vi_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u1*u1*u2/vi_mag_3 +
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u2_3/vi_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u2/vi_mag +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u0*u0*u2/vi_mag_4 +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u1*u1*u2/vi_mag_4 -
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u2_3/vi_mag_4 + 2*sin_alpha_i_2*u2/vi_mag_sqr));

    Vector3D Dphi_pDvi(Dphi_pDu0, Dphi_pDu1, Dphi_pDu2);
    return Dphi_pDvi;
}

/*
*Derivative of phi_p wrt xi
 */
Vector3D OPSBody::Dphi_pDxi(Vector3D vi, Vector3D xi, Vector3D xj){
    double u0, u1, u2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, sin_alpha_i, cos_alpha_i;
    double sin_alpha_i_2, cos_alpha_i_2, phi_p0;
    double Dphi_pDx0, Dphi_pDx1, Dphi_pDx2;

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vi_mag = sqrt(vi_mag_sqr);
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;

    Dphi_pDx0 = (-
                 4*cos_alpha_i*sin_alpha_i*u1/vi_mag -
                 4*sin_alpha_i_2*u0*u2/vi_mag_sqr)*((-x0 +
                                                     y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                                                          2*sin_alpha_i_2*u0*u2/vi_mag_sqr) + (-x1 + y1)*(-
                                                                                                          2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                                                          2*sin_alpha_i_2*u1*u2/vi_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 -
                                                                                                                                                          sin_alpha_i_2*u0*u0/vi_mag_sqr - sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                                                          sin_alpha_i_2*u2*u2/vi_mag_sqr));

    Dphi_pDx1 = (4*cos_alpha_i
                 *sin_alpha_i*u0/vi_mag - 4*sin_alpha_i_2*u1*u2/vi_mag_sqr)*((-x0 +
                                                                              y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                                                                                   2*sin_alpha_i_2*u0*u2/vi_mag_sqr) + (-x1 + y1)*(-
                                                                                                                                   2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                                                                                   2*sin_alpha_i_2*u1*u2/vi_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 -
                                                                                                                                                                                   sin_alpha_i_2*u0*u0/vi_mag_sqr - sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                                                                                   sin_alpha_i_2*u2*u2/vi_mag_sqr));

    Dphi_pDx2 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 -
                                                                                                                       sin_alpha_i_2*u0*u0/vi_mag_sqr - sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                       sin_alpha_i_2*u2*u2/vi_mag_sqr))*(-2*cos_alpha_i_2 +
                                                                                                                                                         2*sin_alpha_i_2*u0*u0/vi_mag_sqr + 2*sin_alpha_i_2*u1*u1/vi_mag_sqr -
                                                                                                                                                         2*sin_alpha_i_2*u2*u2/vi_mag_sqr);

    Vector3D Dphi_pDxi(Dphi_pDx0, Dphi_pDx1, Dphi_pDx2);
    return Dphi_pDxi;
}

/*
*Derivative of phi_p wrt xj
 */
Vector3D OPSBody::Dphi_pDxj(Vector3D vi, Vector3D xi, Vector3D xj){
    double u0, u1, u2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, sin_alpha_i, cos_alpha_i;
    double sin_alpha_i_2, cos_alpha_i_2, phi_p0;
    double Dphi_pDy0, Dphi_pDy1, Dphi_pDy2;

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vi_mag = sqrt(vi_mag_sqr);
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;

    Dphi_pDy0 = (4*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                 4*sin_alpha_i_2*u0*u2/vi_mag_sqr)*((-x0 +
                                                     y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                                                          2*sin_alpha_i_2*u0*u2/vi_mag_sqr) + (-x1 + y1)*(-
                                                                                                          2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                                                          2*sin_alpha_i_2*u1*u2/vi_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 -
                                                                                                                                                          sin_alpha_i_2*u0*u0/vi_mag_sqr - sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                                                          sin_alpha_i_2*u2*u2/vi_mag_sqr));

    Dphi_pDy1 = (-
                 4*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                 4*sin_alpha_i_2*u1*u2/vi_mag_sqr)*((-x0 +
                                                     y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                                                          2*sin_alpha_i_2*u0*u2/vi_mag_sqr) + (-x1 + y1)*(-
                                                                                                          2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                                                          2*sin_alpha_i_2*u1*u2/vi_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 -
                                                                                                                                                          sin_alpha_i_2*u0*u0/vi_mag_sqr - sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                                                          sin_alpha_i_2*u2*u2/vi_mag_sqr));

    Dphi_pDy2 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 -
                                                                                                                       sin_alpha_i_2*u0*u0/vi_mag_sqr - sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                       sin_alpha_i_2*u2*u2/vi_mag_sqr))*(2*cos_alpha_i_2 -
                                                                                                                                                         2*sin_alpha_i_2*u0*u0/vi_mag_sqr - 2*sin_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                                                         2*sin_alpha_i_2*u2*u2/vi_mag_sqr);

    Vector3D Dphi_pDxj(Dphi_pDy0, Dphi_pDy1, Dphi_pDy2);
    return Dphi_pDxj;
}

/*
*Derivative of phi_n wrt vi
 */
Vector3D OPSBody::Dphi_nDvi(Vector3D vi, Vector3D vj){
    double u0, u1, u2, v0, v1, v2;
    double vi_mag_sqr, vi_mag, vj_mag_sqr, vj_mag;
    double sin_alpha_i, cos_alpha_i, sin_alpha_j, cos_alpha_j;
    double sin_alpha_i_2, cos_alpha_i_2, sin_alpha_j_2, cos_alpha_j_2, phi_n0;
    double Dphi_nDu0, Dphi_nDu1, Dphi_nDu2;
    double vi_mag_3, vi_mag_4;
    double u0_3, u1_3, u2_3;

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    v0 = vj[0]; v1 = vj[1]; v2 = vj[2];
    
    u0_3 = u0*u0*u0;
    u1_3 = u1*u1*u1;
    u2_3 = u2*u2*u2;
    
    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vj_mag_sqr = v0*v0 + v1*v1 + v2*v2;
    vi_mag = sqrt(vi_mag_sqr);
    vj_mag = sqrt(vj_mag_sqr);
    vi_mag_3 = vi_mag*vi_mag*vi_mag;
    vi_mag_4 = vi_mag*vi_mag*vi_mag*vi_mag;
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;
    sin_alpha_j = sin( 0.5*vj_mag );
    cos_alpha_j = cos( 0.5*vj_mag );
    sin_alpha_j_2 = sin_alpha_j*sin_alpha_j;
    cos_alpha_j_2 = cos_alpha_j*cos_alpha_j;

    Dphi_nDu0 = (-
                 2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                 2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                 2*sin_alpha_i_2*u1*u2/vi_mag_sqr -
                 2*sin_alpha_j_2*v1*v2/vj_mag_sqr)*(-
                                                    2*cos_alpha_i_2*u0*u0/vi_mag_sqr + 4*
                                                    cos_alpha_i*sin_alpha_i*u0*u0/vi_mag_3 +
                                                    4*cos_alpha_i*sin_alpha_i*u0*u1*u2/vi_mag_3 -
                                                    4*cos_alpha_i*sin_alpha_i/vi_mag + 2*sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                    8*sin_alpha_i_2*u0*u1*u2/vi_mag_4) +
            (2*cos_alpha_i*sin_alpha_i*u1/vi_mag -
             2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
             2*sin_alpha_i_2*u0*u2/vi_mag_sqr -
             2*sin_alpha_j_2*v0*v2/vj_mag_sqr)*(2*cos_alpha_i_2*u0*u1/vi_mag_sqr +
                                                4*cos_alpha_i*sin_alpha_i*u0*u0*u2/vi_mag_3 -
                                                4*cos_alpha_i*sin_alpha_i*u0*u1/vi_mag_3 -
                                                8*sin_alpha_i_2*u0*u0*u2/vi_mag_4 -
                                                2*sin_alpha_i_2*u0*u1/vi_mag_sqr + 4*sin_alpha_i_2*u2/vi_mag_sqr) +
            (cos_alpha_i_2 - cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
             sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr +
             sin_alpha_j_2*v0*v0/vj_mag_sqr + sin_alpha_j_2*v1*v1/vj_mag_sqr -
             sin_alpha_j_2*v2*v2/vj_mag_sqr)*(-
                                              2*cos_alpha_i*sin_alpha_i*u0_3/vi_mag_3 -
                                              2*cos_alpha_i*sin_alpha_i*u0*u1*u1/vi_mag_3 +
                                              2*cos_alpha_i*sin_alpha_i*u0*u2*u2/vi_mag_3 -
                                              2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                              4*sin_alpha_i_2*u0_3/vi_mag_4 +
                                              4*sin_alpha_i_2*u0*u1*u1/vi_mag_4 -
                                              4*sin_alpha_i_2*u0*u2*u2/vi_mag_4 -
                                              4*sin_alpha_i_2*u0/vi_mag_sqr);

    Dphi_nDu1 = (-
                 2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                 2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                 2*sin_alpha_i_2*u1*u2/vi_mag_sqr -
                 2*sin_alpha_j_2*v1*v2/vj_mag_sqr)*(-
                                                    2*cos_alpha_i_2*u0*u1/vi_mag_sqr + 4*
                                                    cos_alpha_i*sin_alpha_i*u0*u1/vi_mag_3 +
                                                    4*cos_alpha_i*sin_alpha_i*u1*u1*u2/vi_mag_3 +
                                                    2*sin_alpha_i_2*u0*u1/vi_mag_sqr -
                                                    8*sin_alpha_i_2*u1*u1*u2/vi_mag_4 +
                                                    4*sin_alpha_i_2*u2/vi_mag_sqr) + (2*
                                                                                      cos_alpha_i*sin_alpha_i*u1/vi_mag -
                                                                                      2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                                                                                      2*sin_alpha_i_2*u0*u2/vi_mag_sqr -
                                                                                      2*sin_alpha_j_2*v0*v2/vj_mag_sqr)*(2*cos_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                         4*cos_alpha_i*sin_alpha_i*u0*u1*u2/vi_mag_3 -
                                                                                                                         4*cos_alpha_i*sin_alpha_i*u1*u1/vi_mag_3 +
                                                                                                                         4*cos_alpha_i*sin_alpha_i/vi_mag -
                                                                                                                         8*sin_alpha_i_2*u0*u1*u2/vi_mag_4 -
                                                                                                                         2*sin_alpha_i_2*u1*u1/vi_mag_sqr) + (cos_alpha_i_2 -
                                                                                                                                                              cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                                                              sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr +
                                                                                                                                                              sin_alpha_j_2*v0*v0/vj_mag_sqr + sin_alpha_j_2*v1*v1/vj_mag_sqr -
                                                                                                                                                              sin_alpha_j_2*v2*v2/vj_mag_sqr)*(-
                                                                                                                                                                                               2*cos_alpha_i*sin_alpha_i*u0*u0*u1/vi_mag_3 -
                                                                                                                                                                                               2*cos_alpha_i*sin_alpha_i*u1_3/vi_mag_3 +
                                                                                                                                                                                               2*cos_alpha_i*sin_alpha_i*u1*u2*u2/vi_mag_3 -
                                                                                                                                                                                               2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                                                                                                                                                                                               4*sin_alpha_i_2*u0*u0*u1/vi_mag_4 +
                                                                                                                                                                                               4*sin_alpha_i_2*u1_3/vi_mag_4 -
                                                                                                                                                                                               4*sin_alpha_i_2*u1*u2*u2/vi_mag_4 -
                                                                                                                                                                                               4*sin_alpha_i_2*u1/vi_mag_sqr);

    Dphi_nDu2 = (-
                 2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                 2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                 2*sin_alpha_i_2*u1*u2/vi_mag_sqr -
                 2*sin_alpha_j_2*v1*v2/vj_mag_sqr)*(-
                                                    2*cos_alpha_i_2*u0*u2/vi_mag_sqr + 4*
                                                    cos_alpha_i*sin_alpha_i*u0*u2/vi_mag_3 +
                                                    4*cos_alpha_i*sin_alpha_i*u1*u2*u2/vi_mag_3 +
                                                    2*sin_alpha_i_2*u0*u2/vi_mag_sqr -
                                                    8*sin_alpha_i_2*u1*u2*u2/vi_mag_4 +
                                                    4*sin_alpha_i_2*u1/vi_mag_sqr) + (2*
                                                                                      cos_alpha_i*sin_alpha_i*u1/vi_mag -
                                                                                      2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                                                                                      2*sin_alpha_i_2*u0*u2/vi_mag_sqr -
                                                                                      2*sin_alpha_j_2*v0*v2/vj_mag_sqr)*(2*cos_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                                                                         4*cos_alpha_i*sin_alpha_i*u0*u2*u2/vi_mag_3 -
                                                                                                                         4*cos_alpha_i*sin_alpha_i*u1*u2/vi_mag_3 -
                                                                                                                         8*sin_alpha_i_2*u0*u2*u2/vi_mag_4 + 4*sin_alpha_i_2*u0/vi_mag_sqr -
                                                                                                                         2*sin_alpha_i_2*u1*u2/vi_mag_sqr) + (cos_alpha_i_2 -
                                                                                                                                                              cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                                                              sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr +
                                                                                                                                                              sin_alpha_j_2*v0*v0/vj_mag_sqr + sin_alpha_j_2*v1*v1/vj_mag_sqr -
                                                                                                                                                              sin_alpha_j_2*v2*v2/vj_mag_sqr)*(-
                                                                                                                                                                                               2*cos_alpha_i*sin_alpha_i*u0*u0*u2/vi_mag_3 -
                                                                                                                                                                                               2*cos_alpha_i*sin_alpha_i*u1*u1*u2/vi_mag_3 +
                                                                                                                                                                                               2*cos_alpha_i*sin_alpha_i*u2_3/vi_mag_3 -
                                                                                                                                                                                               2*cos_alpha_i*sin_alpha_i*u2/vi_mag +
                                                                                                                                                                                               4*sin_alpha_i_2*u0*u0*u2/vi_mag_4 +
                                                                                                                                                                                               4*sin_alpha_i_2*u1*u1*u2/vi_mag_4 -
                                                                                                                                                                                               4*sin_alpha_i_2*u2_3/vi_mag_4 + 4*sin_alpha_i_2*u2/vi_mag_sqr);

    Vector3D Dphi_nDvi( Dphi_nDu0, Dphi_nDu1, Dphi_nDu2 );
    return Dphi_nDvi;
}

/*
*Derivative of phi_n wrt vj
 */
Vector3D OPSBody::Dphi_nDvj(Vector3D vi, Vector3D vj){
    double u0, u1, u2, v0, v1, v2;
    double vi_mag_sqr, vi_mag, vj_mag_sqr, vj_mag;
    double sin_alpha_i, cos_alpha_i, sin_alpha_j, cos_alpha_j;
    double sin_alpha_i_2, cos_alpha_i_2, sin_alpha_j_2, cos_alpha_j_2, phi_n0;
    double Dphi_nDv0, Dphi_nDv1, Dphi_nDv2;
    double vj_mag_3, vj_mag_4;
    double v0_3, v1_3, v2_3;

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    v0 = vj[0]; v1 = vj[1]; v2 = vj[2];

    v0_3 = v0*v0*v0;
    v1_3 = v1*v1*v1;
    v2_3 = v2*v2*v2;
    
    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vj_mag_sqr = v0*v0 + v1*v1 + v2*v2;
    vi_mag = sqrt(vi_mag_sqr);
    vj_mag = sqrt(vj_mag_sqr);
    vj_mag_3 = vj_mag*vj_mag*vj_mag;
    vj_mag_4 = vj_mag*vj_mag*vj_mag*vj_mag;
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;
    sin_alpha_j = sin( 0.5*vj_mag );
    cos_alpha_j = cos( 0.5*vj_mag );
    sin_alpha_j_2 = sin_alpha_j*sin_alpha_j;
    cos_alpha_j_2 = cos_alpha_j*cos_alpha_j;

    Dphi_nDv0 = (-
                 2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                 2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                 2*sin_alpha_i_2*u1*u2/vi_mag_sqr -
                 2*sin_alpha_j_2*v1*v2/vj_mag_sqr)*(2*cos_alpha_j_2*v0*v0/vj_mag_sqr -
                                                    4*cos_alpha_j*sin_alpha_j*v0*v0/vj_mag_3 -
                                                    4*cos_alpha_j*sin_alpha_j*v0*v1*v2/vj_mag_3 +
                                                    4*cos_alpha_j*sin_alpha_j/vj_mag - 2*sin_alpha_j_2*v0*v0/vj_mag_sqr +
                                                    8*sin_alpha_j_2*v0*v1*v2/vj_mag_4) +
            (2*cos_alpha_i*sin_alpha_i*u1/vi_mag -
             2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
             2*sin_alpha_i_2*u0*u2/vi_mag_sqr -
             2*sin_alpha_j_2*v0*v2/vj_mag_sqr)*(-
                                                2*cos_alpha_j_2*v0*v1/vj_mag_sqr - 4*
                                                cos_alpha_j*sin_alpha_j*v0*v0*v2/vj_mag_3 +
                                                4*cos_alpha_j*sin_alpha_j*v0*v1/vj_mag_3 +
                                                8*sin_alpha_j_2*v0*v0*v2/vj_mag_4 +
                                                2*sin_alpha_j_2*v0*v1/vj_mag_sqr - 4*sin_alpha_j_2*v2/vj_mag_sqr) +
            (cos_alpha_i_2 - cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
             sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr +
             sin_alpha_j_2*v0*v0/vj_mag_sqr + sin_alpha_j_2*v1*v1/vj_mag_sqr -
             sin_alpha_j_2*v2*v2/vj_mag_sqr)*(2*cos_alpha_j*sin_alpha_j*v0_3/
                                              vj_mag_3 + 2*cos_alpha_j*sin_alpha_j*v0*v1*v1/vj_mag_3 -
                                              2*cos_alpha_j*sin_alpha_j*v0*v2*v2/vj_mag_3 +
                                              2*cos_alpha_j*sin_alpha_j*v0/vj_mag -
                                              4*sin_alpha_j_2*v0_3/vj_mag_4 -
                                              4*sin_alpha_j_2*v0*v1*v1/vj_mag_4 +
                                              4*sin_alpha_j_2*v0*v2*v2/vj_mag_4 +
                                              4*sin_alpha_j_2*v0/vj_mag_sqr);

    Dphi_nDv1 = (-
                 2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                 2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                 2*sin_alpha_i_2*u1*u2/vi_mag_sqr -
                 2*sin_alpha_j_2*v1*v2/vj_mag_sqr)*(2*cos_alpha_j_2*v0*v1/vj_mag_sqr -
                                                    4*cos_alpha_j*sin_alpha_j*v0*v1/vj_mag_3 -
                                                    4*cos_alpha_j*sin_alpha_j*v1*v1*v2/vj_mag_3 -
                                                    2*sin_alpha_j_2*v0*v1/vj_mag_sqr +
                                                    8*sin_alpha_j_2*v1*v1*v2/vj_mag_4 -
                                                    4*sin_alpha_j_2*v2/vj_mag_sqr) + (2*
                                                                                      cos_alpha_i*sin_alpha_i*u1/vi_mag -
                                                                                      2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                                                                                      2*sin_alpha_i_2*u0*u2/vi_mag_sqr -
                                                                                      2*sin_alpha_j_2*v0*v2/vj_mag_sqr)*(-
                                                                                                                         2*cos_alpha_j_2*v1*v1/vj_mag_sqr - 4*
                                                                                                                         cos_alpha_j*sin_alpha_j*v0*v1*v2/vj_mag_3 +
                                                                                                                         4*cos_alpha_j*sin_alpha_j*v1*v1/vj_mag_3 -
                                                                                                                         4*cos_alpha_j*sin_alpha_j/vj_mag +
                                                                                                                         8*sin_alpha_j_2*v0*v1*v2/vj_mag_4 +
                                                                                                                         2*sin_alpha_j_2*v1*v1/vj_mag_sqr) + (cos_alpha_i_2 -
                                                                                                                                                              cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                                                              sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr +
                                                                                                                                                              sin_alpha_j_2*v0*v0/vj_mag_sqr + sin_alpha_j_2*v1*v1/vj_mag_sqr -
                                                                                                                                                              sin_alpha_j_2*v2*v2/vj_mag_sqr)*(2*cos_alpha_j*sin_alpha_j*v0*v0*
                                                                                                                                                                                               v1/vj_mag_3 + 2*cos_alpha_j*sin_alpha_j*v1_3/vj_mag_3 -
                                                                                                                                                                                               2*cos_alpha_j*sin_alpha_j*v1*v2*v2/vj_mag_3 +
                                                                                                                                                                                               2*cos_alpha_j*sin_alpha_j*v1/vj_mag -
                                                                                                                                                                                               4*sin_alpha_j_2*v0*v0*v1/vj_mag_4 -
                                                                                                                                                                                               4*sin_alpha_j_2*v1_3/vj_mag_4 +
                                                                                                                                                                                               4*sin_alpha_j_2*v1*v2*v2/vj_mag_4 +
                                                                                                                                                                                               4*sin_alpha_j_2*v1/vj_mag_sqr);

    Dphi_nDv2 = (-
                 2*cos_alpha_i*sin_alpha_i*u0/vi_mag +
                 2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                 2*sin_alpha_i_2*u1*u2/vi_mag_sqr -
                 2*sin_alpha_j_2*v1*v2/vj_mag_sqr)*(2*cos_alpha_j_2*v0*v2/vj_mag_sqr -
                                                    4*cos_alpha_j*sin_alpha_j*v0*v2/vj_mag_3 -
                                                    4*cos_alpha_j*sin_alpha_j*v1*v2*v2/vj_mag_3 -
                                                    2*sin_alpha_j_2*v0*v2/vj_mag_sqr +
                                                    8*sin_alpha_j_2*v1*v2*v2/vj_mag_4 -
                                                    4*sin_alpha_j_2*v1/vj_mag_sqr) + (2*
                                                                                      cos_alpha_i*sin_alpha_i*u1/vi_mag -
                                                                                      2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                                                                                      2*sin_alpha_i_2*u0*u2/vi_mag_sqr -
                                                                                      2*sin_alpha_j_2*v0*v2/vj_mag_sqr)*(-
                                                                                                                         2*cos_alpha_j_2*v1*v2/vj_mag_sqr - 4*
                                                                                                                         cos_alpha_j*sin_alpha_j*v0*v2*v2/vj_mag_3 +
                                                                                                                         4*cos_alpha_j*sin_alpha_j*v1*v2/vj_mag_3 +
                                                                                                                         8*sin_alpha_j_2*v0*v2*v2/vj_mag_4 - 4*sin_alpha_j_2*v0/vj_mag_sqr +
                                                                                                                         2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (cos_alpha_i_2 -
                                                                                                                                                              cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                                                              sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr +
                                                                                                                                                              sin_alpha_j_2*v0*v0/vj_mag_sqr + sin_alpha_j_2*v1*v1/vj_mag_sqr -
                                                                                                                                                              sin_alpha_j_2*v2*v2/vj_mag_sqr)*(2*cos_alpha_j*sin_alpha_j*v0*v0*
                                                                                                                                                                                               v2/vj_mag_3 + 2*cos_alpha_j*sin_alpha_j*v1*v1*v2/vj_mag_3 -
                                                                                                                                                                                               2*cos_alpha_j*sin_alpha_j*v2_3/vj_mag_3 +
                                                                                                                                                                                               2*cos_alpha_j*sin_alpha_j*v2/vj_mag -
                                                                                                                                                                                               4*sin_alpha_j_2*v0*v0*v2/vj_mag_4 -
                                                                                                                                                                                               4*sin_alpha_j_2*v1*v1*v2/vj_mag_4 +
                                                                                                                                                                                               4*sin_alpha_j_2*v2_3/vj_mag_4 - 4*sin_alpha_j_2*v2/vj_mag_sqr);

    Vector3D Dphi_nDvj( Dphi_nDv0, Dphi_nDv1, Dphi_nDv2 );
    return Dphi_nDvj;
}

/*
*Derivative of phi_c wrt vi
 */
Vector3D OPSBody::Dphi_cDvi(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj){
    double u0, u1, u2, v0, v1, v2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, vj_mag_sqr, vj_mag;
    double sin_alpha_i, cos_alpha_i, sin_alpha_j, cos_alpha_j;
    double sin_alpha_i_2, cos_alpha_i_2, sin_alpha_j_2, cos_alpha_j_2, phi_c0;
    double Dphi_cDu0, Dphi_cDu1, Dphi_cDu2;
    double vi_mag_3, vi_mag_4;
    double u0_3, u1_3, u2_3;

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    v0 = vj[0]; v1 = vj[1]; v2 = vj[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    u0_3 = u0*u0*u0;
    u1_3 = u1*u1*u1;
    u2_3 = u2*u2*u2;
    
    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vj_mag_sqr = v0*v0 + v1*v1 + v2*v2;
    vi_mag = sqrt(vi_mag_sqr);
    vj_mag = sqrt(vj_mag_sqr);
    vi_mag_3 = vi_mag*vi_mag*vi_mag;
    vi_mag_4 = vi_mag*vi_mag*vi_mag*vi_mag;
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;
    sin_alpha_j = sin( 0.5*vj_mag );
    cos_alpha_j = cos( 0.5*vj_mag );
    sin_alpha_j_2 = sin_alpha_j*sin_alpha_j;
    cos_alpha_j_2 = cos_alpha_j*cos_alpha_j;

    Dphi_cDu0 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(2*(-x0 +
                                                                                                                                                            y0)*(cos_alpha_i_2*u0*u1/vi_mag_sqr +
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u0*u2/vi_mag_3 -
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u1/vi_mag_3 -
                                                                                                                                                                 4*sin_alpha_i_2*u0*u0*u2/vi_mag_4 -
                                                                                                                                                                 sin_alpha_i_2*u0*u1/vi_mag_sqr + 2*sin_alpha_i_2*u2/vi_mag_sqr) +
                                                                                                                                                         2*(-x1 + y1)*(-cos_alpha_i_2*u0*u0/vi_mag_sqr +
                                                                                                                                                                       2*cos_alpha_i*sin_alpha_i*u0*u0/vi_mag_3 +
                                                                                                                                                                       2*cos_alpha_i*sin_alpha_i*u0*u1*u2/vi_mag_3 -
                                                                                                                                                                       2*cos_alpha_i*sin_alpha_i/vi_mag + sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                                                                       4*sin_alpha_i_2*u0*u1*u2/vi_mag_4) + 2*(-x2 + y2)*(-
                                                                                                                                                                                                                          cos_alpha_i*sin_alpha_i*u0_3/vi_mag_3 -
                                                                                                                                                                                                                          cos_alpha_i*sin_alpha_i*u0*u1*u1/vi_mag_3 +
                                                                                                                                                                                                                          cos_alpha_i*sin_alpha_i*u0*u2*u2/vi_mag_3 -
                                                                                                                                                                                                                          cos_alpha_i*sin_alpha_i*u0/vi_mag +
                                                                                                                                                                                                                          2*sin_alpha_i_2*u0_3/vi_mag_4 +
                                                                                                                                                                                                                          2*sin_alpha_i_2*u0*u1*u1/vi_mag_4 -
                                                                                                                                                                                                                          2*sin_alpha_i_2*u0*u2*u2/vi_mag_4 -
                                                                                                                                                                                                                          2*sin_alpha_i_2*u0/vi_mag_sqr));

    Dphi_cDu1 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(2*(-x0 +
                                                                                                                                                            y0)*(cos_alpha_i_2*u1*u1/vi_mag_sqr +
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u1*u2/vi_mag_3 -
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u1*u1/vi_mag_3 +
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i/vi_mag -
                                                                                                                                                                 4*sin_alpha_i_2*u0*u1*u2/vi_mag_4 -
                                                                                                                                                                 sin_alpha_i_2*u1*u1/vi_mag_sqr) + 2*(-x1 + y1)*(-
                                                                                                                                                                                                                 cos_alpha_i_2*u0*u1/vi_mag_sqr +
                                                                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u1/vi_mag_3 +
                                                                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u1*u1*u2/vi_mag_3 +
                                                                                                                                                                                                                 sin_alpha_i_2*u0*u1/vi_mag_sqr -
                                                                                                                                                                                                                 4*sin_alpha_i_2*u1*u1*u2/vi_mag_4 + 2*sin_alpha_i_2*u2/vi_mag_sqr) + 2*(-x2 + y2)*(-
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u0*u0*u1/vi_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u1_3/vi_mag_3 +
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u1*u2*u2/vi_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u1/vi_mag +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u0*u0*u1/vi_mag_4 +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u1_3/vi_mag_4 -
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u1*u2*u2/vi_mag_4 -
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u1/vi_mag_sqr));

    Dphi_cDu2 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(2*(-x0 +
                                                                                                                                                            y0)*(cos_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u2*u2/vi_mag_3 -
                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u1*u2/vi_mag_3 -
                                                                                                                                                                 4*sin_alpha_i_2*u0*u2*u2/vi_mag_4 + 2*sin_alpha_i_2*u0/vi_mag_sqr -
                                                                                                                                                                 sin_alpha_i_2*u1*u2/vi_mag_sqr) + 2*(-x1 + y1)*(-
                                                                                                                                                                                                                 cos_alpha_i_2*u0*u2/vi_mag_sqr +
                                                                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u0*u2/vi_mag_3 +
                                                                                                                                                                                                                 2*cos_alpha_i*sin_alpha_i*u1*u2*u2/vi_mag_3 +
                                                                                                                                                                                                                 sin_alpha_i_2*u0*u2/vi_mag_sqr -
                                                                                                                                                                                                                 4*sin_alpha_i_2*u1*u2*u2/vi_mag_4 + 2*sin_alpha_i_2*u1/vi_mag_sqr) + 2*(-x2 + y2)*(-
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u0*u0*u2/vi_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u1*u1*u2/vi_mag_3 +
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u2_3/vi_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_i*sin_alpha_i*u2/vi_mag +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u0*u0*u2/vi_mag_4 +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u1*u1*u2/vi_mag_4 -
                                                                                                                                                                                                                                                                                                    2*sin_alpha_i_2*u2_3/vi_mag_4 + 2*sin_alpha_i_2*u2/vi_mag_sqr));

    Vector3D Dphi_cDvi( Dphi_cDu0, Dphi_cDu1, Dphi_cDu2 );
    return Dphi_cDvi;
}

/*
*Derivative of phi_c wrt vj
 */
Vector3D OPSBody::Dphi_cDvj(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj){
    double u0, u1, u2, v0, v1, v2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, vj_mag_sqr, vj_mag;
    double sin_alpha_i, cos_alpha_i, sin_alpha_j, cos_alpha_j;
    double sin_alpha_i_2, cos_alpha_i_2, sin_alpha_j_2, cos_alpha_j_2, phi_c0;
    double Dphi_cDv0, Dphi_cDv1, Dphi_cDv2;
    double vj_mag_3, vj_mag_4;
    double v0_3, v1_3, v2_3;

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    v0 = vj[0]; v1 = vj[1]; v2 = vj[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    v0_3 = v0*v0*v0;
    v1_3 = v1*v1*v1;
    v2_3 = v2*v2*v2;
    
    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vj_mag_sqr = v0*v0 + v1*v1 + v2*v2;
    vi_mag = sqrt(vi_mag_sqr);
    vj_mag = sqrt(vj_mag_sqr);
    vj_mag_3 = vj_mag*vj_mag*vj_mag;
    vj_mag_4 = vj_mag*vj_mag*vj_mag*vj_mag;
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;
    sin_alpha_j = sin( 0.5*vj_mag );
    cos_alpha_j = cos( 0.5*vj_mag );
    sin_alpha_j_2 = sin_alpha_j*sin_alpha_j;
    cos_alpha_j_2 = cos_alpha_j*cos_alpha_j;

    Dphi_cDv0 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(2*(-x0 +
                                                                                                                                                            y0)*(cos_alpha_j_2*v0*v1/vj_mag_sqr +
                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j*v0*v0*v2/vj_mag_3 -
                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j*v0*v1/vj_mag_3 -
                                                                                                                                                                 4*sin_alpha_j_2*v0*v0*v2/vj_mag_4 -
                                                                                                                                                                 sin_alpha_j_2*v0*v1/vj_mag_sqr + 2*sin_alpha_j_2*v2/vj_mag_sqr) +
                                                                                                                                                         2*(-x1 + y1)*(-cos_alpha_j_2*v0*v0/vj_mag_sqr +
                                                                                                                                                                       2*cos_alpha_j*sin_alpha_j*v0*v0/vj_mag_3 +
                                                                                                                                                                       2*cos_alpha_j*sin_alpha_j*v0*v1*v2/vj_mag_3 -
                                                                                                                                                                       2*cos_alpha_j*sin_alpha_j/vj_mag + sin_alpha_j_2*v0*v0/vj_mag_sqr -
                                                                                                                                                                       4*sin_alpha_j_2*v0*v1*v2/vj_mag_4) + 2*(-x2 + y2)*(-
                                                                                                                                                                                                                          cos_alpha_j*sin_alpha_j*v0_3/vj_mag_3 -
                                                                                                                                                                                                                          cos_alpha_j*sin_alpha_j*v0*v1*v1/vj_mag_3 +
                                                                                                                                                                                                                          cos_alpha_j*sin_alpha_j*v0*v2*v2/vj_mag_3 -
                                                                                                                                                                                                                          cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                                                                                                                                                                          2*sin_alpha_j_2*v0_3/vj_mag_4 +
                                                                                                                                                                                                                          2*sin_alpha_j_2*v0*v1*v1/vj_mag_4 -
                                                                                                                                                                                                                          2*sin_alpha_j_2*v0*v2*v2/vj_mag_4 -
                                                                                                                                                                                                                          2*sin_alpha_j_2*v0/vj_mag_sqr));

    Dphi_cDv1 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(2*(-x0 +
                                                                                                                                                            y0)*(cos_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j*v0*v1*v2/vj_mag_3 -
                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j*v1*v1/vj_mag_3 +
                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j/vj_mag -
                                                                                                                                                                 4*sin_alpha_j_2*v0*v1*v2/vj_mag_4 -
                                                                                                                                                                 sin_alpha_j_2*v1*v1/vj_mag_sqr) + 2*(-x1 + y1)*(-
                                                                                                                                                                                                                 cos_alpha_j_2*v0*v1/vj_mag_sqr +
                                                                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j*v0*v1/vj_mag_3 +
                                                                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j*v1*v1*v2/vj_mag_3 +
                                                                                                                                                                                                                 sin_alpha_j_2*v0*v1/vj_mag_sqr -
                                                                                                                                                                                                                 4*sin_alpha_j_2*v1*v1*v2/vj_mag_4 + 2*sin_alpha_j_2*v2/vj_mag_sqr) + 2*(-x2 + y2)*(-
                                                                                                                                                                                                                                                                                                    cos_alpha_j*sin_alpha_j*v0*v0*v1/vj_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_j*sin_alpha_j*v1_3/vj_mag_3 +
                                                                                                                                                                                                                                                                                                    cos_alpha_j*sin_alpha_j*v1*v2*v2/vj_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_j*sin_alpha_j*v1/vj_mag +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_j_2*v0*v0*v1/vj_mag_4 +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_j_2*v1_3/vj_mag_4 -
                                                                                                                                                                                                                                                                                                    2*sin_alpha_j_2*v1*v2*v2/vj_mag_4 -
                                                                                                                                                                                                                                                                                                    2*sin_alpha_j_2*v1/vj_mag_sqr));

    Dphi_cDv2 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(2*(-x0 +
                                                                                                                                                            y0)*(cos_alpha_j_2*v1*v2/vj_mag_sqr +
                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j*v0*v2*v2/vj_mag_3 -
                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j*v1*v2/vj_mag_3 -
                                                                                                                                                                 4*sin_alpha_j_2*v0*v2*v2/vj_mag_4 + 2*sin_alpha_j_2*v0/vj_mag_sqr -
                                                                                                                                                                 sin_alpha_j_2*v1*v2/vj_mag_sqr) + 2*(-x1 + y1)*(-
                                                                                                                                                                                                                 cos_alpha_j_2*v0*v2/vj_mag_sqr +
                                                                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j*v0*v2/vj_mag_3 +
                                                                                                                                                                                                                 2*cos_alpha_j*sin_alpha_j*v1*v2*v2/vj_mag_3 +
                                                                                                                                                                                                                 sin_alpha_j_2*v0*v2/vj_mag_sqr -
                                                                                                                                                                                                                 4*sin_alpha_j_2*v1*v2*v2/vj_mag_4 + 2*sin_alpha_j_2*v1/vj_mag_sqr) + 2*(-x2 + y2)*(-
                                                                                                                                                                                                                                                                                                    cos_alpha_j*sin_alpha_j*v0*v0*v2/vj_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_j*sin_alpha_j*v1*v1*v2/vj_mag_3 +
                                                                                                                                                                                                                                                                                                    cos_alpha_j*sin_alpha_j*v2_3/vj_mag_3 -
                                                                                                                                                                                                                                                                                                    cos_alpha_j*sin_alpha_j*v2/vj_mag +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_j_2*v0*v0*v2/vj_mag_4 +
                                                                                                                                                                                                                                                                                                    2*sin_alpha_j_2*v1*v1*v2/vj_mag_4 -
                                                                                                                                                                                                                                                                                                    2*sin_alpha_j_2*v2_3/vj_mag_4 + 2*sin_alpha_j_2*v2/vj_mag_sqr));

    Vector3D Dphi_cDvj( Dphi_cDv0, Dphi_cDv1, Dphi_cDv2 );
    return Dphi_cDvj;
}

/*
*Derivative of phi_c wrt xi
 */
Vector3D OPSBody::Dphi_cDxi(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj){
    double u0, u1, u2, v0, v1, v2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, vj_mag_sqr, vj_mag;
    double sin_alpha_i, cos_alpha_i, sin_alpha_j, cos_alpha_j;
    double sin_alpha_i_2, cos_alpha_i_2, sin_alpha_j_2, cos_alpha_j_2, phi_c0;
    double Dphi_cDx0, Dphi_cDx1, Dphi_cDx2;

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    v0 = vj[0]; v1 = vj[1]; v2 = vj[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vj_mag_sqr = v0*v0 + v1*v1 + v2*v2;
    vi_mag = sqrt(vi_mag_sqr);
    vj_mag = sqrt(vj_mag_sqr);
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;
    sin_alpha_j = sin( 0.5*vj_mag );
    cos_alpha_j = cos( 0.5*vj_mag );
    sin_alpha_j_2 = sin_alpha_j*sin_alpha_j;
    cos_alpha_j_2 = cos_alpha_j*cos_alpha_j;

    Dphi_cDx0 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(-
                                                                                                                                                         4*cos_alpha_i*sin_alpha_i*u1/vi_mag -
                                                                                                                                                         4*cos_alpha_j*sin_alpha_j*v1/vj_mag -
                                                                                                                                                         4*sin_alpha_i_2*u0*u2/vi_mag_sqr -
                                                                                                                                                         4*sin_alpha_j_2*v0*v2/vj_mag_sqr);

    Dphi_cDx1 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(4*cos_alpha_i*sin_alpha_i*
                                                                                                                                                         u0/vi_mag + 4*cos_alpha_j*sin_alpha_j*v0/vj_mag -
                                                                                                                                                         4*sin_alpha_i_2*u1*u2/vi_mag_sqr -
                                                                                                                                                         4*sin_alpha_j_2*v1*v2/vj_mag_sqr);

    Dphi_cDx2 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(-2*cos_alpha_i_2 -
                                                                                                                                                         2*cos_alpha_j_2 + 2*sin_alpha_i_2*u0*u0/vi_mag_sqr +
                                                                                                                                                         2*sin_alpha_i_2*u1*u1/vi_mag_sqr - 2*sin_alpha_i_2*u2*u2/vi_mag_sqr +
                                                                                                                                                         2*sin_alpha_j_2*v0*v0/vj_mag_sqr + 2*sin_alpha_j_2*v1*v1/vj_mag_sqr -
                                                                                                                                                         2*sin_alpha_j_2*v2*v2/vj_mag_sqr);

    Vector3D Dphi_cDxi( Dphi_cDx0, Dphi_cDx1, Dphi_cDx2 );
    return Dphi_cDxi;
}

/*
*Derivative of phi_c wrt xj
 */
Vector3D OPSBody::Dphi_cDxj(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj){
    double u0, u1, u2, v0, v1, v2, x0, x1, x2, y0, y1, y2;
    double vi_mag_sqr, vi_mag, vj_mag_sqr, vj_mag;
    double sin_alpha_i, cos_alpha_i, sin_alpha_j, cos_alpha_j;
    double sin_alpha_i_2, cos_alpha_i_2, sin_alpha_j_2, cos_alpha_j_2, phi_c0;
    double Dphi_cDy0, Dphi_cDy1, Dphi_cDy2;

    u0 = vi[0]; u1 = vi[1]; u2 = vi[2];
    v0 = vj[0]; v1 = vj[1]; v2 = vj[2];
    x0 = xi[0]; x1 = xi[1]; x2 = xi[2];
    y0 = xj[0]; y1 = xj[1]; y2 = xj[2];

    vi_mag_sqr = u0*u0 + u1*u1 + u2*u2;
    vj_mag_sqr = v0*v0 + v1*v1 + v2*v2;
    vi_mag = sqrt(vi_mag_sqr);
    vj_mag = sqrt(vj_mag_sqr);
    sin_alpha_i = sin( 0.5*vi_mag );
    cos_alpha_i = cos( 0.5*vi_mag );
    sin_alpha_i_2 = sin_alpha_i*sin_alpha_i;
    cos_alpha_i_2 = cos_alpha_i*cos_alpha_i;
    sin_alpha_j = sin( 0.5*vj_mag );
    cos_alpha_j = cos( 0.5*vj_mag );
    sin_alpha_j_2 = sin_alpha_j*sin_alpha_j;
    cos_alpha_j_2 = cos_alpha_j*cos_alpha_j;

    Dphi_cDy0 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(4*cos_alpha_i*sin_alpha_i*
                                                                                                                                                         u1/vi_mag + 4*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                                                                                                                                                         4*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                                                                                                                                                         4*sin_alpha_j_2*v0*v2/vj_mag_sqr);

    Dphi_cDy1 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(-
                                                                                                                                                         4*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                                                                                                         4*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                                                                                                         4*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                                                                                                         4*sin_alpha_j_2*v1*v2/vj_mag_sqr);

    Dphi_cDy2 = ((-x0 +
                  y0)*(2*cos_alpha_i*sin_alpha_i*u1/vi_mag +
                       2*cos_alpha_j*sin_alpha_j*v1/vj_mag +
                       2*sin_alpha_i_2*u0*u2/vi_mag_sqr +
                       2*sin_alpha_j_2*v0*v2/vj_mag_sqr) + (-x1 + y1)*(-
                                                                       2*cos_alpha_i*sin_alpha_i*u0/vi_mag -
                                                                       2*cos_alpha_j*sin_alpha_j*v0/vj_mag +
                                                                       2*sin_alpha_i_2*u1*u2/vi_mag_sqr +
                                                                       2*sin_alpha_j_2*v1*v2/vj_mag_sqr) + (-x2 + y2)*(cos_alpha_i_2 +
                                                                                                                       cos_alpha_j_2 - sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                       sin_alpha_i_2*u1*u1/vi_mag_sqr + sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                       sin_alpha_j_2*v0*v0/vj_mag_sqr - sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                       sin_alpha_j_2*v2*v2/vj_mag_sqr))*(2*cos_alpha_i_2 +
                                                                                                                                                         2*cos_alpha_j_2 - 2*sin_alpha_i_2*u0*u0/vi_mag_sqr -
                                                                                                                                                         2*sin_alpha_i_2*u1*u1/vi_mag_sqr + 2*sin_alpha_i_2*u2*u2/vi_mag_sqr -
                                                                                                                                                         2*sin_alpha_j_2*v0*v0/vj_mag_sqr - 2*sin_alpha_j_2*v1*v1/vj_mag_sqr +
                                                                                                                                                         2*sin_alpha_j_2*v2*v2/vj_mag_sqr);

    Vector3D Dphi_cDxi( Dphi_cDy0, Dphi_cDy1, Dphi_cDy2 );
    return Dphi_cDxi;
}


}

