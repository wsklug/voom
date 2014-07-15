//----------------------------------------------------------------------
//
//                    William S. Klug, Feng Feng
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//


namespace voom {

	
  //! convert elements' connectivities into pairs <int, int>
  template< class Material_t >
  void TwoPhaseBody< Material_t >::_convertConnectivityToPairs
  ( const tvmet::Vector<int, 3>& v, IntPairsContainer& ipc) {
    //
    //  Connectivity:   2, 3, 4    To pairs: <2,3>, <3,4>, <4, 2>
    //
    //
    ipc.clear();
//     std::pair<int,int> p;

//     p.first = v(0); p.second = v(1); ipc.push_back(p);
//     p.first = v(1); p.second = v(2); ipc.push_back(p);
//     p.first = v(2); p.second = v(0); ipc.push_back(p);
    ipc.push_back(std::make_pair(v(0),v(1)));
    ipc.push_back(std::make_pair(v(1),v(2)));
    ipc.push_back(std::make_pair(v(2),v(0)));
    return;
  }// end converting


  //! reverse <int, int> pairs in a vecotr <1,2> -> <2,1>
  template< class Material_t >
  void TwoPhaseBody< Material_t >::_reversePairs(IntPairsContainer& ipc)
  {
    IntPairsContainer::iterator p = ipc.begin();
    int u;
    for( ; p != ipc.end(); p++ ){
      u = p->first;
      p->first = p->second;
      p->second = u;		
    }
    return;
  } // end reverseParis()


  //! initialize halfedge map
  template< class Material_t >
  void TwoPhaseBody< Material_t >::_initializeMap( MapHalfedge& mh, 
						    MapVertex& mv, 
						    FaceHandleContainer& ufs )
  {
    //
    // initialize MAP
    std::vector<Halfedge_handle> hh(3);         // Triangular element
    std::vector<Vertex_handle> vh(3);
	
    Face_handle fh0;
  
    //ItrVectorInt itrUncheckedElementIndex = ufs.begin();
    IntPairsContainer ipc;
    FaceHandleIterator pfh = ufs.begin();
    //
    // get connectivity for the currently considered face
    const tvmet::Vector<int,3>& c = (*pfh)->getConnectivity();

    _convertConnectivityToPairs( c, ipc);
    //
    //  initialize the halfedge connection in the first element
    //
    //  The Halfedge_handle in HDSFaceConnectivity is pointed to 
    //  the first halfedge_handle in hds.
    //
    //  each halfedge has:
    //               opposite, prev, next, vertex, face
    //  each face has:
    //               halfedge
    //  each vertex has:
    //               halfedge 
    //
    //      IMPORTANT:  vertex->_halfedge()-> vertex() = vertex;
    //                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    assert(_hds.is_valid());		
    //
    // make sure that the vertex handles exist
    assert( mv.find(c[0]) != mv.end() );
    assert( mv.find(c[1]) != mv.end() );
    assert( mv.find(c[2]) != mv.end() );

    // set vertex handle
    vh[0] = mv.find(c[0]) -> second;
    vh[1] = mv.find(c[1]) -> second;
    vh[2] = mv.find(c[2]) -> second;
	
    _hds.insert_halfedge();
    _hds.insert_halfedge();
    _hds.insert_halfedge();

    hh[0] = _hds.halfedges_begin();
    hh[1] = hh[0] + 2;    
    hh[2] = hh[1] + 2;    

    fh0 = _hds.faces_begin();
  
    vh[0] -> halfedge() = hh[0];
    vh[1] -> halfedge() = hh[1];
    vh[2] -> halfedge() = hh[2];

    hh[0] -> prev() = hh[2];
    hh[0] -> next() = hh[1];
    hh[0] -> vertex() = vh[0];
    hh[0] -> face() = fh0;

    hh[1] -> prev() = hh[0];
    hh[1] -> next() = hh[2];
    hh[1] -> vertex() = vh[1];
    hh[1] -> face() = fh0;

    hh[2] -> prev() = hh[1];
    hh[2] -> next() = hh[0];
    hh[2] -> vertex() = vh[2];
    hh[2] -> face() = fh0;

    fh0->halfedge() = hh[0];

    assert(vh[0] -> halfedge() -> vertex() == vh[0]);
    assert(vh[1] -> halfedge() -> vertex() == vh[1]);
    assert(vh[2] -> halfedge() -> vertex() == vh[2]);
    assert(ipc.size() == 3);   // triangle only

    mh.insert( make_pair(ipc[0], hh[0]) );
    mh.insert( make_pair(ipc[1], hh[1]) );
    mh.insert( make_pair(ipc[2], hh[2]) );		
		
    ufs.erase(pfh);
    assert ( _hds.is_valid() );
  }


  //! creating HDS
  template< class Material_t >
  void TwoPhaseBody< Material_t >::_createHDS
  ( ConnectivityContainer & connectivities,
    const int numberOfBoundaries
   ) {
    // ++++++++++++++++++++++++++++++ initialize all necessage information +++++++++++++++++++++++++
    // 
    MapHalfedge mh;
    MapVertex   mv;
    FaceHandleContainer ufs;
    ///
    _hds.clear();
    //
    // reserve the space of HDS 
    //
    // total number of nodes in the current body
    const int numberOfNodes =  _shellNodes.size();
    //
    // one connectivity gives one face
    const int numberOfFaces =  connectivities.size();

    //
    // For triangles connected with each other at least by one edge
    const int numberOfHalfedges = 2 * ( numberOfFaces + numberOfNodes - 2 + numberOfBoundaries);
    //
    _hds.reserve(numberOfNodes, numberOfHalfedges, numberOfFaces);

    assert(_hds.is_valid());

    // insert vertices for all nodes
    for (FeNodeIterator pn=_shellNodes.begin(); pn!=_shellNodes.end(); pn ++)
      _hds.insert_vertex();
    //
    // connect vertices with nodes
    Vertex_handle vh = _hds.vertices_begin();
    //
    for(FeNodeIterator pn=_shellNodes.begin(); 
	vh!=_hds.vertices_end(); vh++,pn++){
      vh -> setNodePointer( *pn );
      // make sure no same nodal indices exist
      assert( mv.find((*pn)->id()) == mv.end() );
      mv.insert(make_pair( (*pn)->id(), vh ));
    }

    assert(_hds.is_valid());
    // initialize faces in HDS
    //
    // allocate memory for face
    ConnectivityContainer::iterator c = connectivities.begin();
    for ( ; c != connectivities.end(); c++)
      _hds.insert_face();
    //
    // connect faces with connectivities
    c = connectivities.begin();
    for (Face_handle h = _hds.faces_begin(); h != _hds.faces_end(); h++, c++){
      h -> setConnectivity( *c );
    }

    //
    // push all faces in the ufs Vector
    
    for (Face_handle h = _hds.faces_begin(); h != _hds.faces_end(); h++)
      ufs.push_back(h);

    _initializeMap(mh, mv, ufs);

    // +++++++++++++++ End initializing all necessage information ++++++++++
    //
    unsigned previousSize = 0;
    std::vector<MapHalfedge::iterator> vctItrMapHalfedge;
    Face_handle fh = _hds.faces_begin();
    MapHalfedge::iterator itrMapHandle;
    FaceHandleIterator itrUncheckedFace;
    //
    //  searching the vector storing the elemnet index, find at least
    //  one then erase that index in each loop;
    //
    //  when the size of this vector is zero, then the halfedge
    //  searching is finished.
    //
    //  in order to avoid the infinite loop, compare the current size
    //  of vector with previous one.
    //
    while ( ufs.size() != 0 ) {
      if( ufs.size() == previousSize ) {
	std::cout << std::endl;
	std::cout << " There is(are) face(s) left, but program can not check its(their) connectivities with the others." << std::endl;
	std::cout << " Possibly, the whole surface is divided into several individual pieces...." << std::endl;
	std::cout << std::endl;
	exit(0);
      }
      
      previousSize = ufs.size();  // avoid infinite loop
      //
	itrUncheckedFace = ufs.begin();
	//
	//  Searching all elements once 
	//
	while (itrUncheckedFace != ufs.end()) {
	  const tvmet::Vector<int, 3>& c = (*itrUncheckedFace)->getConnectivity();
	  IntPairsContainer ipc;
	  _convertConnectivityToPairs(c, ipc); 
	  _reversePairs(ipc);
	  vctItrMapHalfedge.clear();
	  //
	  //   in each face, need to search three pairs <int, int>
	  //
	  IntPairsIterator itr = ipc.begin();
	  for( ; itr != ipc.end(); itr++){
	    itrMapHandle = mh.find( *itr );
	    if(itrMapHandle != mh.end())
	      vctItrMapHalfedge.push_back(itrMapHandle);
	  }
	  //
	  //  check the halfedge connection in another direction if nothing is found in previous searching.
	  //
	    if ( vctItrMapHalfedge.size() == 0 ) {
	      _reversePairs(ipc);
	      IntPairsIterator itr = ipc.begin();
	      for( ; itr != ipc.end(); itr++){
		itrMapHandle = mh.find( *itr );
		if(itrMapHandle != mh.end())
		  vctItrMapHalfedge.push_back(itrMapHandle);
	      }
	    }

	    _reversePairs(ipc);
	    //  
	    //  Finish searching shared edges. results are stored in vector vctItrMapHalfedge.
	    //
	    std::vector<Halfedge_handle> hhdl;
	    if (vctItrMapHalfedge.size() == 0) {
	      itrUncheckedFace ++ ;
	      continue;
	    }
	    else if (vctItrMapHalfedge.size() == 1) {
	      hhdl.clear();
	      hhdl.resize(2);         // Triangular element
	      _hds.insert_halfedge();
	      _hds.insert_halfedge();
	      
	      hhdl[0] = _hds.halfedges_end() - 2;
	      hhdl[1] = _hds.halfedges_end() - 4;
	      
	      const int nodeIndex = _newNodalIndex(vctItrMapHalfedge[0], c);
	      if (nodeIndex == -1 ) {
		std::cout << "Check fuction -- newNodalIndex()." <<std::endl;
		exit(0);
	      }
	      //
	      assert(mv.find(nodeIndex) != mv.end());					
	      Vertex_handle   v0 = mv.find(nodeIndex) -> second;
	      
	      
	      Halfedge_handle hh = (vctItrMapHalfedge[0] -> second ) -> opposite();
	      
	      Face_handle f0 = *itrUncheckedFace;
	      //
	      //  hh     // the existing halfedge
	      hh -> prev() = hhdl[0];
	      hh -> next() = hhdl[1];
	      hh -> vertex() = ((hh -> opposite())  -> next()) -> vertex();
	      hh -> face() = f0;
	      //
	      //  h0     // the new halfedge:first one
	      hhdl[0] -> prev() = hhdl[1];
	      hhdl[0] -> next() = hh;
	      hhdl[0] -> vertex() = v0;
	      hhdl[0] -> face() = f0;
	      //
	      //  
	      v0 -> halfedge() = hhdl[0];
	      //
	      //  h1    // another new halfedge:second one
	      hhdl[1] -> prev() = hh;
	      hhdl[1] -> next() = hhdl[0];
	      hhdl[1] -> vertex() = (hh -> opposite()) -> vertex();
	      hhdl[1] -> face() = f0;
	      
	      // getNthElementInIndex(*itrUncheckedElementIndex).setHalfedgeHandle(hh);
	      f0 -> halfedge() = hh;
	      //
	      //  make sure that the vertex of the halfedge of current vertex
	      //  is pointed to current vertex
	      assert( v0 -> halfedge() -> vertex() == v0);
	      assert( _hds.is_valid() );
	    }
	    else if (vctItrMapHalfedge.size() == 2) {
	      Halfedge_handle hh1, hh2;
	      hhdl.clear();
	      hhdl.resize(1);  
	      _hds.insert_halfedge();
	      hhdl[0] = _hds.halfedges_end() - 2;
	      
	      Face_handle f1 = *itrUncheckedFace;
	      
	      hh1 = (vctItrMapHalfedge[0] -> second ) -> opposite();
	      hh1->vertex() = hh1->opposite()->next() -> vertex();
	      
	      hh2 = (vctItrMapHalfedge[1] -> second ) -> opposite();
	      hh2->vertex() = hh2-> opposite() -> next() -> vertex();
	      
	      if( hh1 -> opposite()-> vertex() == hh2 -> vertex()) {
		//
		//  hh1   //  existing halfedge 1
		hh1->next() = hh2; 
		hh1->prev() = hhdl[0]; 
		hh1->face() = f1;
		//
		//  hh2   //  existing halfedge 2
		hh2->prev() = hh1; 
		hh2->next() = hhdl[0]; 
		hh2->face() = f1;
		//
		//  hhdl[0]
		hhdl[0]->prev() = hh2;
		hhdl[0]->next() = hh1;
		hhdl[0]->vertex() = hh2 -> opposite() -> vertex();
		hhdl[0]->face() = f1;
	      }
	      else {
		//  hh1
		hh1->next() = hhdl[0];  hh1->prev() = hh2; hh1->face() = f1;
		//  hh2
		hh2->prev() = hhdl[0];  hh2->next() = hh1; hh2->face() = f1;
		// hhdl[0]
		hhdl[0]->prev() = hh1;
		hhdl[0]->next() = hh2;
		hhdl[0]->vertex() = hh1 -> opposite() -> vertex();
		hhdl[0]->face() = f1;
	      }
	      f1->halfedge() = hhdl[0];
	      assert( _hds.is_valid() );
	    }
	    else if (vctItrMapHalfedge.size() == 3) {
	      Halfedge_handle hh0, hh1, hh2;
	      Face_handle f0 = *itrUncheckedFace;
	      
	      hh0 = vctItrMapHalfedge[0] -> second -> opposite();
	      hh0->vertex() = hh0 -> opposite() -> next() -> vertex();
	      
	      hh1 = vctItrMapHalfedge[1] -> second -> opposite();
	      hh1->vertex() = hh1 -> opposite() -> next() -> vertex();
	      
	      hh2 = vctItrMapHalfedge[2] -> second -> opposite();
	      hh2->vertex() = hh2 -> opposite() -> next() -> vertex();
	      
	      if( hh0 -> opposite() -> vertex() == hh1 -> vertex()) {
		hh0 -> prev() = hh2;  hh0 -> next() = hh1; hh0 -> face() = f0;
		hh1 -> prev() = hh0;  hh1 -> next() = hh2; hh1 -> face() = f0;
		hh2 -> prev() = hh1;  hh2 -> next() = hh0; hh2 -> face() = f0;
	      }
	      else {
		hh0 -> prev() = hh1;  hh0 -> next() = hh2; hh0 -> face() = f0;
		hh1 -> prev() = hh2;  hh1 -> next() = hh0; hh1 -> face() = f0;
		hh2 -> prev() = hh0;  hh2 -> next() = hh1; hh2 -> face() = f0;
	      }
	      f0 -> halfedge() = hh0;
	      assert( _hds.is_valid() );
	    }
	    else
	      std::cout << "Error in halfedge assignment." << std::endl;

	    //
	    //  if find one element, erase its index from the vector.
	    if (vctItrMapHalfedge.size() != 0 ) ufs.erase(itrUncheckedFace);
	    
	    //
	    //  reassign the MAP structure with new pairs
	    //
	    //  add new Pairs first.
	    std::vector<Halfedge_handle>::iterator hhi = hhdl.begin();
	    for ( ; hhi != hhdl.end(); hhi++) {
	      int iFirst, iSecond;
	      iFirst  = (*hhi) -> vertex() -> getNodePointer() -> id();
	      iSecond = (*hhi) -> next() -> vertex() -> getNodePointer() -> id();
	      std::pair<int, int> ip(iFirst, iSecond);
	      mh.insert(make_pair(ip, *hhi));
	    }
	    //
	    //  delete used pairs
	    std::vector<MapHalfedge::iterator>::iterator imhi 
	      = vctItrMapHalfedge.begin();
	    for ( ; imhi != vctItrMapHalfedge.end(); imhi++) {
	      mh.erase( *imhi );
	    }
	    
	}
	
    }
  } // end creating HDS
  

  /*! find the third nodal index in a triangle element if have already
   *  known others
  */
  template< class Material_t >
  int TwoPhaseBody< Material_t >::_newNodalIndex
  ( const MapHalfedge::iterator itrMapHandle,
    const tvmet::Vector<int, 3>& c
    ) {
    for (int i = 0; i < int( c.size()); i ++) {
      if ( c[i] != (itrMapHandle -> first).first &&
	   c[i] != (itrMapHandle -> first).second)
	return c[i];
    }
    return -1;
  }


//   //! initialize neighbouring nodes for each node using HDS
//   template < class Material_t >
//   void LoopShellBody< Material_t >::_initNodeNeighbors()
//   {
//     std::vector<DeformationNode<3>*> neighbors;
//     Vertex_handle vh = _hds.vertices_begin();
//     for ( ; vh != _hds.vertices_end(); vh ++) {
//       assert( vh -> halfedge() -> vertex() == vh );
//       //
//       Halfedge_handle h0 = vh->halfedge()-> opposite();
//       neighbors.push_back( h0 -> vertex() -> getNodePointer() );
//       Halfedge_handle h1 = h0 -> next() -> opposite();
//       while (h1 != h0 ) {
// 	neighbors.push_back( h1 -> vertex() -> getNodePointer() );
// 	h1 = h1 -> next() -> opposite();
//       } // end while
//       vh -> getNodePointer() -> setNeighbouringNodes(neighbors);
//       vh -> getNodePointer() -> initLimitSurfaceReferencePosition();
//     } // end for
//   } // end _initNodeNeighbors
	

  //! generate necessary informations for creating a new element 
  template< class Material_t >
  void TwoPhaseBody< Material_t >::_createElement
  ( Material_t material, 
    const Face_handle fh, 
    const unsigned quadOrder
    )
  {
    typename TwoPhaseElement<Material_t>::NodeContainer nds;
    CornerValences v;
    
    nds.clear();

    // for initialize elements
    // int v(0), v(1), v(2);
    // total vertices around one vertex of the currently considered element. 
      
    v(0) = 1;
    v(1) = 1;
    v(2) = 1;
      
    Halfedge_handle h0 = fh -> halfedge();
      
    Halfedge_handle h1 = h0 -> opposite() -> next();
    while ( h0 != h1 ) {
      v(0) ++;
      h1 = h1 -> opposite() -> next();
    }
      
    h0 = h0 -> next();
    h1 = h0 -> opposite() -> next();
    while ( h0 != h1 ) {
      v(1) ++;
      h1 = h1 -> opposite() -> next();
    }
    
      
    h0 = h0 -> next();
    h1 = h0 -> opposite() -> next();
    while ( h0 != h1 ) {
      v(2) ++;
      h1 = h1 -> opposite() -> next();
    }
      
    h0 = fh -> halfedge();
    //
    // pointer version
    XCNode<3>* 
    pnode = (FeNode_t*) (h0 -> vertex() -> getNodePointer() );
    nds.push_back(pnode);
    pnode = (FeNode_t*) (h0 -> next() -> vertex() -> getNodePointer() );
    nds.push_back(pnode);
    pnode = (FeNode_t*) (h0 -> prev() -> vertex() -> getNodePointer() );
    nds.push_back(pnode);
    //
    //  pointer version
    //
    h1 = h0 -> opposite() ->prev();
    pnode = (FeNode_t*) ( h1 -> vertex() -> getNodePointer() );
    if ( find(nds.begin(), nds.end(), pnode) == nds.end() )
      nds.push_back(pnode);
      
    for (unsigned i = 0; i < v(1)-4; i ++) {
      h1 = h1 ->opposite() -> prev();
      pnode = (FeNode_t*) ( h1 -> vertex() -> getNodePointer() );
      nds.push_back(pnode);
    }
      
    h0 = h0 -> next();
    h1 = h0 -> opposite() ->prev();
    pnode = (FeNode_t*) ( h1 -> vertex() -> getNodePointer() );
    if ( find(nds.begin(), nds.end(), pnode) == nds.end() )
      nds.push_back(pnode);
      
    for (unsigned i = 0; i < v(2)-4; i ++) {
      h1 = h1 ->opposite() -> prev();
      pnode = (FeNode_t*) ( h1 -> vertex() -> getNodePointer() );
      nds.push_back(pnode);
    }
      
    h0 = h0 -> next();
    h1 = h0 -> opposite() ->prev();
    pnode = (FeNode_t*) ( h1 -> vertex() -> getNodePointer() );
    if ( find(nds.begin(), nds.end(), pnode) == nds.end() )
      nds.push_back(pnode);
      
    for (unsigned i = 0; i < v(0)-4; i ++) {
      h1 = h1 ->opposite() -> prev();
      pnode = (FeNode_t*) ( h1 -> vertex() -> getNodePointer() );
      nds.push_back(pnode);
    }
            
    if (false) {
      std::cout << std::endl;
      for (int i = 0; i < 3; i++) {
	std::cout << "Vertex ["
		  << i
		  << "] has "
		  << v(i)
		  << " neighbouring vertices."
		  << std::endl;
      }
	  
      for (FeNodeIterator p = nds.begin(); p != nds.end(); p ++) {
	int dist = distance(nds.begin(), p);
	std::cout << "The "
		  << std::setw(5)
		  << dist
		  << "th vertex is " 
		  << std::setw(5)
		  << (*p)->id()
		  << "th Node"
		  << std::endl;
      }
      std::cout << "------------------------- end ----------------------------"
		<< std::endl;
    }

    const unsigned npn = v(0) + v(1) + v(2) - 6;		
    FeElement_t * elem 
      = new FeElement_t( TriangleQuadrature(quadOrder), 
			 material,
			 nds,
			 v,
			 _pressureNode,
			 _tensionNode,
			 _chemicalTensionNode,
			 _volumeConstraint,
			 _areaConstraint,
			 _areaOneConstraint
			 );
    _shells.push_back(elem);

    return;
  } // end initializing element's information
	
	
}; // namespace voom
