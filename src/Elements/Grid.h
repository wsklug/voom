// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Andrew R. Missel
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__Grid_h__)
#define __Grid_h__

#include "VoomMath.h"
#include "Node.h"
#include "PeriodicBox.h"

using namespace std;
using namespace voom;

namespace voom
{
  template<class Elem,class PosClass,int N>
  class Grid {
  public:
    
    typedef tvmet::Vector<double,N> VectorND;
    typedef tvmet::Vector<int,N> VectorNI;
    typedef typename std::vector<int> IndexList;
    typedef typename std::vector<Elem*> ElemList;
    typedef typename ElemList::iterator ElemListIterator;
    typedef typename std::set<Elem*> ElemContainer;
    typedef typename ElemContainer::iterator ElemIterator;
    typedef typename std::vector<ElemContainer> ElemBoxes;
    typedef typename ElemBoxes::iterator ElemBoxesIterator;
    typedef typename std::map<Elem*,ElemContainer> ElemPairContainer;
    typedef typename ElemPairContainer::iterator ElemPairIterator;

    typedef const VectorND & (PosClass::*PosFunc)();

    Grid() {;}

    // must specify a (proposed) grid spacing and pass in a box //
    Grid(PeriodicBox* box, VectorND & gridSpace, PosFunc pf, bool computeNeighbors=false) : _box(box), getElemPos(pf), _computeNeighbors(computeNeighbors) {
      getGridSpacing(gridSpace);
      int nBoxes = 1;
      for(int i=0; i<N; i++) {
	nBoxes *= _nBoxes[i];
      }
      _elemBoxes.resize(nBoxes);
    }

    // can also pass in vector of elements upon grid creation //
    Grid(PeriodicBox* box, VectorND & gridSpace, ElemList & elems, PosFunc pf, bool computeNeighbors=false) : _box(box), getElemPos(pf), _computeNeighbors(computeNeighbors) {
      getGridSpacing(gridSpace);
      int nBoxes = 1;
      for(int i=0; i<N; i++) {
	nBoxes *= _nBoxes[i];
      }
      _elemBoxes.resize(nBoxes);
      addElems(elems);
    }
    
    ~Grid() {}
    
    void setPosFunc(PosFunc pf) {
      getElemPos = pf;
    }

    void setBox(PeriodicBox* box) { _box = box; }

    void setGridSpace(VectorND & gridSpace) {
      getGridSpacing(gridSpace);
      int nBoxes = 1;
      for(int i=0; i<N; i++) {
	nBoxes *= _nBoxes[i];
      }
      _elemBoxes.resize(nBoxes);
    }

    void setComputeNeighbors(bool computeNeighbors) {
      _computeNeighbors = computeNeighbors;
    }

    // adds an element to the grid //
    void addElem(Elem* e) {
      VectorND epos;
      //getElemPos(e,epos);
      epos = (e->*getElemPos)();
      VectorNI ecoords;
      getGridCoords(epos,ecoords);
      int idx = getBoxIndex(ecoords);

      if(_computeNeighbors) {
	// update double lookup table //
	IndexList idxlst;
	getNeighborIndices(ecoords,idxlst);
	idxlst.push_back(idx);
	for(int i=0; i<idxlst.size(); i++) {
	  ElemContainer & ec = _elemBoxes[idxlst[i]];
	  for(ElemIterator ei=ec.begin(); ei!=ec.end(); ei++) {
	    _elemPairs[*ei].insert(e);
	    _elemPairs[e].insert(*ei);
	  }
	}
      }

      // add element to grid //
      _elemBoxes[idx].insert(e);
    }

    // adds a set of elements //
    void addElems(ElemContainer & es) {
      for(ElemIterator ei=es.begin(); ei!=es.end(); ei++) {
	addElem(*ei);
      }
    }
    
    // adds a vector of elements //
    void addElems(ElemList & es) {
      for(ElemListIterator ei=es.begin(); ei!=es.end(); ei++) {
	addElem(*ei);
      }
    }

    // removes an element //
    void removeElem(Elem* e) {
      VectorND epos;
      //getElemPos(e,epos);
      epos = (e->*getElemPos)();
      VectorNI ecoords;
      getGridCoords(epos,ecoords);
      int idx = getBoxIndex(ecoords);

      // look in the grid square it should be in //
      if(_elemBoxes[idx].find(e) != _elemBoxes[idx].end()) {
	_elemBoxes[idx].erase(e);
      }   
      // if it's not there, just search through all other boxes to find it //
      else {
	bool erased = false;
	ElemBoxesIterator ebi = _elemBoxes.begin();
	while(!erased && ebi!=_elemBoxes.end()) {
	  if(ebi->find(e) != ebi->end()) {
	    ebi->erase(e);
	    erased = true;
	  }
	}
	if(!erased) std::cout << "Grid: sorry, could not find that element in this grid!" << std::endl;
      }

      // now erase neighbors from list //
      for(ElemIterator ei=_elemPairs[e].begin(); ei!=_elemPairs[e].end(); ei++) {
	_elemPairs[*ei].erase(e);
      }
      _elemPairs[e].clear();
    }
    
    // clears out the grid and neighbor lists and resets everything //
    void resetGrid() {
      ElemList tmpEC;
      for(ElemBoxesIterator ebi=_elemBoxes.begin(); ebi!=_elemBoxes.end(); ebi++) {
	for(ElemIterator ei=ebi->begin(); ei!=ebi->end(); ei++) {
	  tmpEC.push_back(*ei);
	}
	ebi->clear();
      }
      _elemBoxes.clear();
      for(ElemPairIterator epi=_elemPairs.begin(); epi!=_elemPairs.end(); epi++) {
	epi->second.clear();
      }
      _elemPairs.clear();
      addElems(tmpEC);
    }

    void resetGridSpace(VectorND & gridSpace) {
      getGridSpacing(gridSpace);
      int nBoxes = 1;
      for(int i=0; i<N; i++) {
	nBoxes *= _nBoxes[i];
      }
      ElemContainer tmpEC;
      for(ElemBoxesIterator ebi=_elemBoxes.begin(); ebi!=_elemBoxes.end(); ebi++) {
	tmpEC.insert(ebi->begin(),ebi->end());
	ebi->clear();
      }
      _elemBoxes.resize(nBoxes);

      for(ElemPairIterator epi=_elemPairs.begin(); epi!=_elemPairs.end(); epi++) {
	epi->second.clear();
      }
      _elemPairs.clear();

      addElems(tmpEC);
    }

    void resetGridSpace(int changeNBoxes) {
      
      
    }

    // resets the grid's box and re-inserts the elements correctly // 
    void resetBox(PeriodicBox* box) {
      _box = box;
      getGridSpacing(_gridSpace);
      ElemContainer tmpEC;
      for(ElemBoxesIterator ebi=_elemBoxes.begin(); ebi!=_elemBoxes.end(); ebi++) {
	tmpEC.insert(ebi->begin(),ebi->end());
	ebi->clear();
      }
      _elemBoxes.clear();
      int nBoxes = 1;
      for(int i=0; i<N; i++) {
	nBoxes *= _nBoxes[i];
      }
      _elemBoxes.resize(nBoxes);
      for(ElemPairIterator epi=_elemPairs.begin(); epi!=_elemPairs.end(); epi++) {
	epi->second.clear();
      }
      _elemPairs.clear();
      addElems(tmpEC);
    }

    // clears out all elements from grid //
    void clearGrid() {
      for(ElemBoxesIterator ebi=_elemBoxes.begin(); ebi!=_elemBoxes.end(); ebi++) {
	ebi->clear();
      }
      _elemBoxes.clear();
      for(ElemPairIterator epi=_elemPairs.begin(); epi!=_elemPairs.end(); epi++) {
	epi->second.clear();
      }
      _elemPairs.clear();
    }

    
    // this is the most important function: checks to make sure elements are in the right place and corrects errors //
    void correctGrid() {
      // first go through and make sure everything is in the right place //
      for(int ebi=0; ebi<_elemBoxes.size(); ebi++) {
	ElemIterator ei=_elemBoxes[ebi].begin();
	while(ei!=_elemBoxes[ebi].end()) {
	  VectorND epos;
	  //getElemPos(*ei,epos);
	  epos = ((*ei)->*getElemPos)();
	  VectorNI ecoords;
	  getGridCoords(epos,ecoords);
	  int idx = getBoxIndex(ecoords);
	  // if element is in the wrong box, fix it //
	  if(idx != ebi) {
	    Elem* tmpElem = *ei;
	    ei++;
	    _elemBoxes[ebi].erase(tmpElem);
	    _elemBoxes[idx].insert(tmpElem);

	    // now remove pairs that involve this element //
	    if(_elemPairs.find(tmpElem) != _elemPairs.end()) {
	      for(ElemIterator eci=_elemPairs[tmpElem].begin(); eci!=_elemPairs[tmpElem].end(); eci++) {
		_elemPairs[*eci].erase(tmpElem);
	      }
	      _elemPairs[tmpElem].clear();
	    }

	    // now add pairs //
	    IndexList nidxs;
	    getNeighborIndices(ecoords,nidxs);
	    for(ElemIterator nei=_elemBoxes[idx].begin(); nei!=_elemBoxes[idx].end(); nei++) {
	      if(*nei != tmpElem) {
		_elemPairs[*nei].insert(tmpElem);
		_elemPairs[tmpElem].insert(*nei);
	      }
	    }
	    for(int j=0; j<nidxs.size(); j++) {
	      for(ElemIterator nei=_elemBoxes[nidxs[j]].begin(); nei!=_elemBoxes[nidxs[j]].end(); nei++) {
		_elemPairs[*nei].insert(tmpElem);
		_elemPairs[tmpElem].insert(*nei);
	      }
	    }
	  }
	  else ei++;
	}
      }      
    }

    // this is the money function: returns a reference to a container full of all pairs of neighbor elements //
    ElemPairContainer & getNeighbors() {
      correctGrid();
      return _elemPairs;
    }

    // constructs and returns an element pair container with all pairs containing an element from a particular box //
    ElemPairContainer getNeighbors(int i) {
      VectorNI coords;
      getGridCoords(i,coords);
      IndexList idx;
      getNeighborIndices(coords,idx);
      ElemPairContainer tmpEC;
      ElemContainer & ec = _elemBoxes[i];
      for(ElemIterator ei=ec.begin(); ei!=ec.end(); ei++) {
	for(int j=0; j<idx.size(); j++) {
	  ElemContainer & ec2 = _elemBoxes[idx[j]];
	  for(ElemIterator ei2=ec2.begin(); ei2!=ec2.end(); ei2++) {
	    tmpEC[*ei].insert(*ei2);
	  }
	}
	for(ElemIterator ei3=ec.begin(); ei3!=ec.end(); ei3++) {
	  if(*ei != *ei3) tmpEC[*ei].insert(*ei3);
	}
      }

      return tmpEC;
    }
    
    ElemList getBoxNeighbors(int i) {
      VectorNI coords;
      getGridCoords(i,coords);
      IndexList idx;
      getNeighborIndices(coords,idx);
      //idx.push_back(i);
      ElemList boxNeighbs;
      for(int j=0; j<idx.size(); j++) {
	ElemContainer & ec = _elemBoxes[idx[j]];
	for(ElemIterator ei=ec.begin(); ei!=ec.end(); ei++) {
	  boxNeighbs.push_back(*ei);
	}
      }
      
      return boxNeighbs;
    }

    ElemContainer & getBoxElems(int i) {
      return _elemBoxes[i];
      
    }

    // function to get element's position //
    //void getElemPos(Elem* e, VectorND & epos) {
    //  epos = e->point();
    //}

    // given an element's position, stores its grid coordinates in the second argument //
    void getGridCoords(VectorND & epos, VectorNI & ecoords) {
      _box->mapPoint(epos);
      for(int i=0; i<N; i++) {
	ecoords[i] = (int)(floor(epos[i]/_gridSpace[i])+.5);
      }
    }

    void getGridCoords(int indx, VectorNI & ecoords) {
      int tmpindex = indx;
      for(int i=0; i<N; i++) {
	int curboxmults = _boxMults[i];
 	int curind = indx/curboxmults;
	//std::cout << "index = " << indx << ", curboxmults = " << curboxmults << ", division = " << indx/curboxmults << std::endl;
 	indx = indx - curind*curboxmults;
 	ecoords[i] = curind;
      }

      //std::cout << "boxMults = " << _boxMults << ", index = " << indx << ", ecoords = " << ecoords << std::endl;

      // TEST //
      if(getBoxIndex(ecoords) != tmpindex) {
	std::cout << "Error: box indices not mapped to position correctly." << std::endl;
	std::cout << "Index = " << tmpindex << ", ecoords = " << ecoords << ", mapped index = " << getBoxIndex(ecoords) << std::endl;
      }
    }

    // given an element's grid coordinates, returns the index at which its square is found in the container //
    int getBoxIndex(VectorNI & gcs) {
      int indx = 0;
      for(int i=0; i<N; i++) {
	indx += gcs[i]*_boxMults[i];
      }
      return indx;
    }

    // given an element's grid coordinates, stores a list of indices of adjacent grid squares in the second argument //
    void getNeighborIndices(VectorNI & gcs, IndexList & idxs) {
      if(N==2) {
	idxs.clear();
	idxs.reserve(8);
	for(int ix=-1;ix<2;ix++) {
	  for(int iy=-1;iy<2;iy++) {
	    if(ix!=0 || iy!=0) {
	      VectorNI tmpInd;
	      tmpInd[0] = (gcs[0]+ix+_nBoxes[0])%_nBoxes[0];
	      tmpInd[1] = (gcs[1]+iy+_nBoxes[1])%_nBoxes[1];
	      idxs.push_back(getBoxIndex(tmpInd));
	    }
	  }
	}
      }
      
      else if(N==3) {

      }

      else {
	std::cerr << "We are not interested in your physically irrelevant dimensions, mathematician." << std::endl;
      }
    }

    void getNeighborIndices(int indx, IndexList & idxs) {
      if(N==2) {
	idxs.clear();
	idxs.reserve(8);
	VectorNI gcs;
	getGridCoords(indx,gcs);
	for(int ix=-1;ix<2;ix++) {
	  for(int iy=-1;iy<2;iy++) {
	    if(ix!=0 || iy!=0) {
	      VectorNI tmpInd;
	      tmpInd[0] = (gcs[0]+ix+_nBoxes[0])%_nBoxes[0];
	      tmpInd[1] = (gcs[1]+iy+_nBoxes[1])%_nBoxes[1];
	      idxs.push_back(getBoxIndex(tmpInd));
	    }
	  }
	}
      }
      
      else if(N==3) {

      }

      else {
	std::cerr << "We are not interested in your physically irrelevant dimensions, mathematician." << std::endl;
      }
    }

    // given a proposed grid spacing, finds a close number that allows for perfect tiling of the system //
    void getGridSpacing(VectorND & gridspace) {
      const VectorND & size = _box->size();
      int nbxs = 1;
      for(int i=0;i<N;i++) {
	int nb = (int)(floor(size[i]/gridspace[i])+.5);
	if(nb%2 !=0) nb--;
	_nBoxes[i] = nb;
	_gridSpace[i] = size[i]/nb;
	nbxs *= nb;
      }

      for(int i=0;i<N;i++) {
	nbxs /= _nBoxes[i];
	_boxMults[i] = nbxs;
	//std::cout << "_boxMults[" << i << "] = " << _boxMults[i] << std::endl;
      }
    }

    VectorND & gridSpace() { return _gridSpace; }

    std::vector< VectorND > gridCorners(int k) {
      std::vector< VectorND > corners;
      VectorNI ecoords;
      getGridCoords(k,ecoords);
      
      if(N==2) {
	for(int i=0; i<2; i++) {
	  for(int j=0; j<2; j++) {
	    VectorND edgeCoords;
	    edgeCoords[0] = (ecoords[0]+i)*_gridSpace[0];
	    edgeCoords[1] = (ecoords[1]+j)*_gridSpace[1];
	    corners.push_back(edgeCoords);
	  }
	}
      }

      return corners;
      
    }

    int elementCount() {
      int elems = 0;
      for(ElemBoxesIterator ebi=_elemBoxes.begin(); ebi!=_elemBoxes.end(); ebi++) {
	elems += ebi->size();
      }
      return elems;
    }

    ElemBoxes & elemBoxes() {
      correctGrid();
      return _elemBoxes; 
    }

    ElemList allElems() {
      ElemList alles;
      for(ElemBoxesIterator ebi=_elemBoxes.begin(); ebi!=_elemBoxes.end(); ebi++) {
	for(ElemIterator ei=ebi->begin(); ei!=ebi->end(); ei++) {
	  alles.push_back(*ei);
	}
      }
      return alles;
    }

    int nBoxes() { return _elemBoxes.size(); }

    VectorNI & nBoxesDim() { return _nBoxes; }

  private:
    VectorND _gridSpace;
    VectorNI _nBoxes;
    VectorNI _boxMults;
    PeriodicBox* _box;
    ElemBoxes _elemBoxes;

    bool _computeNeighbors;
    ElemPairContainer _elemPairs;

    PosFunc getElemPos;

  };

//   template<int N>
//   void Grid<Node,N>::getElemPos(Node * n, tvmet::Vector<double,N> & epos) {
//     int nNodes = f->nodes.size();
//     epos = f->nodes[0]->point() + f->nodes[nNodes-1]->point();
//     epos /= 2.0;
//   }

//   template<>
//   void Grid<SemiflexibleGel<3>::Filament,3>::getElemPos(SemiflexibleGel<3>::Filament * f, tvmet::Vector<double,3> & epos) {
//     int nNodes = f->nodes.size();
//     epos = f->nodes[0]->point() + f->nodes[nNodes-1]->point();
//     epos /= 2.0;
//   }

//   template<int N>
//   void Grid<SemiflexibleGel<N>::Filament,N>::getElemPos(SemiflexibleGel<N>::Filament * f, VectorND & epos) {
//     int nNodes = f->nodes.size();
//     epos = f->nodes[0]->point() + f->nodes[nNodes-1]->point();
//     epos /= 2.0;
//   }
  
}

#endif //__Grid_h__
