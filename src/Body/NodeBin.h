#ifndef __NodeBin__
#define __NodeBin__

#include<vector>
#include<iostream>
//#include"NodeBin.cc"
#include<limits>
namespace voom
{
  
  class NodeBin
  {
  public:
    //Constructor/destructors:
    
    //!Default constructor
    NodeBin(){}
    
    //!Construct from a single grid size and a set of data
// wrong implementation    NodeBin(double size, std::vector<std::vector<double> >& in_points){NodeBin(size,size,size,in_points);}
    
    //!Construct from three grid sizes in three directions
//    NodeBin(double dx, double dy, double dz, std::vector<std::vector<double> >& in_points);

NodeBin(double dx, double dy, double dz, const std::vector< DeformationNode<3>* >& Nodes)
  {  //std::cout<<"in the function"<<dx<<"\t"<<dy<<"\t"<<dz<<std::endl;
    
    _dx=dx; _dy=dy; _dz=dz;

    //create the vector of points from nodes
    int nPoints=Nodes.size();
    std::vector<std::vector<double> > in_points;
    in_points.reserve(nPoints);

    for(int i=0;i<nPoints;i++){
      std::vector<double> temp(3);
      temp[0]=Nodes[i]->getPosition(0);
      temp[1]=Nodes[i]->getPosition(1);
      temp[2]=Nodes[i]->getPosition(2);
      in_points.push_back(temp);
    }

    //find minimum and maximum of in_points
    _minx = _miny = _minz = std::numeric_limits<double>::max();
    _maxx = _maxy = _maxz = std::numeric_limits<double>::min();
    for(int i=0;i<in_points.size();i++){
      
      if(_minx > in_points[i][0]) _minx = in_points[i][0];
      if(_miny > in_points[i][1]) _miny = in_points[i][1];
      if(_minz > in_points[i][2]) _minz = in_points[i][2];
      
      if(_maxx < in_points[i][0]) _maxx = in_points[i][0];
      if(_maxy < in_points[i][1]) _maxy = in_points[i][1];
      if(_maxz < in_points[i][2]) _maxz = in_points[i][2];
      
    }
    //increase the size by dx,dy and dz at every corner
    _minx=_minx-dx; _miny=_miny-dy; _minz=_minz-dz;
    _maxx=_maxx+dx; _maxy=_maxy+dy; _maxz=_maxz+dz;
    
    //calculate number of grids
    _nx = int((_maxx-_minx)/dx);
    _ny = int((_maxy-_miny)/dy);
    _nz = int((_maxz-_minz)/dz);
    
    //reserve the size
    _grid.resize(_nx*_ny*_nz);
    
    for(int i=0;i<in_points.size();i++){
      int ix = int((in_points[i][0] - _minx)/dx);
      int iy = int((in_points[i][1] - _miny)/dy);
      int iz = int((in_points[i][2] - _minz)/dz);
      _grid[ix+iy*_nx+iz*_nx*_ny].push_back(i);
    }
//std::cout<<"values are "<<_nx<<"\t"<<_ny<<"\t"<<_nz<<"\n"<<_dx<<"\t"<<_dy<<"\t"<<_dz<<std::endl;
  }
    
    //!Construct from three number of cells in three directions
  //  NodeBin(int nx, int ny, int nz, std::vector<std::vector<double> >& in_points);

NodeBin(int nx, int ny, int nz, std::vector<std::vector<double> >& in_points)
  {
    
    _nx=nx; _ny=ny; _nz=nz;
    //find minimum and maximum of in_points
    _minx = _miny = _minz = std::numeric_limits<double>::max();
    _maxx = _maxy = _maxz = std::numeric_limits<double>::min();
    for(int i=0;i<in_points.size();i++){
      
      if(_minx > in_points[i][0]) _minx = in_points[i][0];
      if(_miny > in_points[i][1]) _miny = in_points[i][1];
      if(_minz > in_points[i][2]) _minz = in_points[i][2];
      
      if(_maxx < in_points[i][0]) _maxx = in_points[i][0];
      if(_maxy < in_points[i][1]) _maxy = in_points[i][1];
      if(_maxz < in_points[i][2]) _maxz = in_points[i][2];
      
    }
    
    //calculate spacing of grids
    _dx = (_maxx-_minx)/_nx;
    _dy = (_maxy-_miny)/_ny;
    _dz = (_maxz-_minz)/_nz;
    
    //add 2 to nx, ny and nz
    _nx+=2; _ny+=2; _nz+=2; 
    //and increase the size by dx,dy and dz at every corner
    _minx=_minx-_dx; _miny=_miny-_dy; _minz=_minz-_dz;
    _maxx=_maxx+_dx; _maxy=_maxy+_dy; _maxz=_maxz+_dz;

    //reserve the size
    _grid.resize(_nx*_ny*_nz);
    
    for(int i=0;i<in_points.size();i++){
      int ix = int((in_points[i][0] - _minx)/_dx);
      int iy = int((in_points[i][1] - _miny)/_dy);
      int iz = int((in_points[i][2] - _minz)/_dz);
      _grid[ix+iy*_nx+iz*_nx*_ny].push_back(i);
    }
  }
    
    //!Get back the list of int for a given point based on it's location
  //  std::vector<int>& GetList(std::vector<double>& point);

//const std::vector<int>& GetList(const std::vector<double>& point){
const std::vector<int>& GetList(const MFShape::CoordinateArray & point){
    //calculate the box number for point
    int ix=int((point(0)-_minx)/_dx);
    int iy=int((point(1)-_miny)/_dy);
    int iz=int((point(2)-_minz)/_dz);
//std::cout<<"values are "<<_nx<<"\t"<<_ny<<"\t"<<_nz<<"\n"<<ix<<"\t"<<iy<<"\t"<<iz<<std::endl;
    _returnlist.clear();
    bool flag=false;
    for(int ii=std::max(0,ix-1);ii<=std::min(ix+1,_nx-1);ii++)
      for(int jj=std::max(0,iy-1);jj<=std::min(iy+1,_ny-1);jj++)
	for(int kk=std::max(0,iz-1);kk<=std::min(iz+1,_nz-1);kk++){
          flag=true;
	  for(int mm=0;mm<_grid[ii+jj*_nx+kk*_nx*_ny].size();mm++)
	    _returnlist.push_back(_grid[ii+jj*_nx+kk*_nx*_ny][mm]);
	}
    //if(flag==false) std::cout<<"Point does not lie in any of the grid cells"<<std::endl;
    return _returnlist;
  }
    
  private:
    //!Minimum and maximum positions
    double _minx, _maxx, _miny, _maxy, _minz, _maxz;
    
    //!Size of the grid in x, y and z direction
    int _nx, _ny, _nz;
    
    //!Spacing of the grid in x, y and z direction
    double _dx, _dy, _dz;
    
    //!std::vector of lists, where each element contains list of points
    std::vector< std::vector<int> > _grid;

    //!std::vector of list of points to be returned
    std::vector<int> _returnlist;
  };
}
#endif //NodeBin


