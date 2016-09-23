#include "ShapeTet4.h"

namespace voom {  
  void ShapeTet4::compute(const CoordinateArray & s) {

    // feed in Na^

    _functions[0] = s(0);
    _functions[1] = s(1);
    _functions[2] = s(2);
    _functions[3] = 1.0-s(0)-s(1)-s(2);

    //    for (int q=0; q<4; q++)  
    // std::cout << "functions[" << q << "] = " << _functions[q] << std::endl;

    /*
    // we can output the barycentric curvilinear coordinates as a check,  
    // they should be what we input initially, but will change during tests.  
    
    std::cout << "s1 =" << std::endl 
    << s(0) << std::endl;
    std::cout << "s2 =" << std::endl 
    << s(1) << std::endl;
    std::cout << "s3 =" << std::endl 
    << s(2) << std::endl;
    */
    
    //  feed in the Na^,j matrix and call it "dnds"
    
    tvmet::Matrix<double,4,3> dnds;
    dnds =  1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
	    0.0, 0.0, 1.0,
	   -1.0,-1.0,-1.0;

    //  output "dnds" as a check, which is the matrix of Na^,j components.
    
    /*
      std::cout << "dnds =" << std::endl 
      << dnds << std::endl;
      
      
    */
    
    // create a 3x3 "Jacobian" Jij matrix called "dxds" and make it equal 
    //to positions * [Na^,j]
    
    tvmet::Matrix<double,3,3> dxds(0.0);
    for(int i=0; i<3; i++) {
      for(int alpha=0; alpha<3; alpha++) {
	for(int a=0; a<4; a++) {
	  dxds(i,alpha) += _positions[a](i)*dnds(a,alpha);
	  //   std::cout << "dxds(" << i << "," << alpha << ") = " << dxds(i,alpha) << std::endl;
	}
      }
    }
    
    // create a numeric scalar determinant called "xsj" from the 
    //Jij matrix, which is called "dxds"
   
   double xsj = dxds(0,0)*(dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1))
               -dxds(1,0)*(dxds(0,1)*dxds(2,2)-dxds(0,2)*dxds(2,1))
               +dxds(2,0)*(dxds(0,1)*dxds(1,2)-dxds(0,2)*dxds(1,1));
   
   if (xsj <= 0.0) {
     std::cerr << "ShapeTet4 negative jacobian!" << std::endl;
     exit(0);
   }

   tvmet::Matrix<double,3,3> invJac(0.0);   
   invJac(0,0) = (dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1))/xsj;
   invJac(0,1) =-(dxds(0,1)*dxds(2,2)-dxds(0,2)*dxds(2,1))/xsj;
   invJac(0,2) = (dxds(0,1)*dxds(1,2)-dxds(0,2)*dxds(1,1))/xsj;
   invJac(1,0) =-(dxds(1,0)*dxds(2,2)-dxds(1,2)*dxds(2,0))/xsj;
   invJac(1,1) = (dxds(0,0)*dxds(2,2)-dxds(0,2)*dxds(2,0))/xsj;
   invJac(1,2) =-(dxds(0,0)*dxds(1,2)-dxds(0,2)*dxds(1,0))/xsj;
   invJac(2,0) = (dxds(1,0)*dxds(2,1)-dxds(1,1)*dxds(2,0))/xsj;
   invJac(2,1) =-(dxds(0,0)*dxds(2,1)-dxds(0,1)*dxds(2,0))/xsj;
   invJac(2,2) = (dxds(0,0)*dxds(1,1)-dxds(0,1)*dxds(1,0))/xsj;
   
#if 0   
   double xsj =   dxds(0,0)*(dxds(1,1)*dxds(2,2) - dxds(1,2)*dxds(2,1))
                - dxds(0,1)*(dxds(1,0)*dxds(2,2) - dxds(1,2)*dxds(2,0))
                + dxds(0,2)*(dxds(1,0)*dxds(2,1) - dxds(1,1)*dxds(2,0));

   // make a cofactor matrix 
   tvmet::Matrix<double,3,3> Cofactors(0.0);
   for (int y=0; y<3; y++) {
     if (y == 0) {int a = 1; int b = 2; }
     if (y == 1) {int a = 0; int b = 2; }
     if (y == 2) {int a = 0; int b = 1; }

     for(int z=0; z<3; z++) {
       if (z == 0) {int c = 1; int d = 2; }
       if (z == 1) {int c = 0; int d = 2; }
       if (z == 2) {int c = 0; int d = 1; }
       
       // here the actual formula would be Cofactors(i,j)=-1^(i+j) where i,j are 123, but since we have
       // i,j = 012, it is actually -1^(i+1+j+1). But can simplify this to -1^2 * -1^(i+j) = -1^(i+j)
       // so we can use the same notation even though ij is 012 as opposed to 123
       
       Cofactors(y,z) = std::pow(-1,y+z)*(dxds(a,c)*dxds(b,d) - dxds(a,d)*dxds(b,c));
     }
   }
   
   // making the adjoint matrix from the cofactor matrix
   tvmet::Matrix<double,3,3> Adjoint(0.0);
   Adjoint =  trans(Cofactors);
   
   // creating the inverse jacobian from the Adjoint matrix
   tvmet::Matrix<double,3,3> invJac(0.0);
   invJac = Adjoint / xsj;
#endif
   
   // volume normalize the scalar jacobian by multiplying by the volume = 1/6 for this tetrahedron
   _jacobian = xsj/6.0;
   
   // lets check the value of xsj and show the jacobian matrix
   /*
     std::cout << "Jacobian Matrix =:" << dxds << "with determinant of: " << xsj << "\n";
     std::cout << "Jacobian Inverse =:" << invJac ;
     std::cout << "this should give the identity matrix when multiplied:" << dxds*invJac ;
   */
   
   // set derivatives Na,i = Na^,j * Jji^-1
   //     _derivatives = dnds*invJac;
   tvmet::Matrix<double,4,3> dndx;
   dndx = dnds*invJac;
   for(int a=0; a<4; a++)
     for(int i=0; i<3; i++) _derivatives[a](i) = dndx(a,i);
   
  }
};


//   xsj = dxds[0]*(dxds[4]*dxds[8]-dxds[5]*dxds[7])
//     -dxds[3]*(dxds[1]*dxds[8]-dxds[2]*dxds[7])
//     +dxds[6]*(dxds[1]*dxds[5]-dxds[2]*dxds[4]);
  
//   if (xsj <= 0.0)
//     ErrorFunction("MyShapeTNT::ShapeFunctions",
//                   "Negative Jacobian");
  
//   xjac[0] = (dxds[4]*dxds[8]-dxds[5]*dxds[7])/xsj;
//   xjac[1] =-(dxds[1]*dxds[8]-dxds[2]*dxds[7])/xsj;
//   xjac[2] = (dxds[1]*dxds[5]-dxds[2]*dxds[4])/xsj;
//   xjac[3] =-(dxds[3]*dxds[8]-dxds[5]*dxds[6])/xsj;
//   xjac[4] = (dxds[0]*dxds[8]-dxds[2]*dxds[6])/xsj;
//   xjac[5] =-(dxds[0]*dxds[5]-dxds[2]*dxds[3])/xsj;
//   xjac[6] = (dxds[3]*dxds[7]-dxds[4]*dxds[6])/xsj;
//   xjac[7] =-(dxds[0]*dxds[7]-dxds[1]*dxds[6])/xsj;
//   xjac[8] = (dxds[0]*dxds[4]-dxds[1]*dxds[3])/xsj;
  
//   Jacobian = xsj/24.0;
