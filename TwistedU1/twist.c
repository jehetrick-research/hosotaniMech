////////////////////////////////////////////////////
// Various utility functions for Twisted B.C.s
////////////////////////////////////////////////////

#include "su3_trP_includes.h"



// set w-phases, if cc=-1 complex conjugate
void set_wphase(double w, int gen, int cc, int y, su3_matrix *wphase) {
   int a,b;
   double wcs, wsn, w2cs, w2sn;

   //   printf("set_wphase: w= %f gen= %d cc= %d y= %d\n", w,gen,cc,y);

   switch(gen) {
   case 3:
      wcs = cos(w*y);
      wsn = cc*sin(w*y);
      
      //      printf("set_wphase: wcs= %f wsn= %f\n", wcs, wsn);
      
      for(a=0; a<3; a++)  {
	 for(b=0; b<3; b++)  {
	    if (a != b)  {
	       wphase->e[a][b] = cmplx(0.0,0.0);
	    }
	 }
      }
      wphase->e[0][0] = cmplx(wcs, wsn);	
      wphase->e[1][1] = cmplx(wcs,-wsn);
      wphase->e[2][2] = cmplx(1.0, 0.0);
      break;
      
   case 8:
      wcs = cos(w*y);
      wsn = cc*sin(w*y);
      w2cs = cos(2*w*y);
      w2sn = cc*sin(2*w*y);
      
      for(a=0; a<3; a++)  {
	 for(b=0; b<3; b++)  {
	    if (a != b)  {
	       wphase->e[a][b] = cmplx(0.0, 0.0);
	       wphase->e[a][b] = cmplx(0.0,0.0);
	    }
	 }
      }
      wphase->e[0][0] = cmplx(wcs,  wsn);	
      wphase->e[1][1] = cmplx(wcs,  wsn);
      wphase->e[2][2] = cmplx(w2cs,-w2sn);
      break;
   } 
}
