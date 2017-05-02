/************************** flux.c *******************************/
//
//  fluxplane(w,gen):
//  Create a uniform flux through xy-planes @ z=0 on all t-slices
//  in a given Cartan sub-group, 3 or 8.
//
//  load_plaq(ssplaq, stplaq)
//  compute and store 6 local plaquettes at sites
//
//  get_flux(zplane)
//  compute the average flux across an XY plane at z=zplane
//  average over time-slices
//   
/////////////////////////////////////////////////////////////////////

#include "su3_trP_includes.h"


void fluxplane(double w, int gen) {
   int i, a, b;
   register site *s;
   float wcs, wsn;
   float wcs2, wsn2;
   
   //////////////////////////////////////////////
   // Assumes COLD Lattice has been initialized
   //////////////////////////////////////////////

   FORALLSITES(i,s) {
      if(s->z == 0) {
	 switch(gen) {
	 case 3:
	    wcs = cos(w*(s->y));
	    wsn = -sin(w*(s->y));
	    
	    //	    for(a=0; a<3; a++)  {
	    //	       for(b=0; b<3; b++)  {
	    //		  if (a != b)  {
	    //		     s->link[0].e[a][b] = cmplx(0.0,0.0);
	    //		     s->link[1].e[a][b] = cmplx(0.0,0.0);
	    //		  }
		     s->link[0].e[0][0] = cmplx(wcs, wsn);	
		     s->link[0].e[1][1] = cmplx(wcs,-wsn);
		     s->link[0].e[2][2] = cmplx(1.0, 0.0);
	    break;

	 case 8:
	    wcs = cos(w*(s->y));
	    wsn = -sin(w*(s->y));
	    wcs2 = cos(2*w*(s->y));
	    wsn2 = -sin(2*w*(s->y));
	    
	    for(a=0; a<3; a++)  {
	       for(b=0; b<3; b++)  {
		  if (a != b)  {
		     s->link[0].e[a][b] = cmplx(0.0, 0.0);
		     s->link[1].e[a][b] = cmplx(0.0,0.0);
		  }
		  else  {
		     s->link[0].e[0][0] = cmplx(wcs,  wsn);	
		     s->link[0].e[1][1] = cmplx(wcs,  wsn);
		     s->link[0].e[2][2] = cmplx(wcs2,-wsn2);
		     // needed?
		     s->link[1].e[0][0] = cmplx(1.0, 0.0);	
		     s->link[1].e[1][1] = cmplx(1.0, 0.0);
		     s->link[1].e[2][2] = cmplx(1.0, 0.0);
		  }
	       }
	    }
	    break;
	 } 
      }
   }

} // flux(int tslice, int fluxdir, int q, int gen)




void load_plaq_w(double w, double *ss_plaq,double *st_plaq) {
   su3_matrix *su3mat;
   register int i,dir1,dir2;
   site *s;
   register su3_matrix *m1,*m4;
   su3_matrix mtmp, mtmp2;
   double ss_sum,st_sum;
   msg_tag *mtag0,*mtag1;
   int pln; // tracks x-y, x-z, y-z, x-t, y-t, z-t planes for s->plaq[pln]
   //   double wcs,wsn,w2cs,w2sn; // flux amount
   su3_matrix wphase;
   int munu;



   pln = 0;
   ss_sum = st_sum = 0.0;


   set_wphase(w, 3, -1, ny, &wphase);
   //   printf("load_plaq: set_wphase:\n");
   //   dumpmat(&wphase);

   su3mat = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
   if(su3mat == NULL) {
      printf("plaquette: can't malloc su3mat\n");
      fflush(stdout); terminate(1);
   }
   

   munu = 0;
   for(dir1=YUP;dir1<=TUP;dir1++){
      for(dir2=XUP;dir2<dir1;dir2++){

	 mtag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
				    dir1, EVENANDODD, gen_pt[0] );
	 mtag1 = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
				    dir2, EVENANDODD, gen_pt[1] );
	 
	 FORALLSITES(i,s){
	    m1 = &(s->link[dir1]);
	    m4 = &(s->link[dir2]);
	    mult_su3_an(m4,m1,&su3mat[i]);
	 }
	 
	 wait_gather(mtag0);
	 wait_gather(mtag1);
	 
	 FORALLSITES(i,s){
	    // add w-phase to xy B.C.
	    if((s->y==ny-1) && (dir1==1) && (dir2==0)) {
	       mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
			 &mtmp2);
	       mult_su3_nn(&mtmp2, &wphase, &mtmp);
	    } else {
	       mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
			    &mtmp);
	    }
	    mult_su3_an((su3_matrix *)(gen_pt[1][i]), &mtmp, 
			&(s->plaq[munu]));

	    if(dir1==TUP ) { st_sum += trace_su3(&(s->plaq[munu])).real; }
	    else { ss_sum += (trace_su3(&(s->plaq[munu])).real); }
	 }
	       
	 cleanup_gather(mtag0);
	 cleanup_gather(mtag1);
	 munu++;
      }
   }
   g_doublesum( &ss_sum );
   g_doublesum( &st_sum );
   *ss_plaq = ss_sum /((Real)(3*nx*ny*nz*nt));
   *st_plaq = st_sum /((double)(3*nx*ny*nz*nt));
   
   free(su3mat);
}


void load_plaq(double *ss_plaq,double *st_plaq) {
   su3_matrix *su3mat;
   register int i,dir1,dir2;
   site *s;
   register su3_matrix *m1,*m4;
   su3_matrix mtmp, mtmp2;
   double ss_sum,st_sum;
   msg_tag *mtag0,*mtag1;
   int pln; // tracks x-y, x-z, y-z, x-t, y-t, z-t planes for s->plaq[pln]
   //   double wcs,wsn,w2cs,w2sn; // flux amount
   su3_matrix wphase;
   int munu;


   printf("FIX THIS:  load_plaq()\n");
   exit(0);


   pln = 0;
   ss_sum = st_sum = 0.0;


   //   set_wphase(w, 3, -1, ny, &wphase);
   //   printf("load_plaq: set_wphase:\n");
   //   dumpmat(&wphase);

   su3mat = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
   if(su3mat == NULL) {
      printf("plaquette: can't malloc su3mat\n");
      fflush(stdout); terminate(1);
   }
   

   munu = 0;
   for(dir1=YUP;dir1<=TUP;dir1++){
      for(dir2=XUP;dir2<dir1;dir2++){

	 mtag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
				    dir1, EVENANDODD, gen_pt[0] );
	 mtag1 = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
				    dir2, EVENANDODD, gen_pt[1] );
	 
	 FORALLSITES(i,s){
	    m1 = &(s->link[dir1]);
	    m4 = &(s->link[dir2]);
	    mult_su3_an(m4,m1,&su3mat[i]);
	 }
	 
	 wait_gather(mtag0);
	 wait_gather(mtag1);
	 
	 FORALLSITES(i,s){
	    // add w-phase to xy B.C.
	    if((s->y==ny-1) && (dir1==1) && (dir2==0)) {
	       mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
			 &mtmp2);
	       mult_su3_nn(&mtmp2, &wphase, &mtmp);
	    } else {
	       mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
			    &mtmp);
	    }
	    mult_su3_an((su3_matrix *)(gen_pt[1][i]), &mtmp, 
			&(s->plaq[munu]));

	    if(dir1==TUP ) { st_sum += trace_su3(&(s->plaq[munu])).real; }
	    else { ss_sum += (trace_su3(&(s->plaq[munu])).real); }
	 }
	       
	 cleanup_gather(mtag0);
	 cleanup_gather(mtag1);
	 munu++;
      }
   }
   g_doublesum( &ss_sum );
   g_doublesum( &st_sum );
   *ss_plaq = ss_sum /((Real)(3*nx*ny*nz*nt));
   *st_plaq = st_sum /((double)(3*nx*ny*nz*nt));
   
   free(su3mat);
}





///////////////////////////////////////////
// if loadPlaq = 1, will (re)load_plaq() //
/////////////////////////////////////////// 

void get_flux_w(double w, int zplane, int loadPlaq) {
   register int i,a,b;
   register site *s;
   complex lam[3];
   double theta[3];
   double flux[3];
   //   su3_matrix Pave;
   double ssplaq, stplaq;


   for(a=0; a<3; a++) {
      flux[a] = theta[a] = 0.0;
      //      for(j=0; j<3; j++) { 
      //	 Pave.e[a][j] = cmplx(0.0, 0.0);
      //      }
   }

   if(loadPlaq==1) {
      load_plaq_w(w, &ssplaq, &stplaq);
   }

   FORALLSITES(i,s){
      if((s->z)==zplane) { 
	 for(a=0; a<3; a++) {
	    CSUM(lam[a], s->plaq[0].e[a][a]);
	    theta[a] = atan2(s->plaq[0].e[a][a].imag, s->plaq[0].e[a][a].real);
	    flux[a] += theta[a];
	    //	    for(b=0; b<3; b++) {
	    //	       CSUM(Pave.e[a][b], (s->plaq[0].e[a][b]));
	    //	    }
	 }
      }
   }
   
   for(a=0; a<3; a++) { 
      g_floatsum(&(flux[a])); 
      flux[a] /= nx * ny * nt;
   }
   
   
   node0_printf("FLUX b= %f h= %f w= %f zplane= %d flux_th= %f %f %f\n", beta, h, w, zplane, 
		flux[0], flux[1], flux[2]);
}


