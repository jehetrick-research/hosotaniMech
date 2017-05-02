/****** dsdu_qhb.c  -- compute the staple ******************/
/* MIMD version 6 */
/* UMH: Combined with Schroedinger functional version, Jan 2000 */

//#include "../generic_pg/generic_pg_includes.h"
//#include "../include/loopend.h"

#include "su3_trP_includes.h"



void dsdu_twist(float w, int dir1,int parity, int gen) {
   register int i,dir2,otherparity=0;
   site *st;
   msg_tag *tag0,*tag1,*tag2,*tag3;
   int start;
   su3_matrix tmat1,tmat2,tmat3,tfin;
   su3_matrix wphase;
   
   switch(parity) {
     case EVEN:		otherparity=ODD;	break;
     case ODD:		otherparity=EVEN;	break;
     case EVENANDODD:	otherparity=EVENANDODD;	break;
   }
   
   /* Loop over other directions, computing force from plaquettes in
      the dir1,dir2 plane */

   start=1; /* indicates staple sum not initialized */
   
   set_wphase(w, 3, 1, ny, &wphase);  // wphase*
   //	 printf("set_w(-1):\n");
   //	 dumpmat(&wphase);
   
   for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1) {
	 /* get link[dir2] from direction dir1 on other parity */
	 tag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
				   dir1, otherparity, gen_pt[0] );
	 
	 /* get link[dir2] from direction dir1 */
	 tag1 = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
				   dir1, parity, gen_pt[1] );
	 
	 /* get link[dir1] from direction dir2 */
	 tag2 = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
				   dir2, parity, gen_pt[2] );
	 

	 /////////////////////////////////////////////////
	 /* Construct Lower staple (computed at backward site) */
	 wait_gather(tag0);

	 FORSOMEPARITY(i,st,otherparity){
	    su3mat_copy((su3_matrix *)gen_pt[0][i], &tfin);

	    if((st->y==ny-1) && (dir1==0) && (dir2==1)) {
		  mult_su3_nn( &wphase, &(st->link[dir1]), &tmat3 );
		  mult_su3_an( &(st->link[dir2]), &tmat3, &tmat1); 
	    } else 
	    if((st->y==ny-1) && (dir1==1) && (dir2==0)) {
		  // mult: "an" here because, wphase set with -1 above
                  mult_su3_an( &wphase, (su3_matrix *)gen_pt[0][i], &tfin );   // nn -> an
		  mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1);
	    } else {
		  mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1);
	    }
	    mult_su3_nn( &tmat1, &tfin, &(st->tempmat1));
	 }
	 cleanup_gather(tag0);

	 /* get tempmat1 from direction -dir2 */
	 tag3 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
				   OPP_DIR(dir2), parity, gen_pt[3] );
	 
	 //////////////////////////////////////////////////
	 /* Upper staple */
	 wait_gather(tag1);
	 wait_gather(tag2);

	 set_wphase(w, 3, -1, ny, &wphase);
	 //	 printf("set_w(1):\n");
	 //	 dumpmat(&wphase);
	 if(start){  /* this is the first contribution to staple */
	    FORSOMEPARITY(i,st,parity){
	       su3mat_copy((su3_matrix *)gen_pt[1][i], &tfin);

	       if((st->y==ny-1) && (dir1==0) && (dir2==1)) {
		     mult_su3_nn( &wphase, (su3_matrix *)gen_pt[2][i], &tmat3 );
		     mult_su3_nn( &(st->link[dir2]), &tmat3, &tmat1 );
	       } else 
	       if((st->y==ny-1) && (dir1==1) && (dir2==0)) {
		  mult_su3_nn( &wphase, (su3_matrix *)gen_pt[1][i], &tfin );
		  mult_su3_nn( &(st->link[dir2]), (su3_matrix *)gen_pt[2][i], &tmat1);
	       } else {
		  mult_su3_nn( &(st->link[dir2]), (su3_matrix *)gen_pt[2][i],
			       &tmat1 );
	       }
	       mult_su3_na( &tmat1, &tfin, &(st->staple) );
	       /*
	       printf("stp[%d]: x= %d y= %d z= %d t= %d dir2= %d ",
		      dir1, st->x,st->y,st->z,st->t,dir2);   
	       printf("up.staple= %f %f\n",
		   st->staple.e[0][0].real,st->staple.e[0][0].imag);

	       printf("up y= %d: tmat1\n", st->y);
	       dumpmat(&tmat1);
	       printf("link[dir2=%d]\n",dir2);
	       dumpmat(&(st->link[dir2]));
	       printf("gen2\n");
	       dumpmat((su3_matrix *)gen_pt[2][i]);
	       */

	    }
	    start=0; 
	 }
	 else {
	    FORSOMEPARITY(i,st,parity){
	       su3mat_copy((su3_matrix *)gen_pt[1][i], &tfin);

	       if((st->y==ny-1) && (dir1==0) && (dir2==1)) {
		  mult_su3_nn( &wphase, (su3_matrix *)gen_pt[2][i], &tmat3 );
		  mult_su3_nn( &(st->link[dir2]), &tmat3, &tmat1 );
	       } else 
	       if((st->y==ny-1) && (dir1==1) && (dir2==0)) {
		  mult_su3_nn( &wphase, (su3_matrix *)gen_pt[1][i], &tfin );
		  mult_su3_nn( &(st->link[dir2]), (su3_matrix *)gen_pt[2][i], &tmat1);
	       } else {
		  mult_su3_nn( &(st->link[dir2]), (su3_matrix *)gen_pt[2][i],
			       &tmat1 );
	       }
	       mult_su3_na( &tmat1, &tfin, &tmat2 );
	       add_su3_matrix( &(st->staple), &tmat2, &(st->staple));
	       /*
	       printf("stp[%d]: x= %d y= %d z= %d t= %d dir2= %d ",
		      dir1, st->x,st->y,st->z,st->t,dir2);   
	       printf("up+staple= %f %f\n",
		      tmat2.e[0][0].real, tmat2.e[0][0].imag);
	       */
	    }
	 } /* upper staple */
	 cleanup_gather(tag1);
	 cleanup_gather(tag2);
	
	/* Add lower staple */
	wait_gather(tag3);
	FORSOMEPARITY(i,st,parity){
	   /*
		 printf("stp[%d]: x= %d y= %d z= %d t= %d dir2= %d ",
			dir1, st->x,st->y,st->z,st->t,dir2);   
		 printf("dn+staple= %f %f\n",
			((su3_matrix *)gen_pt[3][i])->e[0][0].real, 
			((su3_matrix *)gen_pt[3][i])->e[0][0].imag );
	   */
	   add_su3_matrix( &(st->staple), (su3_matrix *)gen_pt[3][i],
			   &(st->staple));
	}	/* lower staple */
	      cleanup_gather(tag3);
    }
}

