/************************ control.c ******************************/
/* MIMD version 7 */
/* Main procedure for pure gauge SU3 */

/* This version combines code for the refreshed molecular dynamics
   algorithm with the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM 
   should be defined.  (Actually, the changes to control.c are minimal
   and the real differences will appear in update.c */

#define CONTROL
//#include "pure_gauge_includes.h"
#include "su3_trP_includes.h"


int main(int argc, char *argv[])  {
   int meascount,todo;
   int prompt;
   double dssplaq,dstplaq;
   complex plp;
   double dtime;
   int i;
   site *s;
   //   int x,y,z,t;
   int dir;

   double w;

initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

 g_sync();
    /* set up */
    prompt = setup();

    /* loop over input sets */
    while( readin(prompt) == 0){

        /* perform warmup trajectories */
        dtime = -dclock();


	w = 0.0;


	d_plaquette(&dssplaq,&dstplaq);
	plp = ploop();
	if(this_node==0)printf("INITIAL %f %f %f %f %f %f %f\n",
			       beta, h, (double)(plp.real-1)/2.0,
			       (double)plp.imag,(dssplaq-1)/2,(dstplaq-1)/2,
			       (dssplaq+dstplaq-2)/4);


	//	fluxplane(w, 3);
	// load_plaq_w(w, &dssplaq,&dstplaq);

	/**
	if(this_node==0)printf("load_plaq %f %f %f\n",
	       (dssplaq-1)/2,(dstplaq-1)/2,(dssplaq+dstplaq-2)/4);

	plaq_twist(w, 3,&dssplaq,&dstplaq);
	if(this_node==0)printf("plaq_twist %f %f %f\n",
	       (dssplaq-1)/2,(dstplaq-1)/2,(dssplaq+dstplaq-2)/4);
	**/

	/**
	FORSOMEPARITY(i, s, EVEN) { 
	   printf("x= %d y= %d z= %d t= %d:\n",
		  s->x,s->y,s->z,s->t);
	   for(dir=0; dir<4; dir++) {
	      printf("PLAQ[%d]: %f %f\n", dir, 
		     s->plaq[dir].e[0][0].real, s->plaq[dir].e[0][0].imag);
	      dumpmat(&(s->plaq[dir]));
	      //printf("link[%d]:\n", dir);
	      //dumpmat(&(s->link[dir]));
	   }
	}
	FORSOMEPARITY(i, s, ODD) { 
	   printf("x= %d y= %d z= %d t= %d:\n",
		  s->x,s->y,s->z,s->t);
	   for(dir=0; dir<4; dir++) {
	      printf("PLAQ[%d]: %f %f\n", dir,
		     s->plaq[dir].e[0][0].real, s->plaq[dir].e[0][0].imag);
	      dumpmat(&(s->plaq[dir]));
	      //	      printf("link[%d]:\n", dir);
	      //	      dumpmat(&(s->link[dir]));
	   }
	}
	**/
	
	//	get_flux_w(w, nz/2, 0);

	/**
	printf("BEGIN staples\n");
	for(dir=0; dir<4; dir++) {
	   dsdu_twist(dir, EVEN, 3);
	   FORSOMEPARITY(i, s, EVEN) { 
	      printf("x= %d y= %d z= %d t= %d dir= %d ",
		     s->x,s->y,s->z,s->t,dir);   
	      printf("link[%d]:\n", dir);
	      dumpmat(&(s->link[dir]));
	      printf("staple= %f %f\n",
		     s->staple.e[0][0].real,s->staple.e[0][0].imag);
	   }
	   dsdu_twist(dir, ODD, 3);
	   FORSOMEPARITY(i, s, ODD) { 
	      printf("x= %d y= %d z= %d t= %d dir= %d ",
		     s->x,s->y,s->z,s->t,dir);   
	      printf("link[%d]:\n", dir);
	      dumpmat(&(s->link[dir]));
	      printf("staple= %f %f\n",
		     s->staple.e[0][0].real,s->staple.e[0][0].imag);
	   }
	}
	exit(0);
	**/

	

        for(todo=warms; todo > 0; --todo ){
	   //	   update_flux();
	   update_u1();
        }
        if(this_node==0)printf("WARMUPS COMPLETED\n");




	// install flux plane
	// timeslice=0, fluxdir=0, q=1, gen=3
	//	flux(0, 0, 1, 3);
	//	load_plaq(&dssplaq, &dstplaq);
	//	get_flux(0);


        // measure
	d_plaquette(&dssplaq,&dstplaq);

	/**
	plp = ploop();
	if(this_node==0)printf("MEASU1TWIST %f %f %f %f %f %f %f\n",
			       beta, h, (double)(plp.real-1)/2.0,
			       (double)plp.imag,(dssplaq-1)/2,(dstplaq-1)/2,
			       (dssplaq+dstplaq-2)/4);
	plaq_twist(w, 3, &dssplaq,&dstplaq);
	if(this_node==0)printf("MEASU1TWIST %f %f %f %f %f %f %f\n",
			       beta, h, (double)(plp.real-1)/2.0,
			       (double)plp.imag,(dssplaq-1)/2,(dstplaq-1)/2,
			       (dssplaq+dstplaq-2)/4);
	**/



        for(todo=trajecs; todo > 0; --todo ){ 

            /* do the trajectories */
	   //   	   update_flux();
	   update_u1();
	  

	   /* measure every "propinterval" trajectories */
	   if((todo%propinterval) == 0){
	      
	       /**/
	       printf("BEGIN staples\n");
	       for(dir=0; dir<4; dir++) {
		  dsdu_twist(dir, EVEN, 3);
		  FORSOMEPARITY(i, s, EVEN) { 
		     printf("x= %d y= %d z= %d t= %d dir= %d ",
			    s->x,s->y,s->z,s->t,dir);   
		     printf("link[%d]:\n", dir);
		     dumpmat(&(s->link[dir]));
		     printf("staple= %f %f\n",
			    s->staple.e[0][0].real,s->staple.e[0][0].imag);
		  }
		  dsdu_twist(dir, ODD, 3);
		  FORSOMEPARITY(i, s, ODD) { 
		     printf("x= %d y= %d z= %d t= %d dir= %d ",
			    s->x,s->y,s->z,s->t,dir);   
		     printf("link[%d]:\n", dir);
		     dumpmat(&(s->link[dir]));
		     printf("staple= %f %f\n",
			    s->staple.e[0][0].real,s->staple.e[0][0].imag);
		  }
	       }
	       /**/



	      //	      load_plaq_w(w, &dssplaq, &dstplaq);
	      //	      get_flux_w(w, nz/2, 0);
	      
	      d_plaquette(&dssplaq,&dstplaq);
	      plp = ploop();
	      if(this_node==0)printf("MEASU1 %f %f %f %f %f %f %f\n",
				     beta, h, (double)(plp.real-1)/2.0,
				     (double)plp.imag,(dssplaq-1)/2,(dstplaq-1)/2,
				     (dssplaq+dstplaq-2)/4);

	      /**
	      plaq_twist(w, 3, &dssplaq,&dstplaq);
	      if(this_node==0)printf("MEASU1TWIST %f %f %f %f %f %f %f\n",
				     beta, h, (double)(plp.real-1)/2.0,
				     (double)plp.imag,(dssplaq-1)/2,(dstplaq-1)/2,
				     (dssplaq+dstplaq-2)/4);
	      
	      **/

	      /* Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq */
	      ++meascount;
	      
	      
	      /**
		 FORALLSITES(i,s) {
		    if(s->t == 0) {
		       node0_printf("x= %d y= %d z= %d t= %d\n", s->x, s->y, s->z, s->t);
		       dumpmat(&(s->link[0]));
		    }
		 }
		**/


                fflush(stdout);
            }
        }       /* end loop over trajectories */
	/*
#ifdef ORA_ALGORITHM
       // gaugefix if requested
       if( fixflag == COULOMB_GAUGE_FIX){
	 gaugefix(TUP,(Real)1.8,600,(Real)GAUGE_FIX_TOL);
           if(this_node==0)printf("FIXED TO COULOMB GAUGE\n");
           fflush(stdout);
       }
       else if( fixflag == LANDAU_GAUGE_FIX){
	 gaugefix(8,(Real)1.8,600,(Real)GAUGE_FIX_TOL);
           if(this_node==0)printf("FIXED TO LANDAU GAUGE\n");
           fflush(stdout);
       }
#endif
	*/


        if(this_node==0)printf("RUNNING COMPLETED\n");

        dtime += dclock();
        if(this_node==0){
            printf("Time = %e seconds\n",dtime);
        }
        fflush(stdout);
	dtime = -dclock();

        /* save lattice if requested */
        if( saveflag != FORGET ){
	  save_lattice( saveflag, savefile, stringLFN );
        }
    }
    return 0;
}
