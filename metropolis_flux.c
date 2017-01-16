/**************************** metropolis_trP.c ****************************/
/* MIMD version 6 */

/************************** metropolis_dense.c *******************************/
/* Metropolis updating for SU3 pure gauge */
/* MIMD version 4 */
/* J. Hetrick and D. Toussaint June 1995 */
/* update with "almost quenched" approximation for nonzero density 6/1/95 */
/* monte_space() does spatial links, monte_time() does temporal */
//

// Now: metropolis_trP.c    MILC v.7.x
//
// The update here incorporates the (h Tr_A P) term in the action, as 
// as in Ogilvie and Myers arXiv:0707.1869
//
/////////////////////////////////////////////////////////////////////

#include "su3_trP_includes.h"

/* Generic definitions - could be useful elsewhere */
#define FORSPACEUPDIR(dir) for(dir=XUP; dir<=ZUP; dir++)
#define OPP_PAR(parity) (0x03 ^ parity) /* Switches EVEN and ODD. 
                                           Nulls EVENANDODD*/
#define CMULINC(a,b,c) { (c).real += (a).real*(b).real - (a).imag*(b).imag; \
                      (c).imag += (a).real*(b).imag + (a).imag*(b).real; }
#define CMULDEC(a,b,c) { (c).real -= (a).real*(b).real - (a).imag*(b).imag; \
                      (c).imag -= (a).real*(b).imag + (a).imag*(b).real; }

float make_change(su3_matrix *mat, site *st, float scale);
//void base_measurements(int time);
//complex numerator(float c, complex T);
//complex det_special(float c, complex T);

/* det. phase, P. loop, and total density before this processor starts
   updating this parity and slice */
float base_phase;
complex base_ploop, base_density;




// this monte_space() does not update the 'fluxdir' links on t-slice 'tslice' 

void monte_space_flux(double w, int NumStp, int tslice, int fluxdir) {
  int Nhit;
  int parity;
  float scale;		/* limits size of change matrix */
  float theta,theta_rms;	/* rms angle of change matrix */
  su3_matrix change;	/* proposed change in link */
  su3_matrix newlink;	/* change * oldlink */
  int dir, i;
  register site *st;
  int accept, reject;	/* number of accepts and rejects */
  float oldaction, newaction;
  su3_matrix oldUstaple, newUstaple;
  

  accept = reject = 0;
  theta_rms = 0.0;
  scale = 1.4;

  for(parity=ODD;parity<=EVEN;parity++) {
    FORSPACEUPDIR(dir) {
      /* compute the gauge force */
       dsdu_twist(w,dir,parity,3); 

      /* now for the Metropolis updating */
      FORSOMEPARITY(i,st,parity) {

	 // leave s->z==0 xy-plane alone:
	 //////////////////////////////////////
	 if((st->z == 0) && ((dir == 0) || (dir == 1))) { continue; }

	    for( Nhit = 0 ; Nhit < NumStp; Nhit++) { 
	       theta = make_change(&change,st,scale); 
	       theta_rms += theta*theta;
	       mult_su3_nn( &change, &(st->link[dir]), &newlink );
	       
	       /* compute old and new SU(3) action */
	       /*
	       oldaction=(0.333333*beta)*realtrace_su3( &(st->link[dir]), 
							&(st->staple) );
	       newaction=(0.333333*beta)*realtrace_su3( &newlink, 
	       					&(st->staple) );
	       */

	       /* compute old action and new action */
	       mult_su3_na(&(st->link[dir]), &(st->staple), &oldUstaple);
	       mult_su3_na(&newlink, &(st->staple), &newUstaple);

	       //newaction = beta*cos(theta_new - theta_staple);
	       oldaction = beta * oldUstaple.e[0][0].real;
	       newaction = beta * newUstaple.e[0][0].real;


	       /* accept or reject */
	       if( newaction > oldaction ){
		  st->link[dir]=newlink;
		  accept++;
	       }
	       else{ /* trace decreased */
		  printf("CALL myrand() in monte_space_flux() accept/reject\n");
		  if( myrand(&(st->site_prn)) < exp( newaction-oldaction ) ){
		     st->link[dir]=newlink;
		     accept++;
		  }
		  else{
		     reject++;
		  }
	       }
	    } /* Nhit */
      } /*   st */
    } /*  direction */
  } /* parity */

  /* diagnostics */
/**
  printf("monte_space: accept= %d reject= %d fraction= %.2f  theta_rms= %.3f\n",
    accept,reject,accept/(float)(accept+reject),
    sqrt(theta_rms/(accept+reject)) );
**/

} /* monte_space_flux */


/* time direction update includes h tr_A P in its action */

void monte_time_twist(double w, int NumStp) {
  int Nhit;
  int time,parity;
  float scale;		/* limits size of change matrix */
  float theta,theta_rms;	/* rms angle of change matrix */
  su3_matrix change;	/* proposed change in link */
  su3_matrix newlink;	/* change * oldlink */
  su3_matrix tmat;	/* scratch */
  register int i;
  register site *st;
  int accept, reject;	/* number of accepts and rejects */
  int mcount;		/* number of measurements */
  float oldaction,newaction;
  complex trPf;
  float trPA;
  complex phase_aver,ploop_aver,mag_ploop_aver,density_aver;

  scale = 1.4;
  mcount = accept = reject = 0;
  theta_rms = 0.0;
  phase_aver = ploop_aver = mag_ploop_aver = density_aver = cmplx(0.0,0.0);

  for(parity=ODD;parity<=EVEN;parity++) {
    /* compute the gauge force */
     dsdu_twist(w,TUP,parity,3); 


    /* update time slices separately, because they are linked by P. loops */
    for(time=0;time<nt;time++){
      /* evaluate "almost Polyakov loops", total P. loop, total phase */
      /* see extended comment at bottom */

      ploop_less_slice(time,parity);
      ploop_less_slice(time,OPP_PAR(parity)); /* temporary? */




      /* now for the Metropolis updating */
      FORSOMEPARITY(i,st,parity)if(st->t==time) {
	/* generate random SU(3) matrix */
	for( Nhit = 0 ; Nhit < NumStp; Nhit++) {
	   printf("CALL make_change() in monte_time_twist()\n");
	  theta = make_change(&change,st,scale); 
	  theta_rms += theta*theta;
	  mult_su3_nn( &change, &(st->link[TUP]), &newlink );

	  // compute old action
	  // Wilson part S_w
	  oldaction=(0.333333*beta)*realtrace_su3(&(st->link[TUP]),
		  			    	  &(st->staple));
	  // add (-h tr_A P)
	  mult_su3_nn( &(st->link[TUP]), &(st->ploop_t), &tmat );
	  trPf = trace_su3(&tmat);
	  trPA = cabs_sq(&trPf) - 1.0;

	  oldaction -= h*trPA;

	  //////////////////////
	  // compute new action
	  newaction=(0.333333*beta)*realtrace_su3( &newlink, &(st->staple) );
	  mult_su3_nn( &newlink, &(st->ploop_t), &tmat );
	  trPf = trace_su3(&tmat);
	  trPA = cabs_sq(&trPf) - 1.0;

	  newaction -= h*trPA;

	  /* accept or reject */
	  if( newaction > oldaction ){
	    st->link[TUP]=newlink;
	    accept++;
	  }
	  else{ /* trace decreased */
	    if( myrand(&(st->site_prn)) < exp( newaction-oldaction ) ){
	      st->link[TUP]=newlink;
	      accept++;
	    }
	    else{
	      reject++;
	    }
	  }

	} /* Nhit: number of Metropolis updates */
      } /* st: timeslices sites */
    } /* time slice */
  } /* parity */

  // diagnostics
  /**
  printf("monte_time:  accept= %d reject= %d fraction= %.2f  theta_rms= %.3f\n",
    accept,reject,accept/(float)(accept+reject),
    sqrt(theta_rms/(accept+reject)) );
  **/

} /* monte_time */



/* Calculate product of timelike links for all time slices except "time",
   for sites where parity at "time" is "parity".
   Put the result in the matrix "ploop_t"
*/
void ploop_less_slice(int time,int parity){
  int l_time;	/* time at which we are multiplying */
  int c_time;   /* l_time%nt - use this one in accessing sites */
  int l_parity; 
	/* parity of sites at l_time where parity at "time" is "parity" */
  register int i;
  register site *s;
  msg_tag *tag;

  if(parity==EVENANDODD){
    if(this_node==0)printf("Bad parity in ploop_less_slice()\n");
    terminate(0);
  }

  FORSOMEPARITY( i,s,OPP_PAR(parity) )if(s->t==(time-1+nt)%nt){
    s->ploop_t = s->link[TUP];
  }

  for( l_time = time + nt-2; l_time > time; l_time--){
    c_time = l_time >= nt ? l_time-nt : l_time;
    if( (l_time-time)%2 ==0){l_parity=parity;}
    else {l_parity = OPP_PAR(parity);}

    /* gather current product from slice above */
    tag = start_gather_site( F_OFFSET(ploop_t), sizeof(su3_matrix), TUP,
      l_parity, gen_pt[0]);
    wait_gather(tag);
    FORSOMEPARITY( i,s,l_parity )if(s->t==c_time){
      mult_su3_nn( &(s->link[TUP]), (su3_matrix *)gen_pt[0][i],
	&(s->ploop_t) );
	/* since only on one time slice, don't have to worry if
	   gen_pt points to ploop_t */
    } /* end loop over sites at time l_time */
    cleanup_gather(tag);
  } /* end loop over l_times */

  /* finish by bringing result to desired time slice */
  tag = start_gather_site( F_OFFSET(ploop_t), sizeof(su3_matrix), TUP,
    parity, gen_pt[0]);
  wait_gather(tag);
  FORSOMEPARITY( i,s,parity )if(s->t==time)
    s->ploop_t = *(su3_matrix *)gen_pt[0][i];
  cleanup_gather(tag);
}


/* Make ** U(1) ** change matrix for Metropolis step. */

float make_change(su3_matrix *mat, site *st, float scale){
  register int ia,ib;
  register float theta,c,s;

  /* load identity into "mat" */
  for(ia=0; ia<3; ia++) {
    for(ib=0; ib<3; ib++){
      mat->e[ia][ib].real = mat->e[ia][ib].imag = 0.0;
    }
    mat->e[ia][ia].real = 1.0;
  }

  /* random (symmetric) change */
  /* angles range from +- theta/4 to theta/2 */
  //  do{ theta = scale * (myrand(&(st->site_prn)) - 0.5); }
  //	while( fabs(theta) < 0.25*scale );

  printf("CALL make_change()\n");
  theta = scale * 2*(myrand(&(st->site_prn)) - 0.5);
  c = cos(theta); s = sin(theta);

  mat->e[0][0].real = mat->e[1][1].real = c;
  mat->e[0][0].imag = s;
  mat->e[1][1].imag = -s;
  
// printf("make_change: theta = %e\n",theta);

  return(theta);
}


