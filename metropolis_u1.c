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
// This code does:
//   update_u1: a U(1) update ONLY in the e[1][1] element of the links
//   update_cartan: U(1)^2 update of the diagonal elements of the links  
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


float make_change_u1(su3_matrix *mat, site *st, float scale);
//float make_change_cartan(su3_matrix *mat, site *st, float scale);
//float make_change_su3(su3_matrix *mat, site *st, float scale);
void base_measurements(int time);
complex numerator(float c, complex T);
complex det_special(float c, complex T);





void monte_u1(int NumStp) {
  int Nhit;
  int parity;
  float scale;		/* limits size of change matrix */
  float theta,theta_rms;	/* rms angle of change matrix */
  float theta_old, theta_new, dtheta, theta_staple;
  su3_matrix change, Unew;	/* proposed change in link */
  su3_matrix newlink;	/* change * oldlink */
  su3_matrix tmpsu3;
  su3_matrix oldUstaple, newUstaple;
  int dir, i, a,b;
  register site *st;
  int accept, reject;	/* number of accepts and rejects */
  float oldaction, newaction, r;

  accept = reject = 0;
  theta_rms = 0.0;
  scale = 1.4;

  for(parity=ODD;parity<=EVEN;parity++) {
     for(dir=0; dir<4; dir++) {
      /* compute the gauge force--i.e. load the staples */
      dsdu_qhb(dir,parity); 


      /* now for the Metropolis updating */
      FORSOMEPARITY(i,st,parity) {

	/* generate random U(1) matrix */
	 for( Nhit = 0 ; Nhit < NumStp; Nhit++) { 
	    //	    printf("update %d %d %d %d -> %d p %d drand: %f\n", 
	    //		   st->x,st->y,st->z,st->t,dir, parity,drand48());

	    printf("CALL make_change_u1() in monte_u1()\n");

	    dtheta = make_change_u1( &change, st, scale); 
	    mult_su3_nn(&change, &(st->link[dir]), &Unew);
	    theta_rms += dtheta*dtheta;
	    
	  /* compute old action and new action */
	    mult_su3_na(&(st->link[dir]), &(st->staple), &oldUstaple);
	    mult_su3_na(&Unew, &(st->staple), &newUstaple);
	    
	    //newaction = beta*cos(theta_new - theta_staple);
	    oldaction = beta * oldUstaple.e[0][0].real;
	    newaction = beta * newUstaple.e[0][0].real;
	    
	    /* accept or reject */
	    if( newaction > oldaction ){
	       st->link[dir] = Unew;
	       accept++;
	    }
	    else{ /* trace decreased */
	       printf("CALL myrand() in accept/reject update\n");
	       r = myrand(&(st->site_prn)); 
	       //r = drand48(); 
	       if( r < exp( newaction-oldaction ) ){
		  st->link[dir] = Unew;
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
  printf("monte_u1: accept= %d reject= %d fraction= %.2f  theta_rms= %.3f\n",
    accept,reject,accept/(float)(accept+reject),
    sqrt(theta_rms/(accept+reject)) );
  **/

} /* monte_u1 */




void monte_space_u1(int NumStp) {
  int Nhit;
  int parity;
  float scale;		/* limits size of change matrix */
  float theta,theta_rms;	/* rms angle of change matrix */
  float theta_old, theta_new, theta_staple;
  su3_matrix change;	/* proposed change in link */
  su3_matrix newlink;	/* change * oldlink */
  su3_matrix tmpsu3;
  int dir, i;
  register site *st;
  int accept, reject;	/* number of accepts and rejects */
  float oldaction, newaction;

  accept = reject = 0;
  theta_rms = 0.0;
  scale = 1.0;

  for(parity=ODD;parity<=EVEN;parity++) {
    FORSPACEUPDIR(dir) {
      /* compute the gauge force--i.e. load the staples */
      dsdu_qhb(dir,parity); 
      theta_staple = atan2((st->staple.e[0][0].imag), (st->staple.e[0][0].real));

      /* now for the Metropolis updating */
      FORSOMEPARITY(i,st,parity) {
	/* generate random U(1) matrix */

        for( Nhit = 0 ; Nhit < NumStp; Nhit++) { 
	  theta = make_change_u1(&change,st,scale); theta_rms += theta*theta;

	  // could optimize this
	  //	  mult_su3_nn( &change, &(st->link[dir]), &newlink, &tmpsu3 );

	  theta_old = atan2((st->link[dir].e[0][0].imag), (st->link[dir].e[0][0].real));
	  /* compute old action and new action */
	  oldaction=(beta)*(1-cos(theta_old));
	  newaction=(beta)*(1-cos(theta));

	  /* accept or reject */
	  if( newaction > oldaction ){
	     st->link[dir].e[0][0] = cmplx(cos(theta), sin(theta));
	     st->link[dir].e[1][1] = cmplx(cos(theta), -sin(theta));
	    accept++;
	  }
	  else{ /* trace decreased */
	    if( myrand(&(st->site_prn)) < exp( newaction-oldaction ) ){
	       st->link[dir].e[0][0] = cmplx(cos(theta), sin(theta));
	       st->link[dir].e[1][1] = cmplx(cos(theta), -sin(theta));
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

} /* monte_space */




// For updates preserving a flux winding installed on time slice 'tslice', on the 'fluxdir' links
// monte_space_flux() does not update the 'fluxdir' links on t-slice 'tslice' 

void monte_space_flux_u1(int NumStp, int tslice, int fluxdir) {
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

  accept = reject = 0;
  theta_rms = 0.0;
  scale = 1.1;
  for(parity=ODD;parity<=EVEN;parity++) {
    FORSPACEUPDIR(dir) {
      /* compute the gauge force */
      dsdu_qhb(dir,parity); 
      /* now for the Metropolis updating */
      FORSOMEPARITY(i,st,parity) {
	/* generate random SU(3) matrix */
	/* scale < 2/sqrt(3), so vector magnitude < 1 */

	 // leave tslice link[fluxdir]'s alone:

	 if((st->t != tslice) && (dir != 0) && (dir != 1)) {

	    for( Nhit = 0 ; Nhit < NumStp; Nhit++) { 
	       theta = make_change_u1(&change,st,scale); 
	       theta_rms += theta*theta;
	       mult_su3_nn( &change, &(st->link[dir]), &newlink );
	       
	       /* compute old action and new action */
	       oldaction=(0.333333*beta)*realtrace_su3( &(st->link[dir]), 
							&(st->staple) );
	       newaction=(0.333333*beta)*realtrace_su3( &newlink, 
							&(st->staple) );
	       /* accept or reject */
	       if( newaction > oldaction ){
		  st->link[dir]=newlink;
		  accept++;
	       }
	       else{ /* trace decreased */
		  if( myrand(&(st->site_prn)) < exp( newaction-oldaction ) ){
		     st->link[dir]=newlink;
		     accept++;
		  }
		  else{
		     reject++;
		  }
	       }
	    } /* Nhit */
	 } // tslice
      } /*   st */
    } /*  direction */
  } /* parity */


  /* diagnostics */
  /**/
  printf("monte_space: accept= %d reject= %d fraction= %.2f  theta_rms= %.3f\n",
    accept,reject,accept/(float)(accept+reject),
    sqrt(theta_rms/(accept+reject)) );
  /**/

} /* monte_space_flux() */






/* time direction update includes h tr_A P in its action */

void monte_time(int NumStp) {
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

  scale = 1.1;
  mcount = accept = reject = 0;
  theta_rms = 0.0;
  phase_aver = ploop_aver = mag_ploop_aver = density_aver = cmplx(0.0,0.0);

  for(parity=ODD;parity<=EVEN;parity++) {
    /* compute the gauge force */
    dsdu_qhb(TUP,parity); 


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
	  theta = make_change_u1(&change,st,scale); 
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


/* Make change matrix for U(1) Metropolis step.  */

float make_change_u1(su3_matrix *mat, site *st, float scale){
  register int ia,ib,gen;
  register float theta,c,s;
  int p;

  /* load identity into "mat" */
  for(ia=0; ia<3; ia++) {
    for(ib=0; ib<3; ib++){
      mat->e[ia][ib].real = mat->e[ia][ib].imag = 0.0;
    }
    mat->e[ia][ia].real = 1.0;
  }

  printf("CALL myrand() in make_change_u1()\n");

  theta = scale * 2 * (myrand(&(st->site_prn)) - 0.5);
  //theta = scale * 2 * (drand48() - 0.5);
  c = cos(theta); s = sin(theta);
  mat->e[0][0].real = mat->e[1][1].real = c;
  mat->e[0][0].imag = s;
  mat->e[1][1].imag = -s;
  

  //  printf("make_change: theta = %f\n",theta);
  return(theta);
}








/* Make change matrix for Metropolis step.  Choose random generator
  and rotation angle */

float make_change_cartan(su3_matrix *mat, site *st, float scale){
  register int ia,ib,gen;
  register float theta,c,s;

  /* load identity into "mat" */
  for(ia=0; ia<3; ia++) {
    for(ib=0; ib<3; ib++){
      mat->e[ia][ib].real = mat->e[ia][ib].imag = 0.0;
    }
    mat->e[ia][ia].real = 1.0;
  }

  /* random (symmetric) change */
  do{ theta = scale * (myrand(&(st->site_prn)) - 0.5); }
	while( fabs(theta) < 0.25*scale );
  c = cos(theta); s = sin(theta);

  /* random generator, range 1 (gen 3) or 2 (gen 8) */
  gen = ((int)(2.0*myrand(&(st->site_prn))))+1;

  switch(gen){
  case 1:  // gen=3
      mat->e[0][0].real = mat->e[1][1].real = c;
      mat->e[0][0].imag = s;
      mat->e[1][1].imag = -s;
      break;
  case 2:  // gen=8
      mat->e[0][0].real = mat->e[1][1].real = cos(theta/sqrt(3.));
      mat->e[0][0].imag = mat->e[1][1].imag = sin(theta/sqrt(3.));
      mat->e[2][2].real = cos(-2.*theta/sqrt(3.));
      mat->e[2][2].imag = sin(-2.*theta/sqrt(3.));
      break;
  }
/*printf("make_change: theta = %e\n",theta);*/
  return(theta);
}


////////////////////////////////////////////
// this is a placeholder
////////////////////////////////////////////

float make_change_su3(su3_matrix *mat, site *st, float scale){
  register int ia,ib,gen;
  register float theta,c,s;

  /* load identity into "mat" */
  for(ia=0; ia<3; ia++) {
    for(ib=0; ib<3; ib++){
      mat->e[ia][ib].real = mat->e[ia][ib].imag = 0.0;
    }
    mat->e[ia][ia].real = 1.0;
  }

  /* random (symmetric) change */
  //  do{ theta = scale * (myrand(&(st->site_prn)) - 0.5); }
  //	while( fabs(theta) < 0.25*scale );
  //  c = cos(theta); s = sin(theta);

  /* random generator, range 1 (gen 3) or 2 (gen 8) */
  gen = ((int)(2.0*myrand(&(st->site_prn))))+1;

  switch(gen){
  case 1:  // gen=3
      mat->e[0][0].real = mat->e[1][1].real = c;
      mat->e[0][0].imag = s;
      mat->e[1][1].imag = -s;
      break;
  case 2:  // gen=8
      mat->e[0][0].real = mat->e[1][1].real = cos(theta/sqrt(3.));
      mat->e[0][0].imag = mat->e[1][1].imag = sin(theta/sqrt(3.));
      mat->e[2][2].real = cos(-2.*theta/sqrt(3.));
      mat->e[2][2].imag = sin(-2.*theta/sqrt(3.));
      break;
  }
}


