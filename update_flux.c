/********** update_ora.c ****************************************************/
/* MIMD version 6 */

/*
 Update lattice with microcanonical overrelaxed algorithm
*/

#include "../generic_pg/generic_pg_includes.h"
#include "su3_trP_includes.h"

int update_flux(double w) {
   int iters=0;

   /* check unitarity before doing anything */
   check_unitarity();
   
   /* do "stepsQ" metropolis updates, preserving links[0,1] at z=0 on tslice=0 */
   // flux tslice=0, fluxdir=0
   monte_space_flux(w, stepsQ, 0, 0);
   monte_time_twist(w, stepsQ);
   
   /* reunitarize the gauge field */
   

   reunitarize();         
   
   if(steps > 0) return (iters/steps);
   else return(-99);
}

