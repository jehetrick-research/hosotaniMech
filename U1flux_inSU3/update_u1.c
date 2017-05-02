/********** update_ora.c ****************************************************/
/* MIMD version 6 */

/*
 Update lattice with microcanonical overrelaxed algorithm
*/

#include "../generic_pg/generic_pg_includes.h"
#include "su3_trP_includes.h"

int update_u1()  {
   int iters=0;

   /* check unitarity before doing anything */
   check_unitarity();
   
   /* do "stepsQ" metropolis updates */
   monte_u1(stepsQ);
   //   monte_space_u1(stepsQ);
   //   monte_time_u1(stepsQ);
   
   /* reunitarize the gauge field */
   
   reunitarize_u1();         
   
   if(steps > 0) return (iters/steps);
   else return(-99);
}


#ifdef COMMENT

int update_cartan()  {
   int iters=0;

   /* check unitarity before doing anything */
   check_unitarity();
   
   /* do "stepsQ" metropolis updates */
   monte_space_cartan(stepsQ);
   monte_time_cartan(stepsQ);
   
   /* reunitarize the gauge field */
   
   reunitarize();         
   
   if(steps > 0) return (iters/steps);
   else return(-99);
}

#endif
