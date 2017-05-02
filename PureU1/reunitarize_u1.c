/********** update_ora.c ****************************************************/
/* MIMD version 6 */

/*
Reunitarize diagonal elements
*/

#include "../generic_pg/generic_pg_includes.h"
#include "su3_trP_includes.h"

void reunitarize_u1() {
   int i, a, dir;
   site *s;
   double mag;
   complex phase;

   FORALLSITES(i, s) {
      for(dir=0; dir<4; dir++) {
	 for(a=0; a<2; a++) {
	    phase = s->link[dir].e[a][a];
	    mag = sqrt(phase.real*phase.real + phase.imag*phase.imag);
	    phase.real /= mag; phase.imag /= mag;
	    s->link[dir].e[a][a] = phase;
	 }
      }
   }
}
