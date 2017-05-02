#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"

/* prototypes for functions in high level code */
int setup();
int readin(int prompt);
int congrad( int niter, float rsqmin, int parity, float *rsq );
void scalar_mult_latvec(field_offset src, float scalar,
			field_offset dest, int parity);
void scalar_mult_add_latvec(field_offset src1, field_offset src2,
			    float scalar, field_offset dest, int parity);
void grsource(int parity);
void reunitarize();

int update();
double d_action();

void dslash( field_offset src, field_offset dest, int parity );
void dslash_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );

void rephase( int flag );
int update_dense();
void measure();
void dsdu_qhb(int dir1,int parity);
void ploop_less_slice(int time,int parity);

void plaq_twist(double w, int gen, double *ss_plaq,double *st_plaq);
void fluxplane(double w, int gen);
void dsdu_twist(int dir1,int parity, int gen);
void get_flux();
