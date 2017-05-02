#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defines.h"

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
void measure();
void dsdu_qhb(int dir1,int parity);

// needed for trP action
void ploop_less_slice(int time,int parity);


// U(1) stuff
float make_change_u1(su3_matrix *change, site *s, float scale);
void monte_u1(int NumStp);
int update_u1();
void reunitarize_u1();


// twisted BCs, and flux stuff
void flux(int tslice, int fluxdir, int q, int gen);
void monte_flux(double w, int NumStp, int tslice, int fluxdir);
void monte_space_flux(double w, int NumStp, int tslice, int fluxdir);
void monte_time_twist(double w, int NumStp);
int update_except_fluxplane();
void plaq_twist(double w, int gen, double *ss_plaq,double *st_plaq);
void load_plaq(double *ss_plaq,double *st_plaq);
void load_plaq_w(double w, double *ss_plaq,double *st_plaq);
void fluxplane(double w, int gen);
void dsdu_twist(float w, int dir1, int parity, int gen);
void get_flux(int zplane, int loadPlaq);
void get_flux_w(double w, int zplane, int loadPlaq);
void set_wphase(double w, int gen, int cc, int y, su3_matrix *wphase);
int update_flux(double w);
void meas_flux(double **mflux, double w);
