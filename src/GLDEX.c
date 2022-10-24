/* Part of Steve Su's GLDEX package for the R statistical
 * language.
 *
 *
 * This program is free software; you can redistribute it and/or modify it 
 * under the terms of the GNU General Public License as published by the Free 
 * Software Foundation; either version 3 of the License, or (at your option) 
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA  02111-1307  USA
 *
 */

/* Checking GLD function, single GLD only */
#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

void check_gld(double *a,double *b,double *c,double *d,int *param,int *out){

if(*param==1){

if(*b>0){
 *out=1;
}
else{
 *out=0;}
}

else if (*param==2) {

if((*c < -1) && (*d > 1) && (*b<0)){
 
 *out=1;}
 
else if((*c > 1) && (*d < -1) && (*b<0)){
 
 *out=1;}
 
else if((*d > 1) && (*c > -1) && (*c < 0) && ((pow((1 - *c),
(1 -*c)) * pow((*d - 1),(*d - 1)))/pow((*d - *c),
(*d - *c)) < -*c/ *d) && (*b<0)){
 
 *out=1;}

else if((*c < 0) && (*d <= 0) && *b<0){
 
 *out=1;}
 
else if((*c == 0) && (*d < 0)  && *b<0){
 
 *out=1;}
 
else if((*c > 1) && (*d > -1) && (*d < 0) && ((pow((1 - *d),
(1 -*d)) * pow((*c - 1),(*c - 1)))/pow((*c - *d),
(*c - *d)) < -*d/ *c) && *b<0){
 
 *out=1;}
 
else if((*c > 0) && (*d >= 0) && *b>0){
 
 *out=1;}
 
else if((*c == 0) && (*d > 0) && *b>0){
 
 *out=1;}
 
else{
 
 *out=0;}

}}

/* Checking GLD function, multiple GLDs */
#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

void mult_check_gld(double *a,double *b,double *c,double *d,int *param,int *n,
int *out){
int i;

for(i=0; i<*n; i++){

if(*param==1){

if(b[i]>0){
 out[i]=1;
}
else{
 out[i]=0;}
}

else if (*param==2) {

if((c[i] < -1) && (d[i] > 1) && (b[i]<0)){
 
 out[i]=1;}
 
else if((c[i] > 1) && (d[i] < -1) && (b[i]<0)){
 
 out[i]=1;}
 
else if((d[i] > 1) && (c[i] > -1) && (c[i] < 0) && ((pow((1 - c[i]),
(1 -c[i])) * pow((d[i] - 1),(d[i] - 1)))/pow((d[i] - c[i]),
(d[i] - c[i])) < -c[i]/ d[i]) && (b[i]<0)){
 
 out[i]=1;}

else if((c[i] < 0) && (d[i] <= 0) && b[i]<0){
 
 out[i]=1;}
 
else if((c[i] == 0) && (d[i] < 0)  && b[i]<0){
 
 out[i]=1;}
 
else if((c[i] > 1) && (d[i] > -1) && (d[i] < 0) && ((pow((1 - d[i]),
(1 -d[i])) * pow((c[i] - 1),(c[i] - 1)))/pow((c[i] - d[i]),
(c[i] - d[i])) < -d[i]/ c[i]) && b[i]<0){
 
 out[i]=1;}
 
else if((c[i] > 0) && (d[i] >= 0) && b[i]>0){
 
 out[i]=1;}
 
else if((c[i] == 0) && (d[i] > 0) && b[i]>0){
 
 out[i]=1;}
 
else{
 
 out[i]=0;}

}}}


/* Create a qgl function in C, not really faster so abandoned this but useful training exercise */ 
#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

void q_rs_gld(double *u,double *a,double *b,double *c,double *d,int *n, 
double *out){

int i;

for(i=0; i<*n; i++){
out[i]=*a + (pow(u[i],*c) - pow(1 - u[i],*d))/(*b);
}
  
}

#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

void q_fmkl_gld(double *u,double *a,double *b,double *c,double *d,int *n, 
double *tol, double *out){

int i;

/* c and d are non zeros */
if(fabs(*c)>*tol && fabs(*d)>*tol){
for(i=0; i<*n; i++){
out[i]=*a + ((pow(u[i],*c) - 1)/(*c) - (pow((1 - u[i]),*d) - 1)/(*d))/(*b);}
}

/* c and d are zeros */
if(fabs(*c)<=*tol && fabs(*d)<=*tol){
for(i=0; i<*n; i++){
out[i]=*a + (log(u[i]) - log(1 - u[i]))/(*b);}
}

/* c is zero and d is non zero */
if(fabs(*c)<=*tol && fabs(*d)>*tol){
for(i=0; i<*n; i++){
out[i]=*a + (log(u[i]) - (pow((1 - u[i]),*d) - 1)/(*d))/(*b);}
}

/* c is non zero and d is zero */
if(fabs(*c)>*tol && fabs(*d)<=*tol){
for(i=0; i<*n; i++){
out[i]=*a + ((pow(u[i],*c) - 1)/(*c) - log(1 - u[i]))/(*b);}
}}

 /* Do the min max check for optimisation functions */ 

#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

void q_fmkl_gld_minmax_check(double *mindata, double *maxdata, 
int *lessequalmin, int *greaterequalmax, double *a, 
double *b, double *c,double *d,int *n, double *tol, double *out){

int i;

for(i=0; i<*n; i++){
/* c and d are non zeros */
if (fabs(c[i])>*tol && fabs(d[i])>*tol){

if ((*lessequalmin==1) && (*greaterequalmax==1)){
out[i]= (int) ((a[i] + ((pow(0,c[i]) - 1)/(c[i]) - 
(pow(1,d[i]) - 1)/(d[i]))/(b[i]))<=*mindata) &
(int) ((a[i] + ((pow(1,c[i]) - 1)/(c[i]) - 
(pow(0,d[i]) - 1)/(d[i]))/(b[i]))>=*maxdata);
}

if ((*lessequalmin==1) && (*greaterequalmax==0)){
out[i]= (int) ((a[i] + ((pow(0,c[i]) - 1)/(c[i]) - 
(pow(1,d[i]) - 1)/(d[i]))/(b[i]))<=*mindata) &
(int) ((a[i] + ((pow(1,c[i]) - 1)/(c[i]) - 
(pow(0,d[i]) - 1)/(d[i]))/(b[i]))>*maxdata);
}

if ((*lessequalmin==0) && (*greaterequalmax==1)){
out[i]= (int) ((a[i] + ((pow(0,c[i]) - 1)/(c[i]) - 
(pow(1,d[i]) - 1)/(d[i]))/(b[i]))<*mindata) &
(int) ((a[i] + ((pow(1,c[i]) - 1)/(c[i]) - 
(pow(0,d[i]) - 1)/(d[i]))/(b[i]))>=*maxdata);
}

if ((*lessequalmin==0) && (*greaterequalmax==0)){
out[i]= (int) ((a[i] + ((pow(0,c[i]) - 1)/(c[i]) - 
(pow(1,d[i]) - 1)/(d[i]))/(b[i]))<*mindata) &
(int) ((a[i] + ((pow(1,c[i]) - 1)/(c[i]) - 
(pow(0,d[i]) - 1)/(d[i]))/(b[i]))>*maxdata);
}}

/* c and d are zeros */
if (fabs(c[i])<=*tol && fabs(d[i])<=*tol){
out[i]=1;}

/* c is zero and d is non zero */
if (fabs(c[i])<=*tol && fabs(d[i])>*tol){

if (((*lessequalmin==1) || (*lessequalmin==0)) && (*greaterequalmax==1)){
out[i]=(a[i] + (- (pow(0,d[i]) - 1)/(d[i]))/(b[i]))>=*maxdata;}

if (((*lessequalmin==1) || (*lessequalmin==0)) && (*greaterequalmax==0)){
out[i]=(a[i] + (- (pow(0,d[i]) - 1)/(d[i]))/(b[i]))>*maxdata;}}

/* c is non zero and d is zero */
if (fabs(c[i])>*tol && fabs(d[i])<=*tol){

if (((*greaterequalmax==1) || (*greaterequalmax==0)) && (*lessequalmin==1)){
out[i]=(a[i] + ((pow(0,c[i]) - 1)/(c[i]))/(b[i]))<=*mindata;}

if (((*greaterequalmax==1) || (*greaterequalmax==0)) && (*lessequalmin==0)){
out[i]=(a[i] + ((pow(0,c[i]) - 1)/(c[i]))/(b[i]))<*mindata;}}

}}

#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

void q_rs_gld_minmax_check(double *mindata, double *maxdata, 
int *lessequalmin, int *greaterequalmax, double *a, 
double *b, double *c,double *d,int *n, double *tol, double *out){  

int i;

for(i=0; i<*n; i++){

if ((*lessequalmin==1) && (*greaterequalmax==1)){
out[i]= (int) ((a[i] + (pow(0,c[i]) - 1)/(b[i])) <=*mindata) &
(int) ((a[i] + (1 - pow(0,d[i]))/(b[i]))>=*maxdata);
}

if ((*lessequalmin==1) && (*greaterequalmax==0)){
out[i]= (int) ((a[i] + (pow(0,c[i]) - 1)/(b[i])) <=*mindata) &
 (int) ((a[i] + (1 - pow(0,d[i]))/(b[i]))>*maxdata);
}

if ((*lessequalmin==0) && (*greaterequalmax==1)){
out[i]= (int) ((a[i] + (pow(0,c[i]) - 1)/(b[i])) <*mindata) &
(int) ((a[i] + (1 - pow(0,d[i]))/(b[i]))>=*maxdata);
}


if ((*lessequalmin==0) && (*greaterequalmax==0)){
out[i]= (int) ((a[i] + (pow(0,c[i]) - 1)/(b[i])) <*mindata) &
(int) ((a[i] + (1 - pow(0,d[i]))/(b[i]))>*maxdata);
}}


}



/* pgl functions */

/* gld.rs.fx.c - Part of Robert King's gld package for R, modified by
Steve Su, the instability of using variable==0 when variable is double precision is fixed here. */

#include <stdio.h> 
#include <math.h> 
#include <stdlib.h>

static double la, lb, lc, ld, x, tol;

void funcd(double u,double *a,double *b);

/* the function that finds the root */

void gl_rs_distfunc_p( double *pa,double *pb,double *pc,double *pd, 
double  *px1,double *px2,double *pxacc, int *max_it,
double **ecks, double *u, int *pl,double *tolR)
{
/* pa to pd: 	pointers to the values of the parameters of the gld (rs param)
 * px1:		minimum value of u, should be zero
 * px2:		maximum value of u, should be 1
 * pxacc:	desired accuracy of the calculation
 * max_it: 	maximum iterations for N-R root finder 
 * ecks:	the quantiles of the gld given
 * u:		array to put the calculated depths
 * pl:		length of the data
 */

		
	double  x1, x2, xacc; 		
	double  a, b, c, d;
	int l;

	int i,j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;

	x1 = *px1; x2 = *px2; xacc = *pxacc;
	a = *pa; b = *pb; c = *pc; d = *pd;
	l = *pl, tol = *tolR;

	la = a; lb = b; lc = c; ld = d;

/* Robert King's comment: The C version force the limits to be xacc and 1-xacc
rather than 0 and 1 if lambda3 and lambda4 are negative. */

for (i=0;i<l;i++)
{
        x = *ecks[i];
	u[i] = 0.0;
	funcd(x1,&fl,&df);
	funcd(x2,&fh,&df);


	if (fl*fh >= 0.0) {
 
error("C code numerical failure");

/* 	fprintf(stderr,"Program aborted during calculation of F(x)");
		fprintf(stderr,"at parameter values %e, %e, %e, %e\n", *pa, *pb, *pc, *pd);
		fprintf(stderr,"The x value being investigated was index: %d",i);
		fprintf(stderr," value: %f\n",x);
		exit(1);
*/

		}
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
		}
	else {
		xh = x1;
		xl = x2;
		}
	rts = 0.5*(x1+x2);
	dxold = fabs(x2-x1);
	dx = dxold;
	funcd(rts,&f,&df);
	for (j=1;j<= *max_it;j++) {
		if ((((rts - xh)*df - f)* ( (rts-xl)*df - f) >= 0.0 ) ||
( fabs(2.0*f) > fabs (dxold*df))) {
			dxold = dx;
			dx = .5* (xh - xl);
			rts = xl +dx;
			if (xl == rts ) { u[i] = rts; break; }
			}
		else {
			dxold = dx;
			dx = f/df;
			temp = rts;
			rts -= dx;
			if (temp == rts) { u[i] = rts; break; }
			}
		if (fabs(dx) < xacc) { u[i] = rts; break; }
		funcd(rts,&f,&df);
		if (f < 0.0)
			xl =rts;
		else 
			xh =rts;
		}
}
}


void funcd(double u,double *a,double *b)

{
/* function is the gld F-1(u)  */
/* df_du is dF-1(u)/du  */
double function, df_du;
function=0, df_du=0;

/*  l3 non-zero, l4 non-zero    */
if (fabs(lc)>tol && fabs(ld) >tol){
function = (pow((u),lc) - pow((1.0-u),ld) )/ lb + la - x;
df_du  = (lc * (pow((u),(lc-1.0))) + ld * ( pow((1.0-u),(ld-1.0)) ))/ lb;}
/*      Both l3 and l4 zero     */
if (fabs(lc) <=tol && fabs(ld) <=tol){
function = la - x;
df_du  = 0;}
/*      l3 zero, l4 non-zero    */
if (fabs(lc) <=tol && fabs(ld) >tol){
function = la + ( 1 - pow((1-u),ld) )/ lb - x;
df_du  = (ld) * ( pow((1.0-u),(ld - 1.0)) / lb);}
/*      l4 zero, l3 non-zero     */
if (fabs(lc)>tol && fabs(ld) <=tol){
function = la + ( pow(u,lc) - 1 ) / lb - x;
df_du  = lc * ( pow(u,(lc- 1.0)) ) / lb;}

*a = function;
*b = df_du;

}



//double u,*a,*b;
void funcd2(double u,double *a,double *b)
{

/* function is the gld F-1(u)  */
/* df_du is dF-1(u)/du  */

double function, df_du;

if ( lc == 0 ) {
        if ( ld == 0 ) {
                /*      Both l3 and l4 zero     */
                function = la - x;
                df_du = 0;
                }
        else {
                /*      l3 zero, l4 non-zero    */
                function = la + ( 1 - pow((1-u),ld) )/ lb - x;
                df_du = (ld) * ( pow((1.0-u),(ld - 1.0)) / lb);
                }
        }
else {
        if ( ld == 0 ) {
                /*  l4 zero, l3 non-zero    */
                function = la + ( pow(u,lc) - 1 ) / lb - x;
                df_du = lc * ( pow(u,(lc- 1.0)) ) / lb;
                }
        else {
                /*  l3 non-zero, l4 non-zero    */
                function = ( pow((u),lc) - pow((1.0-u),ld) )/ lb + la - x;
                df_du = (lc * (pow((u),(lc-1.0))) + ld * ( pow((1.0-u),(ld-1.0)) ))/ lb;
                }
        }


*a = function;
*b = df_du;

}




/* gld.fmkl.fx.c Part of Robert King's gld package for R, modified by
Steve Su, the instability of using variable==0 when variable is double precision is fixed here. */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* proto */

void fmkl_funcd( double , double , double *, double *, double *, double *, double *, double *, double *);


/* the function that finds the root */

void gl_fmkl_distfunc_p( double *pa,double *pb,double *pc,double *pd, 
double  *pu1,double *pu2,double *pxacc, int *max_it,
double **ecks, double *u, int *lengthofdata,double *tolr)
{

/* pa to pd:    pointers to the values of the parameters of the gld (rs param)
 * pu1:         minimum value of u, should be zero
 * pu2:         maximum value of u, should be 1
 * pxacc:       desired accuracy of the calculation
 * max_it:      maximum iterations for n-r root finder
 * ecks:        the quantiles of the gld given
 * u:           array to put the calculated depths
 * pl:          length of the data
 */


double  u1, u2, xacc; 		

int i,j;
double df,dx,dxold,f,fh,fl;
double x;
double temp,xh,xl,rts;
/* trying initialising things */
i=0; j=0;
df=0.0;dx=0.0;dxold=0.0;f=0.0;fh=0.0;fl=0.0;
x=0.0;
temp=0.0;xh=0.0;xl=0.0;rts=0.0;

u1 = *pu1; u2 = *pu2; xacc = *pxacc;

if (*pc < 0) {
	if (u1 == 0) {
		u1 = xacc;
		}
	if (u2 == 1) {
		u2 = 1-xacc;
		}	
	}
if (*pd < 0) {
	if (u1 == 0) {
		u1 = xacc;
		}
	if (u2 == 1) {
		u2 = 1-xacc;
		}	
	}


for (i=0;i<*lengthofdata;i++)
{
    x = *ecks[i];
	u[i] = 0.0;
	fmkl_funcd(u1,x,&fl,&df,pa,pb,pc,pd,tolr);
	fmkl_funcd(u2,x,&fh,&df,pa,pb,pc,pd,tolr);
	if (fl*fh >= 0.0) 
	{   error("C code numerical failure");
		/* this is suggested in writing r extensions, but still gives the warning */
		/* error("program aborted at parameter values %f, %f, %f, %f\n the data value being investigated was index %d, value: %f\n", *pa, *pb, *pc, *pd, i, x); */
		/* fprintf(stderr,"program aborted at parameter values %f, %f, %f, %f\n", *pa, *pb, *pc, *pd);
		fprintf(stderr,"the data value being investigated was index %d",i);
		fprintf(stderr," value: %f\n",x);
		exit(1);
		*/
	}
	if (fl < 0.0) {
		xl = u1;
		xh = u2;
		}
	else {
		xh = u1;
		xl = u2;
		}
	rts = 0.5*(u1+u2);
	dxold = fabs(u2-u1);
	dx = dxold;
	fmkl_funcd(rts,x,&f,&df,pa,pb,pc,pd,tolr);
	for (j=1;j<=*max_it;j++) {
		if ((((rts - xh)*df - f)* ( (rts-xl)*df - f) >= 0.0 ) ||
( fabs(2.0*f) > fabs (dxold*df))) {
			dxold = dx;
			dx = .5* (xh - xl);
			rts = xl +dx;
			if (xl == rts ) { 
				u[i] = rts; 
				break; }
			}
		else {
			dxold = dx;
			dx = f/df;
			temp = rts;
			rts -= dx;
			if (temp == rts) { 
				u[i] = rts; 
				break; }
			}
		if (fabs(dx) < xacc) { 
			u[i] = rts; 
			break; }
		fmkl_funcd(rts,x,&f,&df,pa,pb,pc,pd,tolr);
		if (f < 0.0)
			xl =rts;
		else 
			xh =rts;
		}
}
}

void fmkl_funcd( double u, double x, double *F, double *dFdu, double *pa, double *pb, double *pc, double *pd, double *tolR)
{

/* *F is the gld F-1(u)  */            
/* *dFdu is dF-1(u)/du	*/
/* NOTE: *dFdu is 1/f(x) a.k.a. 1/f(F-1(u)) - the opposite of what's
	required for the pdf */

    if (fabs(*pc)>*tolR && fabs(*pd) >*tolR){
		/*  l3 non-zero, l4 non-zero    */
		*F = ( ( pow((u),*pc) -1 )/ *pc  - (pow((1.0-u),*pd) -1 )/ *pd )/ *pb + *pa - x;
		*dFdu = ( ( pow((u),(*pc-1.0)) ) + ( pow( (1.0-u),(*pd-1.0)) ) )/ *pb ;
    }
	
	  if (fabs(*pc) <=*tolR && fabs(*pd) <=*tolR){	
		/*	Both l3 and l4 zero 	*/
		*F = *pa + (log(u) - log(1.0-u))/ *pb - x;
		*dFdu = ( 1 /(u * (1-u))) / *pb  ;  /* correct but confusing, should be 1/u + 1/(1-u) */
		}
	
    if (fabs(*pc) <=*tolR && fabs(*pd) >*tolR){
		/*	l3 zero, l4 non-zero	*/
		*F = *pa + ( log(u) - (( pow((1-u),*pd)-1)/ *pd ) ) / *pb - x;
		*dFdu = (1/u + pow((1-u), *pd-1) )/ *pb ;
		}

    if (fabs(*pc)>*tolR && fabs(*pd) <=*tolR){
		/*  l4 zero, l3 non-zero    */
		*F = *pa + ( (( pow(u,*pc)- 1)/ *pc) - log(1-u) ) / *pb - x;
		*dFdu = (pow(u, *pc-1) + 1/(1-u)  ) / *pb;
		}

}

int cmp_dblp(const void *x1, const void *x2)
{
    double **dx1 = (double **) x1;
    double **dx2 = (double **) x2;

    if (**dx1 < **dx2)
      return -1;
    else if (**dx1 > **dx2)
      return 1;
    else
      return 0;
}


int cmp_dbl_p(const void *x1, const void *x2)
{
    double **dx1 = (double **) x1;
    double **dx2 = (double **) x2;

    if (**dx1 < **dx2)
      return -1;
    else if (**dx1 > **dx2)
      return 1;
    else
      return 0;
}


void pgl (char **param, double *lambda1, double *lambda2, double *lambda3, double *lambda4,
    double *inverse_eps, int *max_iterations, double *q, double *out, int *n, double *tolR)
{
  // the array in between min and max
  double **q_cond = malloc (sizeof(double *) * *n);
  int q_cond_length = 0;

  double civps[2];
  civps[0] = *inverse_eps;
  civps[1] = 1 - *inverse_eps;

  double q_rs_gld_out[2];
  int qdgl_n = 2;
  
  if ( !(strcmp("fmkl", *param)) || !(strcmp("fkml", *param)))
    q_fmkl_gld (civps, lambda1, lambda2, lambda3, lambda4, &qdgl_n, tolR, q_rs_gld_out);
  else if (strcmp("rs", *param) == 0)
    q_rs_gld (civps, lambda1, lambda2, lambda3, lambda4, &qdgl_n, q_rs_gld_out);
  else
  { 
        error("Error: Parameterisation must be fmkl or rs");
		
		/* fprintf (stderr, "Error: Parameterisation must be fmkl or rs");
        exit(EXIT_FAILURE); */
  }
  double min_e = q_rs_gld_out[0];
  double max_e = q_rs_gld_out[1];


  int i = 0; 
  for (i = 0; i < *n; i++)
  {
    out[i] = q[i];
    if (q[i] <= min_e)
      out[i] = 0;
    else if (q[i] >= max_e)
      out[i] = 1;
    // only use pointers to values so that ordering is retained in original
    // after sorting and use
    else 
      q_cond[q_cond_length++] = &out[i];
  }


  // sort q values that lie inbetween 
  qsort(q_cond, q_cond_length, sizeof(double *), cmp_dbl_p);

  double *u = malloc (sizeof(double) * q_cond_length);
  for (i = 0; i < q_cond_length; i++)
    u[i] = 0;

  double pu1 = 0;
  double pu2 = 1;
	    
  if ( !(strcmp("fmkl", *param)) || !(strcmp("fkml", *param)))
    gl_fmkl_distfunc_p(lambda1, lambda2, lambda3, lambda4, &pu1, &pu2,
        inverse_eps, max_iterations, q_cond, u, &q_cond_length, tolR);
  else if (strcmp("rs", *param) == 0)
    gl_rs_distfunc_p(lambda1, lambda2, lambda3, lambda4, &pu1, &pu2,
        inverse_eps, max_iterations, q_cond, u, &q_cond_length, tolR);
  else
  {  
        error("Error: Parameterisation must be fmkl or rs");
        /* fprintf (stderr, "Error: Parameterisation must be fmkl or rs");
        exit(EXIT_FAILURE); */
  }


  for (i = 0; i < q_cond_length; i++)
    *q_cond[i] = u[i];

  free(q_cond);

  free(u);

}


void qdgl_rs (double *lambda1, double *lambda2, double *lambda3, double
    *lambda4, double *u, int *n, double *dens)
{
  for (int i = 0; i < *n; i++)
  {
    dens[i] = *lambda2 / 
      (*lambda3 * pow(u[i], *lambda3-1) + *lambda4 * pow(1-u[i], *lambda4-1));
  }
}

void qdgl_fmkl (double *lambda1, double *lambda2, double *lambda3, double
    *lambda4, double *u, int *n, double *dens)
{
  for (int i = 0; i < *n; i++)
  {
    dens[i] = *lambda2/(pow(u[i], *lambda3-1) + pow(1-u[i], *lambda4-1));
  }
}

void dgl (char **param, double *lambda1, double *lambda2, double *lambda3, 
    double *lambda4, double *inverse_eps, int *max_iterations, double *x, double *out, int *n, double *tolR)
{
  double p[2];
  int p_n = 2;
  p[0] = 0;
  p[1] = 1;

  // 3. Calculate the theoretical minimum and maximum values based on the
  // parameters using qgl.
  double extreme[2];
  if ( !(strcmp("fmkl", *param)) || !(strcmp("fkml", *param)))
    q_fmkl_gld (p, lambda1, lambda2, lambda3, lambda4, &p_n, tolR, extreme);
  else if (strcmp("rs", *param) == 0)
    q_rs_gld  (p, lambda1, lambda2, lambda3, lambda4, &p_n, extreme);
  else
  {  
  error("Error: Parameterisation must be fmkl or rs");
   /* fprintf (stderr, "Error: Parameterisation must be fmkl or rs");
    exit(EXIT_FAILURE); */
  }
 
  double *u = malloc (sizeof(double) * *n);
  // Plug into PGL
  pgl (param, lambda1, lambda2, lambda3, lambda4, inverse_eps,
      max_iterations, x, u, n, tolR);

  if ( !(strcmp("fmkl", *param)) || !(strcmp("fkml", *param)))
    qdgl_fmkl (lambda1, lambda2, lambda3, lambda4, u, n, out);
  else if (strcmp("rs", *param) == 0)
    qdgl_rs (lambda1, lambda2, lambda3, lambda4, u, n, out);
  else
  { error("Error: Parameterisation must be fmkl or rs");
    /* fprintf (stderr, "Error: Parameterisation must be fmkl or rs");
    exit(EXIT_FAILURE); */
  }

  for (int i = 0; i < *n; i++)
  {
    if (x[i] < extreme[0] || x[i] > extreme[1])
      out[i] = 0;
  }
  free(u);
}


void check_gld_m(double *a,double *b,double *c,double *d,char **param,int *out){

 if ( (strcmp("fmkl", *param) == 0) || (strcmp("fkml", *param) == 0)){

if(*b>0){
*out=1;
}
else{
*out=0;}
}

else if (strcmp("rs", *param) == 0){

if((*c < -1) && (*d > 1) && (*b<0)){

*out=1;}

else if((*c > 1) && (*d < -1) && (*b<0)){

*out=1;}

else if((*d > 1) && (*c > -1) && (*c < 0) && ((pow((1 - *c),
(1 -*c)) * pow((*d - 1),(*d - 1)))/pow((*d - *c),
(*d - *c)) < -*c/ *d) && (*b<0)){

*out=1;}

else if((*c < 0) && (*d <= 0) && *b<0){

*out=1;}

else if((*c == 0) && (*d < 0)  && *b<0){

*out=1;}

else if((*c > 1) && (*d > -1) && (*d < 0) && ((pow((1 - *d),
(1 -*d)) * pow((*c - 1),(*c - 1)))/pow((*c - *d),
(*c - *d)) < -*d/ *c) && *b<0){

*out=1;}

else if((*c > 0) && (*d >= 0) && *b>0){

*out=1;}

else if((*c == 0) && (*d > 0) && *b>0){

*out=1;}

else{

*out=0;}

}}


void optim_fun3 (char **param, double *lambda1, double *lambda2, double *lambda3,
  double *lambda4, double *inverse_eps, int *max_iterations, double *data, double *out, int *n, double *tolR)
{
 double p[2];
 int p_n = 2;
 p[0] = 0;
 p[1] = 1;


 double extreme[2];
 int check;
 
 check_gld_m(lambda1, lambda2, lambda3, lambda4,param,&check);

 // if check = 0 set to NaN. Will require checking for NaN in R
 if (check == 0)
 {
   *out = NAN;
   // set to NA
   return;
 }

 // 3. Calculate the theoretical mininum and maximum values based on the
 // parameters using qgl.

 if ( !(strcmp("fmkl", *param)) || !(strcmp("fkml", *param)))
  q_fmkl_gld (p, lambda1, lambda2, lambda3, lambda4, &p_n, tolR, extreme);
 else if (strcmp("rs", *param) == 0)
  q_rs_gld  (p, lambda1, lambda2, lambda3, lambda4, &p_n, extreme);
 else
 { error("Error: Parameterisation must be fmkl or rs");
  /* fprintf (stderr, "Error: Parameterisation must be fmkl or rs");
  exit(EXIT_FAILURE); */
 }

 
 // otherwise if not within bounds, set to NA
 int withinbounds = 1;
 for (int i = 0; i < *n && withinbounds; i++)
 {
   if (extreme[0] > data[i] || extreme[1] < data[i])
     withinbounds = 0;
 }

 if (!withinbounds)
 {
   // set to NA
   *out = NAN;
   return;
 }


   // otherwise must be within bounds and check must not be 1, so sum etc
   double *dgl_out = malloc (sizeof(double) * *n);
   dgl (param, lambda1, lambda2, lambda3, lambda4, inverse_eps, max_iterations, data, dgl_out, n, tolR);
   double sum = 0;
   for (int i = 0; i < *n; i++)
   {
     sum += log(dgl_out[i]);
   }
   free(dgl_out);

   *out = -sum;
}

void optim_fun3_v (char **param, double *lambda1, double *lambda2, double *lambda3,
  double *lambda4, double *inverse_eps, int *max_iterations, double *data, double *out, int *n, double *tolR, int *npar)
{
 for (int i = 0; i < *npar; i++)
  {
    optim_fun3(param,&lambda1[i],&lambda2[i],&lambda3[i],&lambda4[i],inverse_eps,max_iterations,data,&out[i],n,tolR);
  }
}

