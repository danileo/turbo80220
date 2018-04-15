/* turbo_code.cpp
 * turbo encoder
 * MEX-function, use with:
 * y=turbo_code(u);
 * u: input bits
 * y: output coded bits
 */
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include "mex.h"

#include "utility.cpp"
#include "turbo_code.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *c_serial;
    double *u;
    int c_l,u_l;
    
    // get input
    u=mxGetPr(prhs[0]);
    u_l=mxGetN(prhs[0]);
    
    // prepare output
	c_l=u_l*5+18;
    plhs[0]=mxCreateDoubleMatrix(1, c_l, mxREAL);
    c_serial=mxGetPr(plhs[0]);
    
    // call the main algorithm
    code(u,u_l,c_serial);
	
    return;
}

/* turbo coding
 * u[]: input bits
 * mu: length of u[]
 * y_serial: coded bits, length mu*5+18
 */
void code(double u[],int mu,double y_serial[])
{
	unsigned char *y1=new unsigned char[mu+3];
	unsigned char *y2=new unsigned char[mu+3];
	double *u_interleaved=new double[mu];
	
	buildHelpers();
	
	conv(u,y1,mu); // constituent encoder 1
	interleave(u,u_interleaved,mu); //interleaving
	conv(u_interleaved,y2,mu); // constituent encoder 2
	
	// puncturing and multiplexing
	int k=0;
	unsigned char bits[3];
	for(int l=0;l<mu;l++)
	{
		dec2bin(y1[l],bits);
		y_serial[k++]=bits[0];
		y_serial[k++]=bits[1];
		y_serial[k++]=bits[2];
		dec2bin(y2[l],bits);
		y_serial[k++]=bits[1];
		y_serial[k++]=bits[2];
	}
	
	for(int l=mu;l<mu+3;l++)
	{
		dec2bin(y1[l],bits);
		y_serial[k++]=bits[0];
		y_serial[k++]=bits[1];
		y_serial[k++]=bits[2];
	}
	
	for(int l=mu;l<mu+3;l++)
	{
		dec2bin(y1[l],bits);
		y_serial[k++]=bits[0];
		y_serial[k++]=bits[1];
		y_serial[k++]=bits[2];
	}
    
    delete[] y1;
    delete[] y2;
    delete[] u_interleaved;
}

/* decimal to binary conversion
 * a: decimal input value
 * bits[]: binary output value, one bit per element, length 3
 */
void dec2bin(unsigned char a,unsigned char bits[])
{
	bits[2]=a&1;
	bits[1]=(a&2)>>1;
	bits[0]=a>>2;
}

static const unsigned char helper[4]={0,1,1,0};
/* convolutional coding
 * u[]: input bits
 * y[]: output bits, length mu+3
 * mu[]: length of u[]
 */
void conv(double u[],unsigned char y[],int mu)
{
	unsigned char state=0;
	for(int k=0;k<mu;k++)
	{
		y[k]=O(u[k],state);
		state=S(u[k],state);
	}
	unsigned char sym;
	for(int k=mu;k<mu+3;k++)
	{
		sym=helper[state&3];
		y[k]=O(sym,state);
		state=S(sym,state);
	}
}

/* interleaver
 * u[]: input bits
 * u_interleaved[]: interleaved bits
 * mu: length of u[] and u_interleaved[]
 */
void interleave(double u[],double u_interleaved[],int mu)
{
	int *interleaver_array=new int[mu];
	interleaving_procedure(interleaver_array,mu);
	for(int i=0;i<mu;i++)
	{
		u_interleaved[i]=u[interleaver_array[i]];
	}
    delete[] interleaver_array;
}
