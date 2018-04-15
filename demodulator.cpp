/* demodulator.cpp
 * 4-QAM demodulator and channel deinterleaver
 * MEX-function, use with:
 * r=demodulator(rI,rQ,y_l,sig_sq);
 * rI: input modulated symbols, in-phase component
 * rQ: input modulated symbols, in-quadrature component
 * y_l: non-modulated bit sequence length
 * sig_sq: noise variance
 * r: demodulated and deinterleaved bit sequence (the output is 2-PAM modulated)
 */
#include <iostream>
#include <cmath>
#include "mex.h"
#include "common.cpp"
#include "deinterleaver.cpp"

using namespace std;

/* approx. of log(1+e^-x)
 * parameter
 * return approx. of the function
 */
double log_approx(double x)
{
	if(x<2)
        return -0.2831095848*x+0.6931471806;
	if(x<4.3337056)
		return -0.05438904156*x+0.2357060942;
	return 0;
}

/* return the maximum between x and y
 * (slightly better performing than library function) */
double mymax(double x,double y)
{
    return (x > y) ? x : y;
}

/* return log(e^a+e^b) */
double log_sum(double a,double b)
{
	if(a==-infinity||b==-infinity)
		return mymax(a,b);
    return mymax(a,b)+log_approx(abs(a-b));
}

/* get squared norm of the difference between two QAM symbols
 * sigma_sq: noise variance
 * a0: symbol a, real component
 * a1: symbol a, imaginary component
 * b0: symbol b, real component
 * b1: symbol b, imaginary component
 * return (-1/(2*sigma^2))*||modulate(b0+j*b1)-(a0+j*a1)||^2
 */
double distance(double sigma_sq, double a0, double a1, unsigned char b0, unsigned char b1)
{
	unsigned int index=b0|(b1<<1);
	double c0=symbol_table_I[index];
	double c1=symbol_table_Q[index];
	
	double dist=(c0-a0)*(c0-a0)+(c1-a1)*(c1-a1);
	return -dist/(2*sigma_sq);
}

/* demodulate 4-QAM symbols arrays
 * inI[]: in-phase symbols components
 * inQ[]: in-quadrature symbols components
 * symbolLength: length of inI[] and inQ[]
 * usefulLength: demodulated length (it's different from 2*symbolLength if usefulLength is odd)
 * out[]: demodulated symbols mapped using 2-PAM, length usefulLength
 * sigma_sq: noise variance
 */
void demodulator(double inI[], double inQ[], int symbolLength, int usefulLength, double out[], double sigma_sq)
{
	int totBit=0;
	for(int sym=0;sym<symbolLength && totBit<usefulLength;sym++) // symbols index
		for(int bit=0;bit<2 && totBit<usefulLength;bit++) // bits in symbol index
		{
			double loglike0=-infinity;
			double loglike1=-infinity;
			
			// span every possible symbol
			for(unsigned char i=0;i<2;i++)
				for(unsigned char j=0;j<2;j++)
				{
					if(bit==0&&i==0 || bit==1&&j==0)
						loglike0=log_sum(loglike0,distance(sigma_sq,inI[sym],inQ[sym],i,j));
					else if(bit==0&&i==1 || bit==1&&j==1)
						loglike1=log_sum(loglike1,distance(sigma_sq,inI[sym],inQ[sym],i,j));
				}
			
			out[totBit++]=loglike1-loglike0;
		}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *inI;
	double *inQ;
	double *out;
	int inLength,usefulLength;
	double sigma_sq;
	
	// get input
	inI=mxGetPr(prhs[0]);
	inQ=mxGetPr(prhs[1]);
	inLength=mxGetN(prhs[0]);
	usefulLength=mxGetScalar(prhs[2]);
	sigma_sq=mxGetScalar(prhs[3]);
	
	// prepare output
	plhs[0]=mxCreateDoubleMatrix(1, usefulLength, mxREAL);
	out=mxGetPr(plhs[0]);
	
	// call the main algorithm
	double *demodulated=new double[usefulLength];
	demodulator(inI,inQ,inLength,usefulLength,demodulated,sigma_sq);
	bicm_deinterleaver(demodulated,usefulLength,out);
	delete[] demodulated;
	
	return;
}
