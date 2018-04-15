/* modulator.cpp
 * channel interleaver and 4-QAM modulator
 * MEX-function, use with:
 * [sI,sQ]=modulator(y);
 * y: input coded bits
 * sI: output modulated symbols, in-phase component
 * sQ: output modulated symbols, in-quadrature component
 */
#include <iostream>
#include <cmath>
#include "mex.h"
#include "common.cpp"
#include "interleaver.cpp"

using namespace std;

/* modulate bit sequence using QAM
 * in[]: input array of length length, one bit per entry
 * length: length of in[]
 * outI[]: in-phase component of modulated symbols, must have length ceil(length/2.0)
 * outQ[]: in-quadrature component of modulated symbols, same length as outI[]
 */
void modulate(unsigned char in[], int length, double outI[], double outQ[])
{
	for(unsigned int i=0;i<length/2;i++)
	{
		unsigned int index=in[2*i]|(in[2*i+1]<<1);
		outI[i]=symbol_table_I[index];
		outQ[i]=symbol_table_Q[index];
	}
	if(length%2>0)
	{
		unsigned int index=in[2*(length/2)];
		outI[length/2]=symbol_table_I[index];
		outQ[length/2]=symbol_table_Q[index];
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *in;
	double *outI;
	double *outQ;
	int inLength,outLength;
	
	// get input
	in=mxGetPr(prhs[0]);
	inLength=mxGetN(prhs[0]);
	
	// prepare output
	outLength=inLength/2+inLength%2;
	plhs[0]=mxCreateDoubleMatrix(1, outLength, mxREAL);
	outI=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(1, outLength, mxREAL);
	outQ=mxGetPr(plhs[1]);
	
	// call the main algorithm
	unsigned char *interleaved=new unsigned char[inLength];
	bicm_interleaver(in,inLength,interleaved);
	modulate(interleaved,inLength,outI,outQ);
	delete[] interleaved;
	
	return;
}
