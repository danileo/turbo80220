/* turbo.cpp
 * turbo decoder
 * MEX-function, use with:
 * u_hat=turbo(r,sig_sq);
 * r: input encoded bits
 * sig_sq: noise variance
 * u_hat: decoded bits
 */
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include "mex.h"

using namespace std;

#include "utility.cpp"
#include "turbo.h"

#define MAX_ITERATIONS 20

double conversion_matrix[][N]={
	{-1,-1,-1},{-1,-1,1},{-1,1,-1},{-1,1,1},
	{1,-1,-1},{1,-1,1},{1,1,-1},{1,1,1}};

const double infinity=std::numeric_limits<double>::infinity();
int mu;
double sigma_w2=1;
int nIT=0;
double (*B_matrix)[NUM_STATES];
double (*F_matrix)[NUM_STATES];
double (*q_array)[2];

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *r_serial;
    double *u_hat;
    int r_l,u_l;
	double sigma2;
    
    // get input
    r_serial=mxGetPr(prhs[0]);
    r_l=mxGetN(prhs[0]);
	sigma2=mxGetScalar(prhs[1]);
    
    // prepare output
	u_l=(r_l-18)/5;
    plhs[0]=mxCreateDoubleMatrix(u_l, MAX_ITERATIONS, mxREAL);
    u_hat=mxGetPr(plhs[0]);
    
    // call the main algorithm
    turbo(r_serial, r_l, u_hat, u_l, sigma2);
	
    return;
}

/* turbo decoder
 * r_serial[]: input, 2-PAM modulated, values
 * r_l: length of r_serial
 * u_hat[]: output values, length (r_l-18)/5
 * sigma2: noise variance
 */
void turbo(double r_serial[],int r_l,double u_hat[],int u_l,double sigma2)
{
	mu=u_l;
	sigma_w2=sigma2;
	double (*r1)[N]=new double[mu+NU][N];
	double (*r2)[N]=new double[mu+NU][N];
	double (*q1)[2]=new double[mu+NU][2];
	double (*q2)[2]=new double[mu+NU][2];
	double (*E_matrix1)[2]=new double[mu][2];
	double (*E_matrix2)[2]=new double[mu][2];
	B_matrix=new double[mu+NU][NUM_STATES];
	F_matrix=new double[mu+NU][NUM_STATES];
	
	buildHelpers();
	initialize_q_array(q1);
	initialize_q_array(q2);
	
	int* interleaver_array=new int[mu];
	interleaving_procedure(interleaver_array,mu);
	
	depuncture(r_serial,r_l,r1,r2,u_l);
	
	for(int iteration=0;iteration<MAX_ITERATIONS;iteration++)
	{
		nIT=iteration;
		q_array=q1;
		bcjr(r1,E_matrix1);
		
		marginalization(q1,E_matrix1,u_hat);
		
		propagate(interleaver_array,E_matrix1,q2,true);
		
		q_array=q2;
		bcjr(r2,E_matrix2);
		
		propagate(interleaver_array,q1,E_matrix2,false);
	}
	
	delete[] B_matrix;
	delete[] F_matrix;
	delete[] q1;
	delete[] q2;
	delete[] E_matrix1;
	delete[] E_matrix2;
	delete[] interleaver_array;
    delete[] r1;
    delete[] r2;
}

/* depuncturing/demultiplexing for turbo code
 * serial[]: input values
 * serial_l: length of serial[]
 * r1[][N],r2[][N]: (u_l+NU)xN matrices of demultiplexed bits for const. enc. 1
 *                  and 2, respectively
 * u_l: number of input bits to turbo encoder
 */
void depuncture(double serial[],int serial_l,double r1[][N],double r2[][N],int u_l)
{
	int i;
	int i_serial=0;
	for(i=0;i<u_l;i++)
	{
		r1[i][0]=serial[i_serial++];
		r1[i][1]=serial[i_serial++];
		r1[i][2]=serial[i_serial++];
		r2[i][0]=0;
		r2[i][1]=serial[i_serial++];
		r2[i][2]=serial[i_serial++];
	}
	int i2=i;
	while(i_serial<serial_l-9)
	{
		r1[i][0]=serial[i_serial++];
		r1[i][1]=serial[i_serial++];
		r1[i++][2]=serial[i_serial++];
	}
	while(i_serial<serial_l)
	{
		r2[i2][0]=serial[i_serial++];
		r2[i2][1]=serial[i_serial++];
		r2[i2++][2]=serial[i_serial++];
	}
}

/* propagate extrinsic information between constituent encoders
 * interleaver_array[]: array used for interleaving, length mu
 * q1[][2]: (mu)x(2) matrix, extrinsic information for constituent encoder 1
 * q2[][2]: (mu)x(2) matrix, extrinsic information for constituent encoder 2
 * direct: true if the propagation happens from const. enc. 1 to const. enc. 2,
 *         false otherwise
 */
void propagate(int interleaver_array[],double q1[][2],double q2[][2],bool direct)
{
	for(int i=0;i<mu;i++)
	{
		if(direct)
		{
			q2[i][0]=q1[interleaver_array[i]][0];
			q2[i][1]=q1[interleaver_array[i]][1];
		}
		else
		{
			q1[interleaver_array[i]][0]=q2[i][0];
			q1[interleaver_array[i]][1]=q2[i][1];
		}
	}
}

/* initialize a priori probabilities
 * q[][2]: (mu+NU)x2 matrix of a priori probabilities
 */
void initialize_q_array(double q[][2])
{
	for(int i=0;i<mu;i++)
		q[i][0]=q[i][1]=-0.69314718;
	for(int i=mu;i<mu+NU;i++)
		q[i][0]=q[i][1]=0;
}

/* BCJR algorithm
 * r[][N]: (mu+NU)xN matrix of coded bits from the channel
 * E_matrix[][2]: (mu)x2 matrix of output extrinsic information
 */
void bcjr(double r[][N],double E_matrix[][2])
{
	B_matrix[mu+NU-1][0]=0;
	for(int i=1;i<NUM_STATES;i++)
		B_matrix[mu+NU-1][i]=-infinity;
	F_matrix[0][0]=0;
	for(int i=1;i<NUM_STATES;i++)
		F_matrix[0][i]=-infinity;
	
	for(int l=1;l<mu+NU;l++)
    {
        double tot=-infinity;
		for(unsigned char s=0;s<NUM_STATES;s++){
			F_matrix[l][s]=F(l,s,r);
            tot=log_sum(tot,F_matrix[l][s]);
        }
        for(unsigned char s=0;s<NUM_STATES;s++){
			F_matrix[l][s]-=tot;
        }
    }
	
	for(int l=mu+NU-2;l>=0;l--){
        double tot=-infinity;
		for(unsigned char s=0;s<NUM_STATES;s++){
			B_matrix[l][s]=B(l,s,r);
            tot=log_sum(tot,B_matrix[l][s]);
        }
        for(unsigned char s=0;s<NUM_STATES;s++){
			B_matrix[l][s]-=tot;
        }
    }
	
	double epsilon;
	for(int l=0;l<mu;l++)
	{
		E_matrix[l][0]=E(l,0,r);
		E_matrix[l][1]=E(l,1,r);
		epsilon=log_sum(E_matrix[l][0],E_matrix[l][1]);
		E_matrix[l][0]-=epsilon;
		E_matrix[l][1]-=epsilon;
	}
}

/* backward message
 * l: bit index
 * s: symbol for which calculate the backward message
 * r[][N]: (mu+NU)xN matrix of coded bits from the channel
 * return the backward message
 */
double B(int l,unsigned char s,double r[][N])
{
	double res=-infinity;
	res=q(l+1,0)+B_matrix[l+1][S(0,s)]+g(l+1,O(0,s),r);
	res=log_sum(res,q(l+1,1)+B_matrix[l+1][S(1,s)]+g(l+1,O(1,s),r));
	return res;
}

/* forward message
 * l: bit index
 * s: symbol for which calculate the forward message
 * r[][N]: (mu+NU)xN matrix of coded bits from the channel
 * return the forward message
 */
double F(int l,unsigned char s,double r[][N])
{
	double res=-infinity;
	for(unsigned char s1=0;s1<NUM_STATES;s1++)
	{
		if(s==S(0,s1))
			res=log_sum(res,q(l-1,0)+g(l-1,O(0,s1),r)+F_matrix[l-1][s1]);
        else if(s==S(1,s1))
			res=log_sum(res,q(l-1,1)+g(l-1,O(1,s1),r)+F_matrix[l-1][s1]);
	}
	return res;
}

/* extrinsic information
 * l: bit index
 * u: symbol for which calculate the extrinsic information
 * r[][N]: (mu+NU)xN matrix of coded bits from the channel
 * return the extrinsic information
 */
double E(int l,unsigned char u,double r[][N])
{
	double res=-infinity;
	for(unsigned char s=0;s<NUM_STATES;s++)
		res=log_sum(res,F_matrix[l][s]+B_matrix[l][S(u,s)]+g(l,O(u,s),r));
	return res;
}

/* marginalization
 * q_arr[][2]: (mu+NU)x2 matrix of a priori probabilities
 * E_matrix[][2]: (mu)x2 matrix of extrinsic information
 * u_hat[]: array of decoded symbols, length mu
 */
void marginalization(double q_arr[][2],double E_matrix[][2],double u_hat[])
{
	double val0,val1;
	for(int l=0;l<mu;l++)
	{
		val0=q_arr[l][0]+E_matrix[l][0];
		val1=q_arr[l][1]+E_matrix[l][1];
		if(val0>val1) u_hat[l+mu*nIT]=0;
		else u_hat[l+mu*nIT]=1;
	}
}

/* a priori probabilities
 * l: bit index
 * u: symbol
 * return a priori probability for symbol u at index l
 */
double q(int l, unsigned char u)
{
	return q_array[l][u];
}

/* g function
 * l: bit index
 * y: uncoded tentative symbol
 * r[][N]: (mu+NU)xN matrix of coded bits from the channel
 */
double g(int l,unsigned char y,double r[][N])
{
	double* coded=conversion_matrix[y];
	return -1/(2*sigma_w2)*calculate_diff_squared_norm(r[l],coded);
}

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

/* squared norm of r-coded
 * r[]: array of length N
 * coded[]: array of length N
 * return the squared norm of the difference between r and coded
 */
double calculate_diff_squared_norm(double r[],double coded[])
{
	double res=0;
	for(int i=0;i<N;i++)
		res+=(r[i]-coded[i])*(r[i]-coded[i]);
	return res;
}
