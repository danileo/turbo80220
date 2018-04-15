using namespace std;

/* BICM deinterleave using PBRI
 * in[]: symbols to be deinterleaved
 * length: length of in[] and out[]
 * out[]: deinterleaved symbols
 */
void pruned_bit_reversal_deinterleaver(double in[], int length, double out[])
{
	int *interleaver_array=new int[length];
	pruned_bit_reversal_interleaver_core(interleaver_array,length);
	for(int i=0;i<length;i++)
		out[interleaver_array[i]]=in[i];
	delete[] interleaver_array;
}

/* descramble using the Common real scrambler algorithm
 * input[]: symbols to be descrambled
 * output[]: descrambled symbols
 * length[]: length of input[] and output[]
 */
void descrambler(double input[], double output[], int length)
{
	unsigned char *scrambling_array=new unsigned char[length];
	generate_scrambling_random_sequence(scrambling_array,length);
	for(int i=0;i<length;i++)
		if(scrambling_array[i]==0)
			output[i]=input[i];
		else
			output[i]=-input[i];
	delete[] scrambling_array;
}

/* multiplex the five turbo components in a single array
 * output[]: multiplexed array
 * out_length: length of output[]
 * U[]: first component, length (out_length-18)/5+6
 * V0[],V1[],V00[],V11[]: second, third, fourth and fifth components, length (out_length-18)/5+3 each
 */
void multiplexer(double output[], int out_length, double U[], double V0[], double V1[], double V00[], double V11[])
{
	int Nturbo=(out_length-18)/5;
	int k=0;
	int i;
	
	for(i=0;i<Nturbo;i++)
	{
		output[k++]=U[i];
		output[k++]=V0[i];
		output[k++]=V1[i];
		output[k++]=V00[i];
		output[k++]=V11[i];
	}
	
	// take care of tail bits
	
	output[k]=U[i];
	output[k+3]=U[i+1];
	output[k+6]=U[i+2];
	output[k+9]=U[i+3];
	output[k+12]=U[i+4];
	output[k+15]=U[i+5];
	
	output[k+1]=V0[i];
	output[k+4]=V0[i+1];
	output[k+7]=V0[i+2];
	
	output[k+2]=V1[i];
	output[k+5]=V1[i+1];
	output[k+8]=V1[i+2];
	
	output[k+10]=V00[i];
	output[k+13]=V00[i+1];
	output[k+16]=V00[i+2];
	
	output[k+11]=V11[i];
	output[k+14]=V11[i+1];
	output[k+17]=V11[i+2];
}

/* split rate 1/5 array in the three components to be elaborated
 * u[]: first component
 * lengthu: length of u[]
 * v1[],v2[]: second and third components, length lengthv each
 * lengthv: length of v1[] and v2[]
 * in[]: array to split in the three components, length lengthu+2*lengthv
 */
void split_array_1_5(double u[], int lengthu, double v1[], double v2[], int lengthv, double in[])
{
	memcpy(u,in,lengthu*sizeof(double));
	memcpy(v1,&in[lengthu],lengthv*sizeof(double));
	memcpy(v2,&in[lengthu+lengthv],lengthv*sizeof(double));
}

/* main BICM DEINTERLEAVER algorithm
 * in[]: input array to be deinterleaved
 * in_length: length of in[] and out[]
 * out[]: deinterleaved array
 */
void bicm_deinterleaver(double in[], int in_length, double out[])
{
	int out_length=in_length;
	int Nturbo=(out_length-18)/5;
	int U_length=Nturbo+6;
	int V_length=Nturbo+3;
	
	double *U=new double[U_length];
	double *V0=new double[V_length];
	double *V1=new double[V_length];
	double *V00=new double[V_length];
	double *V11=new double[V_length];
	
	double *U_perm=new double[U_length];
	double *V0_perm=new double[V_length];
	double *V1_perm=new double[V_length];
	double *V00_perm=new double[V_length];
	double *V11_perm=new double[V_length];
	double *V0_V00=new double[2*V_length];
	double *V1_V11=new double[2*V_length];
	
	double *descrambled=new double[in_length];
	descrambler(in,descrambled,in_length);
	
	split_array_1_5(U_perm,U_length,V0_V00,V1_V11,2*V_length,descrambled);
	
	for(int i=0;i<V_length;i++)
	{
		V0_perm[i]=V0_V00[2*i];
		V00_perm[i]=V0_V00[2*i+1];
		V1_perm[i]=V1_V11[2*i];
		V11_perm[i]=V1_V11[2*i+1];
	}
	
	pruned_bit_reversal_deinterleaver(U_perm,U_length,U);
	pruned_bit_reversal_deinterleaver(V0_perm,V_length,V0);
	pruned_bit_reversal_deinterleaver(V1_perm,V_length,V1);
	pruned_bit_reversal_deinterleaver(V00_perm,V_length,V00);
	pruned_bit_reversal_deinterleaver(V11_perm,V_length,V11);
	
	multiplexer(out,out_length,U,V0,V1,V00,V11);
	
	delete[] U;
	delete[] V0;
	delete[] V1;
	delete[] V00;
	delete[] V11;
	
	delete[] U_perm;
	delete[] V0_perm;
	delete[] V1_perm;
	delete[] V00_perm;
	delete[] V11_perm;
	delete[] V0_V00;
	delete[] V1_V11;
	
	delete[] descrambled;
}
