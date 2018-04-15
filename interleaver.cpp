using namespace std;

/* interleave data using PBRI algorithm
 * in[]: input data to be interleaved of length length
 * length: length of in[]
 * out[]: interleaved data, same length as in[]
 */
void pruned_bit_reversal_interleaver(unsigned char in[], int length, unsigned char out[])
{
	int *interleaver_array=new int[length];
	pruned_bit_reversal_interleaver_core(interleaver_array,length);
	for(int i=0;i<length;i++)
		out[i]=in[interleaver_array[i]];
	delete[] interleaver_array;
}

/* scramble data using the Common real scrambling algorithm
 * input[]: data to be scrambled of length length
 * output[]: scrambled data, same length as input[]
 * length: length of input[]
 */
void scrambler(unsigned char input[], unsigned char output[], int length)
{
	unsigned char *scrambling_array=new unsigned char[length];
	generate_scrambling_random_sequence(scrambling_array,length);
	for(int i=0;i<length;i++)
			output[i]=input[i]^scrambling_array[i];
	delete[] scrambling_array;
}

/* demultiplex input array into the 5 turbo coder components
 * input[]: input data of length length
 * length: length of input[]
 * U[]: first output component, length (in_length-18)/5+6
 * V0[],V1[],V00[],V11[]: subsequent output components, each of length (in_length-18)/5+3
 */
void demultiplexer(double input[], int length, unsigned char U[], unsigned char V0[], unsigned char V1[], unsigned char V00[], unsigned char V11[])
{
	int Nturbo=(length-18)/5;
	int k=0;
	int i;
	
	for(i=0;i<Nturbo;i++)
	{
		U[i]=input[k++];
		V0[i]=input[k++];
		V1[i]=input[k++];
		V00[i]=input[k++];
		V11[i]=input[k++];
	}
	
	// take care of tail bits
	
	U[i]=input[k];
	U[i+1]=input[k+3];
	U[i+2]=input[k+6];
	U[i+3]=input[k+9];
	U[i+4]=input[k+12];
	U[i+5]=input[k+15];
	
	V0[i]=input[k+1];
	V0[i+1]=input[k+4];
	V0[i+2]=input[k+7];
	
	V1[i]=input[k+2];
	V1[i+1]=input[k+5];
	V1[i+2]=input[k+8];
	
	V00[i]=input[k+10];
	V00[i+1]=input[k+13];
	V00[i+2]=input[k+16];
	
	V11[i]=input[k+11];
	V11[i+1]=input[k+14];
	V11[i+2]=input[k+17];
}

/* join components into a single array to obtain a RATE 1/5 code
 * u[]: first component of length lengthu
 * lengthu: length of u[]
 * v1[],v2[]: second and third components of length lengthv each
 * lengthv: length of v1[] and v2[]
 * result[]: output of joined arrays, it must be of length lengthu+2*lengthv
 */
void join_arrays_1_5(unsigned char u[], int lengthu, unsigned char v1[], unsigned char v2[], int lengthv, unsigned char result[])
{
	memcpy(result,u,lengthu*sizeof(unsigned char));
	memcpy(&result[lengthu],v1,lengthv*sizeof(unsigned char));
	memcpy(&result[lengthu+lengthv],v2,lengthv*sizeof(unsigned char));
}

/* main algorithm for BICM INTERLEAVING
 * in[]: bits to be (bicm)interleaved (= PBRI+scrambling) of length in_length
 * in_length: length of in[]
 * out[]: (bicm)interleaved bits of length in_length (if using join_arrays_1_5)
 */
void bicm_interleaver(double in[], int in_length, unsigned char out[])
{
	int Nturbo=(in_length-18)/5;
	int U_length=Nturbo+6;
	int V_length=Nturbo+3;
	
	unsigned char *U=new unsigned char[U_length];
	unsigned char *V0=new unsigned char[V_length];
	unsigned char *V1=new unsigned char[V_length];
	unsigned char *V00=new unsigned char[V_length];
	unsigned char *V11=new unsigned char[V_length];
	
	demultiplexer(in,in_length,U,V0,V1,V00,V11);
	
	unsigned char *U_perm=new unsigned char[U_length];
	unsigned char *V0_perm=new unsigned char[V_length];
	unsigned char *V1_perm=new unsigned char[V_length];
	unsigned char *V00_perm=new unsigned char[V_length];
	unsigned char *V11_perm=new unsigned char[V_length];
	
	pruned_bit_reversal_interleaver(U,U_length,U_perm);
	pruned_bit_reversal_interleaver(V0,V_length,V0_perm);
	pruned_bit_reversal_interleaver(V1,V_length,V1_perm);
	pruned_bit_reversal_interleaver(V00,V_length,V00_perm);
	pruned_bit_reversal_interleaver(V11,V_length,V11_perm);
	
	delete[] U;
	delete[] V0;
	delete[] V1;
	delete[] V00;
	delete[] V11;
	unsigned char *V0_V00=new unsigned char[2*V_length];
	unsigned char *V1_V11=new unsigned char[2*V_length];
	
	for(int i=0;i<V_length;i++)
	{
		V0_V00[2*i]=V0_perm[i];
		V0_V00[2*i+1]=V00_perm[i];
		V1_V11[2*i]=V1_perm[i];
		V1_V11[2*i+1]=V11_perm[i];
	}
	
	unsigned char *interleaved=new unsigned char[U_length+4*V_length];
	join_arrays_1_5(U_perm,U_length,V0_V00,V1_V11,2*V_length,interleaved);
	
	int out_length=U_length+4*V_length;
	scrambler(interleaved,out,out_length);
	
	delete[] U_perm;
	delete[] V0_perm;
	delete[] V1_perm;
	delete[] V00_perm;
	delete[] V11_perm;
	delete[] V0_V00;
	delete[] V1_V11;
	
	delete[] interleaved;
}
