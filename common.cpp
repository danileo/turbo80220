#include <limits>
#include <climits>
#include <cstring>
#define D 0.7071067812

using namespace std;

const double infinity=std::numeric_limits<double>::infinity();

// 4-QAM symbol table
double symbol_table_I[]={D,-D,D,-D};
double symbol_table_Q[]={D,D,-D,-D};

/* bit-reverse the input, considering only the least significant n bits,
 *    others bit must be 0 and will be 0
 * v: input value to be reversed
 * n: bit to consider for reversing
 * return the bit-reversed input
 */
unsigned int bit_reversal(unsigned int v, unsigned int n)
{
	unsigned int r = v; // r will be reversed bits of v; first get LSB of v
	int s = sizeof(v) * CHAR_BIT - 1; // extra shift needed at end

	for (v >>= 1; v; v >>= 1)
	{   
		r <<= 1;
		r |= v & 1;
		s--;
	}
	r <<= s; // shift when v's highest bits are zero
	
	return r>>(sizeof(v)*CHAR_BIT-n);
}

/* generate the array used for interleaving using PBRI
 * result[]: interleaving array of length length
 * length: length of result[]
 */
void pruned_bit_reversal_interleaver_core(int result[], unsigned int length)
{
	unsigned int n=0;
	while(length>(((unsigned int)1)<<n))
		n++;
	unsigned int i=0;
	unsigned int j=0;
	while(i<length)
	{
		unsigned int x=bit_reversal(j,n);
		if(x<length)
			result[i++]=x;
		j++;
	}
}

/* generate the array used in the Common real scrambling algorithm
 * output[]: array to be used in the scrambling algorithm of length length
 * length: length of output
 */
void generate_scrambling_random_sequence(unsigned char output[], int length)
{
	// seed to initialize the algorithm, must be of length 20
	unsigned char state[]={1,0,1,0,0,1,0,0,1,1,1,1,0,1,0,1,0,1,0,0};
	
	for(int i=0;i<length;i++)
	{
		output[i]=state[19]^state[18]^state[15]^state[13];
		memmove(&state[1],state,19*sizeof(unsigned char));
		state[0]=output[i];
	}
}
