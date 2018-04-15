#include <iostream>
using namespace std;
#define NUM_STATES 8 // number of states = 2^NU
#define NU 3
#define N 3 // output bits for each input bit in a constituent encoder

// matrix definitions
// S function: SA*u+[SM1;SM2;SM3]*s
#define SA 4
#define SM1 3
#define SM2 4
#define SM3 2
// O function: OA*u+[OM1;OM2;OM3]*s
#define OA 7
#define OM1 0
#define OM2 6
#define OM3 4

static const unsigned char O_helper[2]={0,OA};
static unsigned char O_helper2[8];
static const unsigned char S_helper[2]={0,SA};
static unsigned char S_helper2[8];

static const bool ParityTable256[256] = 
{
#   define P2(n) n, n^1, n^1, n
#   define P4(n) P2(n), P2(n^1), P2(n^1), P2(n)
#   define P6(n) P4(n), P4(n^1), P4(n^1), P4(n)
	P6(0), P6(1), P6(1), P6(0)
};

/* get parity bit of an unsigned char */
unsigned char getParity(unsigned char v)
{
	return ParityTable256[v];
}

/* build auxiliary matrices for the calculation of output function
 * and state update function */
void buildHelpers()
{
	for(unsigned char s=0;s<NUM_STATES;s++)
	{
		O_helper2[s]=(getParity(OM1&s)<<2)|(getParity(OM2&s)<<1)|(getParity(OM3&s));
		S_helper2[s]=(getParity(SM1&s)<<2)|(getParity(SM2&s)<<1)|(getParity(SM3&s));
	}
}

/* output function
 * u: input bit
 * s: current state
 * return the coded bit
 */
unsigned char O(unsigned char u,unsigned char s)
{
	return O_helper[u]^O_helper2[s];
}

/* state update function
 * u: input bit
 * s: current state
 * return the new state
 */
unsigned char S(unsigned char u,unsigned char s)
{
	return S_helper[u]^S_helper2[s];
}

static const int interleavingHelperMatrix[32][8]=
{
	{3, 1, 5, 27,  3,  15,  3,  13},
	{3,  1,  15,  3,  27,  127,  1,  335},
	{3,  3,  5,  1,  15,  89,  5,  87},
	{1,  5,  15,  15,  13,  1,  83,  15},
	{3,  1,  1,  13,  29,  31,  19,  15},
	{1,  5,  9,  17,  5,  15,  179,  1},
	{3,  1,  9,  23,  1,  61,  19,  333},
	{1,  5,  15,  13,  31,  47,  99,  11},
	{1,  3,  13,  9,  3,  127,  23,  13},
	{1,  5,  15,  3,  9,  17,  1,  1},
	{3,  3,  7,  15,  15,  119,  3,  121},
	{1,  5,  11,  3,  31,  15,  13,  155},
	{1,  3,  15,  13,  17,  57,  13,  1},
	{1,  5,  3,  1,  5,  123,  3,  175},
	{1,  5,  15,  13,  39,  95,  17,  421},
	{3,  1,  5,  29,  1,  5,  1,  5},
	{3,  3,  13,  21,  19,  85,  63,  509},
	{1,  5,  15,  19,  27,  17,  131,  215},
	{3,  3,  9,  1,  15,  55,  17,  47},
	{3,  5,  3,  3,  13,  57,  131,  425},
	{3,  3,  1,  29,  45,  15,  211,  295},
	{1,  5,  3,  17,  5,  41,  173,  229},
	{3,  5,  15,  25,  33,  93,  231,  427},
	{1,  5,  1,  29,  15,  87,  171,  83},
	{3,  1,  13,  9,  13,  63,  23,  409},
	{1,  5,  1,  13,  9,  15,  147,  387},
	{3,  1,  9,  23,  15,  13,  243,  193},
	{1,  5,  15,  13,  31,  15,  213,  57},
	{3,  3,  11,  13,  17,  81,  189,  501},
	{1,  5,  3,  1,  5,  57,  51,  313},
	{1,  5,  15,  13,  15,  31,  15,  489},
	{3,  3,  5,  13,  33,  69,  67,  391}
};

#define KEEP_ONLY_LSBs(orig,n) ((orig)&(~((~(unsigned int)0)<<n)))

/* create interleaving array
 * interleaver_array[]: interleaving array
 * mu: length of interleaver_array[]
 */
void interleaving_procedure(int interleaver_array[],int mu)
{
	unsigned int n=0;
	while(mu>(((unsigned int)1)<<(n+5)))
		n++;
	unsigned int counter=0;
	int array_position=0;
	while(array_position<mu)
	{
		unsigned int c=KEEP_ONLY_LSBs((counter>>5)+1,n);
		unsigned int d1=counter&(unsigned int)31;
		unsigned int d2=interleavingHelperMatrix[d1][n-2];
		unsigned int e=KEEP_ONLY_LSBs(c*d2,n);
		
		unsigned char f=d1;
		f = (f * 0x0202020202ULL & 0x010884422010ULL) % 1023;
		f=f>>3;
		
		unsigned int g=((unsigned int)f<<n)+e;
		
		if(g<(unsigned int)mu)
		{
			interleaver_array[array_position]=g;
			array_position++;
		}
		
		counter++;
	}
}
