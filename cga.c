#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

uint32_t* freq;
uint32_t taumin, taumax, taumaxminus1K, tauminplus1K, maxminus1K, update;
int *x, *y, *tmp;
unsigned long violations;
int i,j,K,n;
int repeat = 3000;
int maxiter = 250000;
int iter,aver;
int arg1, arg2, arg3, arg4;
int mx,my,cx,cy,r;
int bestfit;





/* Now three random number generators follow */

/* PCG32: very fast */

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}




/* WELL512: not so fast */

#define W 32
#define R 16
#define P 0
#define M1 13
#define M2 9
#define M3 5

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT3NEG(t,v) (v<<(-(t)))
#define MAT4NEG(t,b,v) (v ^ ((v<<(-(t))) & b))

#define V0            STATE[state_i                   ]
#define VM1           STATE[(state_i+M1) & 0x0000000fU]
#define VM2           STATE[(state_i+M2) & 0x0000000fU]
#define VM3           STATE[(state_i+M3) & 0x0000000fU]
#define VRm1          STATE[(state_i+15) & 0x0000000fU]
#define VRm2          STATE[(state_i+14) & 0x0000000fU]
#define newV0         STATE[(state_i+15) & 0x0000000fU]
#define newV1         STATE[state_i                 ]
#define newVRm1       STATE[(state_i+14) & 0x0000000fU]

#define FACT 2.32830643653869628906e-10

static unsigned int state_i = 0;
static unsigned int STATE[R];
static unsigned int z0, z1, z2;

void InitWELLRNG512a (unsigned int *init){
   int j;
   state_i = 0;
   for (j = 0; j < R; j++)
     STATE[j] = init[j];
}

double WELLRNG512a (void){
  z0    = VRm1;
  z1    = MAT0NEG (-16,V0)    ^ MAT0NEG (-15, VM1);
  z2    = MAT0POS (11, VM2)  ;
  newV1 = z1                  ^ z2;
  newV0 = MAT0NEG (-2,z0)     ^ MAT0NEG(-18,z1)    ^ MAT3NEG(-28,z2) ^ MAT4NEG(-5,0xda442d24U,newV1) ;
  state_i = (state_i + 15) & 0x0000000fU;
  return ((double) STATE[state_i]) * FACT;
}





/* Mersenne Twister, similar speed as WELL512*/

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */

void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
        (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */






int main(int argc, char* argsv[])
{

 
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4; // init MT

  unsigned int wellinit[16] = {786059013,787674301,621199908,607659975,461951514,378960130,436303928,498836333,918757912,635468445,26789061,736223013,380898898,471176312,272161869,711242386};

  InitWELLRNG512a(wellinit); // init Well512

  pcg32_random_t mypcg_state = {1786059013,6687694301}; // init PCG
    

  if(argc != 5)
    {
      printf("Usage program n minK maxK repetitions\n");
      exit(1);
    }

  arg1=atoi(argsv[1]);
  arg2=atoi(argsv[2]);
  arg3=atoi(argsv[3]);
  arg4=atoi(argsv[4]);

  printf("Using n=%d minK=%d maxK=%d repeat=%d\n",arg1,arg2,arg3,arg4);

  init_by_array(init, length);

  repeat=arg4;
        
  for(n = arg1; n<=2*arg1; n+=arg1/10)
    {
      taumax = UINT32_MAX - UINT32_MAX / n;
 
      taumin = UINT32_MAX / n;
     
      freq = malloc(sizeof(uint32_t)*n);

      for(K = arg2; K <= arg3; K += 1)
    {

      update = UINT32_MAX / K;

      taumaxminus1K = taumax - update;

      tauminplus1K = taumin + update;
        
      
      x = (int*) malloc(sizeof(int)*n);
      y = (int*) malloc(sizeof(int)*n);

      violations = 0;
      aver=0;
            
      for(r=0; r<repeat; r++)
        {
          for(i=0; i<n; i++)
        {
          freq[i] = UINT32_MAX / 2;
        }
          
          bestfit = 0;
          iter = 0;
          maxiter = n * K  * 20;
          bestfit = -1;
          
          while((bestfit < n) && (iter < maxiter))
        {
          
          cx=cy=0;
          for(i=0; i<n; i++)
            {
              mx = (x[j] = (pcg32_random_r(&mypcg_state) < freq[i] ) ? 1 : 0);
              my = (y[i] = (pcg32_random_r(&mypcg_state) < freq[i] ) ? 1 : 0);
              cx += mx;
              cy += my;
            }
          if(cx > bestfit) bestfit = cx;
          if(cy > bestfit) bestfit = cy;

          
          if(cx < cy) {
            tmp = x;
            x = y;
            y = tmp;
          }
          

                        
          for(i=0;i<n;i++)
            {
              if(x[i] > y[i])
            {
              if(freq[i] < taumaxminus1K)
                {
                  freq[i] += update;
                }
              else
                {
                  freq[i] = taumax;
                }
            }
              else
            {
              if(x[i] < y[i])
                {
                  if(freq[i] > tauminplus1K)
                {
                  freq[i] -= update;
                }
                  else
                {
                  freq[i] = taumin;
                  violations++;
                }
                }
            }
            }
                        
          iter++;
        } // until bestfit=n
          
          aver += iter;
          
        } // repeat x times
      
      
      printf("Average at n=%d, K=%d, repeat=%d is %.2f; violations=%.2f\n", n, K, repeat, ((double) aver / repeat), ((double) violations) / repeat);
      

    } // loop over lambda
      
      free(freq);
      free(x);
      free(y);
    } // loop over n
  
}



