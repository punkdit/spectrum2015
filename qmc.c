
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "gsl_rng.h"

/*

https://www.gnu.org/software/gsl/doc/html/rng.html

    unsigned long int gsl_rng_get (const gsl_rng * r)
This function returns a random integer from the generator
r. The minimum and maximum values depend on the algorithm
used, but all integers in the range [min,max] are equally
likely. The values of min and max can be determined using
the auxiliary functions gsl_rng_max (r) and gsl_rng_min (r).

    double gsl_rng_uniform (const gsl_rng * r)
This function returns a double precision floating point
number uniformly distributed in the range [0,1). The
range includes 0.0 but excludes 1.0. The value is typically
obtained by dividing the result of gsl_rng_get(r) by
gsl_rng_max(r) + 1.0 in double precision. Some generators
compute this ratio internally so that they can provide
floating point numbers with more than 32 bits of randomness
(the maximum number of bits that can be portably represented
in a single unsigned long int).

    double gsl_rng_uniform_pos (const gsl_rng * r)
This function returns a positive double precision floating
point number uniformly distributed in the range (0,1),
excluding both 0.0 and 1.0. The number is obtained by
sampling the generator with the algorithm of gsl_rng_uniform
until a non-zero value is obtained. You can use this
function if you need to avoid a singularity at 0.0.

    unsigned long int gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n)
This function returns a random integer from 0 to n-1
inclusive by scaling down and/or discarding samples from
the generator r. All integers in the range [0,n-1] are
produced with equal probability. For generators with
a non-zero minimum value an offset is applied so that
zero is returned with the correct probability.

Note that this function is designed for sampling from
ranges smaller than the range of the underlying generator.
The parameter n must be less than or equal to the range
of the generator r. If n is larger than the range of
the generator then the function calls the error handler
with an error code of GSL_EINVAL and returns zero.

In particular, this function is not intended for generating
the full range of unsigned integer values [0,2^32-1].
Instead choose a generator with the maximal integer range
and zero minimum value, such as gsl_rng_ranlxd1, gsl_rng_mt19937
or gsl_rng_taus, and sample it directly using gsl_rng_get.
The range of each generator can be found using the auxiliary
functions described in the next section.

The above methods do not expose the random number ‘state’
which changes from call to call. It is often useful to
be able to save and restore the state. To permit these
practices, a few somewhat more advanced functions are
supplied. These include:

    int gsl_rng_memcpy (gsl_rng * dest, const gsl_rng * src)
This function copies the random number generator src
into the pre-existing generator dest, making dest into
an exact copy of src. The two generators must be of the
same type.

    gsl_rng * gsl_rng_clone (const gsl_rng * r)
This function returns a pointer to a newly created generator
which is an exact copy of the generator r.

*/

/*
    const gsl_rng_type **t, **t0;
    
    t0 = gsl_rng_types_setup ();
    
    printf ("Available generators:\n");
    
    for (t = t0; *t != 0; t++)
    {
        printf ("%s ", (*t)->name);
    }
    printf("\n");
*/
    /* Available generators:
        borosh13 cmrg coveyou fishman18 fishman20 fishman2x gfsr4
        knuthran knuthran2 knuthran2002 lecuyer21 minstd mrg
        mt19937 mt19937_1999 mt19937_1998 r250 ran0 ran1 ran2
        ran3 rand rand48 random128-bsd random128-glibc2 random128-libc5
        random256-bsd random256-glibc2 random256-libc5 random32-bsd
        random32-glibc2 random32-libc5 random64-bsd random64-glibc2
        random64-libc5 random8-bsd random8-glibc2 random8-libc5
        random-bsd random-glibc2 random-libc5 randu ranf ranlux
        ranlux389 ranlxd1 ranlxd2 ranlxs0 ranlxs1 ranlxs2 ranmar
        slatec taus taus2 taus113 transputer tt800 uni uni32
        vax waterman14 zuf 
    */



typedef uint64_t spins_t;
//typedef __uint128_t spins_t;


#include "sse.h"
/*

const int n = 4;
const int mx = 4;
const int mz = 4;
const int nletters = 5;
spins_t Gx[5] = {0, 8, 4, 2, 1};
spins_t Gz[4] = {12, 6, 3, 9};
*/


gsl_rng * rng;


void
init_word(int *word, int order)
{
    int i;
    for(i=0; i<order; i++)
    {
        int idx = gsl_rng_uniform_int(rng, nletters);
        word[i] = idx;
        word[i+order] = idx; // XXX choose random place
    }
}


// See also ex3.c for faster version
int 
countbits(spins_t v)
{
    int c;
    for(c = 0; v; v >>= 1)
    {   
      c += v & 1;
    }   
    return c;
}


void
dump_word(int *word, int order)
{
    int i;
    printf("[");
    for(i=0; i<2*order; i++)
        printf("%d", word[i]);
    printf("]");
}

double
eval1(spins_t u)
{
    double w;
    int j;
    w = offset + mz;
    for(j=0; j<mz; j++)
    {
//        printf("\tw = %f\n", w);
//        printf("\tGz[j] = %ld\n", Gz[j]);
        w -= 2*(countbits(Gz[j] & u) & 1);
    }
    assert(w>=0.0);
    return w;
}


double 
eval_word(int *word, int order, spins_t *psi)
{
    spins_t u, v;
    int i;
    double weight;

    u = 0;
    weight = 1.0;
//    printf("eval_word\n");
    for(i=0; i<2*order; i++)
    {
        int idx = word[i];
        v = u ^ Gx[idx];
//        printf("<%ld|H|%ld>\n", u, v);
        if(u==v)
        {
            double w;
            w = eval1(u);
            weight *= w;
        }
        // else: weight *= Jx
        u = v;
        if(i==order)
            *psi = u; // wavefunction sample
    }
    assert(u==0);
    return weight;
}


int 
eq_word(int *word, int *word1, int order)
{
    int idx;
    for(idx=0; idx<order*2; idx++)
        if(word[idx]!=word1[idx])
            return 0;
    return 1;
}


double
move_word(int *word, int *word1, int order)
{
    double ratio;

    memcpy(word1, word, sizeof(int)*2*order);

    if(gsl_rng_uniform(rng) < 0.5)
    {
        // replace a pair
    
        int idx, jdx, delta;
        int counts[nletters];
        double pairs, fwd, rev;
        

        memset(counts, 0, sizeof(counts));
        for(idx=0; idx<2*order; idx++)
            counts[word[idx]] += 1;
//        printf("counts:");
//        for(idx=0; idx<nletters; idx++)
//            printf("%d ", counts[idx]);
//        printf("\n");

        pairs = 0;
        for(idx=0; idx<nletters; idx++)
        {
            int count = counts[idx];
            assert(count%2==0);
            if(count)
                pairs += count * (count-1) / 2;
        }
        fwd = 1./(4*pairs);

//        dump_word(word, order);
//        printf("pairs = %f\n\n", pairs);


        while(1)
        {
            idx = gsl_rng_uniform_int(rng, 2*order);
            jdx = gsl_rng_uniform_int(rng, 2*order);
            if(idx==jdx || word[idx]!=word[jdx])
                continue;

            delta = gsl_rng_uniform_int(rng, nletters-1) + 1;
            word1[idx] = (word1[idx] + delta) % nletters;
            word1[jdx] = word1[idx];
            break;
        }


        memset(counts, 0, sizeof(counts));
        for(idx=0; idx<2*order; idx++)
            counts[word1[idx]] += 1;
        pairs = 0;
        for(idx=0; idx<nletters; idx++)
        {
            int count = counts[idx];
            assert(count%2==0);
            if(count)
                pairs += count * (count-1) / 2;
        }
        rev = 1./(4*pairs);

        ratio = rev/fwd;

//        dump_word(word1, order);
//        printf("pairs = %f\n\n", pairs);
    }
    else
    {
        // swap two letters
        int idx;
        idx = gsl_rng_uniform_int(rng, 2*order-1);
        spins_t tmp;
        tmp = word1[idx];
        word1[idx] = word1[idx+1];
        word1[idx+1] = tmp;

        ratio = 1.0;
    }

    return ratio;
}


int
main(int argc, char *argv[])
{

    rng = gsl_rng_alloc (gsl_rng_gfsr4);   // fastest rng

    int order, seed, trials, counts;
    order = 10;
    trials = 100;
    counts = 1;

    if(argc>=2)
    {
        order = atoi(argv[1]);
    }
    if(argc>=3)
    {
        trials = atoi(argv[2]);
    }
    if(argc>=4)
    {
        counts = atoi(argv[3]);
    }

    if(argc>=5)
    {
        seed = atoi(argv[4]);
        gsl_rng_set(rng, seed);
    }
    else
    {
        seed = (int)time(NULL);
        gsl_rng_set(rng, seed);
    }

    printf("order = %d\n", order);
    printf("trials = %d\n", trials);
    printf("counts = %d\n", counts);
    printf("seed = %d\n", seed);

    int *support;
    assert(n<=10);
    support = (int *)calloc(powl(2, n), sizeof(int));

    double h;
    int trial, count;
    int *word, *word1;
    word = (int *)malloc(sizeof(int)*order*2);
    word1 = (int *)malloc(sizeof(int)*order*2);

    int accept = 0;

    h = 0;
    for(count=0; count<counts; count++)
    {
        init_word(word, order);
    
        spins_t psi=0, psi1=0;
        double weight, weight1, x, ratio;
    
        weight = eval_word(word, order, &psi);

        for(trial=0; trial<trials; trial++)
        {
//            dump_word(word, order); printf("\n");
            
            ratio = move_word(word, word1, order);
            weight1 = eval_word(word1, order, &psi1);
    
//            dump_word(word1, order); printf("\n");
//            printf("weight = %f\n", weight);
//            printf("weight1 = %f\n", weight1);
//            printf("ratio=%f\n", ratio);
//            printf("%f\n", ratio*weight1/weight);
    
            x = gsl_rng_uniform(rng);
            if(weight==0.0 || x <= (ratio * weight1 / weight))
            {
                if(!eq_word(word, word1, order))
                {
                    accept += 1;
//                    printf("ACCEPT\n");
                }
                memcpy(word, word1, 2*order*sizeof(int));
                weight = weight1;
                psi = psi1;
            }

//            printf("\n");
        }
        x = eval1(psi);

//        dump_word(word, order);
//        printf("weight = %f ", weight);
//        printf("psi = %ld ", psi);
//        printf("x = %f\n", x);

        h += x;

        support[psi] += 1;
    }

    h /= counts;
    printf("<h> = %f\n", h);
    printf("support: \n");
    int norm = support[0];

    int idx;
    int N = powl(2, n);
    for(idx=0; idx<N; idx++)
    {
        printf("%.6f  ", (double)support[idx] / norm);
        if((idx+1)%8 == 0)
            printf("\n");
    }
    printf("\n");

    printf("accept: %f\n", (double)accept / (trials*counts));

    gsl_rng_free(rng);

    return 0;
}



