
#include <stdio.h>
#include <malloc.h>
#include <stdint.h>
#include <math.h>

#include "gsl_rng.h"

/*

https://www.gnu.org/software/gsl/doc/html/rng.html

Function: unsigned long int gsl_rng_get (const gsl_rng * r)
This function returns a random integer from the generator
r. The minimum and maximum values depend on the algorithm
used, but all integers in the range [min,max] are equally
likely. The values of min and max can be determined using
the auxiliary functions gsl_rng_max (r) and gsl_rng_min (r).

Function: double gsl_rng_uniform (const gsl_rng * r)
This function returns a double precision floating point
number uniformly distributed in the range [0,1). The
range includes 0.0 but excludes 1.0. The value is typically
obtained by dividing the result of gsl_rng_get(r) by
gsl_rng_max(r) + 1.0 in double precision. Some generators
compute this ratio internally so that they can provide
floating point numbers with more than 32 bits of randomness
(the maximum number of bits that can be portably represented
in a single unsigned long int).

Function: double gsl_rng_uniform_pos (const gsl_rng * r)
This function returns a positive double precision floating
point number uniformly distributed in the range (0,1),
excluding both 0.0 and 1.0. The number is obtained by
sampling the generator with the algorithm of gsl_rng_uniform
until a non-zero value is obtained. You can use this
function if you need to avoid a singularity at 0.0.

Function: unsigned long int gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n)
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

Function: int gsl_rng_memcpy (gsl_rng * dest, const gsl_rng * src)
This function copies the random number generator src
into the pre-existing generator dest, making dest into
an exact copy of src. The two generators must be of the
same type.

Function: gsl_rng * gsl_rng_clone (const gsl_rng * r)
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



#define WIDTH (32)
#define HEIGHT (32)

void
dump_state(uint8_t *state)
{
    int i, j, bit;
    for(j=0; j<HEIGHT; j++) // 
    {
        for(i=0; i<WIDTH; i++) // 
        {
            bit = state[i + (WIDTH)*j];
            if(bit)
                printf("*");
            else
                printf(".");
        }
        printf("\n");
    }
}

void
main()
{

    gsl_rng * rng;
    // timings for different rng's
    //    rng = gsl_rng_alloc (gsl_rng_taus);    // 4 s
    //    rng = gsl_rng_alloc (gsl_rng_mt19937); // 5 s
    //    rng = gsl_rng_alloc (gsl_rng_ranlxs0); // 30 s "luxury"
    //    rng = gsl_rng_alloc (gsl_rng_cmrg);    // 20 s
    //    rng = gsl_rng_alloc (gsl_rng_taus2);   // 4 s
    //    rng = gsl_rng_alloc (gsl_rng_mrg);     // 11 s
    //    rng = gsl_rng_alloc (gsl_rng_gfsr4);   // 3.3 s

    /*
    double val = 0.0;
    long int i;
    for(i=0; i<(1024*1024*1024); i++)
    {
        val += gsl_rng_uniform (rng);
        //printf("val = %f\n", val);
    }
    */

    rng = gsl_rng_alloc (gsl_rng_gfsr4);   // fastest rng

    uint8_t *state;
    int i, j;
    unsigned long int bit;

    state = (uint8_t *)malloc(WIDTH * HEIGHT * sizeof(uint8_t));

    // initialize
    for(j=0; j<HEIGHT; j++) // 
    for(i=0; i<WIDTH; i++) // 
      {
        bit = gsl_rng_uniform_int(rng, 2);
        state[i + WIDTH*j] = bit;
      }

    dump_state(state);

    int trial;
    double beta = 1.0 / (2.0);

#define TRIALS (1024*1024)
#define get_state(i, j)     ((double)state[(i)%WIDTH + WIDTH*((j)%HEIGHT)])

    for(trial=0; trial<TRIALS; trial++)
    {
        i = gsl_rng_uniform_int(rng, WIDTH);
        j = gsl_rng_uniform_int(rng, HEIGHT);
        
        double field = 0.0;
        field += 2*get_state(i+1, j) - 1.0;
        field += 2*get_state(i, j+1) - 1.0;
        field += 2*get_state(i-1, j) - 1.0;
        field += 2*get_state(i, j-1) - 1.0;
        double p = 1.0/(1 + exp(-2*beta*field));
        double val = gsl_rng_uniform(rng);
        //printf("val = %f, p = %f\n", val, p);
        if(val < p)
            state[i+WIDTH*j] = 1;
        else
            state[i+WIDTH*j] = 0;
    }

    printf("\n"); dump_state(state);

    free(state);
    gsl_rng_free(rng);
}



