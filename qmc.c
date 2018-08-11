
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "gsl_rng.h"

double
factorial(int n)
{
    double r;
    int i;
    r = 1.0;
    for(i=1; i<=n; i++)
        r *= i;
    return r;
}


// see also:
// https://www.johndcook.com/blog/2010/08/16/how-to-compute-log-factorial/
#define LOG_FACT_MAX (4096)
double _log_fact_cache[LOG_FACT_MAX];

double
_log_factorial(int n)
{
    double r = 0.0;
    int i;
    for(i=1; i<=n; i++)
        r += log((double)i);

    return r;
}

#define LOG_MAX (1024)
double _log_cache[LOG_MAX];

double
ilog(int n)
{
    assert(n>0);
    assert(n-1<LOG_MAX);
    return _log_cache[n-1];
}

void
init_cache()
{
    int n;
    for(n=0; n<LOG_FACT_MAX; n++)
        _log_fact_cache[n] = _log_factorial(n);
    for(n=0; n<LOG_MAX; n++)
        _log_cache[n] = log(n+1);
}

double
log_factorial(int n)
{
    assert(0<=n && n<LOG_FACT_MAX);
    return _log_fact_cache[n];
}


//
// From:
// https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive
//

// Counting bits set by lookup table
const unsigned char BitsSetTable256[256] = 
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};


#if defined(SMALL_PROB)

typedef uint32_t slab_t;
#define SLABS (1)

#elif defined(MED_PROB)

typedef uint64_t slab_t;
#define SLABS (1)

#elif defined(LARGE_PROB)

typedef uint64_t slab_t;
#define SLABS (2)

//typedef __uint128_t slab_t;

#elif defined(HUGE_PROB)

typedef uint64_t slab_t;
#define SLABS (6)

#else
#error "problem size not defined"
#endif


typedef slab_t spins_t[SLABS];

int 
countbits_fast(spins_t v)
{
    uint32_t c = 0;
    uint8_t *p = (uint8_t *)v;
    int i;
    for(i=0; i<sizeof(spins_t); i++)
        c += BitsSetTable256[p[i]];
    return c;
}




#include "model.h"
/*
const int n = 8;
const int mx = 8;
const int offset = 8;
const int mz = 8;
const int nletters = 9;
char *_Gx[8] = {
  "1.......", 
  ".1......", 
  "..1.....", 
  "...1....", 
  "....1...", 
  ".....1..", 
  "......1.", 
  ".......1"
};
char *_Gz[8] = {
  "11......", 
  ".11.....", 
  "..11....", 
  "...11...", 
  "....11..", 
  ".....11.", 
  "......11", 
  "1......1"
};
spins_t Gx[9];
spins_t Gz[8];
int Jx[8] = {1, 1, 1, 1, 1, 1, 1, 1};
int Jz[8] = {1, 1, 1, 1, 1, 1, 1, 1};
*/



void
SPINS_INIT(spins_t u)
{
    memset(u, 0, sizeof(slab_t)*SLABS);
}

void
SPINS_SET(spins_t u, spins_t v)
{
    memcpy(u, v, sizeof(slab_t)*SLABS);
}

void
SPINS_AND(spins_t u, spins_t v, spins_t w)
{
    int i;
    for(i=0; i<SLABS; i++)
        u[i] = v[i] & w[i];
}

void
SPINS_XOR(spins_t u, spins_t v, spins_t w)
{
    int i;
    for(i=0; i<SLABS; i++)
        u[i] = v[i] ^ w[i];
}

int
SPINS_EQ(spins_t u, spins_t v)
{
    int i;
    for(i=0; i<SLABS; i++)
        if(u[i] != v[i])
            return 0;
    return 1;
}


slab_t
bits2slab(char *bits)
{
    slab_t result = 0;
    int idx;
    for(idx = 0; idx<sizeof(slab_t)*8; idx++)
    {
        if(bits[idx]==0)
            break; // null terminator
        result *= 2;
        if(bits[idx]=='1')
            result += 1;
    }
    return result;
}

void
bits2spins(spins_t u, char *bits)
{
    int idx;
    for(idx=0; idx<SLABS; idx++)
    {
        u[idx] = bits2slab(bits);
        bits += sizeof(slab_t)*8;
    }
}


void
init_ops()
{
    int idx;
    SPINS_INIT(Gx[0]); // first one is identity op
    for(idx=0; idx<mx; idx++)
        bits2spins(Gx[idx+1], _Gx[idx]);
    for(idx=0; idx<mz; idx++)
        bits2spins(Gz[idx], _Gz[idx]);
}


gsl_rng * rng;

#define MAX_WORD (4096)

struct state 
{
    spins_t u;
    int size;
    int word[MAX_WORD];
};

void
check_state(struct state *s)
{
    assert(s);
//    assert(s->u<(((spins_t)1)<<n));
    assert(0<=s->size);
    assert(s->size<=MAX_WORD);
}

int verbose;


void
init_state(struct state *s)
{
    assert(s);
    SPINS_INIT(s->u);
    s->size = 0;
    memset(s->word, 0, sizeof(int)*MAX_WORD);
    check_state(s);
}

void
state_insert(struct state *s, int idx, int i);

void
dump_state(struct state *s)
{
    int i;
    check_state(s);
    //printf("%ld", (long)s->u);
    printf("[");
    for(i=0; i<s->size; i++)
        printf("%d,", s->word[i]);
    printf("]");
}

void
init_state_rand(struct state *s, double beta)
{
    int idx;
    assert(s);
    init_state(s);

    // init one byte at a time
    uint8_t *u = (uint8_t *)s->u;
    int nn = n;
    for(idx=0; idx<sizeof(spins_t) && nn>0; idx++)
    {
        if(nn>=8)
            *u = gsl_rng_uniform_int(rng, 1<<8);
        else
            *u = gsl_rng_uniform_int(rng, 1<<nn);
        u++;
        nn -= 8;
    }
    assert(nn<=0);

/*
    int jdx;
    for(jdx=0; jdx<1; jdx++)
      for(idx=0; idx<offset; idx++)
        {
            int letter, i;
            letter = gsl_rng_uniform_int(rng, nletters);
            i = gsl_rng_uniform_int(rng, s->size+1);
            state_insert(s, i, letter);
            i = gsl_rng_uniform_int(rng, s->size+1);
            state_insert(s, i, letter);
        }

    if(verbose>1)
    { dump_state(s); printf("\n");}
*/

    check_state(s);
}


double
eval1(spins_t u)
{
    spins_t v;
    double w;
    int j;
    w = offset;
    for(j=0; j<mz; j++)
    {
        w += Jz[j];
        SPINS_AND(v, Gz[j], u);
        w -= 2*Jz[j]*(countbits_fast(v) & 1);
    }
    assert(w>=0.0);
    return w;
}


double 
eval_state(struct state *s, double beta)
{
    spins_t u, v;
    int i;
    double weight;

    check_state(s);
    SPINS_SET(u, s->u);
    weight = 1.0;
//    printf("eval_word\n");
    for(i=0; i<s->size; i++)
    {
        int idx = s->word[i];
        //v = u ^ Gx[idx];
        SPINS_XOR(v, u, Gx[idx]);
//        printf("<%ld|H|%ld>\n", u, v);
        if(SPINS_EQ(u, v))
        {
            double w;
            w = eval1(u);
            weight *= w;
        }
        else
        {
            weight *= Jx[idx];
        }
        //u = v;
        SPINS_SET(u, v);
    }
    //assert(u==s->u);
    assert(SPINS_EQ(u, s->u));
    weight *= 1.0 * pow(beta, s->size) / factorial(s->size);
    return weight;
}


int 
log_eval_state(struct state *s, double log_beta, double *result)
{
    spins_t u, v;
    int i;
    double weight, c_weight, w;
    int changed = 1;

    check_state(s);
    *result = 0.;
    SPINS_SET(u, s->u);
    weight = 0.0;
    for(i=0; i<s->size; i++)
    {
        int idx = s->word[i];
        //v = u ^ Gx[idx];
        if(idx)
            SPINS_XOR(v, u, Gx[idx]);
        if(idx==0 || SPINS_EQ(u, v))
        {
            if(changed)
                w = eval1(u);
            else
                w = c_weight;
            assert(w>=0.0);
            if(w==0.0)
                return 0; // result is log(zero)
            weight += ilog(w);
            changed = 0;
            c_weight = w;
        }
        else
        {
            assert(idx);
            weight += ilog(Jx[idx-1]);
            SPINS_SET(u, v);
            changed = 1;
        }
    }
    assert(SPINS_EQ(u, s->u));
    weight += s->size * log_beta - log_factorial(s->size);
    *result = weight;
    return 1;
}


int 
eq_state(struct state *s1, struct state *s2)
{
    int idx;
    check_state(s1);
    check_state(s2);
//    if(s1->u != s2->u)
    if(!SPINS_EQ(s1->u, s2->u))
        return 0;
    if(s1->size != s2->size)
        return 0;
    for(idx=0; idx<s1->size; idx++)
        if(s1->word[idx]!=s2->word[idx])
            return 0;
    return 1;
}


void
state_insert(struct state *s, int idx, int i)
{
    int j;
    check_state(s);
    assert(0<=idx);
    assert(idx<=s->size);
    assert(s->size < MAX_WORD);
    assert(0<=i && i<nletters);

    // use memmove ?
    for(j=s->size; j>idx; j--)
        s->word[j] = s->word[j-1];
    s->word[idx] = i;
    s->size++;

    check_state(s);
}


int
state_pop(struct state *s, int idx)
{
    check_state(s);

    int i, j;
    check_state(s);
    assert(0<=idx);
    assert(idx<s->size);
    assert(s->size > 0);
    i = s->word[idx];

    // use memmove ?
    for(j=idx; j+1<s->size; j++)
        s->word[j] = s->word[j+1];
    s->size--;
    check_state(s);
    return i;
}


int
state_count(struct state *s, int i)
{
    int j, count;
    check_state(s);
    assert(0<=i && i<nletters);

    count = 0;
    for(j=0; j<s->size; j++)
        if(s->word[j]==i)
            count++;
    return count;
}


double
move_state(struct state *s, struct state *s1)
{
    double ratio;
    int i, idx, jdx, total;

    check_state(s);
    assert(s1);

    memcpy(s1, s, sizeof(struct state));

//    printf("move_state ");
//    dump_state(s);
//    printf("\n");

    i = gsl_rng_uniform_int(rng, 5);
//    printf("i=%d ", i);

    ratio = 1.0;

    if(i==0)
    {
        // add I
        total = state_count(s, 0);
        idx = gsl_rng_uniform_int(rng, s->size+1);
        ratio = ((double)(s->size+1)) / (total+1);
//        printf(" idx=%d ", idx);
        state_insert(s1, idx, 0);
    }
    else
    if(i==1)
    {
        // remove I
        total = state_count(s, 0);
        if(total)
        {
            while(1)
            {
                idx = gsl_rng_uniform_int(rng, s->size);
                if(s->word[idx]==0)
                    break;
            }
            ratio = ((double)total) / s->size;
            state_pop(s1, idx);
        }
    }
    else
    if(i==2 && s->size>1)
    {
        // replace a pair
        double fwd, rev;

        int counts[nletters];
        int pairs;

        memset(counts, 0, sizeof(counts));
        for(idx=0; idx<s1->size; idx++)
            counts[s1->word[idx]] += 1;

        pairs = 0;
        for(idx=0; idx<nletters; idx++)
        {
            int count = counts[idx];
            assert(idx==0 || count%2==0);
            if(count)
                pairs += count * (count-1) / 2; // XX check overflow
        }
        fwd = 1./pairs;

        while(1)
        {
            idx = gsl_rng_uniform_int(rng, s->size);
            jdx = gsl_rng_uniform_int(rng, s->size);
            if(idx != jdx && s->word[idx]==s->word[jdx])
            {
                int delta = gsl_rng_uniform_int(rng, nletters-1);
                int letter = (s->word[idx] + delta) % nletters;
                s1->word[idx] = letter;
                s1->word[jdx] = letter;
                break;
            }
        }


        memset(counts, 0, sizeof(counts));
        for(idx=0; idx<s1->size; idx++)
            counts[s1->word[idx]] += 1;

        pairs = 0;
        for(idx=0; idx<nletters; idx++)
        {
            int count = counts[idx];
            assert(idx==0 || count%2==0);
            if(count)
                pairs += count * (count-1) / 2; // XX check overflow
        }
        rev = 1./pairs;

        ratio = rev/fwd;
    }
//
//    if(i==2)
//    {
//        // add a pair
//        idx = gsl_rng_uniform_int(rng, s->size+1);
//        jdx = gsl_rng_uniform_int(rng, s->size+2);
//        i = gsl_rng_uniform_int(rng, nletters);
//        state_insert(s1, idx, i);
//        state_insert(s1, jdx, i);
//
//        int counts[nletters];
//        int pairs;
//
//        memset(counts, 0, sizeof(counts));
//        for(idx=0; idx<s1->size; idx++)
//            counts[s1->word[idx]] += 1;
//
//        pairs = 0;
//        for(idx=0; idx<nletters; idx++)
//        {
//            int count = counts[idx];
//            assert(idx==0 || count%2==0);
//            if(count)
//                pairs += count * (count-1) / 2; // XX check overflow
//        }
//
//        // XXX This calculation is not exactly right
//        ratio = (0.5/pairs) * (s->size+1) * (s->size+2) * nletters;
//    }
//    else
//    if(i==3 && s->size>1)
//    {
//        // remove a pair
//        int counts[nletters];
//        int pairs;
//
//        memset(counts, 0, sizeof(counts));
//        for(idx=0; idx<s1->size; idx++)
//            counts[s1->word[idx]] += 1;
//
//        pairs = 0;
//        for(idx=0; idx<nletters; idx++)
//        {
//            int count = counts[idx];
//            assert(idx==0 || count%2==0);
//            if(count)
//                pairs += count * (count-1) / 2; // XX check overflow
//        }
//
//        while(1)
//        {
//            idx = gsl_rng_uniform_int(rng, s->size);
//            jdx = gsl_rng_uniform_int(rng, s->size);
//            if(idx != jdx && s->word[idx]==s->word[jdx])
//            {
//                if(idx<jdx)
//                {
//                    state_pop(s1, jdx);
//                    state_pop(s1, idx);
//                }
//                else
//                {
//                    state_pop(s1, idx);
//                    state_pop(s1, jdx);
//                }
//                break;
//            }
//        }
//
//        // XXX This calculation is not exactly right
//        ratio = 2.*pairs/s->size/(s->size-1)/nletters;
//    }
    else
    if(i==3 && s->size>1)
    {
        // swap
        int tmp;
        idx = gsl_rng_uniform_int(rng, s1->size-1);
        tmp = s1->word[idx];
        s1->word[idx] = s1->word[idx+1];
        s1->word[idx+1] = tmp;
    }
    else
    if(i==4)
    {
        // flip a spin
        idx = gsl_rng_uniform_int(rng, SLABS); // slab
        jdx = gsl_rng_uniform_int(rng, sizeof(slab_t)*8); // bit
        //s1->u = (s1->u) ^ (((spins_t)1)<<idx);
        s1->u[idx] = (s1->u[idx]) ^ (((slab_t)1)<<jdx);
    }

    check_state(s1);

//    printf(" --> "); dump_state(s1); printf("\n");

    return ratio;
}

int
main(int argc, char *argv[])
{
    init_cache();
    init_ops();

    rng = gsl_rng_alloc (gsl_rng_gfsr4);   // fastest rng

    int chains, length, burnin, period;
    int seed;
    double beta;
    chains = 1;
    period = 1000; // unit
    burnin = 100; // warmup length in period units
    length = 1000;  // total length in period units
    verbose = 0;
    beta = 1.0;
    FILE *output = NULL;

    seed = (int)time(NULL); 
    gsl_rng_set(rng, seed);

    int idx;
    for(idx=1; idx<argc; idx++)
    {
//        if(idx==1) { chains = atoi(argv[idx]); }
//        if(idx==2) { length = atoi(argv[idx]); }
//        if(idx==3) { burnin = atoi(argv[idx]); }
//        if(idx==4) { period = atoi(argv[idx]); }
//        if(idx==5) { beta = atof(argv[idx]); }
//        if(idx==6) { seed = atoi(argv[idx]); gsl_rng_set(rng, seed); }

        char *arg = argv[idx];
        char *value;
        for(value = arg; *value; value++)
            if(*value=='=')
                break;
        if(*value == 0)
            continue;
        *value = 0;
        value++;
        
        if(strcmp("chains", arg)==0)      { chains = atoi(value); }
        else if(strcmp("length", arg)==0) { length = atoi(value); }
        else if(strcmp("burnin", arg)==0) { burnin = atoi(value); }
        else if(strcmp("period", arg)==0) { period = atoi(value); }
        else if(strcmp("verbose", arg)==0) { verbose = atoi(value); }
        else if(strcmp("beta", arg)==0)   { beta = atof(value); }
        else if(strcmp("seed", arg)==0)   
        { seed = atoi(value); gsl_rng_set(rng, seed); }
        else if(strcmp("output", arg)==0)   
        {
            printf("writing to file '%s'\n", value);
            output = fopen(value, "w+");
            if(output==NULL)
                printf("failed to open file '%s'\n", value);
        }
        else
        {
            printf("unrecognized arg '%s'\n", arg);
            return 1;
        }
    }

    printf("chains=%d ", chains); // number of chains
    printf("length=%d ", length); // length of each chain
    printf("burnin=%d ", burnin); // before first sample
    printf("period=%d ", period); // sample period
    printf("beta=%f ", beta);
    printf("seed=%d\n", seed);
    assert(length > burnin);

    printf("n = %d\n", n);
    printf("offset = %d\n", offset);

    burnin *= period;
    length *= period;
//    assert(burnin % period == 0);
//    assert(length % period == 0);

    assert(sizeof(spins_t) >= (n/8));

    int trial, chain;
    int accept = 0;
    int samples = 0;
    double total = 0.0;
    double log_beta = log(beta);

    for(chain=0; chain<chains; chain++)
    {
        struct state s, s1;
        double weight, weight1, x, ratio;
        int nz, nz1; // non-zero
    
        init_state_rand(&s, beta);
//        init_state(&s);
//        weight = eval_state(&s, beta);
        nz = log_eval_state(&s, beta, &weight);

        for(trial=0; trial<length; trial++)
        {
            ratio = move_state(&s, &s1);
            nz1 = log_eval_state(&s1, log_beta, &weight1);

            x = gsl_rng_uniform(rng);
            if(!nz || (nz1 && x <= (ratio * exp(weight1 - weight))))
            {
                if(!eq_state(&s, &s1))
                    accept += 1;
                memcpy(&s, &s1, sizeof(struct state));
                weight = weight1;
            }

            if(trial < burnin)
                continue;

            if((trial+1) % period == 0)
            {
                total += s.size;
                samples += 1;
                if(verbose>0) { printf("."); fflush(stdout); }
                if(output)
                    fprintf(output, "%d\n", s.size);
            }
        }
        if(verbose>0) { printf("/"); fflush(stdout); }
    }
    if(verbose>0) printf("\n");

    printf("accept: %f, ", (double)accept / ((length+1)*chains));
    printf("samples: %d\n", samples);
    printf("<H>: %f\n", total / samples / beta);
    printf("<H> - offset: %f\n", total / samples / beta - offset);

    gsl_rng_free(rng);
    if(output)
        fclose(output);

    return 0;
}



