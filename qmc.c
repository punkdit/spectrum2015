
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "gsl_rng.h"


typedef uint64_t spins_t;
//typedef __uint128_t spins_t;

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

#define MAX_WORD (256)

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
    assert(s->u<(1<<n));
    assert(0<=s->size);
    assert(s->size<=MAX_WORD);
}


void
init_state(struct state *s)
{
    assert(s);
    s->u = 0;
    s->size = 0;
    memset(s->word, 0, sizeof(int)*MAX_WORD);
    check_state(s);
}

void
init_state_rand(struct state *s, double beta)
{
    assert(s);
    // ...
    check_state(s);
}


void
dump_state(struct state *s)
{
    int i;
    check_state(s);
    printf("%ld", s->u);
    printf("[");
    for(i=0; i<s->size; i++)
        printf("%d", s->word[i]);
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
factorial(int n)
{
    double r;
    int i;
    r = 1.0;
    for(i=1; i<=n; i++)
        r *= i;
    return r;
}

double 
eval_state(struct state *s, double beta)
{
    spins_t u, v;
    int i;
    double weight;

    check_state(s);
    u = s->u;
    weight = 1.0;
//    printf("eval_word\n");
    for(i=0; i<s->size; i++)
    {
        int idx = s->word[i];
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
    }
    assert(u==s->u);
    weight *= 1.0 * pow(beta, s->size) / factorial(s->size);
    return weight;
}


int 
eq_state(struct state *s1, struct state *s2)
{
    int idx;
    check_state(s1);
    check_state(s2);
    if(s1->u != s2->u)
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

    i = gsl_rng_uniform_int(rng, 6);
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
    if(i==2)
    {
        // add a pair
        idx = gsl_rng_uniform_int(rng, s->size+1);
        jdx = gsl_rng_uniform_int(rng, s->size+2);
        i = gsl_rng_uniform_int(rng, nletters);
        state_insert(s1, idx, i);
        state_insert(s1, jdx, i);

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

        ratio = (0.5/pairs) * (s->size+1) * (s->size+2) * nletters;
    }
    else
    if(i==3 && s->size>1)
    {
        // remove a pair
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

        while(1)
        {
            idx = gsl_rng_uniform_int(rng, s->size);
            jdx = gsl_rng_uniform_int(rng, s->size);
            if(idx != jdx && s->word[idx]==s->word[jdx])
            {
                if(idx<jdx)
                {
                    state_pop(s1, jdx);
                    state_pop(s1, idx);
                }
                else
                {
                    state_pop(s1, idx);
                    state_pop(s1, jdx);
                }
                break;
            }
        }

        ratio = 2.*pairs/s->size/(s->size-1)/nletters;
    }
    else
    if(i==4 && s->size>1)
    {
        // swap
        spins_t tmp;
        idx = gsl_rng_uniform_int(rng, s1->size-1);
        tmp = s1->word[idx];
        s1->word[idx] = s1->word[idx+1];
        s1->word[idx+1] = tmp;
    }
    else
    if(i==5)
    {
        // flip a spin
        idx = gsl_rng_uniform_int(rng, n);
        s1->u = (s1->u) ^ (1<<idx);
    }

    check_state(s1);

//    printf(" --> "); dump_state(s1); printf("\n");

    return ratio;
}


int
main(int argc, char *argv[])
{

    rng = gsl_rng_alloc (gsl_rng_gfsr4);   // fastest rng

    int seed, trials, counts;
    double beta;
    trials = 100;
    counts = 1;
    beta = 1.0;

    seed = (int)time(NULL); 
    gsl_rng_set(rng, seed);

    if(argc>=2) { trials = atoi(argv[1]); }
    if(argc>=3) { counts = atoi(argv[2]); }
    if(argc>=4) { beta = atof(argv[3]); }
    if(argc>=5) { seed = atoi(argv[4]); gsl_rng_set(rng, seed); }

    printf("trials = %d\n", trials);
    printf("counts = %d\n", counts);
    printf("beta = %f\n", beta);
    printf("seed = %d\n", seed);

    assert(n<=10);

    int trial, count;
    int accept = 0;
    double total = 0.0;

    for(count=0; count<counts; count++)
    {
        struct state s, s1;
        double weight, weight1, x, ratio;
    
        init_state(&s);
        weight = eval_state(&s, beta);

        for(trial=0; trial<trials; trial++)
        {
            ratio = move_state(&s, &s1);
            weight1 = eval_state(&s1, beta);

//            dump_state(&s);
//            printf(" --> ");
//            dump_state(&s1);
//            printf(" weight=%f, weight1=%f, ratio=%f", weight, weight1, ratio);
//            printf("\n");
    
            x = gsl_rng_uniform(rng);
            if(weight==0.0 || x <= (ratio * weight1 / weight))
            {
                if(!eq_state(&s, &s1))
                {
                    accept += 1;
//                    printf("ACCEPT\n");
                }
                memcpy(&s, &s1, sizeof(struct state));
                weight = weight1;
            }

        }

        total += s.size;
    }

    printf("accept: %f\n", (double)accept / (trials*counts));
    printf("<n>: %f\n", total / counts);

    gsl_rng_free(rng);

    return 0;
}



