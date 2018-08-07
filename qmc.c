
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
    printf("%d", s->u);
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
eval_state(struct state *s)
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
    for(j=idx+1; j<s->size; j++)
        s->word[j] = s->word[j-1];
    s->word[idx] = i;
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
    return i;
}


void
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

    memcpy(s1, s, sizeof(struct state*));

    i = gsl_rng_uniform_int(rng, 6);
    ratio = 1.0;

    if(i==0)
    {
        // add I
        total = state_count(s, 0);
        idx = gsl_rng_uniform(rng, s->size+1);
        ratio = ((double)(s->size)) / (total+1);
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
                idx = gsl_rng_uniform(rng, s->size);
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
        idx = gsl_rng_uniform(rng, s->size);
        jdx = gsl_rng_uniform(rng, s->size+1);
        i = gsl_rng_uniform(rng, nletters);
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
            idx = gsl_rng_uniform(rng, s->size);
            jdx = gsl_rng_uniform(rng, s->size);
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
            }
        }

        ratio = 2.*pairs/s->size/(s->size-1)/nletters;
    }
    else
    if(i==4 && s->size>1)
    {
        // swap
        spins_t tmp;
        idx = gsl_rng_uniform(rng, s1->size-1);
        tmp = s1->word[idx];
        s1->word[idx] = s1->word[idx+1];
        s1->word[idx+1] = tmp;
    }
    else
    if(i==5)
    {
        // flip a spin
        idx = gsl_rng_uniform(rng, n);
        s1->u[idx] = 1-s1->u[idx];
    }

    check_state(s1);

    return ratio;
}


int
main(int argc, char *argv[])
{

    rng = gsl_rng_alloc (gsl_rng_gfsr4);   // fastest rng

    int seed, trials, counts;
    trials = 100;
    counts = 1;
    beta = 1.0;

    if(argc>=2)
    {
        trials = atoi(argv[1]);
    }
    if(argc>=3)
    {
        counts = atoi(argv[2]);
    }

    if(argc>=4)
    {
        seed = atoi(argv[3]);
        gsl_rng_set(rng, seed);
    }
    else
    {
        seed = (int)time(NULL);
        gsl_rng_set(rng, seed);
    }

    printf("trials = %d\n", trials);
    printf("counts = %d\n", counts);
    printf("seed = %d\n", seed);

    assert(n<=10);

    int trial, count;

    int accept = 0;

    for(count=0; count<counts; count++)
    {
        init_state(word, order);
    
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



