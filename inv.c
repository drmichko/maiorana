#include <stdio.h>
#include <stdlib.h>
#include "boolean.h"
#include "agltools.h"
#include "sstools.h"
#include "code.h"
#include "mapping.h"
#include "distrib.h"
#include "basistools.h"
#include "orbitools.h"
#include "code.h"
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
int verb = 0;

#define NBALS 10795
#define NBIG 5622036

boole QuadRk2[NBALS];
boole QuadRk4[NBIG];

void *rootj = NULL;
void *rootjs = NULL;
int countj = 0;
int countjs= 0;


void prepare(void)
{
    basis_t base = monomialBasis(2, 2, 8);
    aglGroup GL = mkglGroup();
    aglVectorGroup ldg = aglVectorGroupAction(GL, &base);
    initBrowse(&base);

    boole fct = strtoboole("anf=ab");
    vector v = booletovector(fct, &base);
    browse(v, ldg);
    int count = 0;
    for (v = 0; v < base.size; v++)
	if (base.table[v]) {
	    QuadRk2[count] = vectortoboole(v, &base);
	    count++;
	}
    printf("\ncount=%d\n", count);

    free(base.table);
    aglVectorGroupFree(ldg);
}

void Prepare(void)
{
    basis_t base = monomialBasis(2, 2, 8);
    aglGroup GL = mkglGroup();
    aglVectorGroup ldg = aglVectorGroupAction(GL, &base);
    initBrowse(&base);

    boole fct = strtoboole("anf=ab+cd");
    vector v = booletovector(fct, &base);
    browse(v, ldg);
    int count = 0;
    for (v = 0; v < base.size; v++)
	if (base.table[v]) {
	    QuadRk4[count] = vectortoboole(v, &base);
	    count++;
	}
    printf("\ncount=%d\n", count);

    free(base.table);
    aglVectorGroupFree(ldg);
}
void *rootJ = NULL;
int countJ = 0;

void *rootJs = NULL;
int countJs = 0;
#define MAXP 32
#define MAXB 200000

int cptK[256] = { 0 };

int countK = 0;

int invK(boole f, int r)
{
    code cc = rmcode(0, r, ffdimen);
    code prd = multicode(cc, f);
    int k;
    for (k = 0; k < prd.nbl; k++)
	projboole(r + 2, r + 4, prd.fct[k]);
    int res = pivotage(prd);
    if (!cptK[res]) {
	countK++;
	cptK[res] = 1;
    }
    freecode(cc);
    freecode(prd);
    return res;
}

void *rootX = NULL;
void *rootXs = NULL;
int countX = 0, countXs = 0;

int invX(boole f, int r)
{
    code cc = rmcode(0, r, ffdimen);
    code prd = multicode(cc, f);
    int k;
    for (k = 0; k < prd.nbl; k++)
        projboole(r + 2, r + 4, prd.fct[k]);
    int res = pivotage(prd);

    boole fct = getboole();
    int ts[ffsize];
    int limite = 1 << res;
    int iter = 1;
    int* tp = calloc( limite, sizeof(int) );
    while (iter < limite ) {
            int i = __builtin_ctz(iter);
	    int x;
            for( x = 0; x < ffsize; x++ )
		    fct[x] ^= prd.fct[i][x]; 
            for( x = 0; x < ffsize; x++ )
		    ts[x] = fct[x] ? -1 : +1;
	    int val = findtable( ts, ffsize, &rootXs, &countXs, 1);
            tp[ iter ] = val;
            iter++;
        }

     //res = findtable(tp, limite, &rootX, &countX, 0);
     free( fct );
     free( tp );
    freecode(cc);
    freecode(prd);
    return res;
}


void *rootb = NULL;
int countb = 0;
void *rootbs = NULL;
int countbs = 0;

int invb(boole f)
{
    int q;
    int res, val;
    int ts[ffsize];
    uchar tmp[ffsize];
    int tp[NBALS] = { 0 };
    for (q = 0; q < NBALS; q++) {
	int x;
	for (x = 0; x < ffsize; x++)
	    tmp[x] = f[x] && QuadRk2[q][x];
	projboole(6, 6, tmp);
	for (x = 0; x < ffsize; x++)
	    ts[x] = tmp[x] ? -1 : 1;
	Fourier(ts, ffsize);
	for (x = 0; x < ffsize; x++)
	    ts[x] = abs(ts[x]);
	val = findtable(ts, ffsize, &rootbs, &countbs, 0);
	if (verb > 1 && nouvelle)
	    printf("\ncountJs : %d\n", countbs);
	tp[q] = val;

    }
    res = findtable(tp, NBALS, &rootb, &countb, 0);

    if (verb && nouvelle)
	printf("\ncountb : %d\n", countb);
    return res;
}

void *rootB = NULL;
int countB = 0;
void *rootBs = NULL;
int countBs = 0;

int invB4(boole f)
{
    int q;
    int res, val;
    int ts[ffsize];
    uchar tmp[ffsize];
    int tp[NBALS] = { 0 };
    for (q = 0; q < NBALS; q++) {
        int x;
        for (x = 0; x < ffsize; x++)
            tmp[x] = f[x] && QuadRk2[q][x];
        projboole(6, 6, tmp);
        for (x = 0; x < ffsize; x++)
            ts[x] = tmp[x] ? -1 : 1;
        Fourier(ts, ffsize);
        for (x = 0; x < ffsize; x++)
            ts[x] = abs(ts[x]);
        val = findtable(ts, ffsize, &rootBs, &countBs, 0);
        if (verb > 1 && nouvelle)
            printf("\ncountJs : %d\n", countBs);
        tp[q] = val;

    }
    res = findtable(tp, NBALS, &rootB, &countB, 0);

    if (verb && nouvelle)
        printf("\ncountb : %d\n", countB);
    return res;
}


int invj( boole f)
{
    int q;
    int res, val;
    int ts[ffsize];
    int tp[MAXP] = { 0 };
    for (q = 0; q < NBALS; q++) {
	int x;
	for (x = 0; x < ffsize; x++)
	    ts[x] = f[x] ^ QuadRk2[q][x] ? -1 : +1;
	Fourier(ts, ffsize);
	for (x = 0; x < ffsize; x++)
	    ts[x] = abs(ts[x]);
	val = findtable(ts, ffsize, &rootj, &countjs, 0);
	if (verb > 1 && nouvelle)
	    printf("\ncountjs : %d\n", countjs);
	assert(val < MAXP);
	tp[val]++;

    }
    res = findspltable(tp, MAXP, &rootj, &countj);

    if (verb && nouvelle)
	printf("\ncountj : %d\n", countj);
    return res;
}




int invJ(boole f)
{
    int q;
    int res, val;
    int ts[ffsize];
    int tp[MAXB] = { 0 };
    for (q = 0; q < NBIG; q++) {
	int x;
	for (x = 0; x < ffsize; x++)
	    ts[x] = f[x] ^ QuadRk4[q][x] ? -1 : +1;
	Fourier(ts, ffsize);
	for (x = 0; x < ffsize; x++)
	    ts[x] = abs(ts[x]);
	val = findtable(ts, ffsize, &rootJs, &countJs, 0);
	if (verb > 1 && nouvelle)
	    printf("\ncountJs : %d\n", countJs);
	assert(val < MAXB);
	tp[val]++;

    }
    res = findspltable(tp, MAXP, &rootJ, &countJ);

    if (verb && nouvelle)
	printf("\ncountJ : %d\n", countJ);
    return res;
}
int invjold(boole f)
{
    int q;
    int res = 0;
    uchar tmp[ffsize];
    for (q = 0; q < NBALS; q++) {
	int x;
	for (x = 0; x < ffsize; x++)
	    tmp[x] = f[x] ^ QuadRk2[q][x];
	if (isbent(tmp))
	    res++;
    }
    return res;
}

int countD = 0;
void *rootD = NULL;

int countDs = 0;
void *rootDs = NULL;

int invD(boole f)
{
    shortvec u, x;
    int tmp;
#define MAXD 64
    int t[ffsize];
    int ts[ffsize];
    boole der;

    der = getboole();
    for (u = 0; u < ffsize; u++) {
	for (x = 0; x < ffsize; x++)
	    ts[x] = f[x] ^ f[x ^ u] ^ f[u] ^ f[0] ? -1 : +1;;
	Fourier(ts, ffsize);
	for (x = 0; x < ffsize; x++)
	    ts[x] = abs(ts[x]);

	int val = findtable(ts, ffsize, &rootDs, &countDs, 0);
	t[u] = val;
    }
    free(der);
    tmp = findtable(t, ffsize, &rootD, &countD, 0);
    if (verb && nouvelle)
	printf("\ncountD : %d (%d)\n", countD, countDs);
    return tmp;
}


int countR = 0;
void *rootR = NULL;

int countRs = 0;
void *rootRs = NULL;

void *rootRt = NULL;
int countRt = 0;


int mydegree(boole f, int ffsize)
{
    int x, p, d;
    boole aux;
    d = -1;
    aux = getboolecpy(f);

    xform(aux, ffsize);

    for (x = 0; x < ffsize; x++)
	if (aux[x] == 1) {
	    p = weight(x);
	    if (p > d)
		d = p;
	}
    free(aux);
    return (d);
}

boole duale( boole f )
{ int tfr[ ffsize ] ;
	int x;
	for( x = 0 ; x < ffsize; x++ )
		tfr[x] = f[x] ? -1: +1;
	Fourier( tfr, ffsize );
  boole g = getboole();
  for( x = 0 ; x < ffsize; x++ )
	  g[x] = ( tfr[x] > 0 );
  return g;
}

int L( boole f, int q )
{ int tfr[ ffsize ] ;
	int res = 0;
	int x;
	for( x = 0 ; x < q; x++ )
		tfr[x] = f[x] ? -1: +1;
	Fourier( tfr, q );
  for( x = 0 ; x < ffsize; x++ )
	  res += f[x]*f[x]*f[x]*f[x];
  return res;
}


int invR(boole f)
{
    shortvec u, p, msk, x, y, w;
    int tmp;
    int *tp = calloc(ffsize, sizeof(int));
    int *ts = calloc(ffsize, sizeof(int));
    int t[2];

    boole F = getboole(), G = getboole();
    for (u = 1; u < ffsize; u++) {
	p = 1;
	while ((p & u) == 0)
	    p <<= 1;
	msk = p ^ (ffsize / 2);
	for (x = 0; x < ffsize / 2; x++) {
	    y = x;
	    if (y & p)
		y ^= msk;
	    w = weight(y & u) & 1;
	    if (w)
		y ^= p;
	    F[x] = f[y];
	    G[x] = f[y ^ p];
	}
	
	for( x = 0 ; x < ffsize/2  ; x++ )
	   ts[x] = F[ x ] ? -1 : +1;
	//Fourier( ts, ffsize/2 );
	t[0]  = findtable( ts, ffsize/2,  &rootRs, &countRs, 1 );

	for( x = 0 ; x < ffsize/2  ; x++ )
	   ts[x] = G[ x ] ? -1 : +1;
	 //Fourier( ts, ffsize/2 );
	 t[1]  = findtable( ts, ffsize/2,  &rootRs, &countRs, 1 );

	if (t[0] < t[1]) {
	    tmp = t[0];
	    t[0] = t[1];
	    t[1] = tmp;
	}
	tp[u] = findspltable(t, 2, &rootRt, &countRt);
	if (nouvelle)
	    printf("\ncountRt: %d -> %d  %d\n", countRt, t[0], t[1]);
    }
    free(F);
    free(G);

    tmp = findtable(tp, ffsize, &rootR, &countR, 1);
    if (verb && nouvelle)
	printf("\ncountR : %d (%d)\n", countR, countRs);
    free(ts);
    free(tp);
    return tmp;
}
int line[400] = { 0 };

void numline(char *s)
{
    printf("\nusing line:");
    while (*s) {
	if (isdigit(*s)) {
	    int v;
	    if (sscanf(s, "%d", &v)) {
		line[v] = 1;
		printf(" %d", v);
	    }
	    while (isdigit(*s))
		s++;
	    s--;
	}
	s++;
    }
}

int main(int argc, char *argv[])
{


    int opt, optdeg = 0, optj = 0, optJ = 0, optB = 0, optb = 0, \
    optK = 0, optD = 0, optR = 0, optn = 0, optX = 0, optQ =  0, optdual = 0;
    char *file = NULL;

void *root = NULL;
int   count = 0;

    while ((opt = getopt(argc, argv, "tdf:vhjJbBDRK:X:n:Q")) != -1) {
	switch (opt) {
	case 'd':
	    optdeg = 1;
	    break;
	case 'j':
	    optj = 1;
	    break;
	case 't':
	    optdual = 1;
	    break;
	case 'J':
	    optJ = 1;
	    break;
	case 'Q':
	    optQ = 1;
	    break;
	case 'B':
	    optB = 1;
	    break;
	case 'b':
	    optb = 1;
	    break;
	case 'D':
	    optD = 1;
	    break;
	case 'R':
	    optR = 1;
	    break;
	case 'n':
	    optn = 1;
	    numline(optarg);
	    break;
	case 'K':
	    optK = atoi(optarg);
	    break;
	case 'X':
	    optX = atoi(optarg);
	    break;
	case 'f':
	    file = strdup(optarg);
	    break;
	case 'v':
	    verb++;
	    break;
	case 'h':
	default:		/* '?' */
	    fprintf(stderr, "Usage: %s -h \n", argv[0]);
	    exit(EXIT_FAILURE);
	}
    }


    initboole(8);
    initagldim(8);
    if ( optJ || optj || optb ) prepare();
    if ( optQ || optJ || optB ) Prepare();
    int R[12];
    boole f, g;
    int total = 0;
    FILE *src = fopen(file, "r");
    int nbi;
    while ((f = loadBoole(src))) {
	total++;
	projboole(2, ffdimen, f );
	g = duale( f );
	projboole(2, ffdimen, g );
	if (! optn  || line[total] == 1) {
	    nbi = 0;
	    if (verb > 1)
		panf(stdout, f);
	    if (optdeg)
		R[nbi++] = degree(f);
	    if (optdeg &&  optdual )
		R[nbi++] = degree(g);
	    
	    if (optj)
		R[nbi++] = invj(f);
	    if (optj && optdual )
		R[nbi++] = invj(g);

	    if (optJ)
		R[nbi++] = invJ(f);
	    if (optJ && optdual )
		R[nbi++] = invJ(g);
	    
	    if (optb)
		R[nbi++] = invb(f);

	    if (optb && optdual )
		R[nbi++] = invb( g );

	    if (optK)
		R[nbi++] = invK(f, optK);
	    if (optK && optdual )
		R[nbi++] = invK(g, optK);

	    if (optX)
		R[nbi++] = invX(f, optX);
	    if (optK && optdual )
		R[nbi++] = invX(g, optX);

	    if (optD)
		R[nbi++] = invD(f);
	    if (optD && optdual )
		R[nbi++] = invD(g);
	    
	    if (optR)
		R[nbi++] = invR(f);
	    if (optR && optdual )
		R[nbi++] = invR(g);

	    int val = findspltable(R, nbi, &root, &count);
	    if ( nouvelle ){
		printf("\nnew countj=%d %d", count, total);
	    }
	    if ( verb ){
		printf("\nitem : %d invariant : %d ", total, val );
		panf( stdout, f );
	    }
	}
	    free(f);
	    free(g);
    }
    fclose(src);
    printf("\n#class number : %d  using %d invariants\n", count, nbi );
    if (optj)
	printf("\ncountj : %d ( %d ) \n", countj, countjs);
    if (optJ)
	printf("\ncountJ : %d ( %d ) \n", countJ, countJs);
    if (optB)
	printf("\ncountB : %d ( %d ) \n", countB, countBs);
    if (optD)
	printf("\ncountD : %d ( %d ) \n", countD, countDs);
    if (optK)
	printf("\ncountK: %d ( %d ) \n", countK, optK);
    if (optX)
	printf("\ncountX: %d ( %d ) \n", countX, optX);
    if (optR)
	printf("\ncountR: %d ( %d ) \n", countR, countRs);
    return 0;
}
