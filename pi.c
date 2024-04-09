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
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
int verb = 0;
int limite = 0;
int target = 0;
int  NBLINEAR;
int optaff = 0;

vector*  linear;

int  dmin = 1;
int  dmax = 3;

basis_t base;

int count = 0;
int64_t total = 0;
int cpt = 0;

typedef struct l {
    boole pi;
    int   J;
    struct l *next;
} enrpi, *liste;

liste all = NULL;
liste rep[256] = {NULL};
int   counter[256] = {0};
liste grp = NULL;
int nbclass = 0;

int countpi = 0;

int q;

code V;

typedef struct quad {
    shortvec q[4];
} quad;

quad *W;

uint64_t size;
agl *table = NULL;

int affinity( mapping f )
{ int res = 0;
	shortvec x, y, z;
  for( x = 0; x < ffsize; x++ )
  for( y = x+1; y < ffsize; y++ )
  for( z = y+1; z < ffsize; z++ ) 
	  if ( z < (x^y^z)  && 0 == (f[x]^f[y]^f[z]^f[x^y^z] ) )
		  res++;
  return res;
}


int isaffine( mapping t )
{ int s = 0;
	int k;
	for (k = 0; k < V.nbl && s == 0; k++)
	    s = t[W[k].q[0]] ^ t[W[k].q[1]] ^ t[W[k].q[2]] ^ t[W[k].q[3]];
	return (  k == V.nbl && s == 0  );
}

int equivalent(mapping f, mapping g)
{
    uchar t[16];
    uchar inv[ 16 ] ;

    int x;
    for (x = 0; x < 16; x++)
	inv[ f[x] ] = x;
    liste aux = grp;
    while (aux ) {
	for (x = 0; x < 16; x++)
	    t[x] = inv[ aux->pi[ g[x] ] ];
	if( isaffine( t ) ) return 1;
	aux = aux->next;
    }
    return 0;
}

code flagCode(void)
{
    code R = getcode(16, 16 );
    code T = getcode(16, 16 );
    int k = 0;
    shortvec x, y, z, t;
    for (t = 0; t < R.nbl; t++)
	R.fct[t] = getboole();
    for (t = 0; t < R.nbl; t++)
	T.fct[t] = getboole();
    R.nbl = 0;
    T.nbl = 0;
    for (x = 0; x < ffsize ; x++)
	for (y = x + 1; y < ffsize; y++)
	    for (z = y + 1; z < ffsize; z++)
		if ((x ^ y ^ z) > z) {
		    for (t = 0; t < ffsize; t++)
			R.fct[k][t] = 0;
		    for (t = 0; t < ffsize; t++)
			T.fct[k][t] = 0;
		    R.fct[k][x] = 1;
		    R.fct[k][y] = 1;
		    R.fct[k][z] = 1;
		    R.fct[k][x ^ y ^ z] = 1;
		    R.nbl = k + 1;
		    T.fct[k][x] = 1;
		    T.fct[k][y] = 1;
		    T.fct[k][z] = 1;
		    T.fct[k][x ^ y ^ z] = 1;
		    T.nbl = k + 1;
		    k = pivotage( R );
		    R.nbl = k;
		    T.nbl = k;
		}
    assert(T.nbl == 11);
    // freecode( R );
    pcode("W", T);
    return T;
}

int ispi( mapping g )
{ int tmp[ 16 ] = {0};
	shortvec x;
	for( x = 0; x < 16; x++ ) tmp[g[x]]++;
	for( x = 0; x < 16; x++ ) 
		if ( tmp[x] != 1 ) return 0;
	return 1;
}

void loadpi( char *fn )
{
    FILE *src = fopen( fn , "r" );
    if ( ! src ) {
	    perror( fn );
	    exit(1);
    }
    mapping g;
    int val;

    while ( ( g = loadmap( src ) ) ) {
	    liste aux;
	    aux  = malloc( sizeof(enrpi) );
	    aux->pi = g;
	    aux->next = all;
	    all = aux;
	    countpi++;
    }
  
}

basis_t prepare( int dim )
{
    initboole( dim );
    initagldim( dim );
    basis_t res =  monomialBasis( dmin  , dmax ,  4  );
    int r = orbitBasic( mkaglGroup( 4 ) , & res );
    printf("\norbit : %d\n", r);

  
    vector v;
 
    linear = calloc( res.size  , sizeof(vector ) );
    int k = 1;
    for( v = 0; v < res.size ; v++ ) {
	    boole f = vectortoboole( v, &res );
	    if (  valuation(f) == dmin &&  degree(f) <= dmin ) {
		    linear[k] = v;
		    k++;
	    }
	    free( f );
    }
    NBLINEAR = k;
    return res;
}

void * rootj = NULL;
int   countj = 0;

void * roots = NULL;
int   counts = 0;

int countp = 0;
void * rootp = NULL;

int nouveau;

int invpi( boole g )
{
int b;
vector sp[16];

for( b = 0; b < 16; b++ ){
	 boole f = component( b, g );
	 vector v = booletovector( f, &base );
	 sp[  b ] =   v;
	 free( f );
}

nouveau = 0;
#define MAXP 1024 
#define MAXS 24
    int tp[ MAXP  ] = { 0 };
    int ts[ MAXS ];
    int j;
    for (j = 0; j < NBLINEAR; j++) {
        int k;
        for (k = 0; k < MAXS; k++)
            ts[k] = 0;
	int i;
        for (i = 0; i < 16; i++)
            ts[  base.table[ linear[ j ] ^ sp[i] ] ]++;
        int last = counts;
        int val = findspltable(ts, MAXS , &roots, &counts);
        if (last == val && verb > 2)
            printf("counts=%d\n", counts);
        assert(val < MAXP);
        tp[val]++;

    }
    int last = countp;
    int val  = findspltable(tp , MAXP, &rootp, &countp);
    if ( val == last ) nouveau = 1;
    return val;
}


int test( mapping f )
{ liste aux;
  int v = affinity( f );
  aux = rep[ v ];
  int j  = invpi( f );
  while( aux ) {
	  if ( aux->J == j && equivalent(f, aux->pi) ) 
		  return 0;
	  aux = aux -> next;
  }
   aux  = malloc( sizeof(enrpi) );
            aux->pi = getboolecpy(f);
	    aux->J  = j;
            aux->next = rep[v];
            rep[v] = aux;
	    counter[v]++;
  nbclass++;
  printf("\nclass=%d\n", nbclass);
  pTT( stdout, f );
  return 1;
}


void doit(  uchar pi[]   )
{   liste aux;
	if ( test( pi ) ) {
			pTT( stdout, pi );
		}
}

void mkrandpi( void )
{ uchar f[16];
	int i, j, t;
	for( i = 0; i < 16; i++ )
		f[i] = i;
  while ( nbclass < limite ) {
		  doit( f );
	  i = random() % 16;
	  j = random() % 16;
	  t = f[i];
	  f[i] = f[j];
	  f[j] = t;

  }
}
void permutation(int p, int lib[], uchar pi[], int n )
{
    int i;
    if ( countp >  limite ) 
	    return;
    if (p == n) {
        doit( pi );
        return;
    }
    for (i = 0; i < n; i++)
        if (  lib[i] ) {
            lib[i] = 0;
            pi[ p ] =  i;
            permutation(p + 1, lib, pi, n);
            lib[i] = 1;
        }
}



int mkallpi( void  )
{
    uchar pi[16] = { 0 };
    int lib[16] = { 0 };
    int i;
    for (i = 0; i  < 11; i++) 
	lib[i] = 1;
    for( i = 11; i < 16; i++ )
	    pi[i] = i;
    permutation(0, lib, pi, 11 );

    printf("\ncountp=%d\n", countp);

    return 0;
}


int accept( int t[]  )
{ int x, y;
  for( x = 1; x < 16; x++ )
          t[x]^=t[0];
  t[0] = 0;
  for( x = 0; x < 16; x++ ){
        y = 0;
        int i;
        for( i = 0; i < 4; i++ )
                if ( x & ((1<<i) ) ) y ^= t[ 1 << i ];
        if ( y != t[x] ) return 0;
  }
  return 1;
}

       
 int oldaccept( int t[]  )
{ int x, y;
  for( x = 1; x < 16; x++ )
          t[x]^=t[0];
  t[0] = 0;
  for( x = 0; x < 16; x++ ){
        y = 0;
        int i;
        for( i = 0; i < 4; i++ )
                if ( x & ((1<<i) ) ) y ^= t[ 1 << i ];
        if ( y != t[x] ) return 0;
  }
  return 1;
}


int self( mapping f)
{
        uchar g[ffsize];
        int x;
        for( x = 0; x < ffsize; x++)
                g[ f[x] ] = x;
        return equivalent( f, g) ;
}

       
int stabilizer( mapping f )
{
    int t[16];
    int g[16];
    int i, k, x, s, cpt;

    ssinit();
    for (x = 0; x < 16; x++)
            g[ f[x] ] = x;

    cpt=0;
    for (i = 0; i < size; i++) {
        for (x = 0; x < 16; x++)
            t[x] = f[aglImage(g[x], table[i])];
        s = 0;
        for (k = 0; k < V.nbl && s == 0; k++)
            s = t[W[k].q[0]] ^ t[W[k].q[1]] ^ t[W[k].q[2]] ^ t[W[k].q[3]];
        if ( ( k == V.nbl)  &&  (s == 0 ) ) {
                SchreierSims( table[i] , 0 );
                cpt++;
        }
    }
    uint64_t res, tmp;
    res = ssOrder( );
    assert( res == cpt );
    tmp = aglcard( 4 );
    tmp = tmp * tmp;
    assert( tmp % res == 0 );
    tmp  /= res;
    total += tmp;
    printf("\n");
    pTT( stdout, f);
    printf("\n#affinity = %d\n",  affinity( f )  );
    printf("\n#duality  = %d\n",  self( f )  );
    paglGroup(stdout, ssgen() );
    printf("\nstabSize=%ld\n", res ); 
    ssfree();
    return res;
}
int main(int argc, char*argv[] )
{


    int opt, optstab = 0, opteq =0, optj=0;
    char * file = NULL;
       while ((opt = getopt(argc, argv, "el:vsf:a:j")) != -1) {
	switch (opt) {
	case 'a' :
	   optaff =1 ;
	   target = atoi( optarg );
	   break;
	case 'e' :
	   opteq  =1 ;
	   break;
	case 'l' :
	   limite  = atoi( optarg );
	   break;
	case 'j' :
	   optj  =1 ;
	   break;
	case 'f' :
	   file = strdup( optarg);
	   break;
	case 's' :
	   optstab  =1 ;
	   break;
	case 'v' : verb++;
		   break;
	case 'h':
	default:		/* '?' */
	    fprintf(stderr, "Usage: %s -h \n", argv[0]);
	    exit(EXIT_FAILURE);
	}
    }

    srandom( time(NULL)) ;
    base = prepare( 4 );


    initboole(4);
    initagldim(4);



    V = flagCode();
    W = calloc(11, sizeof(quad));
    int k = 0;
    for (k = 0; k < 11; k++) {
        int r = 0, x;
        for (x = 0; x < ffsize; x++)
            if (V.fct[k][x])
                W[k].q[r++] = x;
        assert( r == 4 );
        }
    
    
    table = devGroup(mkaglGroup(), &size);
    printf("\naglsize=%ld\n", size);
    int i;
    for( i = 0; i < size; i++ ) {
	    liste aux = malloc( sizeof(enrpi) );
	    aux->pi = getboole();
	    shortvec x;
	    for( x = 0; x < 16; x++ )
		    aux->pi[x] = aglImage( x, table[i] );
	    assert( isaffine( aux->pi ) == 1 );
	    //pTT( stdout, aux->pi );
	    aux->next = grp;
	    grp = aux;
    }

    if ( file   ) loadpi( file );
    if ( optj   ) mkallpi(  );
    if ( optaff ) mkrandpi(  );

    printf("\ncountpi=%d\n", countpi);

    if ( optstab ){
	liste aux = all;
    	while( aux ) {
    		stabilizer( aux->pi );
		aux = aux -> next;
    	}    
    	printf("\ntotal=%ld\n", total );
    	uint64_t f = 1;
	int i;
    	for( i=1; i <17; i++)
	    f*=i;
	printf("\nfacto=%ld\n", f);
    }

    if ( opteq ){
	liste aux = all;
	int nbc = 0;
    	while( aux ) {
	        test( aux->pi );
		aux = aux-> next;
    	}    
    }
  
 int v;
  for( v = 0; v < 256; v++ )
	  if ( counter[v] ) printf("%3d", v );
	    printf("\n");
  for( v = 0; v < 256; v++ )
	  if ( counter[v] ) printf("%3d", counter[v] );
    
    return 0;
}
