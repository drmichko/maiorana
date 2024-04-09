
#include <stdio.h>
#include <stdlib.h>
#include "boolean.h"
#include "agltools.h"
#include "sstools.h"
#include "code.h"
#include "mapping.h"
#include "distrib.h"
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>

#define MAX 39916800

int flag = 0;
uint64_t total = 0;
int  *K, *J;
char *R;
mapping *PI;


int count = 0;
int cpt = 0;
typedef struct l {
    int pi[16];
    struct l *next;
} enrpi, *liste;

liste all = NULL;
liste orb = NULL;

int q;

int val[16];

code V;

typedef struct quad {
    shortvec q[4];
} quad;

quad *W;

uint64_t size;
agl *table = NULL;


void saveall(void)
{
    int i, m = 0;
    FILE *dst = fopen("signal.txt", "w");
    for (i = 0; i < count; i++)
	if (K[i] == i) {
	    pTT(dst, PI[i] );
	    m++;
	}
    fclose(dst);
    printf("\n%d maps save\n", m);
}

int edge = 0 , fusion = 0;
void gestionnaire1(int num)
{
    flag = 1 - flag;
    //printf("\nflag=%d   m=%d  p=%d fusion=%d\n", flag, edge, count - fusion, fusion );
}

void gestionnaire2(int num)
{
    //saveall();
}



int affinity;

int flags( mapping f )
{ int res = 0;
	shortvec x, y, z;
  for( x = 0; x < ffsize; x++ )
  for( y = x+1; y < ffsize; y++ )
  for( z = y+1; z < ffsize; z++ ) 
	  if ( z < (x^y^z)  && 0 == (f[x]^f[y]^f[z]^f[x^y^z] ) )
		  res++;
  return res;
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
    printf("\nres=%ld  members=%ld total=%ld\n", res, tmp , total );
    printf("\nres=%ld aff=%d", res, flags( f )  );
    ssfree();
    return res;
}


code flagCode(void)
{
    code R = getcode(16, 18);
    code T = getcode(16, 16);
    int k = 0;
    shortvec x, y, z, t;
    for (t = 0; t < R.nbl; t++)
	R.fct[t] = getboole();
    for (t = 0; t < R.nbl; t++)
	T.fct[t] = getboole();
    R.nbl = 0;
    T.nbl = 0;
    for (x = 0; x < ffsize; x++)
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
		    k = pivotage(R);
		    R.nbl = k;
		    T.nbl = k;
		}
    assert(T.nbl == 11);
    // freecode( R );
    // pcode(T);
    return T;
}

void checkgroup( void )
{ int i, j, k;
  for( i = 0; i < size; i++)
  for( j = i+1; j < size; j++){
	for( k = 0; k < 5 && table[i][k] ==  table[j][k]; k++ );;;
  assert( k < 5 );
  }
  exit(0);
}
int main(int argc, char*argv[] )
{

    signal(SIGUSR1, gestionnaire1);
    signal(SIGUSR2, gestionnaire2);
    printf("\npid =%d\n", getpid());

    initboole(4);
    initagldim( 4 );


    table = devGroup(mkaglGroup(), &size);
    // checkgroup( );

    printf("\naglsize=%ld\n", size);

    V = flagCode();
    W = calloc(11, sizeof(quad));
    int k = 0;
    for (k = 0; k < 11; k++) {
	int r = 0, x;
	for (x = 0; x < ffsize; x++)
	    if (V.fct[k][x])
		W[k].q[r++] = x;
    }


    boole g;
    int num = 0;
    FILE* src = fopen( argv[1] , "r" );
    while ( ( g = loadmap( src ) ) ){ 
		    pTT( stdout , g );
		    stabilizer( g );
		    num++;
    }
    printf("\n#map : %d\n", num );
    printf("\n#agl : %ld\n", size );
    printf("\ntotal=%ld\n", total );
    int i; 
    int64_t f = 1;
    for( i = 1; i <=16; i++ ) f *= i;

    printf("\nfactor=%ld\n", f );
    fclose( src );

    return 0;
}
