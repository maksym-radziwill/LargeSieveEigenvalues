/*
subroutine dsyev_ 	( 	character  	JOBZ,
		character  	UPLO,
		integer  	N,
		double precision, dimension( lda, * )  	A,
		integer  	LDA,
		double precision, dimension( * )  	W,
		double precision, dimension( * )  	WORK,
		integer  	LWORK,
		integer  	INFO 
	) 		

DSYEV_ computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices
Download DSYEV_ + dependencies [TGZ] [ZIP] [TXT]

Purpose:

     DSYEV_ computes all eigenvalues and, optionally, eigenvectors of a
     real symmetric matrix A.

Parameters:
    [in]	JOBZ	

              JOBZ is CHARACTER*1
              = 'N':  Compute eigenvalues only;
              = 'V':  Compute eigenvalues and eigenvectors.

    [in]	UPLO	

              UPLO is CHARACTER*1
              = 'U':  Upper triangle of A is stored;
              = 'L':  Lower triangle of A is stored.

    [in]	N	

              N is INTEGER
              The order of the matrix A.  N >= 0.

    [in,out]	A	

              A is DOUBLE PRECISION array, dimension (LDA, N)
              On entry, the symmetric matrix A.  If UPLO = 'U', the
              leading N-by-N upper triangular part of A contains the
              upper triangular part of the matrix A.  If UPLO = 'L',
              the leading N-by-N lower triangular part of A contains
              the lower triangular part of the matrix A.
              On exit, if JOBZ = 'V', then if INFO = 0, A contains the
              orthonormal eigenvectors of the matrix A.
              If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
              or the upper triangle (if UPLO='U') of A, including the
              diagonal, is destroyed.

    [in]	LDA	

              LDA is INTEGER
              The leading dimension of the array A.  LDA >= max(1,N).

    [out]	W	

              W is DOUBLE PRECISION array, dimension (N)
              If INFO = 0, the eigenvalues in ascending order.

    [out]	WORK	

              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

    [in]	LWORK	

              LWORK is INTEGER
              The length of the array WORK.  LWORK >= max(1,3*N-1).
              For optimal efficiency, LWORK >= (NB+2)*N,
              where NB is the blocksize for DSYTRD returned by ILAENV.

              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA.

    [out]	INFO	

              INFO is INTEGER
              = 0:  successful exit
              < 0:  if INFO = -i, the i-th argument had an illegal value
              > 0:  if INFO = i, the algorithm failed to converge; i
                    off-diagonal elements of an intermediate tridiagonal
                    form did not converge to zero.

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/times.h>
#include <sys/time.h>
#include <getopt.h>
#include "unistd.h"

#define INT_LIMIT 24 /* Always n <= 2^24 */

void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, 
           double *w, double *work, int *lwork, int *info);

void ddisna_(char *job, int *m, int *n, double * d, double *sep, int *info); 

double slamch_(char *c); 

int * sieve (int n) {
  char * table = (char *) calloc(n + 1, sizeof(char));
  if(table == NULL){
    fprintf(stderr, "Could not allocate %d bytes for prime table\n", n + 1);
    exit(-1); 
  }
  
  memset(table, (char) 1, n + 1); 

  for(int i = 2; i <= n; i++)
    if(table[i] == 1)
      for(int j = 2*i; j <= n; j += i)
	table[j] = 0; 

  int num_primes = 0;
  for(int i = 2; i <= n; i++)
    if(table[i] == 1)
      num_primes += 1; 
  
  int * primes = (int * ) calloc(num_primes + 1, sizeof(int));
  if(primes == NULL){
    fprintf(stderr, "Could not allocate %d bytes for prime table\n", (int) (num_primes + 1)*sizeof(int));
    exit(-1); 
  }
  
  int j = 0; 
  
  for(int i = 2; i <= n; i++)
    if(table[i] == 1){
      primes[j] = i;
      j = j + 1; 
    }

  free(table);

  return primes; 
}


int gen_mul(int n, int j, const int * table){
  if (n == 1) return 1;
  if (n <= 0) return 0;
  if (table[j] == 0) return 1; 
  if (table[j] > n) return 1;
  if (table[j] == n) return 0; // This line is special to the Liouville function and should be avoided
  return (gen_mul(n, j + 1,table) - gen_mul(n/table[j], j + 1,table)); 
}

struct moebius_sum {
  int * elem;
  int size; 
};

struct moebius_sum generate_moebius_sum(int n){
  struct moebius_sum mu;
  mu.size = n + 1;
  mu.elem = (int *) calloc((n + 1),sizeof(int));
  mu.elem[0] = 0; 
  
  int * table = (int *) sieve(n+1); 

  for(int i = 1; i <= n; i++)
    mu.elem[i] = gen_mul(i, 0, table); 
  
  free(table); 

  return mu;
}

struct factor {
  int p;
  int alpha; 
};

struct factors {
  struct factor * f;
  int size; 
};

int power(int k, int a){
  int l = 1; 
  for(int i = 0; i < a; i++){
    l = l*k;
  }
  return l;
}

void add_factor(struct factor * f, int * size, int p, int alpha){
  *size = *size + 1;
  f[*size-1].p = p;
  f[*size-1].alpha = alpha;
}

struct factor * factorize(int n, const int * prime_table){
  int m = n; 
  /* Our integers have never more than 16 prime factors, since n < 2^16 always */
  struct factor * f = (struct factor *) calloc(INT_LIMIT, sizeof(int)); 
  int size = 0; 

  /* replace by sqrt(m) afterwards */
  
  for(int i = 0; prime_table[i] <= ((int) m + 1) && prime_table[i] != 0; i++){
    if(m % prime_table[i] == 0){
      
      /* Compute alpha */
      int k = 0;
      int p = prime_table[i];
      int pk = p; 
      while(m % pk == 0){
	pk = pk*p;
	k = k + 1; 
      }
      
      add_factor(f, &size, prime_table[i], k); 
      m = m / power(prime_table[i], k);
    }
  }
    
  if(m != 1)
    add_factor(f, &size, m, 1); 

  return f; 
}

int divisor_sum(int n, struct factor * f, struct moebius_sum mu, int q){
  if(f[0].p == 0) 
    return n*mu.elem[q/n]; /* Sum of divisors */ /* Replace this by n * M(q/n) in application */
  
  int sum = 0; 
  
  for(int i = 0; i <= f[0].alpha; i++)
    sum += divisor_sum(n / (int) power(f[0].p, i), f + 1, mu, q);

  return sum;
}

void print_factors(struct factor * f){
  for(int i = 0; f[i].p != 1; i++)
    printf("(%d, %d)\n", f[i].p, f[i].alpha);  
}

int farey_size(int n, struct moebius_sum mu){
  int sum = 1;
  
  for(int d = 1; d <= n; d++){
    int q = n / d; 
    sum += (mu.elem[d]-mu.elem[d - 1])*q*q; 
  }
  return (sum / 2); 
}

int * generate_values(int n, int q){
  struct moebius_sum mu = generate_moebius_sum(q + 1);
  int * prime_table = (int *) sieve(n + 10); 
  int * values = (int *) calloc(n + 10, sizeof(int));
  values[0] = farey_size(q, mu);
  for(int i = 1; i < n; i++){
    struct factor * f = factorize(i, prime_table); 
    values[i] = divisor_sum(i, f, mu, q);
    free(f); 
  }
  free(prime_table);
  free(mu.elem); 
  return values; 
}

int matrix_entry(int n, int q, int i, int j, int * values){
  return values[abs(i - j)];
}

void line_break(){
  printf("\n---------------------------------------------------\n\n"); 
}

void print_matrix(int n, int q, int * values){
  printf("The large sieve matrix A^{\\star} A for q = %d, n = %d: \n\n",q,n); 
  for(int i = 0; i < n; i++){
    printf("Row %2d : [", i); 
    for(int j = 0; j < n; j++){
      printf("%3d, ", matrix_entry(n,q,i,j,values)); 
    }
    printf("]\n"); 
  }
  line_break(); 
}

void print_row(int * row, int n){
  printf("First row of Toeplitz matrix\n");
  printf("["); 
  for(int i = 0; i < n; i++)
    printf("%3d, ", row[i]);
  printf("]"); 
  line_break(); 
}

void usage (char * s){
  printf("Usage: %s q [-n n] [-v] [-r] [-m]\n"
	 "Computes *all* the eigenvalues of the matrix A^{\\star} A\n"
	 "with A the large sieve matrix of dimension |F_{q}| \\times n.\n"
	 "If only the first argument q is provided then n is set to be |F_{q}|\n"
         "Other parameters are optional: -n sets the values of n\n"
	 "                               -v prints also the eigenvectors (slower!)\n"
	 "                               -r prints the first row of the matrix\n"
	 "                               -m prints the matrix\n"
	 "                               -e prints normalized eigenvalues with no formatting\n", s);
  exit(0); 
}



double machine_precision(){
  char jobz = 'E';
  return slamch_(&jobz); 
}

double max(double a, double b){
  if(a > b) return a;
  return b; 
}

double * eigenvectors_precision(double * w, int n){ 
  char * job = "Eigenvectors";
  int  N = n;
  double * rcondz = calloc(N, sizeof(double));
  double * zerobd = calloc(N, sizeof(double));
  double anorm = max (abs(w[0]), abs(w[n-1])); 
  int info;
  double * W = w; 

  ddisna_(job, &N, &N, W, rcondz, &info);

  for(int i = 0; i < n; i++)
    zerobd[i] = machine_precision()*(anorm / rcondz[i]); 

  free(rcondz);

  return zerobd; 
}


void print_eigenvectors(double * matrix, double * w, int n){

  double * zerobd = eigenvectors_precision(w,n); 
  
  for(int i = 0; i < n; i++){
    printf("Eigenvector with eigenvalue %10f (normalized %10f)\n", w[i], w[i] / (double) n);

    if(zerobd[i] > 0.5){
      printf("Skipping due to error bound exceeding 0.5\n\n");
      continue;
    }else{
      printf("Error bound for acute angle: %10f\n", zerobd[i]); 
      for(int j = 0; j < n; j++){
	printf("[ %10f ]\n", matrix[j + i*n]); 
      }
    }
  }
  
  free(zerobd); 
  line_break(); 
}

double eigenvalues_precision(double * w, int n){
  double anorm  = max (abs(w[0]), abs(w[n-1]));
  return machine_precision()*anorm; 
}

void print_eigenvalues_precision(double * w, int n){
  printf("\nNotes: Eigenvalues computed to precision %10f\n"
	 "Normalized eigenvalues computed to precision %10f\n",
	 eigenvalues_precision(w,n), eigenvalues_precision(w,n) / n + machine_precision()); 
  line_break();
}
			      
void print_eigenvalues(double * w, int n, int q){
  
  printf("Eigenvalues for q = %d, n = %d\n(and their normalized by %d counterparts): \n\n" , q, n, n); 
  for(int i = 0; i < n; i++){
    printf("%6d: %10f (normalized %10f) \n", i + 1, w[i] , w[i] / (double) n); 
  }

  line_break(); 
}

void print_plain_eigenvalues(double * w, int n){
  for(int i = 0; i < n; i++){
    printf("%10f\n", w[i] / (double) n); 
  }
}

int main(int argc, char ** argv)
{
  char        jobz;
  char        uplo;
  int        n = 0;
  double        *a;
  int          lda;
  double        *w;
  double     *work;
  int        lwork;
  double     wwork;
  int         info;
  int            q;
  int printrow = 0; 
  int printmat = 0; 
  char           c;
  int    plain = 0; 
  int      dim = 0; 
  
  jobz = 'N';   /* Change this to V if also want eigenvectors -- these are obtained by printing the matrix */
  uplo = 'U';

  if(argc < 2) usage(argv[0]);
  q = atoi(argv[1]); 
  if(q <= 0) usage(argv[0]);

  if(q >= power(2,INT_LIMIT)){
    printf("Error: q taken to be larger than %d\n", power(2,INT_LIMIT));
    return -1; 
  }
    
  while ((c = getopt (argc, argv, "n:vrme")) != -1)
    switch (c){
      case 'n':	n = atoi(optarg); break;
      case 'v': jobz = 'V'; break;
      case 'r': printrow = 1; break;
      case 'm':	printmat = 1; break;
      case 'e': plain = 1; break;
      default: usage(argv[0]); 
    }  
  
  if(n == 0){
    struct moebius_sum mu = generate_moebius_sum(q);  
    n   = farey_size(q, mu);
    free(mu.elem);
    printf("Setting n = %d\n", n); 
  }

  dim = n; 
  
  lda = n;
  a = malloc(sizeof(double)*lda*n);
  w = malloc(sizeof(double)*n);

  int * values = generate_values(n, q);

  for(int i = 0; i < n; i ++){
    for(int j = 0; j < n; j++){
      a[j + i*n] = matrix_entry(n,q,i,j,values); 
    }
  }

  if(!plain) printf("Matrix initialized\n"); 

  line_break(); 
  
  if(printrow) print_row(values, n); 
  
  if(printmat) print_matrix(n, q, values); 

  free(values); 
  
  lwork = -1;
  // get workspace information:
  dsyev_(&jobz, &uplo, &n, a, &lda, w, &wwork, &lwork, &info);
  lwork = (int) wwork;
  work = malloc(sizeof(double)*lwork);
  // do the real work:
  dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

  if(plain) print_plain_eigenvalues(w, n);
  else      print_eigenvalues(w, n, q); 

  print_eigenvalues_precision(w,n); 
  
  if(jobz == 'V') print_eigenvectors(a, w, dim); 

  if(!plain){
    if(info == 0)
      printf("Computation successful\n");
    else
      printf("Error encountered during computation\n"); 
  }

  free(work);
  free(a);
  free(w); 
  
  return 0;
}
