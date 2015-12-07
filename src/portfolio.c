#include <getopt.h>
#include <stdio.h>            /* C input/output                       */
#include <stdlib.h>           /* C standard library                   */
#include <time.h>
#include <math.h>
#include <glpk.h>             /* GNU GLPK linear/mixed integer solver */
#include <mpi.h>
#include <omp.h>
#include "mt19937p.h"
#include "rnglib.h"
#include "ranlib.h"

/*
    R is the random contraint data in row major memory layout
    ridx is an N array of integers 
    soln is an array of length (n+1) soln[n] is t (objective value)
    active_constr is an N-length 0-1 array
*/
void solve_lp(int N, int n, double* R, int* ridx, double* soln, int* active_constr)
{
    double tol = 1.0e-14;
    int size = (N+1)*(n+1) + 1; // We add one because GLPK indexes arrays
                                // starting at 1 instead of 0.

    glp_prob *lp;
    int* ia = malloc(size * sizeof(int));
    int* ja = malloc(size * sizeof(int));
    double* ar = malloc(size * sizeof(double));

    int i, j;
  
    lp = glp_create_prob();
    glp_set_prob_name(lp, "portfolio");
    glp_set_obj_dir(lp, GLP_MAX);

    glp_add_rows(lp, N+1);

    // Sampled constraints are ">= 0"
    for (i = 1; i <= N; i++) {
        glp_set_row_bnds(lp, i, GLP_LO, 0.0, 0.0);
    }
  
    // Sum = 1 constraint
    glp_set_row_name(lp, N+1, "sum");
    glp_set_row_bnds(lp, N+1, GLP_FX, 1.0, 1.0);
  
    glp_add_cols(lp, n+1);
  
    // Nonnegative variables y
    for (i = 1; i <= n; i++) {
        glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, i, 0.0);
    }
  
    // Free variable t
    glp_set_col_name(lp, n+1, "t");
    glp_set_col_bnds(lp, n+1, GLP_FR, 0.0, 0.0);
    glp_set_obj_coef(lp, n+1, 1.0);
    
    // for (i = 0; i < N*(n-1); i++) {
    //     printf("%d: %g\n", i, R[i]);
    // }

    int idx = 1;
    // Sampled constraints
    for (i = 1; i <= N; i++) {
        // Uncertain assets
        for (j = 1; j < n; j++) {
            ia[idx] = i;
            ja[idx] = j;
            ar[idx] = R[ ridx[(i-1)] * (n-1) + (j-1) ];
            idx += 1;
        }

        // Fixed return asset
        ia[idx] = i;
        ja[idx] = n;
        ar[idx] = 1.05;
        idx += 1;

        // t
        ia[idx] = i;
        ja[idx] = n+1;
        ar[idx] = -1.0;
        idx += 1;
    }

    // Sum = 1 constraint
    for (i = 1; i <= n; i++) {
        ia[idx] = N+1;
        ja[idx] = i;
        ar[idx] = 1.0;
        idx += 1;
    }
    // t
    ia[idx] = N+1;
    ja[idx] = n+1;
    ar[idx] = 0.0;
    idx += 1;

    // for (i = 1; i < size; i++) {
    //     printf("%d %d %g\n", ia[i], ja[i], ar[i]);
    // }

    glp_load_matrix(lp, size-1, ia, ja, ar);

    // glp_scale_prob(lp, GLP_SF_AUTO);

    glp_smcp param;
    glp_init_smcp(&param);
    param.meth = GLP_PRIMAL;
    //glp_std_basis(lp);
    glp_simplex(lp, &param);
  
    double z = glp_get_obj_val(lp);
    // printf("z = %g\n", z);
    if (soln) {
        for (i = 0; i < n; i++) {
            double y = glp_get_col_prim(lp, i+1);
            soln[i] = y;
            // printf("y%d = %g\n", i, y);
        }
        double t = glp_get_col_prim(lp, n+1);
        soln[n] = t;
        // printf("t = %g\n", glp_get_col_prim(lp, n+1));
    }
  

    for (i = 1; i <= N; i++) {
        double slack = glp_get_row_prim(lp, i);
        active_constr[i-1] = fabs(slack) < tol ? 1 : 0;
        // printf("constr%d %d\n", i, active_constr[i-1]);
    }


    glp_delete_prob(lp);
    // glp_free_env();
    free(ia);
    free(ja);
    free(ar);
}

double* generate_random_constraints(char* seed, int N, int n, double* R) 
{
    // struct mt19937p state;
    // sgenrand(seed, &state);

    int i, j;

    initialize();
    int seed1, seed2;
    phrtsd(seed, &seed1, &seed2);
    set_initial_seed(seed1, seed2);
    // for (int i = 0; i < 10; i++)
    //     printf("Random N(0,1): %g\n", gennor(0, 1));

    float meanv[n];
    for (i = 0; i < n; i++) {
        meanv[i] = 0.094278;
    }

    // float varv[n];
    // for (int i = 0; i < n; i++) {
    //     varv[i] = 1;
    // }
    // float* covm = setcov(n, varv, 0);
    float covm[n*n];
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j) {
                covm[i+j*n] = 0.002064;
            }
            else{
                covm[i+j*n] = 0;
            }
        }
    }

    float parm[n*(n+3)/2+1];
    setgmn(meanv, covm, n, parm);

    // Row major layout
    float* mvn;
    for (i = 0; i < N; i++) {
        float* mvn = genmn(parm);
        for (j = 0; j < n; j++) {
            R[i*n + j] = exp(mvn[j]);
            // R[i*n + j] = 1.0 + 0.1*genrand(&state) + 0.01; // [1.01, 1.11]
        }
        free(mvn);
    }

    return R;
}

void write_soln(const char* fname, int n, double* soln)
{
    FILE* fp = fopen(fname, "w+");
    if (fp == NULL) {
        fprintf(stderr, "Could not open output file: %s\n", fname);
        exit(-1);
    }
    int i;
    for (i = 0; i < n; ++i) {
        fprintf(fp, "%g\n", soln[i]);
    }
    fclose(fp);
}

int main(int argc, char** argv)
{
    int rank, num_p;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i;

    int n = 200;
    int N = 500;
    // Option processing
    extern char* optarg;
    const char* optstring = "n:N:";
    int c;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch (c) {
        case 'n': n = atoi(optarg); break;
        case 'N': N = atoi(optarg); break;
        }
    }

    glp_term_out(GLP_OFF);
    
    double* R;
    int* active_constr;
    if (rank == 0) {
        time_t curr_time;
        time(&curr_time);
        char* seed = ctime(&curr_time);
        R = malloc(N*num_p*(n-1)*sizeof(double));
        generate_random_constraints(seed, N*num_p, n-1, R);
        active_constr = (int*) calloc(N*num_p, sizeof(int));
    }

    double t0;
    if (rank == 0) {
        t0 = omp_get_wtime();
    }

    double* R_local = malloc(N*(n-1)*sizeof(double));
    int* active_constr_local = (int*) calloc(N, sizeof(int));

    MPI_Scatter(R, N*(n-1), MPI_DOUBLE, R_local, N*(n-1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int ridx[N];
    for (i = 0; i < N; i++) {
        ridx[i] = i;
    }
    solve_lp(N, n, R_local, ridx, NULL, active_constr_local);

    MPI_Gather(active_constr_local, N, MPI_INT, active_constr, N, MPI_INT, 0, MPI_COMM_WORLD);

    // Re-solve with just the active constraints...
    double soln[n+1];
    if (rank == 0) {

        int num_active_constr = 0;
        for (i = 0; i < N*num_p; i++) {
            num_active_constr += active_constr[i];
        }

        int ridx[num_active_constr];
        int count = 0;
        int i;
        for (i = 0; i < N*num_p; i++) {
            if (active_constr[i] == 1) {
                ridx[count] = i;
                count += 1;
            }
        }
        
        solve_lp(num_active_constr, n, R, ridx, soln, active_constr);
        // write_soln("soln.txt", n+1, soln);
    }

    if (rank == 0) {
        double t1 = omp_get_wtime();
        printf("%d %g %g", num_p, t1-t0, soln[n]);
        for (i = 0; i < n; i++) {
            printf(" %g", soln[i]);
        }
        printf("\n");
    }

    if (rank == 0) {
        free(R);
        free(active_constr);
    }

    free(R_local);
    free(active_constr_local);
    MPI_Finalize();
    return 0;
}
