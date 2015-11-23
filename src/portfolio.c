#include <stdio.h>            /* C input/output                       */
#include <stdlib.h>           /* C standard library                   */
#include <time.h>
#include <glpk.h>             /* GNU GLPK linear/mixed integer solver */
#include "mt19937p.h"

int main(int argc, char** argv)
{
    const int n = 10;
    const int N = 10;
    const int size = (N+1)*(n+1);
  
    /* declare variables */
    glp_prob *lp;
    int ia[size], ja[size];
    double ar[size];
  
    /* create problem */
    lp = glp_create_prob();
    glp_set_prob_name(lp, "portfolio");
    glp_set_obj_dir(lp, GLP_MAX);
  
    /* fill problem */
    glp_add_rows(lp, N+1);

    for (int i = 1; i <= N; i++) {
        glp_set_row_bnds(lp, i, GLP_LO, 0.0, 0.0);
    }
  
    glp_set_row_name(lp, N+1, "sum");
    glp_set_row_bnds(lp, N+1, GLP_FX, 1.0, 1.0);
  
    glp_add_cols(lp, n+1);
  
    // Nonnegative variables y
    for (int i = 1; i <= n; i++) {
        glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, i, 0.0);
    }
  
    // Free variable t
    glp_set_col_name(lp, n+1, "t");
    glp_set_col_bnds(lp, n+1, GLP_FR, 0.0, 0.0);
    glp_set_obj_coef(lp, n+1, 1.0);

    struct mt19937p state;
    sgenrand(time(NULL), &state);
  
    int idx = 1;
    // Random constraints
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= n; j++) {
            ia[idx] = i;
            ja[idx] = j;
            double r = genrand(&state) + 1;
            // printf("r=%f\n", r);
            ar[idx] = r;
            idx += 1;
        }
        // t
        ia[idx] = i;
        ja[idx] = n+1;
        ar[idx] = -1.0;
        idx += 1;
    }


    // Sum = 1 constraint
    for (int i = 1; i <= n; i++) {
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

    glp_load_matrix(lp, size, ia, ja, ar);
  
    /* solve problem */
    glp_simplex(lp, NULL);
  
    /* recover and display results */
    printf("z = %g\n", glp_get_obj_val(lp));
    for (int i = 1; i <= n; i++) {
        printf("y%d = %g\n", i, glp_get_col_prim(lp, i));
    }
    printf("t = %g\n", glp_get_col_prim(lp, n+1));
  
    /* recover dual variables at optimality
       when y > 0, the constraint is active */
    // y1 = glp_get_row_dual(lp, 1);
    // y2 = glp_get_row_dual(lp, 2);
    // y3 = glp_get_row_dual(lp, 3);
    // printf("y1 = %g; y2 = %g; y3 = %g\n", y1, y2, y3);
  
    /* housekeeping */
    glp_delete_prob(lp);
    glp_free_env();
    return 0;
}