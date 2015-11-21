#include <stdio.h>            /* C input/output                       */
#include <stdlib.h>           /* C standard library                   */
#include <glpk.h>             /* GNU GLPK linear/mixed integer solver */

/*
    This solves the LP
        max      4 x_1  +  3 x_2
        s.t.     3 x_1  +    x_2 <= 9
                 3 x_1  +  2 x_2 <= 10
                   x_1  +    x_2 <= 4
                   x_1 >= 0  x_2 >= 0
*/

int main(void)
{
  /* declare variables */
  glp_prob *lp;
  int ia[1+1000], ja[1+1000];
  double ar[1+1000], z, x1, x2, y1, y2, y3;
  /* create problem */
  lp = glp_create_prob();
  glp_set_prob_name(lp, "short");
  glp_set_obj_dir(lp, GLP_MAX);
  /* fill problem */
  glp_add_rows(lp, 3);

  glp_set_row_name(lp, 1, "y1");
  glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 9.0);

  glp_set_row_name(lp, 2, "y2");
  glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 10.0);

  glp_set_row_name(lp, 3, "y3");
  glp_set_row_bnds(lp, 3, GLP_UP, 0.0, 4.0);

  glp_add_cols(lp, 2);

  glp_set_col_name(lp, 1, "x1");
  glp_set_col_bnds(lp, 1, GLP_LO, 0.0, 0.0);
  glp_set_obj_coef(lp, 1, 4.0);

  glp_set_col_name(lp, 2, "x2");
  glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 0.0);
  glp_set_obj_coef(lp, 2, 3.0);

  ia[1] = 1, ja[1] = 1, ar[1] = 3.0; /* a[1,1] = 3 */
  ia[2] = 1, ja[2] = 2, ar[2] = 1.0; /* a[1,2] = 1 */
  ia[3] = 2, ja[3] = 1, ar[3] = 3.0; /* a[2,1] = 3 */
  ia[4] = 2, ja[4] = 2, ar[4] = 2.0; /* a[2,2] = 2 */
  ia[5] = 3, ja[5] = 1, ar[5] = 1.0; /* a[3,1] = 1 */
  ia[6] = 3, ja[6] = 2, ar[6] = 1.0; /* a[3,2] = 1 */
  glp_load_matrix(lp, 6, ia, ja, ar);

  /* solve problem */
  glp_simplex(lp, NULL);

  /* recover and display results */
  z = glp_get_obj_val(lp);
  x1 = glp_get_col_prim(lp, 1);
  x2 = glp_get_col_prim(lp, 2);
  printf("z = %g; x1 = %g; x2 = %g\n", z, x1, x2);

  /* recover dual variables at optimality
     when y > 0, the constraint is active */
  y1 = glp_get_row_dual(lp, 1);
  y2 = glp_get_row_dual(lp, 2);
  y3 = glp_get_row_dual(lp, 3);
  printf("y1 = %g; y2 = %g; y3 = %g\n", y1, y2, y3);

  /* housekeeping */
  glp_delete_prob(lp);
  glp_free_env();
  return 0;
}