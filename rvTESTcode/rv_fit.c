#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define MAXLINE 80


/* get_line: read a char* line and return length */
int get_line(char *line, int max, FILE *filep)
{
  if (fgets(line, max, filep) == NULL) {
    return 0;  /* return zero if line is empty */
  }
  
  return strlen(line); 
}


/* get_file_length: return number of lines from char* line to EOF */
int get_file_length(char *line, FILE *filep)
{
  int flen = 0;

  while (get_line(line, MAXLINE, filep)) {
    if (isdigit(line[0])) {
      flen++;  /* increment flen iff line contains data */
    }
  }
  rewind(filep);  /* set filep back to SOF */

  return flen;
}


/* get_data: read r, v_lum data from file and write to r[], v[]*/
void get_data(char *line, double r[], double v[], int flen, FILE *filep)
{
  char *rtoken, *vtoken = NULL;  /* r, v strings from data file line */
  int i = 0;  /* iterator index for r[], v[] */

  /* fill r[], v[] with data points from data file */
  while (i < flen && get_line(line, MAXLINE, filep)) {
    rtoken = strtok(line, " \t\n");
    vtoken = strtok(NULL, " \t\n");
    /* ignore non-numeric lines and non-positive values */
    if (atof(rtoken) > 0.0 && atof(vtoken) > 0.0) {
      r[i] = atof(rtoken);  /* convert rtoken to double */ 
      v[i] = atof(vtoken);  /* convert vtoken to double */
      // printf("%f %f\n", r[i], v[i]);
      i++;
    }
  }

  return;
}


/* interpolate_vgal: run gsl spline to calculate and print the interpolated 
 * luminous velocity v of the MW at each associated radius r. */
void interpolate_vgal(double r_mw[], double v_mw[], double r_gal[], 
    int milkyway_flen, int galaxy_flen)
{
  double r, v = 0.0;  /* r, v values used by spline interpolator */
  int i = 0; /* iterator index used to calculate vi_gal from r_gal[] */

  /* initial gsl_interp_accel to acc pointer */
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  /* allocate computer memory to spline pointer */
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, milkyway_flen);

  /* init spline using galaxy r, v_lum values, and galaxy file length */
  gsl_spline_init(spline, r_mw, v_mw, milkyway_flen);

  /* iterate over r values from r_gal[] */
  for (r = r_gal[i]; i < galaxy_flen; r = r_gal[++i]) {
    /* update vi_gal value using new ri_mw value */
    v = gsl_spline_eval(spline, r, acc);
    printf("%f %f\n", r, v);
  }

  /* free memory allocated to spline and acc */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return;
}


/* rv_fit: run interpolator on (r, v_lum) data files from command line */
int main(int argc, char* argv[])
{
  /* return 1 if no filename provided to argv[] */
  if (argc <= 2) { 
    printf("Provide MW and galaxy data filenames as command line arguments.\n");

    return 1;
  }

  FILE *milkyway_fp = fopen(*++argv, "r");  /* pointer to MW data file*/
  FILE *galaxy_fp = fopen(*++argv, "r");  /* pointer to galaxy data file*/
  char *line = (char*) calloc(MAXLINE, sizeof(char));  /* data file line */

  /* return 1 if invalid filename provided to argv[] */
  if (milkyway_fp == NULL || galaxy_fp == NULL) {
    printf("Both '%s' and '%s' must be in the working directory.\n",
        *argv, *(argv - 1));

    return 1;
  }

  int milkyway_flen = get_file_length(line, milkyway_fp);  /* MW file length */ 
  int galaxy_flen = get_file_length(line, galaxy_fp);  /* galaxy file length */ 
  double r_mw[milkyway_flen], v_mw[milkyway_flen];  /* MW r, v arrays */ 
  double r_gal[galaxy_flen], v_gal[galaxy_flen];;  /* galaxy r, v arrays */
  
  /* fill r_mw, v_mw, r_gal, v_gal with values from corresponding data files */
  get_data(line, r_mw, v_mw, milkyway_flen, milkyway_fp);
  get_data(line, r_gal, v_gal, galaxy_flen, galaxy_fp);

  /* unlink data files and free memory allocated to file and line pointers */
  fclose(milkyway_fp);
  fclose(galaxy_fp);
  free(line); 

  interpolate_vgal(r_mw, v_mw, r_gal, milkyway_flen, galaxy_flen);
  
  return 0;
}
