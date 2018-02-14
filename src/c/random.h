
  #ifndef _H_MYRANDOM_
    #define _H_MYRANDOM_

    #include <gsl/gsl_rng.h>
    #include <gsl/gsl_randist.h>


      typedef struct  {

          int state, narray, narrtime;
          double *arrtime, *array;

        } randarray;

  
      typedef struct  {

          const gsl_rng_type *T_uniform, *T_gaussian;
          gsl_rng *r_uniform, *r_gaussian;

        } randpara;



      int random_init(randpara *random);

      double random_uniform(randpara *random, double lower, double upper)  ;

      double random_gauss(randpara *random, double sigma) ;


  #endif
