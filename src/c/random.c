  #include<time.h>
  #include<stdlib.h>
  #include <math.h>

  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_randist.h>


  #include "random.h"
  #include "const.h"


  #define RANDLEN_TIME  10000
  #define RANDLEN       50000


//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------


  int random_init(randpara *random)  {



      gsl_rng_env_setup();


      random->T_uniform  = gsl_rng_default;
      random->T_gaussian = gsl_rng_default;

      random->r_uniform  = gsl_rng_alloc(random->T_uniform);
      random->r_gaussian = gsl_rng_alloc(random->T_gaussian);

      gsl_rng_set(random->r_uniform,  time(NULL) );
      gsl_rng_set(random->r_gaussian, time(NULL) );

      return TRUE;
    }



  double random_uniform(randpara *random, double lower, double upper)  {

      double x = gsl_rng_uniform(random->r_uniform);
//      printf("random:%lg\n", x); fflush(stdout);

      return lower+(upper-lower)*x;
    }



  double random_gauss(randpara *random, double sigma)  {

      double x = gsl_ran_gaussian(random->r_gaussian, sigma);

      return x;
    }


  void randvec_gauss(randpara *random, double *array, int n) {

      int i;

      for(i=0; i<n; i++)
        array[i] = random_gauss(random, 1);

      return;
    }
//-----------------------------------------------------------------------------












//-----------------------------------------------------------------------------

//  #define _RNG_TEST_
  #ifdef _RNG_TEST_

  int main( )  {

      int i;
      double x, y;
//      printf("%d\n\n", RAND_MAX);

      randpara random;
      random_init(&random); 


      for(i=0; i<10; i++)  {
        x = random_uniform(&random, 0, 1);
        y = random_gauss(&random, 1);
        printf("%lg %lg\n", x, y);

        } 


    }

  #endif









//--------------------------------------------------------------------------



/*
  int random_seed_init(random_type *random)  {

      int i;
      double N = 1./(RAND_MAX + 1.);

      random->narrtime = RANDLEN_TIME;

      random->arrtime = dvec(RANDLEN_TIME);

      srand((double)time(NULL));

      for(i=0; i<RANDLEN_TIME; i++)
        random->array[i] = (double) (rand() * N );


      return TRUE;
    }


  int random_init(random_type *random)  {

      random->narray  = RANDLEN;
      random->array   =  dvec(RANDLEN);


      return TRUE;
    }

*/
