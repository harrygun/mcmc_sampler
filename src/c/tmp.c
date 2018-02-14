//-----------------------------------//
//           Do some test.           //
//-----------------------------------//


  #include<time.h>
  #include<stdlib.h>
  #include <math.h>


  #include "matrix.h"
  


  #define RANDLEN_TIME  100


//-----------------------------------//
//           Do some test.           //
//-----------------------------------//

  int main( )  {

      int i, j, k, length;
      double **x, **y;
//      printf("%d\n\n", RAND_MAX);



//      for(i=0; i<10; i++)  {
      for(i=10; i>0; i--)  {

        length = RANDLEN_TIME*(i+1);

        x = dmat(length, length);

        for(j=0; j<length; j++)
          for(k=0; k<length; k++)
            x[j][k] = j*2+k + i*1e6;

        printf("%d  ", length);
//        printf("%lf\n", x[1010][1010]);

        printf("%lf\n", x[3][9]);

        freemat((void **)x, length, length);

        printf("%dth allocation complete\n", i);
        } 


    }
