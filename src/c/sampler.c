#include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_randist.h>
  #include <gsl/gsl_errno.h>
  #include <gsl/gsl_math.h>
  #include <gsl/gsl_eigen.h>


  #include <mpi.h>

  #include "iniparser.h"

  #include "const.h"
  #include "varb.h"
  #include "mymath.h"
  #include "myerr.h"
  #include "matrix.h"
  #include "init.h"
  #include "power.h"
  #include "cospara.h"

  #include "optstruct.h"
  #include "random.h"
  #include "sampler.h"

//--------------------------------------------------
//--------------------------------------------------
// MCMC sampler 


//  #ifdef MPI

//--------------------------------------------------------
  // estimate covariance matrix after each UPDATE time
    #define UPDATE  5
  // begin to estimate covariance matrix at this point
    #define BEGIN_COV_UPDATE  100
  // the maximum number of points that covariance matrix will consider
    #define MAX_COV_UPDATE  2000

    #define MIN_CONVERG  2000

    #define R_CRITERION  3.0

  // if the "time" of a point stayed is greater than MAX_STAY_TIME, decrease the step size
    #define MAX_STAY_TIME  5
  // if there are MAX_ACCEPT_TIME points ( whose "time" is equal to one) 
  //  are accepted in succession,increase the step size
    #define MAX_ACCEPT_TIME  3


  // the weight which increase or decrease the step size
    #define INCREASE_STEP  1.05
    #define DECREASE_STEP  0.952381



    #define ANN_PERIOD  2000
    #define ANN_AMPLITUDE 2

     


//--------------------------------------------------------
    #define GIVEMETASK  10 //MPI flag for slave()'s wish to receive a task from master()
    #define TAKERESULT  20 //MPI flag for slave()'s request to master() to take result
    #define TAKETASK  30   //MPI flag for master()'s request to slave to take task.

    #define TASK_FAIL 77

//-----------------------------------------------------------
//master assign tasks to slaves.

  int master(chainpar *chain,  double (*chi2)(void *inform, void *arg, int *x), void *info )  {

      int i, j, k, fail=TRUE;

      if(chain->debug>DEBUGV)
        printf("enter master()\n");
//------------------------------------------------------------------------------
// initialize chains.

     double ***chains, **currentpoint, **trialpoint, ***covmat,
            ***eigenvec, **point_stepsize , **eigenvalue, 
            *scale, *gauss_x;
     int *point_num, *update_num, *loop_num;


     chains = dmat3(chain->nchain, chain->npoint, chain->size);

     currentpoint = dmat(chain->nchain, chain->size);
     trialpoint   = dmat(chain->nchain, chain->size);

     covmat   = dmat3(chain->nchain, chain->npara, chain->npara);
     eigenvec = dmat3(chain->nchain, chain->npara, chain->npara);
     
     point_stepsize =  dmat(chain->nchain, chain->npara);
     eigenvalue     =  dmat(chain->nchain, chain->npara);

     scale      = dvec(chain->nchain);  

     point_num  = ivec(chain->nchain);
     loop_num   = ivec(chain->nchain);
     update_num = ivec(chain->nchain);

//----------------------------------------------------------------
     dvecvalue(scale , chain->nchain, 1.);
     ivecvalue(update_num, chain->nchain, UPDATE );

//----------------------------------------------------------------
  // initialize radom variables.

     randpara random;
     random_init(&random); 

     gauss_x    = dvec(chain->npara);

//??     randvec_gauss(&random, gauss_x, chain->npara);

//---------------------------------------------------------------
//   pack up all the stuffs. 
     chainvar chain_var;

     chain_var.chains         =  chains; 
     chain_var.currentpoint   =  currentpoint;
     chain_var.trialpoint     =  trialpoint;

     chain_var.covmat         =  covmat;   
     chain_var.eigenvec       =  eigenvec; 

     chain_var.point_stepsize =  point_stepsize;
     chain_var.eigenvalue     =  eigenvalue;    

     chain_var.scale          =  scale;       
     chain_var.point_num      =  point_num ;
     chain_var.loop_num      =   loop_num ;
     chain_var.update_num     =  update_num;

     chain_var.gauss_x        =  gauss_x ;
     chain_var.random         =  &random;

     chain_var.min_chain  = 1e7;
     chain_var.converg = 0;


     if(chain->debug>DEBUGV)
       printf("master(): Initialized chain_var.\n");

//----------------------------------------------------------------
// initialize the starting points.
     if(chain->startrandom != TRUE) {
       for(i=0; i<chain->nchain; i++)
         for(j=0; j<chain->npara; j++)
           currentpoint[i][j] = trialpoint[i][j] = chain->startp[j];
       }

     else  {

       for(j=0; j<chain->npara; j++) {
         currentpoint[0][j]= random_uniform(&random, chain->lowbound[j], chain->upbound[j]) ;
         trialpoint[0][j] = currentpoint[0][j];
         for(i=1; i<chain->nchain; i++)
           trialpoint[i][j] = currentpoint[i][j] = currentpoint[0][j];
//           currentpoint[i][j] = trialpoint[i][j] = chain->startp[j] + 
//                             gauss_x[j]*scale[i]*chain->inisigma[j];
         }
       }

     if(chain->debug>2*DEBUGV) {
       printf("master(): Initialized starting points, nchain=%d.\n", chain->nchain);
       for(j=0; j<chain->npara; j++)
       printf("currentpoint[%d]=%lg\n", j, currentpoint[0][j]);
       }

      currentpoint[0][chain->size-2] = chi2(info, currentpoint[0], &fail); //!!should check
      for(i=0; i<chain->nchain; i++)  {
        currentpoint[i][chain->size-2] = currentpoint[0][chain->size-2];
        printf("master(%d): %lg\n", i, currentpoint[i][chain->size-2]);
        currentpoint[i][chain->size-1] = 1.;
        }

//----------------------------------------------------------------
// initialize running variables.
      int run_bool, break_bool, free_bool, *chain_run_bool, *accept_num, 
          *point_converg_num, chain_n, p_num, pass_length;

      run_bool = free_bool = 1; // TRUE;


     double * result, *next, u; 
     int chi2fail;
   // result & next: MPI transforming variables.

     pass_length = chain->size+1;
     result = dvec(pass_length); 
     next   = dvec(pass_length); 
 
     accept_num = ivec( chain->nchain );
     ivecvalue(accept_num, chain->nchain, 1);    

  //points' number of each chain when reach the convergence
     point_converg_num = ivec(chain->nchain); 
  //go on running for each chain 
     chain_run_bool = ivec(chain->nchain);
     ivecvalue(chain_run_bool, chain->nchain, 1); 

  //-------------------------------------------------
     if(chain->debug>2*DEBUGV) {
       printf("master(): Initialized running variables. \n");

       for(i=0; i<chain->nchain; i++){
         printf("%s\n", chain->fname_chain[i]);
         printf("%s\n", chain->fname_best[i]);
         }
       }
  //-------------------------------------------------

     FILE *fp_chain[chain->nchain], *fp_best[chain->nchain];
//     FILE **fp_chain, **fp_best;
//     fp_chain = (FILE **)anyvec(chain->nchain, sizeof(FILE *));
//     fp_best  = (FILE **)anyvec(chain->nchain, sizeof(FILE *));

     for(i=0; i<chain->nchain; i++)  {
       fp_chain[i] = fopen(chain->fname_chain[i], "wb");
       fp_best[i]  = fopen(chain->fname_best[i],  "wb");
       }

     if(chain->debug>2*DEBUGV) 
       printf("master(): Initialized chain files. \n"); fflush(stdout);
 // End of running variables initialization.
//------------------------------------------

//----------------------------------------------------------
// initialize MPI variables
     MPI_Status status;

//-----------------------------------------------------------
// MCMC sampling begins.

     do   {

       MPI_Recv(result, pass_length, MPI_DOUBLE, MPI_ANY_SOURCE, 
                               MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    //----------------------------------------------------
    // assign task to slaves
       if (status.MPI_TAG == GIVEMETASK) {

          chain_n = status.MPI_SOURCE - 1;

    // make sure the points selected within the boundaries
          for( ; ; )  {
            try_newpoint(chain, &chain_var, chain_n); 

            break_bool = 1;

            //if successful step, break
//            if(break_bool)
            for(i=0; i<chain->npara; i++) {
              if (trialpoint[chain_n][i]<chain->lowbound[i] || 
                               trialpoint[chain_n][i]>chain->upbound[i] )
                break_bool = 0;
              }
  
            if(break_bool) break;
  
            // so we have crossed the boundary. Hence we have to increase the weight of current point
            currentpoint[chain_n][chain->size-1]++;
            }

          for(i=0; i<chain->size; i++) 
            next[i] = trialpoint[chain_n][i];

          MPI_Send(next, pass_length, MPI_DOUBLE, status.MPI_SOURCE, 
                                                TAKETASK, MPI_COMM_WORLD);

          }

    //----------------------------------------------------
    // analyse results from slaves
        if (status.MPI_TAG == TAKERESULT) {

          chain_n = status.MPI_SOURCE - 1;

        // if haven't converge yet.
          if(!chain_var.converg) {

            // first, find out how long the shortest chain is
            for(i=0; i<chain->nchain; i++) {
              if( point_num[i]< chain_var.min_chain )
                chain_var.min_chain = point_num[i];
              }
            
            if(chain_var.min_chain > MIN_CONVERG && random_uniform(&random, 0.0,1.0) > 0.9) 
              convergence_test(chain, &chain_var);
 
            chain_var.min_chain = 1e7;
  
            point_converg_num[chain_n] = point_num[chain_n];

            }

   //-----------------------------------------------------------------
          chain_n = status.MPI_SOURCE - 1;

          loop_num[chain_n]++;


          for(i=0; i <chain->size; i++) 
            trialpoint[chain_n][i] = result[i];

          chi2fail = (int) result[pass_length-1];
  //-------------------------------------------------------------------
      // a new point will be accepted randomly.

          temperature_sched( chain, &chain_var, chain_n);

          u = random_uniform(&random, 0., 1.);

          if( (chi2fail>TASK_FAIL+1||chi2fail<TASK_FAIL-1)&& exp(trialpoint[chain_n][chain->size-2]
                  - currentpoint[chain_n][chain->size-2]) > pow(u, chain->T[chain_n]) )   {
//          if( chi2fail!=TASK_FAIL &&exp(trialpoint[chain_n][chain->size-2]
//                  - currentpoint[chain_n][chain->size-2]) > u)   {

            if(point_num[chain_n] > BEGIN_COV_UPDATE) 
              update_num[chain_n]++;

       //-----------------------------------------------------
         // output the chains
         
            for(i=0; i<chain->size; i++)   
              fprintf(fp_chain[chain_n], "%lf  ", currentpoint[chain_n][i]);
            fprintf(fp_chain[chain_n], "\n");
            fflush(fp_chain[chain_n]);

         // and debug information.
            if( chain->debug > 3*DEBUGV && chain_n==1) {
              printf ("%d  ", point_num[chain_n]+1);
              fflush(stdout);
              for (i=0; i<chain->size; i++)   {
                printf ("%lf  ", currentpoint[chain_n][i]);
                fflush(stdout);
                }
              printf ("\n");
              fflush(stdout);
              }

    //-----------------------------------------------------------------------------
      // if not converge, store in 'chains' for convergence calculation.
            if( !(chain_var.converg) ) {
              p_num = point_num[chain_n];

              for(j=0; j<chain->size; j++) 
                chains[chain_n][p_num][j] = currentpoint[chain_n][j]; 
              }

    //-----------------------------------------------------------------------------
            for(k=0; k <chain->size-1; k++) 
              currentpoint[chain_n][k] = trialpoint[chain_n][k];
        
            currentpoint[chain_n][chain->size-1] = 1.0; // set up "time" 
 
            accept_num[chain_n] ++;
            point_num[chain_n] ++;
            
       // output debug information.
            if( chain->debug > 3*DEBUGV && random_uniform(&random, 0.0,1.0) > 0.998) {
              for(i=0; i <chain->nchain; i++)  
                printf("%d  ",point_num[i]);
              printf("\n");
              }


            }

        //--------------------------------------------------              
        // reject, still old point
          else {
            currentpoint[chain_n][chain->size-1] ++; // old point time ++
            accept_num[chain_n] = 1;  // reset the number of continued accepted points
            }


    //------------------------------------------------------------              
    //------------------------------------------------------------              
        // change the step size.
          if(!chain_var.converg) {
            if(currentpoint[chain_n][chain->size-1] > MAX_STAY_TIME && scale[chain_n] > 0.09) 
              scale[chain_n] *= DECREASE_STEP; 
    	    if(accept_num[chain_n] > MAX_ACCEPT_TIME && scale[chain_n] < 10.0) 
              scale[chain_n] *= INCREASE_STEP;     
            }
          else if(chain->npara>= 7) 
            scale[chain_n] = 2.4/sqrt((double)(chain->npara *( chain->npara -2))); 

          }
  //-----------------------------------------------------------------------------
    // End of "TAKERESULT"


//----------------------------------------------------------------------------------------  
//---------------------------------------------------------------------------------------- 
   // still running? free all chains?(since donot need them for covariance estimating)
        if(chain_var.converg && free_bool)   {
          printf("converge & freeall chains. ");
          if(chain->debug>DEBUGV)
            printf("converge & freeall chains. ");
            for(i=0; i<chain->nchain; i++)
              printf("%d  ", point_num[i]);
            printf("\n");

          freemat3((void ***)chains, chain->nchain, chain->npoint, chain->size);
          free_bool = 0;
          } 

        if(point_num[chain_n]-point_converg_num[chain_n] >= chain->max_point) 
          chain_run_bool[chain_n] = 0;
        run_bool = 0;

        for(i=0; i<chain->nchain; i++) 
          run_bool += chain_run_bool[i]; // if all chains > max_points, break      


        } while(run_bool);
//--------------------------------------------------------------------------------------
  // End of loop.

      printf("master(): Out of the loop\n");

//------------------------------------------
// free all

      fcloseall();

      MPI_Abort(MPI_COMM_WORLD, 110);
 

//      free();

      return TRUE;
    }



//----------------------------------------------------------------------
// slaves calculate chi2/FoM

//  int slave( chainpar *chain, double (*chi2)(void *inform, void *arg),  
//                        double *(constraint)(void *inform, void *arg),  void *info)  {

  int slave( chainpar *chain, double (*chi2)(void *inform, void *arg, int *x),  void *info)  {

//----------------------------------------------
     if(chain->debug>20*DEBUGV)
       printf("enter slave()\n"); fflush(stdout);
//----------------------------------------------
     int i, pass_length, fail=TRUE, rank, loopn=0, slave_out = TRUE;
     double *task, *paraset, outstep=100;

     pass_length = chain->size+1;
     task    = dvec(pass_length);
     paraset = dvec(pass_length);
                    
     MPI_Status status;

     MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  
     for (;;)   { 

       fail=TRUE;
       paraset[chain->size] = 0;

     // slave request a task 
       MPI_Send(task, pass_length, MPI_DOUBLE, 0, GIVEMETASK, MPI_COMM_WORLD);
     // master's answer
       MPI_Recv(task, pass_length, MPI_DOUBLE, 0, TAKETASK, MPI_COMM_WORLD, &status);

       for(i=0; i<chain->size; i++) 
         paraset[i] = task[i]; 


     // calculate chi^2/FoM or any other functions.
//       constraint(info, paraset);
       paraset[chain->size-2] = chi2(info, paraset, &fail);

       loopn++;
       if( slave_out==TRUE && fmod(loopn, outstep)< 1./outstep/2. )
         printf("proc(%d): called Fisher matrix %dth times\n", 
                          rank, (int)( outstep*(int)(loopn/outstep) ) );

       if(fail == TRUE)
         paraset[chain->size] = TASK_FAIL;
       else
         paraset[chain->size] = 0;

       if(chain->debug>20*DEBUGV)   {
         printf("slave_%d:  ", rank );

         for(i=0; i<pass_length; i++)
           printf("%lg  ", paraset[i]);

         printf("\n");  fflush(stdout);
         }

     // sending back the result to master()
       MPI_Send(paraset, pass_length, MPI_DOUBLE, 0, TAKERESULT, MPI_COMM_WORLD);
       }       

//-----------------------------------------
      free(task);         free(paraset);

      return TRUE;
    }



//------------------------------------------------------------------------------
// other used routines.

//-----------------------------------------------------------------------------
// estimate and manage the covariance matrix
  void estimate_cov(chainpar *chain, chainvar *chain_var, int chain_n)   {

      int i, j, k, n;
  
      double ** average, *totime; 

      average = dmat(chain->nchain, chain->npara);
      totime  = dvec(chain->nchain); ;
  
     //---------------------------------------------------------------------
     //set # of components in the covariance matrix we will estimate
      int *cov_component_num = ivec(chain->nchain); 	
  	   
      if(chain_var->point_num[chain_n]<2*BEGIN_COV_UPDATE)    // early on
        cov_component_num[chain_n] = chain_var->point_num[chain_n]/2; 
      else  //all but the first few BEGIN_COV_UPDATE points   
        cov_component_num[chain_n] = chain_var->point_num[chain_n]-BEGIN_COV_UPDATE; 
      

      if(cov_component_num[chain_n]>MAX_COV_UPDATE) 
        cov_component_num[chain_n]=MAX_COV_UPDATE;  // at most MAX_COV_UPDATE

    //-----------------------------------------------------------------------
    // initialization.
//      dmatvalue(average, chain->nchain, chain->npara, 0);
      dmat3value(chain_var->covmat, chain->nchain, chain->npara, chain->npara, 0);


    //-----------------------------------------------------------------------
    // get covariance matrix
      int num = chain_var->point_num[chain_n]-cov_component_num[chain_n];

    //for( j=0; j<chain_var->point_num[chain_n]-cov_component_num[chain_n]; j++ ) 
    //   num++; // ignore the first part of the chain, i.e. the ones to be excluded 
  

    //-----------------------------------------------------------
      for(n=num; n<chain_var->point_num[chain_n]; n++)  {
        for(k=0; k<chain->npara; k++) 
          average[chain_n][k] += chain_var->chains[chain_n][n][k]*
                               chain_var->chains[chain_n][n][chain->size-1];
      //------------------------------------------------
           // total weight of all points
        totime[chain_n] += chain_var->chains[chain_n][n][chain->size-1]; 
        }
    
      for(k=0; k<chain->npara; k++) 
        average[chain_n][k] /= totime[chain_n];  //normalize
      	     
      for(n=num; n<chain_var->point_num[chain_n]; n++)  // run over all points
        for(k=0; k<chain->npara; k++) 
          for(j=0; j<chain->npara; j++)
            chain_var->covmat[chain_n][k][j] += (chain_var->chains[chain_n][n][k] - 
                           average[chain_n][k])*(chain_var->chains[chain_n][n][j] - 
                           average[chain_n][j])*chain_var->chains[chain_n][n][chain->size-1];
  	      
  //---------------------------------------------------------------
  // get covariance matrix  
      for( k=0; k<chain->npara; k++)  
        for( j=0; j<chain->npara; j++) 
          chain_var->covmat[chain_n][k][j] /= (totime[chain_n]-1.0); 
      

  //-------------------------------------------------------------------
  //diagonalization 
  // reformat "[][]" into "[]" for creating matrix
  	     
      double *tempArray = dvec(chain->npara*chain->npara);

      k = 0;
      for(j=0; j<chain->npara; j++)  // copy into tempArray
        for(i=0; i<chain->npara; i++) {           
          tempArray[k] = chain_var->covmat[chain_n][i][j];
          k++;
          }

  //----------------------------------------------------------------
    // get eigenvalues and eigenvectors of the new covariance matrix
      gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(chain->npara);

    // inverse the covariance Matrix
      gsl_matrix_view m = gsl_matrix_view_array(tempArray, chain->npara, chain->npara);

     //(but the covariance Matrix is diagonal)
      gsl_vector *eval = gsl_vector_alloc(chain->npara);
      gsl_matrix *evec = gsl_matrix_alloc(chain->npara, chain->npara);
      gsl_eigen_symmv(&m.matrix, eval, evec, w); 

    //obtain eigenvalues and eigenvectors and save in eval, evec
    //now evec is transformation matrix
      gsl_eigen_symmv_free(w);

    //sort in ascending order, (not strictly neccessary)
      gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);  
                                                                 

    // get eigenvalues 
      for(i=0; i<chain->npara; i++) {   
        chain_var->eigenvalue[chain_n][i] = gsl_vector_get(eval, i);

      //get eigenvectors 
        for(j=0; j<chain->npara; j++) 
          chain_var->eigenvec[chain_n][i][j]=gsl_matrix_get(evec,i,j); 
        }
        
    //------------------------------------
    // scale eigenvalues      
      for(i=0; i<chain->npara; i++) 
        chain_var->eigenvalue[chain_n][i] *= squar(chain_var->scale[chain_n]); 
       

    //------------------------------------
    // free memory
      gsl_vector_free(eval);   gsl_matrix_free(evec);	
        
      freemat((void **)average, chain->nchain, chain->npara);
      free(totime);   free(cov_component_num);	free(tempArray);


      chain_var->update_num[chain_n] = 0;    

      return;
    }



//---------------------------------------------------------------------------
         
  void estimate_point_step(chainpar *chain, chainvar *chain_var, int chain_n)  {        

    int i,j;

    double **y = dmat(chain->nchain, chain->npara);

  //-------------------------------------------------------------------
  // generate y random numbers with variances given in estimate_cov()
    for( i=0; i<chain->npara; i++)  
      y[chain_n][i] = chain_var->gauss_x[i] * sqrt(chain_var->eigenvalue[chain_n][i]);  

    for( i=0; i<chain->npara; i++)  {

      chain_var->point_stepsize[chain_n][i] = 0;

    //get step size to the next point
      for( j=0; j<chain->npara; j++) { 
        chain_var->point_stepsize[chain_n][i] += chain_var->eigenvec[chain_n][i][j]*y[chain_n][j];
       // printf("eigenvec=%lg\n", chain_var->eigenvec[chain_n][i][j]);
        }
       // printf("eigenvale=%lg, gauss=%lg\n", chain_var->eigenvalue[chain_n][i], chain_var->gauss_x[i]);
      }

    //----------------------------------
    freemat((void **)y, chain->nchain, chain->npara);

    return;
    } 



//---------------------------------------------------------------------------

  void try_newpoint(chainpar *chain, chainvar *chain_var, int chain_n) {

     int i;

     randvec_gauss(chain_var->random, chain_var->gauss_x, chain->npara);
  
     if(chain_var->point_num[chain_n] >= BEGIN_COV_UPDATE)  {

       if(!chain_var->converg && chain_var->update_num[chain_n] >= UPDATE) 
         estimate_cov(chain, chain_var, chain_n);
  
       estimate_point_step(chain, chain_var, chain_n);
  
       for(i=0; i<chain->npara; i++ )  {
         chain_var->trialpoint[chain_n][i] = chain_var->currentpoint[chain_n][i] + 
                                 chain_var->point_stepsize[chain_n][i];

         //printf("stepsize=%lg, eigenvec=%lg, y=%lg\n", 
         //  chain_var->point_stepsize[chain_n][i], chain_var->eigenvec[chain_n][i][i], y[chain_n][i]);
         }

       }
  
     else  {
       for(i=0; i<chain->npara; i++ ) 
         chain_var->trialpoint[chain_n][i]=chain_var->currentpoint[chain_n][i] + 
                       chain_var->gauss_x[i]*chain_var->scale[chain_n]*
                        chain->inisigma[i]; 
       }

      return;
    }



//----------------------------------------------------------------------------------

  void convergence_test( chainpar *chain, chainvar *chain_var )  {

     //printf("\n\n!!!!Doing the Convergence check.!!!!\n\n");

     if(chain->checkconverg==FALSE) {   // don't chech the convergence, always convergent.
       chain_var->converg= 1;
       }
   // or adopted simple convergence condition.
     else if(chain->checkconverg==TRUE && chain->simpleconverg==TRUE) {

         chain_var->converg= 1;

       }
   // or ...
     else {
       int i, j, k, m, n;
  
       int cut_at = (chain_var->min_chain)/2; 
      
       double N = 0.;
          
       // calculate the mean of each parameter
       double ** y = dmat(chain->nchain, chain->npara);
       
//--  ---------------------------------------------------
       double * dist_y = dvec(chain->npara);
       double * B =  dvec(chain->npara);    // distribution mean	
       double * W =  dvec(chain->npara);    // variance within chains
       double * R =  dvec(chain->npara);    // monitoring parameter
       
       double M = (double)chain->nchain; 

//--  -------------------------------------------------------------

       for (i=0; i<chain->nchain; i++) { // for all chains

         for (n=cut_at; n<chain_var->min_chain; n++) {  // traverse through the chains
           for(k=0; k<chain->npara; k++)  // for all parameters
             y[i][k] += chain_var->chains[i][n][k]*chain_var->chains[i][n][chain->size-1];
           
           N += chain_var->chains[i][n][chain->size-1];
           }
         
         for (k=0; k<chain->npara; k++) {
           y[i][k] /= N;
           dist_y[k] += y[i][k]/M;  // mean is sum/chain->nchain
           }

         }
            
//---------------------------------------------------------------------------
       // variance between chains
       for (k=0; k <chain->npara; k++) 
         for (i=0; i<chain->nchain; i++) 
            B[k] += squar(y[i][k] - dist_y[k])/(M-1.); 
          
       // get W, the variance within chain_var->chains
       // run over all points
       for (i=0; i<chain->nchain; i++)   // run over all chain_var->chains
         for (n=cut_at; n<chain_var->min_chain; n++)  // run over all points
           for (k=0; k<chain->npara; k++) 
             W[k] += squar(chain_var->chains[i][n][k]-y[i][k])*chain_var->chains[i][n][chain->size-1];   
      
       //printf("%lf  %lf  %lf\n", N*(M-1.0)*B[1], W[1], N*(M-1.0)*B[1]/W[1]);
 
       // normalize W and get R
       for (k=0; k<chain->npara; k++) {
         W[k] /= M*(N-1.0);
         R[k]  = (N-1.0)/N * W[k]  + B[k];
         R[k] /= W[k];
         }
       printf("%d  ",chain_var->min_chain);

       for(k = 0; k < chain->npara; k++)  {
         printf("%lf  ", R[k]); 
         fflush(stdout);  
         }

       printf("\n");  

       int *convergence_para = ivec(chain->npara); //= (int *)malloc(chain->npara*sizeof(int));
            
       for(m=0; m<chain->npara; m++) {
         if(R[m] < R_CRITERION) convergence_para[m] = 1;
         chain_var->converg += convergence_para[m];
         }

       if(chain_var->converg < chain->npara-2) 
          chain_var->converg= 0;

//----------------------------------------

       freemat((void *)y, chain->nchain, chain->npara);
       free(dist_y);   free(B);   free(W);   free(R);  free(convergence_para);

     }

     return;
   }




  void temperature_sched( chainpar* chain, chainvar *chain_var, int chain_n)  {

      int period=ANN_PERIOD;

      chain->T[chain_n] = ANN_AMPLITUDE*(cos( 2.*PI* 
                  (double)(chain_var->loop_num[chain_n]/(double)period ) ) + 1.);

      //printf("T[chain_%d][%d]=%lg\n", chain_n, chain_var->loop_num[chain_n], chain->T[chain_n] );

      return;
    }





//----------------------------------------------------------------------------------------
  

