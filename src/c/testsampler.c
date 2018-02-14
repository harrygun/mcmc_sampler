  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>

  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
  #include <iniparser.h>

  #include "const.h"
  #include "varb.h"
  #include "mymath.h"
  #include "myerr.h"
  #include "matrix.h"
  #include "init.h"
  #include "power.h"
  #include "cospara.h"

  #include "optstruct.h"
  #include "getfom.h"
  #include "random.h"
  #include "sampler.h"

//-------------------------
  #define MPI
//-------------------------
  #ifdef MPI
    #include "mpi.h"
  #endif
//------------------------


  double testchi2( void *inform, void *arg, int *fail  )  {

      int i, npar=4;

      double *para = (double *) arg;
      double chi2, fid0[npar], sig0[npar];

      *fail= FALSE;

      //fid0[0]=4;      sig0[0]=1;       
      fid0[0]=1.5e4;  sig0[0]=1.e3;   
      fid0[1]=4000;   sig0[1]=500;    
      fid0[2]=0.5;    sig0[2]=0.1;    
      fid0[3]=1.0;    sig0[3]=0.1;    

      chi2 = 0;
      for(i=0; i<npar; i++)
        chi2 += squar( (para[i]-fid0[i])/sig0[i] );

      chi2 = -0.5*chi2;

      //printf("chi2(): chi2=%lg\n", chi2); fflush(stdout);

      return chi2;
    }



//--------------------------------------------------
//--------------------------------------------------
// driver 
//-------------------------------------------------


  int main( int argc, char *argv[] ) {

//reading parameter .ini file
      int debug= 10, i, j;
      char ini_name[50];


//      if(argc!=2)  //
//        myerr("Input parameters file is needed.", FALSE);
 
      sscanf( argv[1], "%s", ini_name);

      if(debug>=DEBUGV)
      printf("reading %s\n", ini_name);


//------------------------------------------------
// MPI initialization.
    #ifdef MPI
      int mpi_numtasks, mpi_rank, mpi_rc;
      mpi_rc = MPI_Init(&argc, &argv);

      if (mpi_rc != MPI_SUCCESS) {
           printf ("Error starting MPI program. Terminating.\n");
           MPI_Abort(MPI_COMM_WORLD, mpi_rc);
           }


      MPI_Comm_size(MPI_COMM_WORLD,&mpi_numtasks);
      MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

      printf ("Number of tasks= %d My rank= %d\n", mpi_numtasks, mpi_rank);
// MPI initialization ends.

   //--------------------------------------------
      MPI_Barrier(MPI_COMM_WORLD);

      if(mpi_rank==0) 
        printf("%d Sending parameter filename %s to other processes.\n", mpi_rank, ini_name);

      MPI_Bcast( ini_name, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

      if(mpi_rank!=0) 
        printf("%d Received parameter filename %s.\n", mpi_rank, ini_name);

    #endif




//----------------------------------------------------------------------------------------------

  //----------------------------------------------------
  // Begin to initialize survey parameters.
  //----------------------------------------------------

//----------------------------------------------------------------------------------------------
// Global varibles.
    // global cosmological varibles & ones for FoM calculation.
      Cospar  cp;
      Survpar  sp;

      Interpar pk;
      cp.power = &pk;

      Fipar fpara;
        fpara.cp = &cp;

      double ***psarr;
  //-----------------------------------------------
    // universal setting of survey.
      survey_setting  surv;
    // survey configurations.
      survey_config sconf;

    // info for calculating FoM.
      fom_info info;
    // universal chain variable
      chainpar chain;


    #ifdef MPI

//--------------------------------------------------------------------------------
  // Open ini parameter file.
      printf("Opening File %s.\n", ini_name);

      dictionary * dict;
      dict = iniparser_load(ini_name);
//--------------------------------------------------
  // ini file ends.

//----------------------------------
// initialize the cosmological parameters
      init_cospar(&cp, dict); 


// initialize survey parameters
    //------------------------

      surv.name   = iniparser_getstring( dict, "Survey:name", NULL);
      surv.prefix = iniparser_getstring( dict, "Survey:prefix", "test");
      surv.random = iniparser_getboolean( dict, "Survey:random", INITRUE);

      surv.PLANCK = iniparser_getboolean( dict, "Survey:PLANCK", INITRUE);
      surv.SDSS   = iniparser_getboolean( dict, "Survey:SDSS", INITRUE);

      surv.line = iniparser_getboolean( dict, "Survey:Eimission_Line_Galaxies", TRUE);

      surv.cmbfish_file = iniparser_getstring( dict, "Survey:CMB_FISHER", "fmat_cmb.dat");
      surv.dim_cmb = iniparser_getint( dict, "Survey:DIM_FISHER_CMB", 8);

      surv.time   = iniparser_getdouble( dict, "Survey:time", 1500);
      surv.exptime= iniparser_getdouble( dict, "Survey:exposure_time", 90);
      surv.overheadtime= iniparser_getdouble( dict, "Survey:overhead_time", 15);
      surv.minexptime= iniparser_getdouble( dict, "Survey:min_exposure_time", 15);
      surv.maxexptime= iniparser_getdouble( dict, "Survey:max_exposure_time", 120);
      surv.nptmax= iniparser_getint( dict, "Survey:max_pointing_number", 4);
      surv.nfib= iniparser_getdouble( dict, "Survey:number_of_fibers", 4000);
      surv.fov = iniparser_getdouble( dict, "Survey:field_of_view", 20);


      surv.maxarea= iniparser_getdouble( dict, "Survey:maxarea", 10000);
      surv.zmin   = iniparser_getdouble( dict, "Survey:zmin", 0.);
      surv.zmax   = iniparser_getdouble( dict, "Survey:zmin", 1.);
      surv.dzmin  = iniparser_getdouble( dict, "Survey:dzmin", 0.1);
      surv.dzmax  = iniparser_getdouble( dict, "Survey:dzmax", 0.1);

      surv.camblist = iniparser_getstring( dict, "Input:cambpower_list", "power/camblist.dat");
      surv.camblist_np = iniparser_getint( dict, "Input:camblist_powernumber", 13);

      if(strcmp(surv.prefix, "LAMOST" )==0 ) {
        surv.LAMOST = TRUE;
        if(debug>=DEBUGV)
          printf("Doing LAMOST optimization.\n");
        }

//----------------------------------------------------------------
    // end of initialization of survey parameters


//----------------------------------------------------------------
    // initialize chain parameters
   
      chain.npara     = iniparser_getint( dict, "Chain:parameter_number", 5);
      chain.npoint    = iniparser_getint( dict, "Chain:chain_length", 20000);
      chain.max_point = iniparser_getint( dict, "Chain:chain_max_length", 25000);

      chain.debug   = iniparser_getint( dict, "MCMC_Control:debug", DEBUGV);
      surv.debug = chain.debug;

      chain.nchain  = mpi_numtasks - 1;
      chain.size    = chain.npara + 2;
      chain.temptype= iniparser_getdouble( dict, "Chain:temperature_type", 1.);;

      if( debug>=DEBUGV ) 
        printf("No. of chains = %d\n", chain.nchain);
//----------------------------------------------------------------
    // chains & results files
      chain.fname_chain= charmat(chain.nchain, 100 );
      chain.fname_best = charmat(chain.nchain, 100 );

      for(i=0; i<chain.nchain; i++) {
        sprintf(chain.fname_chain[i], "%s_%s_chain_%d.dat", surv.prefix, surv.name, i+1 );
        sprintf(chain.fname_best[i] , "%s_%s_best_%d.dat" , surv.prefix, surv.name, i+1 );
        }
//----------------------------------------------------------------
    // get the starting points, boundary & initial sigma. 

      chain.T        = dvec(chain.nchain);
      dvecvalue(chain.T, chain.nchain, 1.);

      chain.startp   = dvec(chain.npara);
      chain.lowbound = dvec(chain.npara);
      chain.upbound  = dvec(chain.npara);
      chain.inisigma = dvec(chain.npara);
    //-------------------------------------------

      chain.startp[0] = iniparser_getdouble( dict, "Chain:start_param1", 0.5);
      chain.startp[1] = iniparser_getdouble( dict, "Chain:start_param2", 0.5);
      chain.startp[2] = iniparser_getdouble( dict, "Chain:start_param3", 0.5);
      chain.startp[3] = iniparser_getdouble( dict, "Chain:start_param4", 0.5);
      //chain.startp[4] = iniparser_getdouble( dict, "Chain:start_param5", 0.5);

      chain.lowbound[0] = iniparser_getdouble( dict, "Chain:lower_param1", 0.);
      chain.lowbound[1] = iniparser_getdouble( dict, "Chain:lower_param2", 0.);
      chain.lowbound[2] = iniparser_getdouble( dict, "Chain:lower_param3", 0.);
      chain.lowbound[3] = iniparser_getdouble( dict, "Chain:lower_param4", 0.);
      //chain.lowbound[4] = iniparser_getdouble( dict, "Chain:lower_param5", 0.);

      chain.upbound[0] = iniparser_getdouble( dict, "Chain:upper_param1", 1.0);
      chain.upbound[1] = iniparser_getdouble( dict, "Chain:upper_param2", 1.0);
      chain.upbound[2] = iniparser_getdouble( dict, "Chain:upper_param3", 1.0);
      chain.upbound[3] = iniparser_getdouble( dict, "Chain:upper_param4", 1.0);
      //chain.upbound[4] = iniparser_getdouble( dict, "Chain:upper_param5", 1.0);

      chain.inisigma[0] = iniparser_getdouble( dict, "Chain:sigma_param1", 1.);
      chain.inisigma[1] = iniparser_getdouble( dict, "Chain:sigma_param2", 1.);
      chain.inisigma[2] = iniparser_getdouble( dict, "Chain:sigma_param3", 1.);
      chain.inisigma[3] = iniparser_getdouble( dict, "Chain:sigma_param4", 1.);
      //chain.inisigma[4] = iniparser_getdouble( dict, "Chain:sigma_param5", 1.);

//----------------------------------------------------------------
    // chain control 
      chain.startrandom   = iniparser_getboolean( dict, "Chain:start_random", INITRUE);
      chain.checkconverg  = iniparser_getboolean( dict, "Chain:check_convergence", INITRUE);
      chain.simpleconverg = iniparser_getboolean( dict, "Chain:simple_convergence_model", INITRUE);

//----------------------------------------------------------------
    // end of chain parameters initialization.



//-------------------------------------------------------------------------------------
    // initialize info for FoM calculation.

      info.cp    = &cp;
      info.sp    = &sp;
      info.fpara = &fpara;
//      info.scp   = &scp;

      info.surv  = &surv;
      info.s     = &sconf;

      info.psarr= psarr;

//-----------------------------------
    // get CMB Fisher matrix.
      if( debug>=DEBUGV ) 
        printf("proc_%d: Preparing for CMB Fisher Matrix calculation.\n", mpi_rank);

      char term;
      FILE *fpinput= fopen(surv.cmbfish_file, "r");

      info.fmat_cmb = dmat(surv.dim_cmb, surv.dim_cmb);

      for(i=0; i<surv.dim_cmb; i++) {
        for(j=0; j<surv.dim_cmb-1; j++) 
          fscanf(fpinput, "%lg ", &info.fmat_cmb[i][j] );
        fscanf(fpinput, "%lg%c", &info.fmat_cmb[i][surv.dim_cmb-1], &term );
        }

      if( debug>=DEBUGV ) 
        printf("proc_%d: Read the CMB Fisher Matrix.\n", mpi_rank);

      fclose(fpinput);

//-----------------------------------
    // Initialize powerspectrum for Fisher calculation.

      fpinput = fopen(surv.camblist, "r");
      char **cambpname= charmat( surv.camblist_np, 100);
      char cambnametmp[100];
//      char cambpname[surv.camblist_np][100];


      for(i=0; i<surv.camblist_np; i++) {
//        fscanf(fpinput, "%s\n", &cambpname[i]);
        fscanf(fpinput, "%s\n", cambnametmp);
        sprintf(cambpname[i], "%s", cambnametmp);

        if( debug>=DEBUGV*20) 
          printf("%s\n", cambpname[i]); fflush(stdout);
        }

      fclose(fpinput);

    //----------------------------------------
      psarr = dmat3(700, 2, surv.camblist_np);

      fompower_init(cambpname, psarr, surv.camblist_np); 

      if( debug>=DEBUGV ) 
        printf("proc_%d:  Initialized the power files.\n", mpi_rank);

    //----------------------------------------
      fpinput = fopen(cambpname[0], "r");

      init_powerInterp(cp.power, fpinput);

      fclose(fpinput);

    // End of Initializing powerspectrum.

//----------------------------------------------------------------
    // End of FoM calculation initialization.

    //----------------------------
    //debug printing
      if( debug>=20*DEBUGV || chain.debug>=20*DEBUGV ) {
        printf("survey name=%s\n", surv.name);
        printf("omem=%lg\n", cp.omem);
        printf("omek=%lg\n", cp.omek);
        printf("omeb=%lg\n", cp.omeb);
        printf("h0=%lg\n", cp.h0);
        printf("w0=%lg\n", cp.w0);
        printf("w1=%lg\n", cp.w1);
    //----------------------------
        printf("Include Planck %d (1 for ture)\n", surv.PLANCK);
        printf("Include SDSS   %d (1 for ture)\n", surv.SDSS);
        printf("total time = %lg\n", surv.time);
        printf("exposure time = %lg\n", surv.exptime);
        printf("max pointing number= %d\n", surv.nptmax);

        if(debug>=20*DEBUGV || chain.debug>=20*DEBUGV) {
          printf("npara=%d\n", chain.npara);
          printf("starting random=%d\n", chain.startrandom );
          for(i=0; i<chain.npara; i++) {
            printf("proc_%d:  start[%d]=%lg\n", mpi_rank, i, chain.startp[i]);
            printf("proc_%d:  lower[%d]=%lg\n", mpi_rank, i, chain.lowbound[i]);
            printf("proc_%d:  upper[%d]=%lg\n", mpi_rank, i, chain.upbound[i]);
            printf("proc_%d:  inisig[%d]=%lg\n",mpi_rank, i, chain.inisigma[i]);
            }

          }
        }
    //----------------------------
    // end of debug printing
//------------------------------------------------------------

      iniparser_freedict(dict);
//-------------------------------------------------------------------------------
    // end of reading ini initialization from parameter files.


//--------------------------------------------------------------------------------
// end of parameter initialization.




//-------------------------------------------------------------------------------

    //-----------------------------------------------------------
    //     Starting the MCMC sampling.                         //
    //-----------------------------------------------------------

//-------------------------------------------------------------------------------
// Start the sampler.

    // synchronization.

      MPI_Barrier(MPI_COMM_WORLD);

      //printf("\nproc_%d:cccccccccccccccccc\n", mpi_rank); fflush(stdout);

  #define _NO_TEST_
  #ifdef _NO_TEST_
 //--------------------------------------------------------
   // mpi_rank != 0, calculate fom
      if (mpi_rank) 
        slave(&chain, &testchi2, &info );  

 //--------------------------------------------------------
   // mpi_rank == 0, assign tasks to slaves.
      else  
        master(&chain, &testchi2, &info );
   #endif


      MPI_Barrier(MPI_COMM_WORLD);
//      printf("\n%d:ttttttttttttttttttttttttt\n"); fflush(stdout);

//-------------------------------------------
  // free all
      myinterp_free( cp.power);
      //myinterp_free(&eff);

      freemat((void **)info.fmat_cmb, surv.dim_cmb, surv.dim_cmb );
      freemat3((void ***)psarr, 700, 2, surv.camblist_np );

      freemat((void **)cambpname, surv.camblist_np, 100 );
      free(chain.startp);
      MPI_Finalize();

//      printf("\n%d:zzzzzzzzzzzzzzzzzzzzzzzzzz\n"); fflush(stdout);


//--------------------------------------------------
// THE END


    #endif


    }

//--------------------------------------------------
// END of main()



















//--------------------------------------------------
//--------------------------------------------------
//--------------------------------------------------
//--------------------------------------------------
//--------------------------------------------------
//--------------------------------------------------
// Parkinson's sampler.
// Start the sampler.
/*
      for(i=0; i<ctrl.max_step; i++)  {

    //????temperature:???????   line:313
    //????  temp_sched() ?????

         sconfig_cp(&s_cur, &s_test);

      //  geting new survey configurations.
      //  select "nearby"  surveys
         get_survey( &cp , &s_test, &surv );

//
         get_galaxies( &cp, &s_test, &surv );


        }
*/


//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
 // chain files.
/*
      char *chainame, *bestname;

      sprintf(chainame, "%s_%s_chain_%d.dat", surv.prefix, surv.name, mpi_rank+1);
      sprintf(bestname, "%s_%s_best_%d.dat", surv.prefix, surv.name, mpi_rank+1);


      if( debug>=DEBUGV || ctrl.debug>=DEBUGV ) {
        printf("chain name=%s\n", chainame);
        printf("best name=%s\n", bestname);
        }

      for(i=0; i<mpi_numtasks; i++) {
        chain.[i] = ;
        }
*/
//----------------------------------------------------------------









//--------------------------------------------------
// test iniparser
/*
      dictionary * dict;
      dict = iniparser_load(ini_name);

      if(debug>=DEBUGV) {
  
        printf( "survey = %s\n", iniparser_getstring( dict, "Survey:name", NULL) );
        printf( "%d\n",  iniparser_getint( dict, "Survey:area", -1) );
        }
*/
//--------------------------------------------------
