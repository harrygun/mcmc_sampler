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




//--------------------------------------------------
//--------------------------------------------------
// driver 
//-------------------------------------------------


  int main( int argc, char *argv[] ) {

//reading parameter .ini file
      int debug= 10, i, j;
      char *ini_name;


//      if(argc!=2)  //
//        myerr("Input parameters file is needed.", FALSE);
 
      ini_name =argv[1];

      if(debug>=DEBUGV)
      printf("reading %s\n", ini_name);



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

      surv.line = iniparser_getboolean( dict, "Survey:Eimission_Line_Galaxies", TRUE);

      surv.PLANCK = iniparser_getboolean( dict, "Survey:PLANCK", INITRUE);
      surv.SDSS   = iniparser_getboolean( dict, "Survey:SDSS", INITRUE);

      surv.cmbfish_file = iniparser_getstring( dict, "Survey:CMB_FISHER", "fmat_cmb.dat");
      surv.dim_cmb = iniparser_getint( dict, "Survey:DIM_FISHER_CMB", 7);

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
   
      chain.npara     = iniparser_getint( dict, "Chain:parameter_number", 4);
      chain.npara_tot = iniparser_getint( dict, "Chain:parameter_number_total", 5);
      chain.npoint    = iniparser_getint( dict, "Chain:chain_length", 20000);
      chain.max_point = iniparser_getint( dict, "Chain:chain_max_length", 25000);

      chain.debug   = iniparser_getint( dict, "MCMC_Control:debug", DEBUGV);
      surv.debug = chain.debug;


//----------------------------------------------------------------
    // chains & results files
//      chain.fname_chain= charmat(chain.nchain, 100);
//      chain.fname_best = charmat(chain.nchain, 100);
//
//      for(i=0; i<chain.nchain; i++) {
//        sprintf(chain.fname_chain[i], "%s_%s_chain_%d.dat", surv.prefix, surv.name, i+1 );
//        sprintf(chain.fname_best[i] , "%s_%s_best_%d.dat" , surv.prefix, surv.name, i+1 );
//        }
//----------------------------------------------------------------
    // get the starting points, boundary & initial sigma. 

      char initemp[100];
      int *npara_count, nouse=0;
      double *upb, *lowb, *inisig;

      npara_count = ivec(chain.npara_tot);
      lowb  = dvec(chain.npara_tot);
      upb   = dvec(chain.npara_tot);
      inisig= dvec(chain.npara_tot);

      surv.starp = dvec(chain.npara_tot);
      //for(i=0; i<chain.npara; i++) {
      for(i=0; i<chain.npara_tot; i++) {
        sprintf(initemp, "Chain:start_param%d", i+1);
        surv.starp[i] = iniparser_getdouble( dict, initemp, 1.);

        sprintf(initemp, "Chain:lower_param%d", i+1);
        lowb[i] = iniparser_getdouble( dict, initemp, 1.);

        sprintf(initemp, "Chain:upper_param%d", i+1);
        upb[i] = iniparser_getdouble( dict, initemp, 1.);

        sprintf(initemp, "Chain:sigma_param%d", i+1);
        inisig[i] = iniparser_getdouble( dict, initemp, 1.);
       
        if( (upb[i]-lowb[i])<1e-10*surv.starp[i] ) {
          nouse++;
          npara_count[i] = 99;
          }
        else  
          npara_count[i] = i-nouse;

        printf("npara_count[%d]=%d\n", i, npara_count[i]);
        }

    //-------------------------------------------
      chain.npara = chain.npara_tot-nouse;
      surv.npara     = chain.npara;
      surv.npara_tot = chain.npara_tot;
    //-------------------------------------------
      chain.T        = dvec(chain.nchain);
      dvecvalue(chain.T, chain.nchain, 1.);

      chain.startp   = dvec(chain.npara);
      chain.lowbound = dvec(chain.npara);
      chain.upbound  = dvec(chain.npara);
      chain.inisigma = dvec(chain.npara);
      chain.size    = chain.npara + 2;

      surv.paraindex = ivec(chain.npara);

      if(surv.npara!=surv.npara_tot)
        surv.paraoutdex = ivec(chain.npara_tot-chain.npara);
    //-------------------------------------------
      nouse = 0;
      for(i=0; i<chain.npara_tot; i++) {

        if(npara_count[i]<=chain.npara_tot)
          surv.paraindex[npara_count[i]] = i;
        else{
          surv.paraoutdex[nouse] = i;
          nouse++;
          }
        }

      for(i=0; i<chain.npara; i++) {
        chain.startp[i]   = surv.starp[ surv.paraindex[i]];
        chain.lowbound[i] = lowb[surv.paraindex[i]] ;
        chain.upbound[i]  = upb[surv.paraindex[i]]  ;
        chain.inisigma[i] = inisig[surv.paraindex[i]] ;

        printf("surv.paraindex[%d]=%d\n", i, surv.paraindex[i]);
        }

      for(i=0; i<surv.npara_tot-surv.npara; i++)
        printf("surv.paraoutdex[%d]=%d\n", i, surv.paraoutdex[i]);


//----------------------------------------------------------------
//-------------------------------------------------------------------------------------

//-----------------------------------
    // Initializing efficency.

      Interpar eff;
      surv.eff = &eff;

      if(surv.LAMOST==TRUE) 
        eff_init(&surv);

//-----------------------------------
    // get CMB Fisher matrix.
      if( debug>=DEBUGV ) 
        printf("Preparing for CMB Fisher Matrix calculation.\n" );

      char term;
      FILE *fpinput= fopen(surv.cmbfish_file, "r");
      //printf("%s\n", surv.cmbfish_file);

      info.fmat_cmb = dmat(surv.dim_cmb, surv.dim_cmb);

      for(i=0; i<surv.dim_cmb; i++) {
        for(j=0; j<surv.dim_cmb-1; j++) 
          fscanf(fpinput, "%lg ", &info.fmat_cmb[i][j] );
        fscanf(fpinput, "%lg%c", &info.fmat_cmb[i][surv.dim_cmb-1], &term );
        }

//      mat_sqrtdiag_printf(info.fmat_cmb,  7);

      if( debug>=DEBUGV ) 
        printf("Read the CMB Fisher Matrix.\n");

      fclose(fpinput);

//-----------------------------------
    // Initialize powerspectrum for Fisher calculation.

      fpinput = fopen(surv.camblist, "r");
      char **cambpname= charmat( surv.camblist_np, 100 );
      char cambnametmp[100];
//      char cambpname[surv.camblist_np][100];


      for(i=0; i<surv.camblist_np; i++) {
//        fscanf(fpinput, "%s\n", &cambpname[i]);
        fscanf(fpinput, "%s\n", cambnametmp);
        sprintf(cambpname[i], "%s", cambnametmp);

        if( debug>=DEBUGV*20) 
          printf("%s\n", cambpname[i]);
        }

      fclose(fpinput);

    //----------------------------------------
      //printf("camblist_np=%d\n", surv.camblist_np);
      psarr = dmat3(700, 2, surv.camblist_np);

      fompower_init(cambpname, psarr, surv.camblist_np); 

      if( debug>=DEBUGV ) 
        printf("Initialized the power files.\n" );

    //----------------------------------------
      fpinput = fopen(cambpname[0], "r");

      init_powerInterp(cp.power, fpinput);

      fclose(fpinput);

    // End of Initializing powerspectrum.

//----------------------------------------------------------------
    // End of FoM calculation initialization.


    //----------------------------
    //debug printing
      if( debug>=2*DEBUGV || chain.debug>=20*DEBUGV ) {
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

        if(debug>=2*DEBUGV || chain.debug>=20*DEBUGV) {
          printf("npara=%d\n", chain.npara);
          for(i=0; i<chain.npara; i++) {
            printf("start[%d]=%lg\n", i, chain.startp[i]);
            printf("lower[%d]=%lg\n", i, chain.lowbound[i]);
            printf("upper[%d]=%lg\n", i, chain.upbound[i]);
            }

          }
        }
    //----------------------------
    // end of debug printing
//------------------------------------------------------------
    // initialize info for FoM calculation.

      info.cp    = &cp;
      info.sp    = &sp;
      info.fpara = &fpara;
//      info.scp   = &scp;

      info.surv  = &surv;
      info.s     = &sconf;


      info.psarr= psarr;

      iniparser_freedict(dict);
//-------------------------------------------------------------------------------
    // end of reading ini initialization from parameter files.


//--------------------------------------------------------------------------------
// end of parameter initialization.



      int fail;
      double fom;


  #define _GET_FOM_
  #ifdef _GET_FOM_


    #define _FOM_LATTICE_
    #ifdef _FOM_LATTICE_

      FILE *fpfom=fopen("fom_test_chain_1.dat", "w");
      int pindex, step = 10;
      double p;

      //surv.maxexptime= 2000;
      //surv.minexptime= 1;
      pindex = 4;
      for(i=0; i<step; i++) {
        chain.startp[pindex] =  chain.lowbound[pindex]+ i*
                    (chain.upbound[pindex] -chain.lowbound[pindex] )/(double)step;

        fom = get_fom( &info, chain.startp, &fail );

        for(j=0; j<chain.npara; j++)
          fprintf(fpfom, "%lg  ", chain.startp[j]);
        fprintf(fpfom, "%lg  ", fom);
        fprintf(fpfom, "1\n" );

        fflush(fpfom);

        printf("\n!!!!!!!!!!!!!\nfom=%lg\n!!!!!!!!!!!!!\n\n", exp(fom));

        }

    #else


      for(i=0; i<2; i++) {
      //for(i=0; i<100; i++) {
        chain.startp[0] =  chain.startp[0]+i;
        //chain.startp[4] =  chain.startp[4]+i*0.1;

        fom = get_fom( &info, chain.startp, &fail );

        printf("fom=%lg\n", exp(fom));
        }

    #endif
  #endif



     //#define _GET_GALAXIES_N_
      #ifdef _GET_GALAXIES_N_
      int step=1000;
      FILE *fpout = fopen("ngal.dat", "w" );
      double dt;

      surv.nptmax = 1;
      surv.maxexptime= 2000;
      //surv.minexptime= 10;
      dt = chain.startp[1]/surv.fov;  // area/fov

      //timepfov = s->time*surv->fov/s->area;

      for(i=90; i<91; i++)  {
      //for(i=990; i<991; i++)  {
      //for(i=10; i<step; i++)  {
        chain.startp[0] = (i+1)*dt;

        para2sconfig(chain.startp, info.s, info.surv );
        get_galaxies(info.cp, info.sp,  info.s, info.surv, &fail);
        
        fprintf(fpout, "%lg   %lg\n", chain.startp[0]/dt, info.sp->n[sp.zbin-1]);
        fflush(fpout);
        }

      fclose(fpout);

      #endif



//-------------------------------------------------------------------------------

     myinterp_free( cp.power);
     myinterp_free(&eff);

     freemat((void **)info.fmat_cmb, surv.dim_cmb, surv.dim_cmb );
     freemat3((void ***)psarr, 700, 2, surv.camblist_np );

     freemat((void **)cambpname, surv.camblist_np, 100 );
     free(chain.startp);
      //info.fmat_cmb = dmat(surv.dim_cmb, surv.dim_cmb);
      //psarr = dmat3(700, 2, surv.camblist_np);
      //char **cambpname= charmat( surv.camblist_np, 100 );
      //chain.startp   = dvec(chain.npara);
//-------------------------------------------------------------------------------



    }














//--------------------------------------------------

      //chain.startp[0] = iniparser_getdouble( dict, "Chain:start_param1", 0.5);
      //chain.startp[1] = iniparser_getdouble( dict, "Chain:start_param2", 0.5);
      //chain.startp[2] = iniparser_getdouble( dict, "Chain:start_param3", 0.5);
      //chain.startp[3] = iniparser_getdouble( dict, "Chain:start_param4", 0.5);
      //chain.startp[4] = iniparser_getdouble( dict, "Chain:start_param5", 0.5);

/*


      chain.startp   = dvec(chain.npara);
      chain.lowbound = dvec(chain.npara);
      chain.upbound  = dvec(chain.npara);
    //-------------------------------------------


      for(i=0; i<chain.npara; i++) {
        sprintf(initemp, "Chain:start_param%d", i+1);
        chain.startp[i] = iniparser_getdouble( dict, initemp, 1.);

        sprintf(initemp, "Chain:lower_param%d", i+1);
        chain.lowbound[i] = iniparser_getdouble( dict, initemp, 1.);

        sprintf(initemp, "Chain:upper_param%d", i+1);
        chain.upbound[i] = iniparser_getdouble( dict, initemp, 1.);
        }


*/


