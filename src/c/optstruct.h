
  #ifndef _H_OPT_STRUCT_
    #define _H_OPT_STRUCT_

    #define DEBUGV 2

//----------------------------------------


    typedef struct  {

//      char name[100], prefix[100], cmbfish_file[100], camblist[100];
      char *name, *prefix, *cmbfish_file, *camblist;

      int random, dim_cmb, camblist_np, *paraindex, *paraoutdex, 
          npara, npara_tot; 

      int PLANCK, SDSS, LAMOST;   //include as prior.
      int line, debug;


//--------------------------------------
    // instruments
    // Field-of-View, # of fibres.
      int nptmax;
      double fov, nfib,  
             fibdiam, mirdiam, apert;  //fiber, mirror diameter; aperture factor

      double time, minexptime, maxexptime, exptime, overheadtime, maxarea, 
             zmin, zmax, dzmin, dzmax, *starp;

      Interpar *eff;  // effeciency.

      int eff_init;
//--------------------------------------

      } survey_setting;





    typedef struct {

      int bin, npt; // npt: # of pointings

      double time, area, zmin, zmax, sncont;  // time: exposure time


      } survey_config;





//---------------------------------------

    void sconfig2para( survey_config *s, double *para, survey_setting * surv) ;

    void para2sconfig(double *para, survey_config *s, survey_setting * surv) ;


    void sconfig_cp( survey_config *s1, survey_config *s2);


  #endif




/*
    typedef struct {

        int debug;
        double max_step;

      } prog_ctrl; */

