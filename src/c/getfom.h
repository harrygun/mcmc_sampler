
  #ifndef _H_GET_FOM_
    #define _H_GET_FOM_



    double get_fom( void *inform, void *arg, int *fail  );


    void get_galaxies( Cospar *cp, Survpar * sp, survey_config *s, 
              survey_setting *surv, int *fail );

    int fompower_init(char **name, double ***power, int n) ;

    void eff_init(survey_setting *surv) ;

  #endif
