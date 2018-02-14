
  #ifndef _H_SAMPLER_
    #define _H_SAMPLER_




    typedef struct {

        int startrandom, checkconverg, simpleconverg, debug;

        int size, npoint, npara, npara_tot, nchain, 
            max_point, temptype;

        double *T;

        double *startp, *lowbound, *upbound, *inisigma;


        char **fname_chain, **fname_best; 


      } chainpar;


    typedef struct {

      int converg, min_chain;

      double ***chains, **currentpoint, **trialpoint, ***covmat,
             ***eigenvec, **point_stepsize , **eigenvalue, 
             *scale,  *gauss_x;

     int *point_num, *update_num, *loop_num;

      
      randpara *random;

      } chainvar;



    typedef struct {

      Cospar  *cp;
      Survpar * sp;
//      SCpar *scp;

      Fipar *fpara; 

      survey_config *s;
      survey_setting *surv;

      double *** psarr;
      int npower_input;  //, fail;

      double **fmat_cmb, **fmat_cb, **fmat_cb_inv;


      } fom_info;



//----------------------------------

    void try_newpoint(chainpar *chain, chainvar *chain_var, int chain_n) ;

    void convergence_test( chainpar *chain, chainvar *chain_var ) ;

    void temperature_sched( chainpar* chain, chainvar *chain_var, int chain_n);

    int master(chainpar *chain,  double (*chi2)(void *inform, void *arg, int *x), void *info )  ;

    int slave( chainpar *chain, double (*chi2)(void *inform, void *arg, int *x),  void *info)  ;

  #endif
