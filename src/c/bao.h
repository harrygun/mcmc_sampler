
  #ifndef _H_BAO_
    #define _H_BAO_

    void init_bao( Cospar * cp, Survpar * sp );

//    void baofisher(Fipar * fpara,  Survpar * sp);
    void baofisher(Fipar * fpara,  Survpar * sp, double *** Pk);


    void select_subao( Fipar * fpara, Survpar * sp, int *vec , int dim );


    void baods2de( Fipar * fpara, Survpar * sp );

  #endif
