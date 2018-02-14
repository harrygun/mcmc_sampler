
  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>
  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>

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





//--------------------------------------------------------------
  void sconfig2para( survey_config *s, double *para, survey_setting * surv) {

     // para[] = s->bin; 
      //para[0] = s->time;
      //para[1] = s->area;
      //para[2] = s->zmin;
      //para[3] = s->zmax; 

      int i;
      double *paratot=dvec(surv->npara_tot);

      paratot[0] = s->time;
      paratot[1] = s->area;
      paratot[2] = s->zmin;
      paratot[3] = s->zmax; 
      paratot[4] = s->sncont; 

      for(i=0; i<surv->npara; i++)
        para[i] = paratot[ surv->paraindex[i] ];

      free(paratot);
      return;
    }


//--------------------------------------------------------------
  // cp para to sconfig, and initialize surveys.
  void para2sconfig(double *para, survey_config *s, survey_setting * surv) {

      int i;
      double *paratot=dvec(surv->npara_tot);

      //printf("Entered para2sconfig():\n");
      //for(i=0; i<5; i++)
      //  printf("para2config: para[%d]=%lg\n", i, para[i]);

      for(i=0; i<surv->npara; i++)
        paratot[ surv->paraindex[i] ] = para[i];

      if(surv->npara!=surv->npara_tot) {
        for(i=0; i< surv->npara_tot-surv->npara; i++ )
          paratot[ surv->paraoutdex[i] ] = surv->starp[surv->paraoutdex[i]];
        }

      s->time   = paratot[0];
      s->area   = paratot[1];
      s->zmin   = paratot[2];
      s->zmax   = paratot[3]; 
      s->sncont = paratot[4]; 


//-------------------------------------

      free(paratot);
      return;
    }


//--------------------------------------------------------------







//--------------------------------------------------------------
//   s1-->s2
  void sconfig_cp( survey_config *s1, survey_config *s2) {

      s2->bin = s1->bin  ; 
      s2->time = s1->time  ;
      s2->area = s1->area  ;
      s2->zmin = s1->zmin  ;
      s2->zmax = s1->zmax  ; 

//      s2-> = s1->  ;
//      s2-> = s1->  ;

      return;
    }







