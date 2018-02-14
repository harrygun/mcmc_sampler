#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>


#include "const.h"
#include "varb.h"
#include "mymath.h"
#include "matrix.h"
#include "cospara.h"
#include "power.h"
#include "bao.h"
#include "myerr.h"
#include "ps.h"
//#include "corrf.h"


  #define PN_BAO    5   // # of parameters per bin
  #define CosP_BAO   6  
  
  #define PN_BAO_s    2 
  #define CosP_BAO_s  1// 2
  
  #define _BIN_N_5_
  
  
  
  //temp
  
  //#define PN_BAO    2   // # of parameters per bin
  //#define CosP_BAO   6  
  //
  //#define PN_BAO_s    2 
  //#define CosP_BAO_s  1
  ////
  
  
  //#define Powline 500
  #define Powline 690
  
  //#define Knum 100
  #define Knum 427  //427
  
  
  //#define Knum 500
  //#define Knum 393
  
  #define K_min  232 //0

  #define MUnum 50 //100
  
  
  //#define _NL_POWER_
  
  
  #define WIND_SMOOTH 4.



//#define _DOMEA_
#define _DOMEA_PHOTO_

//------------------------------------------

  void init_bao( Cospar * cp, Survpar * sp )
    {
      int i, j, bin ;
      double zmin, zmax;


  #ifdef _DOMEA_
      bin = 10;
      zmin = 0;
      zmax = 2.5;
      sp->area= 5000;

      sp->phot = FALSE;//TRUE;
      sp->sz0  = 1e-3;

  #endif

  #ifdef _DOMEA_PHOTO_
      bin = 10;
      zmin = 0;
      zmax = 4.;
      sp->area= 5000.;

      sp->phot = TRUE;
      sp->sz0  = 5e-3;//1e-2;

  #endif


  #ifdef _TEST_
      bin = 4;//4;//1;//3;
      zmin = 0.5;//0.2;//0.5;//0.75;//0.4;//0.5;//0.;
      zmax = 1.3;//1.3;//1.25 ; //0.55;//1.;//0.38;
      sp->area= 2000;//8000.;//4000.;

      sp->phot = FALSE;//TRUE;
      sp->sz0  = 1e-3;

  #endif


      sp->zbin = bin;

      sp->z    = dvec(bin);
      sp->b    = dvec(bin);
      sp->z1   = dvec(bin);
      sp->z2   = dvec(bin);
      sp->n    = dvec(bin);
      sp->v    = dvec(bin);
      sp->beta = dvec(bin);


//      sp->b[0]= 1.77;       // 1.8;// 1.77;  // sp->b[1]=1.5;
//      sp->n[0]= 0.00078837;// 1e-4;//0.00078837;  //  sp->n[1]=3e-4;
//      sp->b[0]=  2.;  // sp->b[1]=1.5;
//      sp->n[0]=  1.e-3;  //  sp->n[1]=3e-4;


      for(i=0; i<bin; i++)   {
        sp->z1[i] = zmin +  i * (zmax - zmin)/bin ;
        sp->z2[i] = zmin + (i+1) * (zmax - zmin)/bin ;

        sp->z[i] = (sp->z1[i] + sp->z2[i]) /2.; //0.3 ;//(sp->z1[i] + sp->z2[i]) /2.;
        sp->v[i] = vol(cp, sp->z1[i],  sp->z2[i]) * sp->area;

      #ifdef _DOMEA_
        sp->b[i] =  0.7 + 0.09*i;
      #endif 
      #ifdef _DOMEA_PHOTO_
        sp->b[i] = 0.7 + 0.15*i;
      #endif

      #ifdef _TEST_
        sp->b[i]=  2.*sqrt( 1.+sp->z[i]);//1.8;//2.;  // sp->b[1]=1.5;
      #endif

        sp->n[i]=  8.e-4/ sqrt(1.+sp->z[i]);//1.e-4;//1.e-3;  //  sp->n[1]=3e-4;

        sp->beta[i] = Beta(cp,sp->z[i], sp->b[i] );


        }

//        sp->z[0] = 0.3;//0.475;
//        sp->z[0] = 0.475;
//        sp->z[0] = 1.;

  #ifdef _DOMEA_

      sp->n[0]= 0.124347  ;    
      sp->n[1]= 0.0589665 ; 
      sp->n[2]= 0.0299451 ;    
      sp->n[3]= 0.0162629 ; 
      sp->n[4]= 0.00932294;    
      sp->n[5]= 0.00560978; 
      sp->n[6]= 0.00353169;    
      sp->n[7]= 0.00232016; 
      sp->n[8]= 0.00158632;    
      sp->n[9]= 0.00112556; 


  #endif

  #ifdef _DOMEA_PHOTO_

      sp->n[0]= 0.187778  ;    
      sp->n[1]= 0.0962984 ; 
      sp->n[2]= 0.0537991 ;    
      sp->n[3]= 0.0327748 ; 
      sp->n[4]= 0.0213133 ;    
      sp->n[5]= 0.0145934 ; 
      sp->n[6]= 0.0104278 ;    
      sp->n[7]= 0.00772925; 
      sp->n[8]= 0.00591753;    
      sp->n[9]= 0.00466527; 

  #endif


//-------------------------------------------------
      for(i=0; i<bin; i++)  {
        printf("z1=%lg, z2=%lg, v=%lg\n", sp->z1[i], sp->z2[i], sp->v[i]);
        printf("beta[%d] = %lg, b = %lg  n = %lg\n", i, sp->beta[i], sp->b[i], sp->n[i]);
        }


      return;
    }



  void veff_init( SCpar * scp, double *** veff, 
                                     double * k, double *mu ) 

    {                               // veff[k][mu][bin]
      int i, m, n, bin;
      double p, p0, kmax, dk, dmu, kmid, mumid, 
             sz, sr ;

      kmax = 0.5;

      bin = scp->sp->zbin;

//      for(m=0; m<Knum; m++)
//        k[m] = m * kmax/Knum;

      for(n=0; n<MUnum; n++)
        mu[n] = -1. +  n * 2./(double)MUnum;

//      for(i=0; i<bin; i++)
//        printf("P(k=0.2)= %lg\n", Power(0.2) * 
//           squar(scp->sp->b[i] * Growth(scp->cp,scp->sp->z[i] )  ) *scp->sp->n[i]   );

      for(m=0; m<Knum-1; m++)   
        for(n=0; n<MUnum-1; n++)  {

          dk  = k[m+1]- k[m];
          dmu = mu[n+1]- mu[n];

          kmid = ( k[m] + k[m+1] )/2.;
          mumid = ( mu[n] + mu[n+1] )/2.;

//          p = scp->cp->As * Power( kmid ) ;
//          p = Power(kmid) ;

          p0 = powerInterp(scp->cp->power, kmid) ;


          for(i=0; i<bin; i++ )  {
            p = p0 * squar( scp->sp->b[i] * Growth(scp->cp,scp->sp->z[i] )
                  *  (1. + scp->sp->beta[i]* squar(mumid)  )  ) ;

            if( scp->sp->phot == TRUE)   {
              sz = scp->sp->sz0 * (1. + scp->sp->z[i] );
              sr = Vc * sz / Hubble(scp->cp, scp->sp->z[i]);

              p = p0 * exp(- squar(kmid * mumid * sr)  );
              }

          veff[m][n][i] = squar(scp->sp->n[i]*p / ( scp->sp->n[i]*p +1.) )
                     * scp->sp->v[i] ;

            }


          }


      return;
    }


  double aniw2(Cospar  *cp,Survpar * sp, double k, 
                                    double R, double mu)
    {
      double window, omemz;

//      omemz= cp->omem * cub(1.+sp->z[sp->bini]) /
//                 squar( Eubble(cp, sp->z[sp->bini] ) );


//      window=exp(-squar(k* R)
//              * ( (1.- mu*mu) +
//               mu*mu* squar(1.+ pow(omemz, 0.6) ) ) );

//      window=exp(-squar(k* R)
//              * ( (1.- mu*mu) +
//               mu*mu* squar(1.+ 1. ) ) );

      window=exp(-squar(k* R));

     //printf("f = %lg\n", 1.+ pow(omemz, 0.6)  );

      return window;
    }


  double pobs( Fipar * fpara, Survpar * sp, double *spbini, double *dp,
                      double * kk, double *mmu, int m, int n )
    {
      double da, h, g, beta, pst, 
            da0, h0, g0, beta0 ;
      double po, omemz, b , anist, mu, k ;

      da0    =  exp( spbini[0]);
      h0     =  exp( spbini[1]); 
      g0     =  exp( spbini[2]); 
      beta0  =  exp( spbini[3]); 


      da = exp( spbini[0] + dp[0] );
      h  = exp( spbini[1] + dp[1] );
      g  = exp( spbini[2] + dp[2] );
      beta=exp( spbini[3] + dp[3] );
      pst = dp[4];



      mu = mmu[n];
      k=  sqrt( squar( kk[m]*mu*h/h0 ) + (1-mu*mu) *squar( kk[m]*da0/da ));

      omemz=fpara->cp->omem * cub(1.+ sp->z[sp->bini] ) 
                       / squar( Eubble(fpara->cp, sp->z[sp->bini]) );

      b = pow(omemz, 0.6) / beta;

      anist =  squar( mu* h* da)/ ( squar( mu*h*da) + 
                      squar(da0 * h0)*(1.-squar(mu)) );

//      po = squar(da0/da)* (h/h0) *squar( 1. + beta* anist ) * 
//                       squar(b * g) * Power(k) + pst;
      po = squar(da0/da)* (h/h0) *squar( 1. + beta* anist ) * 
                       squar(b * g) * powerInterp( fpara->cp->power, k) + pst;

      return log(po);
    }


  double derp( Fipar * fpara, Survpar * sp, double *spbini,
                   int pt, double * kk, double * mmu, int m, int n)
    {                       // dp/ d( LSS parameters )
      int p, i;
      double dpdp, *dp, *dp1, *dp2, dd;


      p =  pt - CosP_BAO -  sp->bini * PN_BAO ;


      if(p< 0 || p>= PN_BAO)
        myerr("derp() Error.\n",FALSE);

      dp  = dvec(PN_BAO);
      dp1 = dvec(PN_BAO);
      dp2 = dvec(PN_BAO);



      dd = 2.e-2;
      dp[p] =   dd;

      for(i=0; i< PN_BAO -1 ; i++)  {

        dp1[i] = spbini[i]*dp[i];
        dp2[i] = -dp1[i];
        }

      dp1[PN_BAO-1]=  dp[PN_BAO-1];
      dp2[PN_BAO-1]=  -dp[PN_BAO-1];


//

      dpdp = ( pobs( fpara, sp, spbini, dp1, kk, mmu, m, n ) -  
               pobs( fpara, sp, spbini, dp2, kk, mmu, m, n ))/ 2./ dp1[p];

      free(dp);
      free(dp1);
      free(dp2);

      return dpdp;
    }


  void dpdsp( Fipar * fpara, Survpar * sp, double *spbini, double *dpdp,
                           double *kk, double *mmu, int m, int n  )
    {
      int i, sta, end, bin;

      bin = sp->bini;

      sta = PN_BAO * bin + CosP_BAO;
      end = sta + PN_BAO;

      for(i=0; i< fpara->n; i++)
        dpdp[i] = 0;
      

      for(i=sta; i<end; i++)   {
        dpdp[i] = derp(fpara, sp, spbini, i, kk, mmu, m, n ) ;
//        if(bin>0 )
//          printf("bin=%d, dpdp[%d] = %lg\n", bin, i, dpdp[i]);
        }
        

      return;
    }


  void baofish_one(Fipar * fpara, Survpar * sp, double *dpdp, 
                    double ***veff, double *kk, double *mmu, int m, int n )
    {
      int i, j, sta, end, bin;

      double k, mu, dk, dmu, nl, sper, spar, omemz, aniwindow;

      bin = sp->bini;

      sta = PN_BAO *bin + CosP_BAO;
      end = sta + PN_BAO;

      k = kk[m];                  mu = mmu[n];
      dk = kk[m+1] - kk[m];     dmu= mmu[n+1] - mmu[n];

      
      aniwindow =  squar( aniw2 (fpara->cp, sp, k, WIND_SMOOTH, mu) );

      omemz= //fpara->cp->omem;
        fpara->cp->omem * cub(1.+sp->z[sp->bini]) / 
                 squar( Eubble(fpara->cp, sp->z[sp->bini] ) );

      sper = squar(12.4 * Growth(fpara->cp, sp->z[sp->bini])* 0.758);
      spar = sper * squar(1.+ pow(omemz, 0.6) );

    #ifdef _NL_POWER_

      nl = exp( - squar(k ) * (  (1.- squar(mu) )* sper /2. 
                  + squar(mu) * spar/2.)  );

      nl = nl ;//* squar( aniwindow);
    #else
      nl = 1.;//squar(aniwindow);//nl * aniwindow;
    #endif


      for(i=0; i<CosP_BAO; i++)
        for(j=0; j<CosP_BAO; j++)
          fpara->tmpmat[i][j] = dpdp[i]*dpdp[j]*squar(k)* dk* 
                nl * veff[m][n][bin]* dmu /squar(2.* PI) /2.;


      for(i=0; i<CosP_BAO; i++)
        for(j=sta; j<end; j++ )  {
          fpara->tmpmat[i][j] = dpdp[i]*dpdp[j]*squar(k)*dk* 
                nl * veff[m][n][bin]* dmu /squar(2.* PI) /2.;

         fpara->tmpmat[j][i] = fpara->tmpmat[i][j];
         }


      for(i=sta; i<end; i++ )
        for(j=sta; j<end; j++ )
          fpara->tmpmat[i][j] = dpdp[i]*dpdp[j]*squar(k)*dk* 
                nl * veff[m][n][bin]* dmu /squar(2.* PI) /2.;


      return;
    }



  void baofisher(Fipar * fpara,  Survpar * sp, double *** Pk)
    {

      int i, j, l, m, n, s, t,
          bin, dsn, den, subn, pn ;

      //printf("baofisher(): \n");

// initialize survey

      SCpar  scp;

      scp.sp = sp;
      scp.cp = fpara->cp;



      bin = sp->zbin;

      dsn  =  bin * PN_BAO + CosP_BAO;   // dimof(full Fisher matrix)
      subn =  bin * PN_BAO_s + CosP_BAO_s + 1;  // add da_cmb
      den  = 4;//fpara->paran;             

      pn =  2* CosP_BAO + 1;

//--------------------------------------
      fpara->fmat    =dmat(dsn, dsn);
      fpara->fmat_inv=dmat(dsn, dsn);


      fpara->n    = dsn;
      fpara->subn = subn;
      fpara->convn= den;


// reading power

//      int flin, fcol;
//      FILE *fp;
//      char term, * Pfile[pn];
//      double *** Pk;

//      Pk = dmat3(Powline, 2, pn);

/*
//      Pfile[0]="fiducial_matterpower.dat";
      Pfile[0]="fiducial_long_matterpower.dat";

      Pfile[1]="omebh2+_matterpower.dat";
      Pfile[2]="omebh2-_matterpower.dat";
      Pfile[3]="omech2+_matterpower.dat";
      Pfile[4]="omech2-_matterpower.dat";
      Pfile[5]="h+_matterpower.dat";
      Pfile[6]="h-_matterpower.dat";
      Pfile[7]="ns+_matterpower.dat";
      Pfile[8]="ns-_matterpower.dat";
      Pfile[9]="AS+_matterpower.dat";
      Pfile[10]="AS-_matterpower.dat";
      Pfile[11]="tau+_matterpower.dat";
      Pfile[12]="tau-_matterpower.dat";



      for(i=0; i<pn; i++)
            {
            if((fp=fopen(Pfile[i],"r"))==NULL)
                    printf("Can't open the file %s.\n",Pfile[i]);

            flin = Powline;//countFILEline( fp );
            fcol = 2;

            for(j=0;j<flin; j++) {

              for(l=0; l<fcol-1 ;l++)
                fscanf(fp , "%lg ",&Pk[j][l][i]);

              fscanf(fp,"%lg%c",&Pk[j][fcol-1][i],&term);

              }

            fclose(fp);
            } */

// re-initialize the powerspectrum with fiducial parameters.
/*
       fp=fopen(Pfile[0],"r");
       init_power(fp);

       InitAS(fpara->cp);

       fclose(fp);
*/

//  effective volume

      int kn, kmin, mun;
      kn = Knum;            mun = MUnum;
      kmin = K_min;

      double *** veff, *k, *mu;
      veff = dmat3( kn, mun, bin );   // veff[k][mu][bin]

      k  = dvec(kn);   
      mu = dvec(mun);

      for(i=0; i< kn; i++)
        k[i] = Pk[i][0][0];

      veff_init( &scp, veff, k, mu);




//     d Power/d Cospar  without Da, H
      double ** dpdcp, **dtmp, * dcp,  *dpdp, *spbini;

      dpdcp = dmat( kn, CosP_BAO);
      dcp = dvec(CosP_BAO);
      dpdp = dvec( dsn );
      spbini = dvec(PN_BAO);

      dtmp = dmat( kn, CosP_BAO);


      dcp[0] = 2e-2 * 0.0223*2.;    dcp[1]= 2e-2* 0.1057* 2.; // omch2, ombh2
      dcp[2] = 3e-2 * 0.73 *2. ;      //h
      dcp[3] = 1e-3 * 2. ;  dcp[4] = 1e-3 * 2.3e-9 *2. ;// ns, As
      dcp[5] = 1e-2 * 2. ;   // tau

      for(i=0; i< CosP_BAO; i++)  
        for(m=0; m<kn; m++)  {
          dpdcp[m][i] = ( log(Pk[m][1][2*i+1])- 
                          log(Pk[m][1][2*i+2] ) )/ dcp[i];
          dtmp[m][i] = dpdcp[m][i];
          }

        for(m = 0; m< kn ; m++)  {

          dpdcp[m][2] = -0.5*fpara->cp->h0/ fpara->cp->omem * dtmp[m][2]  ;   // omem

          dpdcp[m][0] = dtmp[m][1] -dpdcp[m][2] /squar(fpara->cp->h0) ;   //omemh2

          dpdcp[m][1] = dtmp[m][0] - dpdcp[m][0] - dpdcp[m][2] / squar(fpara->cp->h0);   //omebh2
          }

      freedmat(dtmp, kn, CosP_BAO);



//  cal distance fisher  matrix

      double ** tmpmat, elem; 

      tmpmat = dmat(dsn, dsn);
      fpara->tmpmat = tmpmat;

//      FILE *fp_te = fopen("tebao.dat","wb");


      for(i=0; i< bin; i++)   {
  //      printf("BAO: %dth bin.\n",i+1);
  
        sp->bini = i;
  
        spbini[0]= log( Da( fpara->cp, sp->z[i] ) );
        spbini[1]= log( Hubble( fpara->cp, sp->z[i] ) );

  #ifdef _BIN_N_5_
        spbini[2]= log( Growth( fpara->cp, sp->z[i] ) );
        spbini[3]= log( Beta( fpara->cp, sp->z[i], sp->b[i] ) );
        spbini[4]= 0.;
  #endif


        for(m=kmin; m<kn-1; m++ )
          for(n=0; n<mun-1; n++)   {
  
            dpdsp(  fpara,  sp, spbini, dpdp, k, mu, m, n  );

            for(j=0; j<CosP_BAO; j++) 
              dpdp[j] =  dpdcp[m][j];  

            dmatzero(tmpmat, dsn, dsn);
            baofish_one(fpara, sp, dpdp, veff, k, mu, m, n );
  

//            fprintf(fp_te,"\n\nfpara->tmpmat; bin=%d\n", i);
//            mat_fprintf(fp_te,fpara->tmpmat, dsn, dsn, "%lg  " );

            for(s=0; s<dsn; s++)
              for(t=0; t<dsn; t++)  {

                elem = fpara->fmat[s][t] + fpara->tmpmat[s][t] ;
                fpara->fmat[s][t]  = elem;
                }
  
            }
        }
  

//      fclose(fp_te);


      myInv(fpara->fmat, fpara->fmat_inv, dsn);


      //printf("baofisher(): End\n");

// free memory

        //free(veff);   // free(dpdcp);    
        //veff = dmat3( kn, mun, bin );   // veff[k][mu][bin]
        //dpdcp = dmat( kn, CosP_BAO);
        freedmat3(veff, kn, mun, bin );
        freedmat(dpdcp, kn, CosP_BAO);
        freedmat(tmpmat, dsn, dsn);
 
        free(dcp);     free(dpdp);
        free(k);       free(mu);       free(spbini);


      return;
    }


    void select_subao( Fipar * fpara, Survpar * sp, int *vec, int dim )
      {
//         int * vec, i, j;
         int  i, j;

        fpara->submat    = dmat ( fpara->subn, fpara->subn);
        fpara->submat_inv= dmat ( fpara->subn, fpara->subn);


//         vec  = ivec(fpara->subn);

         for(i=0; i<CosP_BAO_s; i++)
           vec[i]= i;

         for(i=0; i< sp->zbin; i++) 
           for(j=0; j< PN_BAO_s; j++)
             vec[ CosP_BAO_s + i*PN_BAO_s +j ]=  CosP_BAO + i*PN_BAO + j ;

         vec[fpara->subn -1 ] = dim-1;

         return;
      }


    void baods2de( Fipar * fpara, Survpar * sp )
      {
        int i, j, m, n;
        double ** mat, **ddsdde;

        mat = dmat(fpara->convn, fpara->convn);
        fpara->conv_inv  = dmat (fpara->convn, fpara->convn);

        ddsdde = dmat( fpara->subn, fpara->convn);


// d_{sub p}/ d_{DE}

        for(i=0; i< fpara->convn; i++)  {

          for(m=0 ; m< sp->zbin; m++)     {
            ddsdde[CosP_BAO_s + m*PN_BAO_s ][i] = Da_derp(fpara->cp, sp->z[m], i );
            ddsdde[CosP_BAO_s + m*PN_BAO_s +1 ][i] = H_derp(fpara->cp, sp->z[m], i );
            }

          ddsdde[ fpara->subn-1 ][i] = Da_derp(fpara->cp, 1090.68, i );
          ddsdde[0][i] = 0;
          }
        ddsdde[0][0] = 1;




        for(i=0; i<fpara->convn; i++)
          for(j=0; j<fpara->convn; j++)  
            for(m=0; m<fpara->subn; m++)
              for(n=0; n<fpara->subn; n++)
                mat[i][j] = mat[i][j] + fpara->submat[m][n] * ddsdde[m][i]* 
                                              ddsdde[n][j];


        
        myInv( mat, fpara->conv_inv, fpara->convn );


        freedmat(mat, fpara->convn, fpara->convn);
        freedmat(ddsdde, fpara->subn, fpara->convn);

        //mat = dmat(fpara->convn, fpara->convn);
        //ddsdde = dmat( fpara->subn, fpara->convn);
        return;
      }







