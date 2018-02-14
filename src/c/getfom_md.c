  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>

  #include <gsl/gsl_integration.h>

  #include "const.h"
  #include "varb.h"
  #include "mymath.h"
  #include "myerr.h"
  #include "matrix.h"
  #include "init.h"
  #include "readfile.h"
  #include "power.h"
  #include "cospara.h"
  #include "myinterpolate.h"

  #include "bao.h"
  #include "fom.h"

  #include "random.h"
  #include "optstruct.h"
  #include "sampler.h"
  #include "getfom.h"



  #define POWER_LEN 1000

//----------------------------------------------------------------------------------
//
  int fompower_init(char **name, double ***power, int n) {

      int i, j, line;

      FILE *fp;



      for(i=0; i<n; i++) {

        fp = fopen(name[i], "r");
        line = countFILEline(fp);

        for(j=0; j<line-2; j++)
          fscanf(fp, "%lg  %lg\n",&power[j][0][i], &power[j][1][i]);

        fclose(fp);
        }


//      Pk = dmat3(Powline, 2, pn);


      return TRUE;
    }






//----------------------------------------------------------------------------------
// return minus of the logarithm of Figure-of-Merit:  -Ln(FoM)
// Analogy to chi^2, since we use the general sampler.

  double get_fom( void *inform, void *arg, int *fail  )  {

      fom_info *info = (fom_info *) inform;

      double *para = (double *) arg;
      double fom;
      int i, j, cb_dim;
  
//--------------------------------------------------
     // if(info->surv->debug>DEBUGV) {
     //   printf("Entered get_fom()\n");

     //  // for(j=0; j<5; j++)
     //  // printf("para[%d]=%lg\n", j, para[j]);
     //   }

//--------------------------------------------------
  // copy parameter to survey_config & Survpar
      para2sconfig(para, info->s, info->cp, info->sp);

//--------------------------------------------------
  // get the galaxies # distribution and clustering bias, given restrictions.
      get_galaxies(info->cp, info->sp,  info->s, info->surv, fail);


      sconfig2para(info->s, para );

//--------------------------------------------------
  // if this configuration fails, return.
     // printf("fail=%d\n", *fail);
      if(*fail==TRUE || *fail==MYINFINITY)  {

        free(info->sp->z );
        free(info->sp->b );
        free(info->sp->z1);
        free(info->sp->z2);
        free(info->sp->n );
        free(info->sp->beta);
        free(info->sp->v);

        return -1e10; 
        }

//--------------------------------------------------
  // get the Fisher matrix of BAO.
      baofisher(info->fpara, info->sp, info->psarr);


      cb_dim = info->fpara->n +1;
      info->fmat_cb     = dmat(cb_dim,  cb_dim);
      info->fmat_cb_inv = dmat(cb_dim,  cb_dim);

    //  printf("cb_dim=%d\n", cb_dim);
  
// Adding fisher
  
      dmatcp(info->fpara->fmat, info->fmat_cb, cb_dim-1, cb_dim-1);
  
      //printf("BAO Fisher\n");
      //mat_sqrtdiag_printf( info->fpara->fmat,  cb_dim-1);

      for(i=0; i<6; i++)
        for(j=0; j<6; j++)
          info->fmat_cb[i][j] = info->fpara->fmat[i][j] + info->fmat_cmb[i][j];

       
      for(i=0; i<6; i++) {
        info->fmat_cb[i][cb_dim-1] =  info->fmat_cmb[i][6];
        info->fmat_cb[cb_dim-1][i] = info->fmat_cb[i][cb_dim-1] ;
        }

      info->fmat_cb[cb_dim-1][cb_dim-1] =  info->fmat_cmb[6][6];

      myInv(info->fmat_cb, info->fmat_cb_inv, cb_dim);


      //printf("BAO Fisher\n");
      //mat_sqrtdiag_printf( info->fpara->fmat,  cb_dim-1);

      //printf("\nCOMB Fisher\n");
      //mat_sqrtdiag_printf( info->fmat_cb,  cb_dim);

      //printf("\nCMB Fisher\n");
      //mat_sqrtdiag_printf(info->fmat_cmb,  7);

//---------------------------------------------
  // selsect the sub-matrix.

      int * subao = ivec(info->fpara->subn);
      select_subao( info->fpara, info->sp, subao, cb_dim);
  
      mysub(info->fmat_cb_inv, info->fpara->submat_inv, cb_dim,
                                      info->fpara->subn, subao);
  
      myInv(info->fpara->submat_inv, info->fpara->submat, info->fpara->subn);
 

//---------------------------------------------
  // convert to dark energy.
      baods2de( info->fpara, info->sp);

//---------------------------------------------
  // get the FoM
      int  *deindex;
      deindex = ivec(2);
      deindex[0] = 2;       deindex[1] = 3;

      fom = fomt(info->fpara->conv_inv, deindex,  0);


//-----------------------------------------------
// free all.

      free(subao);    free(deindex);
      freemat((void **)info->fmat_cb, cb_dim,  cb_dim);
      freemat((void **)info->fmat_cb_inv, cb_dim,  cb_dim);
    
      freemat( (void **)info->fpara->fmat, info->fpara->n, info->fpara->n);
      freemat( (void **)info->fpara->fmat_inv, info->fpara->n, info->fpara->n);
      //fpara->fmat    =dmat(dsn, dsn);
      //fpara->fmat_inv=dmat(dsn, dsn);

      free(info->sp->z );
      free(info->sp->b );
      free(info->sp->z1);
      free(info->sp->z2);
      free(info->sp->n );
      free(info->sp->beta);
      free(info->sp->v);
      //-----------------------------------
      //sp->z    = dvec(sp->zbin);
      //sp->b    = dvec(sp->zbin);
      //sp->z1   = dvec(sp->zbin);
      //sp->z2   = dvec(sp->zbin);
      //sp->n    = dvec(sp->zbin);
      //sp->beta = dvec(sp->zbin);
      //sp->v    = dvec(sp->zbin);

      freemat( (void **)info->fpara->submat,  info->fpara->subn, info->fpara->subn);
      freemat( (void **)info->fpara->submat_inv, info->fpara->subn, info->fpara->subn);
      freemat( (void **)info->fpara->conv_inv, info->fpara->convn, info->fpara->convn);
      //-----------------------------------
      //fpara->submat    = dmat ( fpara->subn, fpara->subn);
      //fpara->submat_inv= dmat ( fpara->subn, fpara->subn);
      //fpara->conv_inv  = dmat (fpara->convn, fpara->convn);


      //printf("\n!!!!!!!!!!!!!!!\nfom=%lg\n!!!!!!!!!!!!!!!\n", fom);

      return -log(fom);
//      return (fom);
    }



//--------------------------------------------------------------------------------
//----------------------------------------------------------------------//
//                !!!!!!!!    Galaxies Number Count  !!!!!!!            //
//----------------------------------------------------------------------//
// Quoted from Parkinson's code:                                        //
// "A proper photon counts model and wavelength-dependent sky background//
// is now used to produce the flux and magnitude limits for the line    //
// and continuum surveys."                                              //
//----------------------------------------------------------------------//

  void eff_init(survey_setting *surv) {

      if(surv->LAMOST==TRUE ) {
        int n=12;
        double *eff, *lamb;
        eff  = dvec(n);
        lamb = dvec(n);


        lamb[0] = 370; lamb[1] = 400; lamb[2] = 450; lamb[3] = 500;
        lamb[4] = 550; lamb[5] = 600; lamb[6] = 650; lamb[7] = 700;
        lamb[8] = 750; lamb[9] = 800; lamb[10] = 850; lamb[11] = 900;

        eff[0] = 1.82; eff[1] = 4.54; eff[2] = 9.21; eff[3] = 12.14;
        eff[4] = 10.39; eff[5] = 8.82; eff[6] = 10.34; eff[7] = 11.98;
        eff[8] = 13.14; eff[9] = 13.0; eff[10] = 11.03; eff[11] = 8.48;

        myinterp_init(surv->eff, lamb, eff, n );

        free(eff);  free(lamb);

        }

      return;
    }

  double eff_func(survey_setting *surv, double lam) {

      double effunc;

      if(surv->LAMOST==TRUE ) 
        effunc = myinterp(surv->eff, lam)/100.; //0.1;//myinterp(surv->eff, lam);
      else
        effunc = 0.1;

      //printf("effunc[%lg]=%lg\n", lam, effunc);

      return effunc;
    }


  double tempspectro(double lam) {

       double fnu;

       if(lam<=400)
         fnu = 50;
       else if (lam>=800)
         fnu = 250;
       else
         fnu = 50 + (lam-400)*200./400.;


       return fnu;

    }


  void  getfluxlim(  Cospar *cp, Survpar * sp, survey_config *s, survey_setting *surv, 
       double z, double sncont, double snline, double timeppoint, double *maglim, double *linelim )      {

      //if(surv->debug>= 2*DEBUGV)
      //  printf("Entered getfluxlim()\n");

      int i, j;
      double fibarea, mirarea, apert,
             lam0, lamtmp, lamnm, dlam, rn, npix, rdark,
             fzero, nodshuff, maglimmax, teff, 
             eff, eph, absky, fnu, flam, fphsky, rsky, vback,
             nobj_line, robj_line, nobj_cont, robj_cont, fphobj_line,
             fphobj_cont, fact, fnu1, fnu2, lam1, lam2;

      if(surv->LAMOST==TRUE)  {
        //fibarea = PI*squar(1./2.);  // arcsec^2
        //mirarea = PI*squar(8./2.);  // m^2
        //apert = 0.7;

        fibarea = PI*squar(3.3/2.);  // arcsec^2
        mirarea = PI*squar(4./2.);  // m^2
        apert = 0.7;
        }
      else
        myerr("getfluxlim() error.", FALSE);

      dlam = 1.;          // Width of SRE (nm)
      rn = 0.;            // Read-noise (elec)
      npix = 1.;          // Number of pixels extracted
      rdark = 0.;         // Dark current (elec/s/pix)

      fzero = 3631.;
      nodshuff = 1.;
      maglimmax = 24.;
      teff = timeppoint/nodshuff*60;  // sec

//      printf("getfluxlim(): aaaaa\n");

      lam0 = 372.7;    lamtmp = lam0;
    //Wavelength of observation
      lamnm = lam0*(1.+z);
     // throughput at wavelength
      eff = eff_func(surv, lamnm)*0.1;

    //Energy per photon
      eph = 198.64462/lamnm;                  // 10^-18 J

    // AB sky background at this wavelength from model
      //absky = 22.8 - ((lamnm-400.)/150.);
      absky = 20.5 - ((lamnm-400.)/150.);

    // Calculate sky flux
      fnu = fzero*pow(10., (6.-0.4*absky) );      // 10^-32 W/m^2/Hz/arc^2
      flam = (0.3*fnu)/squar(lamnm/100.);     // 10^-18 W/m^2/nm/arc^2
      fphsky = flam/eph;                      // ph/s/m^2/nm/arc^2
      rsky = fphsky*eff*mirarea*dlam*fibarea; // ph/s
    // Noise term from backgrounds
      vback = nodshuff*rsky*teff + npix*squar(rn) + rdark*teff*npix;

    // Solve quadratic equation to determine number of object electrons
      nobj_line = 0.5*squar(snline)*( 1.+sqrt(1.+((4.*vback)/squar(snline) )) );
      robj_line = nobj_line/teff;                         // ph/s

      nobj_cont = 0.5*squar(sncont)*( 1.+sqrt(1.+((4.*vback)/squar(sncont) )) );
      robj_cont = nobj_cont/teff;                         // ph/s

    // For line only.
      fphobj_line = robj_line/(eff*mirarea);          // ph/s/m^2
      *linelim = (fphobj_line*eph)/apert;            // 10^-18 W/m^2


    // for continum.

     fphobj_cont = robj_cont/(eff*mirarea*dlam);     // ph/s/m^2/nm
     flam = eph*fphobj_cont;                             // 10^-18 W/m^2/nm
     // The flux at lamnm = lam0*(1+z)
     fnu = (flam/0.3)*(squar(lamnm/100.));         // 10^-32 W/m^2/Hz
    // The flux at lamobs using the template
     lam1 = lam0;
     lam2 = lamtmp/(1.+z);
     //interp(lamarr,fnuarr,narr, lam1, fnu1);
     //interp(lamarr,fnuarr,narr, lam2, fnu2);
   //-------------------------------
     fnu1 = tempspectro(lam1 );
     fnu2 = tempspectro(lam2 );

     fact = fnu2/fnu1;
   //-------------------------------
     fnu = fnu*fact;
     *maglim = 15.-2.5*log10(fnu/(apert*fzero));
     //*maglim  = maglimmax;


      //printf("getfluxlim(): z=%lg, eph=%lg, tpoint=%lg, robj=%lg,nobj=%lg\n", z, eph, timeppoint, robj,nobj);
      //printf("getfluxlim(): sn=%lg-%lg,  maglim=%lg, linelim=%lg\n", snline, sncont, *maglim, *linelim);

      return;
    }


  double eqwdist(double m, double z, void *param )  {

      LFpar *lf = (LFpar *) param;

      double dist, xx, wmean, wsig, lam0, lamnm, kc,
             mstarev, app, fnu, flam, wmin;

      wmean = 80.;
      wsig = 40.;
    // Rest wavelength of OII (nm)
      lam0 = 372.7;
    // Observed wavelength
      lamnm = lam0*(1.+z); // nm
    // Absolute magnitude corresponding to this luminosity point
//      mstarev = lf->mstar - 2.5*lf->beta*log10((1.+z)/(1.+lf->z0));
//      abs = mstarev - 2.5*log10(uu)      // rest-frame B-band Vega mag
//      abs = abs - 0.098                  // rest-frame AB mag at B
    // [OII] lies in rest-frame U-band, assume average (U-B)_AB = 0.7
    // This number comes from DEEP plot of [OII] against (U-B)
      m = m + 0.7;
    // K-correction is simple cosmological band-stretching
    // Galaxies get brighter in f_nu as photons bunch together in frequency
      kc = -2.5*log10(1.+z);
    // Apparent magnitude at this redshift

      app = absm2appm(lf->cp, m, z) + kc;          // observed AB mag at about I

    // Fluxes corresponding to this magnitude
      fnu = 3631.*pow(10., (6.-0.4*app) );    // fnu in 10^-32 W m^-2 Hz^-1
      flam = (0.3*fnu)/squar(lamnm/100.); // flam in 10^-18 W m^-2 nm^-1
    // Minimum equivalent width observable (A)
      wmin = 10.*(lf->linelim/flam);
    // Fraction of galaxies possessing higher observed equivalent widths
      xx = (wmin-wmean)/(wsig*sqrt(2.));


      dist = 0.5*(1.- erf(xx));

      //printf("m=%lg, z=%lg, app=%lg\n", m, z, app);
      //printf("wmin=%lf, linelim=%lg, fnu=%lg, fact=%lg\n", wmin, lf->linelim, fnu, dist);
      //printf("wmin=%lf, linelim=%lg, flam=%lg, fact=%lg\n", wmin, lf->linelim, flam, dist);

      return dist;
    }


  double maglim_fiber(double mag, void *para)  {

      LFpar *lf = (LFpar *)para;

      double ang;

      ang = angden(lf->cp, lf, lf->op[1], lf->op[2], mag);

      return ang-lf->op[0];
    }


  void ngal( Cospar *cp, Survpar * sp, survey_config *s, 
            survey_setting *surv, double timeppoint, double *ngobs )           {

    //-------------------------------------------------------------//
      if(surv->debug>20*DEBUGV) 
        printf("ngal()\n");
    //-------------------------------------------------------------//

      int i, j;

      double nfib, sn, snline, sncont, mend, maglim, linelim, fibmaglim,
             zcen;

    // initial cp_lfc
      Cospar cp_lf;

        cp_lf.omem = 0.3;  cp_lf.h0 = 0.7;  cp_lf.omek=0;
        cp_lf.omeb =  cp_lf.omem*0.15;  cp_lf.omex = 1.- cp_lf.omem ;
        cp_lf.tau = 0.092;

        cp_lf.w0= -1.;   cp_lf.w1= 0.;    cp_lf.ns= 0.958;
        cp_lf.AScmb= 2.3e-9;    cp_lf.sig8= 0.9;

        cp_lf.omemh2=  cp_lf.omem * squar( cp_lf.h0 );
        cp_lf.omebh2=  cp_lf.omeb * squar( cp_lf.h0 );
        cp_lf.omech2= cp_lf.omemh2 - cp_lf.omebh2;



    //-------------------------------------------------------------//

      LFpar glf;
        glf.cp = cp;
    //-------------------------------------------------------------//
    // Luminosity function for z=1 blue galaxies, B-band Vega mags //
    // from Willmer et al (DEEP2) Table 4, z=1.1                   //
    // "minimal" weights to be conservative                        //
    // Assumes Omega_m = 0.3, Omega_L = 0.7, H_0 = 70              //
    // Assume L* evolves as (1+z)^beta                             //
    //-------------------------------------------------------------//
     if(surv->line==TRUE)  {
        glf.phistar = 2.08e-3/cub(cp_lf.h0);  // Normalization in (h/Mpc)^3
        glf.mstar   = -21.38;
        glf.alpha   = -1.3 ;
        glf.z0      = 1.1;

        glf.evolution = 2;
        glf.beta = 3.0;   // L^star ~ (1+z)^beta

        glf.kcorrection =1;  //kc = (z-0.5)-0.9; 
        }
      else  {
        glf.phistar = 1.07e-3/cub(cp_lf.h0);  // Normalization in (h/Mpc)^3
        glf.mstar   = -21.11;
        glf.alpha   = -0.5 ;  
        glf.z0      = 0.9;

        glf.evolution = 2;
        glf.beta = 0.;        // L^star ~ (1+z)^beta

        glf.kcorrection =2;  //kc = 2*(z-0.5)-0.9; 
        }



    //------------------------------------------------------------------//
    // Correct the luminosity function parameters for changed cosmology.//
    //------------------------------------------------------------------//
      double fact; 
      fact = rcom(&cp_lf, glf.z0)/rcom(cp, glf.z0);

    //phistar corrected by volume ratio which is distance ratio squared
      glf.phistar = glf.phistar * squar(fact);
    //lstar corrected by luminosity ratio which is distance ratio squared
      glf.mstar = glf.mstar + 5*log10(fact);

    //------------------------------------------------------
    // signal-to-noise ratio: S/N
      snline = 7; //line
      sncont = 0.1; // continuum
    //------------------------------------------------------
      mend   = -26.;   //magnitude bright end
      maglim = 0;

    //------------------------------------------------------
    // total fibers
      if(s->npt==0)
        nfib = surv->nfib;
      else
        nfib = surv->nfib*s->npt;
    //------------------------------------------------------
    // End of Initialization.


  //------------------------------------------------------
  // Begin.
  //-------------------------------------------------------------------//
  // Calculate the magnitude limit corresponding to the fibre density. //
  //-------------------------------------------------------------------//
      double zminuse, zmaxuse, volfact, fibdens, angdens, maglimmax ;

      zmaxuse = s->zmax;
      zminuse = s->zmax - 0.2;
      if(zminuse<s->zmin) 
        zminuse = s->zmin;

      maglim = 21.;

      glf.dofact  = FALSE;

      volfact = vol(cp, zminuse, zmaxuse)/vol(cp, s->zmin, s->zmax) ;

      fibdens = volfact*nfib/surv->fov;   //fiber density between end of bins.
      angdens = 0;
      maglimmax = 28.;

      if(surv->debug>=DEBUGV)  {
        printf("Calculating survey magnitude limit.\n");
        printf("zminuse=%lg\n", zminuse);
        printf("zmaxuse=%lg\n", zmaxuse);
        printf("volfact=%lg\n", volfact);
        printf("nfib=%lg\n", nfib);
        printf("s->area=%lg\n", s->area);

        }

      //for ( ; (angdens<fibdens) && (maglim<maglimmax); ) {

      //  maglim = maglim + 0.02;
      //  glf.maglim  = maglim;

 //   //    printf("maglim=%lg, fibdens=%lg\n", maglim, fibdens);
      //  if(dndz(cp, &glf, zminuse, maglim)==0. || dndz(cp, &glf, zmaxuse, maglim)==0.) 
      //    continue;

      //  angdens = angden(cp, &glf, zminuse, zmaxuse, maglim);
//    //    printf("angdens=%lg\n", angdens);
      //  }

      //fibmaglim = maglim;

//-------------------------------------------------
    // get the magnitude limit corresponding to fiber density

      gsl_function F;

        F.function = & maglim_fiber;
        F.params = &glf;

        glf.op[0] = fibdens;
        glf.op[1] = zminuse;
        glf.op[2] = zmaxuse;

      fibmaglim =  myroot( &F, 18., maglimmax, 1e-4) ;


      if(surv->debug>=2*DEBUGV) 
        printf("fibmaglim=%lg\n", fibmaglim);

  //------------------------------------------------------------------//
  // Convert exposure time to flux or magnitude limit, using a proper //
  // photon counts model including the sky background varying with    //
  // wavelength.                                                      //
  //------------------------------------------------------------------//
      zcen = (zminuse+zmaxuse)/2.;
      getfluxlim(cp, sp, s, surv, zcen, sncont, snline, timeppoint, &maglim, &linelim);


      if(surv->debug>= 2*DEBUGV)
        printf("maglim=%lg, fiber maglim=%lg\n", maglim, fibmaglim);

  //------------------------------------------------------------------//
      maglim = min(maglim, fibmaglim);

      if(surv->debug>= 2*DEBUGV)
        printf("maglim: %lg\n", maglim);

  //-----------------------------------------------------------
    // return number density.
      glf.dofact  = TRUE;//TRUE;
      glf.fact    = &eqwdist;
      glf.factpar = &glf;
      glf.maglim  = maglim;
      glf.linelim = linelim;


      //printf("\n!!!!!\nEWD=%lg \n", eqwdist(maglim, sp->z[0], &glf) ) ;

      angdens= angden(cp, &glf, zminuse, zmaxuse, maglim);
      sp->n[sp->zbin-1] = angdens/vol(cp, zminuse, zmaxuse);

      ngobs[sp->zbin-1] = sp->n[sp->zbin-1]* sp->v[sp->zbin-1]*surv->fov/sp->area;

      for(i=0; i<sp->zbin-1; i++) { 
        sp->n[i]= sp->n[sp->zbin-1];
        ngobs[i] = sp->n[i]* sp->v[i]*surv->fov/sp->area;
        }

     // printf("n=%lg\n", sp->n[sp->zbin-1]);

//---------------------------------------------------------
      //for(i=0; i<sp->zbin; i++) { 
      //  printf("dndz:  z[%d]=%lg \n", i, sp->z[i]);
      //  sp->n[i] = dndz(cp, &glf, sp->z[i], maglim);
      //  ngobs[i] = sp->n[i]*sp->v[i]*surv->fov/s->area;

      //  if(surv->debug>2*DEBUGV)
      //    printf("n@sp->z[%lg]=%lg\n", sp->z[i], sp->n[i]);
      //  }
        

      return;
    }


//--------------------------------------------------
  //Get the galaxies # distribution and clustering bias, given the observing time.
  void get_galaxies( Cospar *cp, Survpar * sp, survey_config *s, 
              survey_setting *surv, int *fail )  {

      int i, j;

      *fail = TRUE;

      if(surv->debug>20*DEBUGV) 
        printf("get_galaxies(): \n");

      if(s->zmin>=s->zmax) {
        if(surv->debug>20*DEBUGV)
          printf("zmin>=zmax \n");
        return;
        }

//--------------------------------------
      s->bin = (int)((s->zmax-s->zmin)/0.2);

      if(s->bin<1) {
        if(surv->debug>20*DEBUGV)
          printf("Too small dz.\n");
        return;
        }

        


//--------------------------------------
      sp->zbin = s->bin;
      sp->area = s->area;


      sp->z    = dvec(sp->zbin);
      sp->b    = dvec(sp->zbin);
      sp->z1   = dvec(sp->zbin);
      sp->z2   = dvec(sp->zbin);
      sp->n    = dvec(sp->zbin);
      sp->beta = dvec(sp->zbin);
      sp->v    = dvec(sp->zbin);

//-------------------------------------
      sp->phot = FALSE;
      sp->sz0  = 1e-3; 

  // redshift, volume, and the bias
      for(i=0; i<sp->zbin; i++)  {

        sp->z1[i] = s->zmin +  i * (s->zmax - s->zmin)/(double)s->bin;
        sp->z2[i] = s->zmin + (i+1) * (s->zmax - s->zmin)/(double)s->bin;
        sp->z[i] = (sp->z1[i] + sp->z2[i]) /2.; ;

        sp->v[i] = vol(cp, sp->z1[i],  sp->z2[i])*sp->area; 
       // printf("v%d=%lg, area=%lg \n ", i, sp->v[i], sp->area);

        if(sp->z[i]<=2.) {
          if(surv->line==TRUE) 
            sp->b[i] = 1.+(1.3 - 1.0)*Growth(cp, 0.55)/Growth( cp, sp->z[i]) ;
          else
            sp->b[i] = 1.+(2. - 1.0)*Growth(cp, 0.55)/Growth( cp, sp->z[i]) ;
          }
        else if(sp->z[i]<=4.9)
          sp->b[i] = 0.53+0.289*squar(1.+sp->z[i]);
        else
          myerr("Bias setting: too high redshift.", FALSE);

        sp->beta[i] = Beta(cp, sp->z[i], sp->b[i] );


        }

//----------------------------------------------------
  // estimate number density

      double timepfov, timeppoint,  // time per fov/pointing
             ntot, *ngobs, ngtot=0;

    //-----------------------------
      if(s->time<=0 || s->area<=0)  {
        dvecvalue(sp->n, sp->zbin, 0);
        s->npt = 0;
        *fail== MYINFINITY;
        return;
        }


      ngobs = dvec(sp->zbin);

    //-----------------------------
      //exposure time per field-of-view
      timepfov = s->time*surv->fov/s->area;
      if(surv->debug>=2*DEBUGV)
        printf("get_galaxies(): fov=%lg,  timepfov=%lg\n", surv->fov, timepfov);

      // too short time
      if(timepfov - surv->overheadtime < surv->minexptime) {
      // decreasing area.
        timepfov = surv->minexptime + surv->overheadtime + 1.e-3;
        s->area = s->time/timepfov*surv->fov;
        }

      if(surv->debug>=2*DEBUGV) {
        printf("get_galaxies():  surv->nptmax=%d\n", surv->nptmax);
        printf("get_galaxies():  s->time=%lg\n", s->time);
        printf("get_galaxies():  s->area=%lg\n", s->area);
        printf("get_galaxies():  s->zmin=%lg\n", s->zmin);
        printf("get_galaxies():  s->zmax=%lg\n", s->zmax);
        }

      for(i=surv->nptmax; i>0; i--)  {  // i, run over the No. of pointing.

        timeppoint = timepfov/(double)i - surv->overheadtime;

        if(timeppoint > surv->maxexptime && surv->debug>=2*DEBUGV) {  //time is too long.
          printf("timeppoint=%lg\n", timeppoint); 
          free(ngobs);
          return; 
          }
        if(timeppoint < surv->minexptime) // too short, go on
          continue;

        // number of pointing
        s->npt = i;

        // get the number density
        //printf("getgalaxies: time/point=%lg,t/fov=%lg,mint=%lg\n", 
        //                   timeppoint, timepfov,surv->minexptime);

        ngal( cp, sp, s, surv, timeppoint, ngobs);


        // the number of galaxies observed per pointing.
        for(j=0; j<sp->zbin; j++)
          ngtot += ngobs[j];

        if(ngtot< s->npt*surv->nfib)   {  // not enough galaxies
          if(s->npt==1)  {
            s->npt==0;   
            break;
            }
          else
            continue;
          }

        break;
        }

       if(surv->debug >=2*DEBUGV)
         printf("npt=%d, ngtot=%lg\n", s->npt, ngtot);

    //------------------------------------------------
      if(s->npt > 0)
        ntot = s->npt*surv->nfib*s->area/surv->fov;
      else
        ntot = ngtot*s->area/surv->fov;
      //less than 1 galaxy
      if(ntot-1<0)  {
        free(ngobs);
        return;
        }


  //-----------------------------------------------
      *fail = FALSE;

      free(ngobs);

      return;
    }
//--------------------------------------------------
  // End of getfom.c






//--------------------------------------------------
//--------------------------------------------------
// not useful
//--------------------------------------------------

  #ifdef _FOM_POWER_INTERP_INIT_
  int fompower_interp_init( char **name,  Interpar ** power, int n )  {

      int i, j, line;

      FILE *fp;

      double *p, *k;

      p= dvec(POWER_LEN);
      k= dvec(POWER_LEN);

      for(i=0; i<n; i++) {

        fp = fopen(name[i], "r");
        line = countFILEline(fp);

        for(j=0; j<line; j++)
          fscanf(fp, "%lg  %lg\n",&k[j], &p[j]);

        myinterp_init( power[i], k, p, line);

        fclose(fp);
        }


      free(p);  free(k);

      return TRUE;
    }

  #endif
