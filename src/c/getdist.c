  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>

  #include <iniparser.h>

  #include "const.h"
  #include "varb.h"
  #include "mymath.h"
  #include "myerr.h"
  #include "matrix.h"
  #include "init.h"
  #include "power.h"
  #include "cospara.h"



  int main( int argc, char *argv[] ) {

//reading parameter .ini file
      int debug= 10, i, j, k, m, n;
      char ini_name[50];

      sscanf( argv[1], "%s", ini_name);


//--------------------------------------------------
  // Open ini parameter file.
      printf("Opening File %s.\n", ini_name);

      dictionary * dict;
      dict = iniparser_load(ini_name);
//--------------------------------------------------
  // ini file ends.
      char *chain_prefix, *output_prefix, fname[100], initemp[100], *sprefix;
      int nchain, npara, chainsize, totpoint, *point;
      double ***chainpoint, **pboundary;

      chain_prefix  = iniparser_getstring( dict, "Chain:chain_file_prefix", NULL);
      output_prefix = iniparser_getstring( dict, "Output:output_prefix", NULL);
      sprefix       = iniparser_getstring( dict, "Output:short_prefix", NULL);

      nchain        = iniparser_getint( dict, "Chain:chain_number", 4);
      npara         = iniparser_getint( dict, "Chain:parameter_number", 5);
      chainsize     = iniparser_getint( dict, "Chain:chain_size", 7);


      point = ivec(nchain);

      FILE *fp_chain[nchain], *fp_output;

      for(i=0; i<nchain; i++) {
        sprintf(fname, "%s_chain_%d.dat", chain_prefix, i+1 );
        printf("%s\n", fname);
        fp_chain[i] = fopen(fname, "r");
        point[i]= countFILEline(fp_chain[i]);
        printf("%dth chain: %d lines\n", i+1, point[i]);
        }



      double maxpoint=0; 
      for(i=0; i<nchain; i++)
        maxpoint = max(point[i], maxpoint);
      printf("maxpoint=%d\n", (int)maxpoint);

      chainpoint = dmat3(nchain, (int)maxpoint, chainsize);
      pboundary   = dmat(npara, 2);


      char term;
      for(i=0; i<nchain; i++)
        for(j=0; j<point[i]; j++) {

          for(k=0; k<chainsize-1; k++)
            fscanf(fp_chain[i], "%lg ",  &chainpoint[i][j][k] );

          fscanf(fp_chain[i], "%lg \n", &chainpoint[i][j][chainsize-1]);

          }


      for(i=0; i<npara; i++)
        for(j=0; j<2; j++) {
          sprintf(initemp, "Plotting:parameter_boundary[%d][%d]", i+1, j+1);
          pboundary[i][j]= iniparser_getdouble( dict, initemp, 1.);
          }


  //-------------------------------------------------------------------------
     // End of Initialization.

  //------------------------------------------//
  //-------------- Begin.---------------------//
  //------------------------------------------//
      int sampledist, chi2dist, pstep, *pstepi, index;
       sampledist = iniparser_getboolean( dict, "Plotting:get_sample_distribution", INIFALSE);
       chi2dist   = iniparser_getboolean( dict, "Plotting:get_chi2_distribution", INITRUE);
      
      pstep  = 1000;
      pstep= iniparser_getint( dict, "Plotting:plot_step_number", 1000);

      pstepi   = ivec(npara);
      ivecvalue(pstepi, npara, pstep);

      double para, p1, p2, ***likelihood, *dp;
      dp         = dvec(npara);
      likelihood = dmat3(npara, pstep, 2);


  //------------------------------------------------------------------------
      if(sampledist== TRUE)  {
      //refining
        double refine=iniparser_getdouble( dict, "Plotting:refinement_fact", 0.2);
        int burnin = iniparser_getint( dict, "Plotting:burn-in", 2000);

      //------------------------------------------------------------------------
      // 1-D plotting.
        printf("Doing Sample distribution\n");

        sprintf(fname, "%s_getdist.dat", output_prefix);
        fp_output=fopen(fname, "w");  

        int **count, *countot, npoint, pposit;

        countot = ivec(npara);
        count    = imat(npara, pstep);

        for(i=0; i<npara; i++)  {
          dp[i] = (pboundary[i][1]- pboundary[i][0])/ (double)pstepi[i];

          for(m=0; m<nchain; m++)   {

            npoint = (point[m]-burnin)*refine;
            //printf("burin=%d, refine=%lg, npoint=%d\n", burnin, refine, npoint);
            for(n=0; n<npoint; n++)  {
              pposit = burnin+ n*(1./refine);
              index = (int) ( (chainpoint[m][pposit][i]- pboundary[i][0])/dp[i] );

              //printf("%d\n", index);
              if(index>=0 && index<pstepi[i])
                count[i][index]++;
                //count[i][index] = count[i][index]+ chainpoint[m][pposit][chainsize-1];
              
              }
            }


          for(j=0; j<pstepi[i]; j++)   {
            p1 = pboundary[i][0] + (double)j*dp[i];
            p2 = pboundary[i][0] + (double)(j+1)*dp[i];
            para = (p1+p2)/2.;


            countot[i] += count[i][j];
            likelihood[i][j][0] = para;
            likelihood[i][j][1] = count[i][j];

            }

          for(j=0; j<pstepi[i]; j++)
            likelihood[i][j][1] = likelihood[i][j][1];///(double)countot[i];

          printf("totpoint[%d]=%d\n", i, countot[i]);
          }

        for(i=0; i<pstep; i++)  {
          for(j=0; j<npara; j++)  
            fprintf(fp_output, "%lg  %lg  ", likelihood[j][i][0], likelihood[j][i][1] );
            
          fprintf(fp_output, "\n");
          }
  
        free(countot);  
        freemat((void **)count, npara, pstep); 
        //count    = imat(npara, pstep);
        }

//--------------------------------------------------------------------
    // getting chi2/FoM distribution: 
    // maximalize the chi2/FoM of other parameters
      if( chi2dist == TRUE) {
        printf("Doing Chi-squar distribution\n");

        sprintf(fname, "%s_fomdist.dat", output_prefix);
        fp_output=fopen(fname, "w"); 

        double **best, fom;
        best = dmat(npara, pstep);
        dmatvalue(best, npara, pstep, 0);

    //----------------------------------------------
      // 1D plotting
        for(i=0; i<npara; i++)  {
          dp[i] = (pboundary[i][1]- pboundary[i][0])/ (double)pstepi[i];
         

          for(m=0; m<nchain; m++)   
            for(n=0; n<point[m]; n++)  {
              index = (int) ( (chainpoint[m][n][i]- pboundary[i][0])/dp[i] );
              fom = exp(chainpoint[m][n][chainsize-2]);

              if(index>=0 && index<pstepi[i] && fom >= best[i][index] )
                best[i][index] = fom;
              }


          for(j=0; j<pstepi[i]; j++)   {
            p1 = pboundary[i][0] + (double)j*dp[i];
            p2 = pboundary[i][0] + (double)(j+1)*dp[i];
            para = (p1+p2)/2.;


            likelihood[i][j][0] = para;
            likelihood[i][j][1] = best[i][j];

            }

          }

        for(i=0; i<pstep; i++)  {
          for(j=0; j<npara; j++)  
            fprintf(fp_output, "%lg  %lg  ", likelihood[j][i][0], likelihood[j][i][1] );
            
          fprintf(fp_output, "\n");
          }

        freedmat(best, npara, pstep);
        }



//-------------------------------------------------
//----------------------------------------------------------
  // Generating supermongo file.

    //----------------------------------------------------------------
        sprintf(fname, "%s_plot1d.sm", output_prefix);
        fp_output=fopen(fname, "w"); 

        fprintf(fp_output, "device postencap %s_1d.eps\n ", sprefix);
        fprintf(fp_output, "\n\ndefine TeX_strings  1\n" );

        fprintf(fp_output, "\n\ndata %s_fomdist.dat\n", sprefix);
        fprintf(fp_output, "#data %s_getdist.dat\n", sprefix);

        fprintf(fp_output, "\n\nlweight 5\nctype black\n" );



        fprintf(fp_output, "\n\nread {  " );
        for(i=0; i<npara; i++)
          fprintf(fp_output, "p%d  %d  l%d  %d  ",i+1, 2*i+1, i+1, 2*i+2 );

        fprintf(fp_output, "}\n" );


        int r1;  int r2 ;

        r1 = 2; 
        if( ( npara - r1*(npara/r1) )<0.2)
          r2 = (int)(npara/r1);
        else
          r2 = (int)(npara/r1)+1;

        printf("(%d, %d)\n", (int)r1, (int)r2 );


        for(i=0; i<npara; i++) {
          fprintf(fp_output, "\n\nwindow %d  %d  %d  %d \n",  (int)r1, (int)r2 , 
                     (int)(i/r2+1), (int)(r2- i+r2*(int)(i/r2)) );

          fprintf(fp_output, "limits p%d l%d\n", i+1, i+1 );
          fprintf(fp_output, "\nbox\nxlabel p%d\nylabel FoM\n", i+1 );
          //fprintf(fp_output, "box\nxlabel p%d\nylabel prob\n", i+1 );

          fprintf(fp_output, "\nconnect p%d  l%d\n",  i+1, i+1);

          }



       // fprintf(fp_output, "\n" );
       // fprintf(fp_output, "\n" );




//----------------------------------------------------------
    // Free all.
      fcloseall();


      free(point);   free(pstepi); 
      freedmat3(likelihood, npara, pstep, 2);   
      freedmat3(chainpoint, nchain, (int)maxpoint, chainsize);   
      freedmat(pboundary, npara, 2);     
      free(dp);    

      //likelihood = dmat3(npara, pstep, 2);
      //chainpoint = dmat3(nchain, (int)maxpoint, chainsize);
      //pboundary   = dmat(npara, 2);

      return;
    }



