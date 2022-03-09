//------------------------------------------------------------------------------------------------------------------------------------------------
// 1st version of PACER BLINK fast imager : produces a dirty images for now (should work for both MWA and SKA-Low stations)
//------------------------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include <bg_globals.h>
#include "bg_fits.h"
#include <mystring.h>

#include <vector>
using namespace std;

// FFTW :
#include <fftw3.h>

// defines :
#include "pacer_imager_defs.h"

// Pacer imager class :
#include "pacer_imager.h"

string in_basename="1276619416_20200619164456";
string gPostfix="";
string out_fits="1276619416_20200619164456_dirty_image.fits";
// string out_rms_fits="out_rms.fits";

// REAL and IMAG external UV GRID files to make dirty image:
string g_in_fits_file_uv_re;
string g_in_fits_file_uv_im;

// options to enable particular functions:
bool gDoGridding=true;
bool gDo_dirty_image=true;

int gImageSize = 128; // was 512

// weighting :
string gUniformWeighting="N"; // uniform, means that sum of visibilities in a given UV cell is divided by a number of points added

double gFrequencyMHz = 159.375; // default of the station in channel 204 is 204*(400.00/512.00) = 159.375 MHz 
double gFOV_degrees  = 180.00;
double gMinUV        = -1000;

// input files :
string gAntennaPositionsFile;

void usage()
{
   printf("pacer_dirty_image VISIBILITY_FITS_BASENAME\n\n\n");
   
   printf("\t-p POSTFIX : default is not postfix\n");
   printf("\t-g 0 / 1 for GRIDDING ON / OFF [default GRIDDING = %d]\n",gDoGridding);
   printf("\t-r FITS_FILE_RE : specify external UV grid real fits file [requires gridding=0]\n");
   printf("\t-i FITS_FILE_IM : specify external UV grid imag fits file [requires gridding=0]\n");
   printf("\t-f FREQ_MHz     : frequency in MHz [default 159.375 MHz]\n");
   printf("\t-F FoV[deg]     : field of view in degrees [default %.4f degree]\n",gFOV_degrees);
   printf("\t-w WEIGHTING    : change weighting schema N - natural, U - uniform [default %s]\n",gUniformWeighting.c_str());
   printf("\t-m MIN_UV_DISTANCE : minimum UV distance in wavelengths for a baseline to be included [default %.4f]\n",gMinUV);
   printf("\t-a antenna_positions.txt : text file with antenna positions in a format : AntName X[m] Y[m] Z[m]\n");

   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "hp:n:g:r:i:f:F:w:m:a:";
   int opt;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
//      printf("opt = %c (%s)\n",opt,optarg);   
      switch (opt) {
         case 'h':
            // antenna1 = atol( optarg );
            usage();
            break;

         case 'a':
            if( optarg && strlen(optarg) ){
               gAntennaPositionsFile = optarg;
            }
            break;

         case 'f':
            gFrequencyMHz = atof( optarg );
            break;

         case 'F':
            gFOV_degrees = atof( optarg );
            break;

         case 'p':
            gPostfix = optarg;
            break;

         case 'n':
            gImageSize = atol( optarg );
            break;

         case 'g':
            gDoGridding = atol( optarg );
            break;
            
         case 'i':
            g_in_fits_file_uv_im = optarg;
            break;

         case 'm':
            gMinUV = atof( optarg );
            break;

         case 'r' :
            g_in_fits_file_uv_re = optarg;
            break; 

         case 'w' :
            if( optarg ){
               gUniformWeighting = optarg;
            }
            break; 

         default:   
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }
}

void print_parameters()
{
    printf("############################################################################################\n");
    printf("PARAMETERS :\n");
    printf("############################################################################################\n");
    printf("Base name for FITS  = %s\n",in_basename.c_str());
    printf("out_fits     = %s\n",out_fits.c_str());
    printf("Postfix      = %s\n",gPostfix.c_str());
    printf("Image size   = (%d x %d)\n",gImageSize,gImageSize);
    printf("Do gridding  = %d\n",gDoGridding);
    if( !gDoGridding ){
       printf("External REAL/IMAG FITS FILES : %s / %s\n",g_in_fits_file_uv_re.c_str(),g_in_fits_file_uv_im.c_str());
    }
    printf("Frequency           = %.4f [MHz]\n",gFrequencyMHz);
    printf("Field of view (FoV) = %.4f [deg]\n",gFOV_degrees);
    printf("Weighting           = %s\n",gUniformWeighting.c_str());
    printf("UV range            = %.4f - Infinity\n",gMinUV);   
    printf("Antenna positions file = %s\n",gAntennaPositionsFile.c_str());
    printf("############################################################################################\n");
}

void print_header()
{
   printf("pacer_imager version 0.00 (C++/object), compiled on %s\n",__TIMESTAMP__);  
}

int main(int argc,char* argv[])
{
   print_header();

  // parsing paratemeters with input and output file names:
  if( argc >= 2 ){
     in_basename = argv[1];
  }
  out_fits="out.fits";
  if( argc >= 3 ){
     out_fits = argv[2];
  }  

  // parse and print parameters :
  parse_cmdline( argc , argv );
  print_parameters();
  
  CPacerImager imager;
  imager.m_ImagerParameters.SetParameters( gAntennaPositionsFile.c_str() );
  
  imager.run_imager( in_basename.c_str(),  gPostfix.c_str(),
                     gFrequencyMHz, gImageSize, gFOV_degrees, gMinUV, 
                     gDoGridding, gDo_dirty_image, 
                     gUniformWeighting.c_str(), g_in_fits_file_uv_re.c_str(), g_in_fits_file_uv_im.c_str()
                   );

}

