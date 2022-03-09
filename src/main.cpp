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

   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "hp:n:g:r:i:f:F:w:m:";
   int opt;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
//      printf("opt = %c (%s)\n",opt,optarg);   
      switch (opt) {
         case 'h':
            // antenna1 = atol( optarg );
            usage();
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
    printf("############################################################################################\n");
}

void fft_shift( CBgFits& dirty_image, CBgFits& out_image )
{
   int xSize = dirty_image.GetXSize();
   int ySize = dirty_image.GetYSize();
   CBgFits tmp_image( xSize, ySize );
   
   int center_freq_x = int( xSize/2 );
   int center_freq_y = int( ySize/2 );
   
   int is_odd = 0;
   if ( (xSize%2) == 1 && (ySize%2) == 1 ){
      is_odd = 1;
   }

   
   // X (horizontal FFT shift) :
   for(int y=0;y<ySize;y++){ 
      float* tmp_data = tmp_image.get_line(y);
      float* image_data = dirty_image.get_line(y);
      
      for(int x=0;x<=center_freq_x;x++){ // check <= -> <
         tmp_data[center_freq_x+x] = image_data[x];
      }
      for(int x=(center_freq_x+is_odd);x<xSize;x++){
         tmp_data[x-(center_freq_x+is_odd)] = image_data[x];
      }      
   }

   for(int x=0;x<xSize;x++){ 
      for(int y=0;y<=center_freq_y;y++){ // check <= -> <
         out_image.setXY(x,center_freq_y+y,tmp_image.getXY(x,y));
      }
      for(int y=(center_freq_y+is_odd);y<ySize;y++){
         out_image.setXY( x , y-(center_freq_y+is_odd),tmp_image.getXY(x,y));
      }      
   }
}

// UV data are with DC in the center -> have to be FFTshfted to form input to FFT function :
void fft_unshift( CBgFits& dirty_image, CBgFits& out_image )
{
   int xSize = dirty_image.GetXSize();
   int ySize = dirty_image.GetYSize();
   CBgFits tmp_image( xSize, ySize );
   
   int center_freq_x = int( xSize/2 );
   int center_freq_y = int( ySize/2 );
   
   int is_odd = 0;
   if ( (xSize%2) == 1 && (ySize%2) == 1 ){
      is_odd = 1;
   }
   
   // X (horizontal FFT shift) :
   for(int y=0;y<ySize;y++){ 
      float* tmp_data = tmp_image.get_line(y);
      float* image_data = dirty_image.get_line(y);
      
      for(int x=0;x<center_freq_x;x++){ // check <= -> <
         tmp_data[center_freq_x+x+is_odd] = image_data[x];
      }
      for(int x=center_freq_x;x<xSize;x++){
         tmp_data[x-center_freq_x] = image_data[x];
      }      
   }

   for(int x=0;x<xSize;x++){ 
      for(int y=0;y<center_freq_y;y++){ // check <= -> <
         out_image.setXY( x, center_freq_y+y+is_odd, tmp_image.getXY(x,y));
      }
      for(int y=center_freq_y;y<ySize;y++){
         out_image.setXY(x,y-center_freq_y,tmp_image.getXY(x,y));
      }      
   }
}



// Based on example : https://github.com/AccelerateHS/accelerate-examples/blob/master/examples/fft/src-fftw/FFTW.c
void dirty_image( CBgFits& uv_grid_real_param, CBgFits& uv_grid_imag_param, CBgFits& uv_grid_counter )
{
   CBgFits uv_grid_real( uv_grid_real_param.GetXSize(), uv_grid_real_param.GetYSize() ), uv_grid_imag( uv_grid_imag_param.GetXSize(), uv_grid_imag_param.GetYSize() );
   fft_unshift( uv_grid_real_param, uv_grid_real );
   fft_unshift( uv_grid_imag_param, uv_grid_imag );

   CBgFits out_image_real( uv_grid_real_param.GetXSize(), uv_grid_real_param.GetYSize() ), out_image_imag( uv_grid_imag_param.GetXSize(), uv_grid_imag_param.GetYSize() );
   out_image_real.SetValue( 0.00 );
   out_image_imag.SetValue( 0.00 );

   int width = uv_grid_real.GetXSize();
   int height = uv_grid_real.GetYSize();
   int size = width*height;
   
   // Allocate input and output buffers
   fftw_complex* in_buffer	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
   fftw_complex* out_buffer	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
   
   // 
   float* real_data = uv_grid_real.get_data();
   float* imag_data = uv_grid_imag.get_data();

   // Copy in image data as real values for the transform.
   for(int i = 0; i < size; i++) {
      double re = real_data[i];
      double im = imag_data[i];
   
      in_buffer[i][0]	= re;
      in_buffer[i][1]	= im;
   }

   // should be the same and it is the same !
/*   for(int y=0;y<height;y++){
      for(int x=0;x<width;x++){
         double re = uv_grid_real.getXY(x,y);
         double im = uv_grid_imag.getXY(x,y);
         
         int pos = y*width + x;
         in_buffer[pos][0]   = re;
         in_buffer[pos][1]   = im;
      }
   }*/
   
   // Transform to frequency space.
   // https://www.fftw.org/fftw3_doc/Complex-Multi_002dDimensional-DFTs.html
   fftw_plan pFwd = fftw_plan_dft_2d( width, height, in_buffer, out_buffer, FFTW_FORWARD, FFTW_ESTIMATE); // was FFTW_FORWARD or FFTW_BACKWARD ???
   fftw_execute(pFwd);
   fftw_destroy_plan(pFwd);
   
   
   // copy resulting data from out_buffer to out_image_real :
   // WARNING : this is image in l,m = cos(alpha), cos(beta) coordinates and still needs to go to SKY COORDINATES !!!
   float* out_data_real = out_image_real.get_data();
   float* out_data_imag = out_image_imag.get_data();
//   double fnorm = 1.00/sqrt(size); //normalisation see : /home/msok/Desktop/PAWSEY/PaCER/logbook/20220119_testing_new_versions_dirty_image_polishing.odt
   double fnorm = 1.00/uv_grid_counter.Sum(); // see RTS : /home/msok/mwa_software/RTS_128t/src/newgridder.cu SumVisibilityWeights and gridKernel.c:650 
                                              // also read TMS (Thomson, Moran, Swenson) about this 
   printf("DEBUG : size = %d (%d x %d), fnorm = %e\n",size,width,height,fnorm);
   for(int i = 0; i < size; i++) {
//      out_data[i] = out_buffer[i][0]*out_buffer[i][0] + out_buffer[i][1]*out_buffer[i][1]; // amplitude ?
     out_data_real[i] = out_buffer[i][0]*fnorm; // real 
     out_data_imag[i] = out_buffer[i][1]*fnorm; // imag
   }   

   char outDirtyImageReal[1024],outDirtyImageImag[1024];   
   sprintf(outDirtyImageReal,"dirty_test_real_%dx%d.fits",width,height);
   sprintf(outDirtyImageImag,"dirty_test_imag_%dx%d.fits",width,height);
   
   out_image_real.WriteFits( outDirtyImageReal );
   out_image_imag.WriteFits( outDirtyImageImag );
   
   // calculate and save FFT-shifted image :
   CBgFits out_image_real2( out_image_real.GetXSize(), out_image_real.GetYSize() ), out_image_imag2( out_image_real.GetXSize(), out_image_real.GetYSize() );
   sprintf(outDirtyImageReal,"dirty_test_real_fftshift_%dx%d.fits",width,height);
   sprintf(outDirtyImageImag,"dirty_test_imag_fftshift_%dx%d.fits",width,height);
   fft_shift( out_image_real, out_image_real2 );
   fft_shift( out_image_imag, out_image_imag2 );
   out_image_real2.WriteFits( outDirtyImageReal );
   out_image_imag2.WriteFits( outDirtyImageImag );

   
   // TODO : re-grid to SKY COORDINATES !!!
   // convert cos(alpha) to alpha - see notes !!!
   // how to do it ???
     
}

bool read_corr_matrix( const char* basename, CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, const char* szPostfix )
{
  // creating FITS file names for REAL, IMAG and U,V,W input FITS files :
  string fits_file_real = basename;
  fits_file_real += "_vis_real";
  if( strlen( szPostfix ) ){
     fits_file_real += szPostfix;
  }
  fits_file_real += ".fits";
 
  string fits_file_imag = basename;
  fits_file_imag += "_vis_imag";
  if( strlen( szPostfix ) ){
     fits_file_imag += szPostfix;
  }
  fits_file_imag += ".fits";
 
  string fits_file_u = basename;
  fits_file_u += "_u";
  if( strlen( szPostfix ) ){
     fits_file_u += szPostfix;
  }
  fits_file_u += ".fits";
 
  string fits_file_v = basename;
  fits_file_v += "_v";
  if( strlen( szPostfix ) ){
     fits_file_v += szPostfix;
  }
  fits_file_v += ".fits";
 
  string fits_file_w = basename;
  fits_file_w += "_w";
  if( strlen( szPostfix ) ){
     fits_file_w += szPostfix;
  }
  fits_file_w += ".fits";


  printf("Expecting the following files to exist:\n");
  printf("\t%s\n",fits_file_real.c_str()); 
  printf("\t%s\n",fits_file_imag.c_str()); 
  printf("\t%s\n",fits_file_u.c_str()); 
  printf("\t%s\n",fits_file_v.c_str()); 
  printf("\t%s\n",fits_file_w.c_str()); 
  
  // REAL(VIS)
  printf("Reading fits file %s ...\n",fits_file_real.c_str());
  if( fits_vis_real.ReadFits( fits_file_real.c_str(), 0, 1, 1 ) ){
     printf("ERROR : could not read visibility FITS file %s\n",fits_file_real.c_str());
     return false;
  }else{
     printf("OK : fits file %s read ok\n",fits_file_real.c_str());
  }

  // IMAG(VIS)
  printf("Reading fits file %s ...\n",fits_file_imag.c_str());
  if( fits_vis_imag.ReadFits( fits_file_imag.c_str(), 0, 1, 1 ) ){
     printf("ERROR : could not read visibility FITS file %s\n",fits_file_imag.c_str());
     return false;
  }else{
     printf("OK : fits file %s read ok\n",fits_file_imag.c_str());
  }

  // U : 
  printf("Reading fits file %s ...\n",fits_file_u.c_str());
  if( fits_vis_u.ReadFits( fits_file_u.c_str(), 0, 1, 1 ) ){
     printf("ERROR : could not read U FITS file %s\n",fits_file_u.c_str());
     return false;
  }else{
     printf("OK : fits file %s read ok\n",fits_file_u.c_str());
  }
  
  // V : 
  printf("Reading fits file %s ...\n",fits_file_v.c_str());
  if( fits_vis_v.ReadFits( fits_file_v.c_str(), 0, 1, 1 ) ){
     printf("ERROR : could not read V FITS file %s\n",fits_file_v.c_str());
     return false;
  }else{
     printf("OK : fits file %s read ok\n",fits_file_v.c_str());
  }
  
  // W : 
  printf("Reading fits file %s ...\n",fits_file_w.c_str());
  if( fits_vis_w.ReadFits( fits_file_w.c_str(), 0, 1, 1 ) ){
     printf("ERROR : could not read W FITS file %s\n",fits_file_w.c_str());
     return false;
  }else{
     printf("OK : fits file %s read ok\n",fits_file_w.c_str());
  }


  return true;
}

void gridding( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
               CBgFits& uv_grid_real, CBgFits& uv_grid_imag, CBgFits& uv_grid_counter, 
               double delta_u, double delta_v, double frequency_mhz, int n_pixels,  double min_uv, const char* weighting )
{
  double u_mean, u_rms, u_min, u_max;
  fits_vis_u.GetStat( u_mean, u_rms, u_min, u_max );
  
  // V : 
  double v_mean, v_rms, v_min, v_max;
  fits_vis_v.GetStat( v_mean, v_rms, v_min, v_max );
  
  // W : 
  double w_mean, w_rms, w_min, w_max;
  fits_vis_w.GetStat( w_mean, w_rms, w_min, w_max );

  // Bacause we are also including conjugates at (-u,-v) UV point in gridding u_min = -u_max and v_min = -v_max :
  // was -35 / +35 
  u_min = -u_max;
//  u_max = +35;  
  v_min = -v_max;
//  v_max = +35;


// calculate using CASA formula from image_tile_auto.py :
// synthesized_beam=(lambda_m/max_baseline)*(180.00/math.pi)*60.00 # in arcmin
// lower=synthesized_beam/5.00
// higher=synthesized_beam/3.00
// cellside_float=(lower+higher)*0.5
// NEW :
//   double alpha = 4.00/15.00; // from (1/5+1/3)*0.5 = 4/15
//   double wrong_factor = 1.00; // 4.00 factor for due to U range -35 to +35 instead of ~-17.5/2 to +17.5/2 (factor 4 )
//   double delta_u = wrong_factor*alpha*(u_max-u_min)/(n_pixels); // factor for due to U range -35 to +35 instead of ~-17.5/2 to +17.5/2 (factor 4 )
//   double delta_v = wrong_factor*alpha*(v_max-v_min)/(n_pixels);
//   int freq_channel = 204;
//   double frequency_mhz = freq_channel*(400.00/512.00);
//   if( gFrequencyMHz > 0 ){
//      frequency_mhz = gFrequencyMHz;
//   }
   double frequency_hz = frequency_mhz*1e6;
   double wavelength_m = VEL_LIGHT / frequency_hz;
   // UV pixel size as function FOVtoGridsize in  /home/msok/mwa_software/RTS_128t/src/gridder.c
   // delta_u = (VEL_LIGHT/frequency_Hz)/(gFOV_degrees*M_PI/180.);
   // delta_v = (VEL_LIGHT/frequency_Hz)/(gFOV_degrees*M_PI/180.);
   printf("DEBUG : wavelength = %.4f [m] , frequency = %.4f [MHz]\n",wavelength_m,frequency_mhz);
// OLD :
//  double delta_u = (u_max - u_min) / n_pixels;  
//  double delta_v = (v_max - v_min) / n_pixels;
  
  printf("DEBUG : U limits %.8f - %.8f , delta_u = %.8f\n",u_min, u_max , delta_u );
  printf("DEBUG : V limits %.8f - %.8f , delta_v = %.8f\n", v_min, v_max , delta_v );
  printf("DEBUG : W limits %.8f - %.8f\n", w_min, w_max );

  // is it ok to chose the UV plane center based on this:  
//  double u_center = (u_min + u_max)/2.00;  
//  double v_center = (v_min + v_max)/2.00;  

  // Limits of UVW :
  // double GetStat( double& mean, double& rms, double& minval, double& maxval, 


  // simple gridding :
  uv_grid_real.SetValue( 0.00 );
  uv_grid_imag.SetValue( 0.00 );
  uv_grid_counter.SetValue( 0.00 );
      
  int added=0, high_value=0;
  for( int ant1 = 0 ; ant1 < fits_vis_real.GetXSize(); ant1++ ){
     for( int ant2 = 0 ; ant2 < fits_vis_real.GetXSize(); ant2++ ){
        if( ant1 > ant2 ){ // was ant1 > ant2 
           double re = fits_vis_real.getXY(ant1,ant2);
           double im = fits_vis_imag.getXY(ant1,ant2);
           
           if( !isnan(re) && !isnan(im) ){
              if ( fabs(re) < MAX_VIS && fabs(im) < MAX_VIS ){
                 // TODO convert [m] -> wavelength 
                 double u = fits_vis_u.getXY(ant1,ant2) / wavelength_m;
                 double v = fits_vis_v.getXY(ant1,ant2) / wavelength_m;
                 double w = fits_vis_w.getXY(ant1,ant2) / wavelength_m;
                 double uv_distance = sqrt(u*u + v*v);
              
                 if( uv_distance > min_uv ){ // check if larger than minimum UV distance 
//                 int u_index = round( (u - u_min)/delta_u );
//                 int v_index = round( (v - v_min)/delta_v );
                    double u_pix = round( u/delta_u );
                    double v_pix = round( v/delta_v );
                    int u_index = u_pix + n_pixels/2; // was u - u_center
                    int v_index = v_pix + n_pixels/2; // was v - v_center
                 
//                 printf("DEBUG : u_index %d vs. %d ( u = %.2f , u_min = %.2f , delta_u = %.2f , u_center =%.2f)\n",u,u_index1,u_index,u_min,delta_u,u_center);
              
                 // Using CELL averaging method or setXY ?
                    uv_grid_real.addXY( u_index, v_index, re );
                    uv_grid_imag.addXY( u_index, v_index, im );
                    int count = uv_grid_counter.getXY( u_index, v_index );
                    uv_grid_counter.setXY( u_index, v_index , count + 1 );
              
                 // add conjugates :
//                 u_index = round( (-u - u_min)/delta_u );
//                 v_index = round( (-v - v_min)/delta_v );
                    u_index = -u_pix + n_pixels/2; // was round( (-u - u_center)/delta_u ) + ...
                    v_index = -v_pix + n_pixels/2; // was round( (-v - v_center)/delta_v ) + ...
                    uv_grid_real.addXY( u_index, v_index, re );
                    uv_grid_imag.addXY( u_index, v_index, -im );
                    count = uv_grid_counter.getXY( u_index, v_index );
                    uv_grid_counter.setXY( u_index, v_index , count + 1 );
                                  
                           
                    added++;
                 }
              }else{
                 printf("DEBUG : visibility value %e +j%e higher than limit %e -> skipped\n",re,im,MAX_VIS);
                 high_value++;
              }
           }
        }
     }
  }  
  printf("DEBUG : added %d UV points to the grid, %d too high values skipped\n",added,high_value);    
  
  // This division is in fact UNIFORM weighting !!!! Not CELL-avareging 
  // normalisation to make it indeed CELL-averaging :
  if( strcmp(weighting, "U" ) == 0 ){
     uv_grid_real.Divide( uv_grid_counter );
     uv_grid_imag.Divide( uv_grid_counter );
  }
  
  char uv_grid_re_name[1024],uv_grid_im_name[1024],uv_grid_counter_name[1024];
  sprintf(uv_grid_re_name,"uv_grid_real_%dx%d.fits",n_pixels,n_pixels);
  sprintf(uv_grid_im_name,"uv_grid_imag_%dx%d.fits",n_pixels,n_pixels);
  sprintf(uv_grid_counter_name,"uv_grid_counter_%dx%d.fits",n_pixels,n_pixels);
  
  if( uv_grid_real.WriteFits( uv_grid_re_name ) ){
     printf("ERROR : could not write output file %s\n",uv_grid_re_name);
  }

  if( uv_grid_imag.WriteFits( uv_grid_im_name ) ){
     printf("ERROR : could not write output file %s\n",uv_grid_im_name);
  }
  
  if( uv_grid_counter.WriteFits( uv_grid_counter_name ) ){
     printf("ERROR : could not write output file %s\n",uv_grid_counter_name);
  }

}

bool run_imager( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, 
                 double frequency_mhz,
                 int    n_pixels,
                 double FOV_degrees,
                 double min_uv,
                 bool   do_gridding,    // excute gridding  (?)                  
                 bool   do_dirty_image, // form dirty image (?)
                 const char* weighting="",    // weighting : U for uniform (others not implemented)
                 const char* in_fits_file_uv_re="", // gridded visibilities can be provided externally
                 const char* in_fits_file_uv_im=""  // gridded visibilities can be provided externally
               )  
{
  // based on RTS : UV pixel size as function FOVtoGridsize in  /home/msok/mwa_software/RTS_128t/src/gridder.c  
  double frequency_hz = frequency_mhz*1e6;
  double wavelength_m = VEL_LIGHT / frequency_hz;
  double delta_u = ( (VEL_LIGHT/frequency_hz)/(FOV_degrees*M_PI/180.) ) /wavelength_m; // in meters (NOT WAVELENGHTS)
  double delta_v = ( (VEL_LIGHT/frequency_hz)/(FOV_degrees*M_PI/180.) ) / wavelength_m; // in meters (NOT WAVELENGHTS)

  CBgFits uv_grid_counter( n_pixels, n_pixels ),uv_grid_real( n_pixels, n_pixels ) , uv_grid_imag( n_pixels, n_pixels );  
  if( do_gridding ){
     gridding( fits_vis_real, fits_vis_imag, fits_vis_u, fits_vis_v, fits_vis_w, uv_grid_real, uv_grid_imag, uv_grid_counter, delta_u, delta_v, frequency_mhz, n_pixels, min_uv, weighting );
  }else{
     if( strlen(in_fits_file_uv_re) && strlen(in_fits_file_uv_im) ){
        uv_grid_counter.SetValue(1.00);
        
        printf("Reading fits file %s ...\n",in_fits_file_uv_re);
        if( uv_grid_real.ReadFits( in_fits_file_uv_re, 0, 1, 1 ) ){
           printf("ERROR : could not read visibility FITS file %s\n",in_fits_file_uv_re);
           exit(-1); 
        }else{
           printf("OK : fits file %s read ok\n",in_fits_file_uv_re);
        }

        printf("Reading fits file %s ...\n",in_fits_file_uv_im);
        if( uv_grid_real.ReadFits( in_fits_file_uv_im, 0, 1, 1 ) ){
           printf("ERROR : could not read visibility FITS file %s\n",in_fits_file_uv_im);
           exit(-1); 
        }else{
           printf("OK : fits file %s read ok\n",in_fits_file_uv_im);
        }

     }else{
        printf("ERROR : when gridding is disabled (-g 0) options -r and -i with REAL and IMAG FITS file names must be specified -> cannot continue !\n");
        exit(-1);
     }
  }

  if( do_dirty_image ){
     // dirty image :  
     printf("PROGRESS : executing dirty image\n");
     dirty_image( uv_grid_real, uv_grid_imag, uv_grid_counter );
  }

  return  true;   
}

int main(int argc,char* argv[])
{
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
  

  // read input data (correlation matrix and UVW) :
  CBgFits fits_vis_real, fits_vis_imag, fits_vis_u, fits_vis_v, fits_vis_w;
  if( read_corr_matrix( in_basename.c_str(), fits_vis_real, fits_vis_imag, fits_vis_u, fits_vis_v, fits_vis_w, gPostfix.c_str() ) ){
     printf("OK : input files read ok\n");
  }else{
     printf("ERROR : could not read one of the input files\n");
  }
  
  // run imager:
  run_imager( fits_vis_real, fits_vis_imag, fits_vis_u, fits_vis_v, fits_vis_w, 
              gFrequencyMHz, gImageSize, gFOV_degrees, gMinUV, 
              gDoGridding, gDo_dirty_image, 
              gUniformWeighting.c_str(), g_in_fits_file_uv_re.c_str(), g_in_fits_file_uv_im.c_str() );
}

