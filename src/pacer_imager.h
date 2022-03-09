#ifndef _PACER_IMAGER_H__
#define _PACER_IMAGER_H__

#include "antenna_positions.h"

#include <string>
using namespace std;

class CBgFits;

class CImagerParameters
{
public :
   string m_AntennaPositionsFile;

   CImagerParameters();
   void SetParameters( const char* szAntennaPositionsFile );
   
//   void SetAntPositionsFile( const char* szAntennaPositionsFile );
   
};

class CPacerImager
{
protected :
   // FFT shift to convert from DC in bin 0 of the FFTW output array to DC in the center bin :
   void fft_shift( CBgFits& dirty_image, CBgFits& out_image );

   // FFT unshift converts from DC in the center bin to DC in bin 0 (as expected input to complex FFTW)
   void fft_unshift( CBgFits& dirty_image, CBgFits& out_image );



public :
   // parameters :
   CImagerParameters m_ImagerParameters;

   // Antenna positions :   
   CAntennaPositions m_AntennaPositions;


   CPacerImager();
   ~CPacerImager();
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // 1st version producing a dirty image (tested on both MWA and SKA-Low).
   // TODO : Test cases can be found in PaCER documentation 
   //-----------------------------------------------------------------------------------------------------------------------------
   void dirty_image( CBgFits& uv_grid_real_param, CBgFits& uv_grid_imag_param, CBgFits& uv_grid_counter, 
                     bool bSaveIntermediate=false, const char* szBaseOutFitsName=NULL, bool bSaveImaginary=true );

   //-----------------------------------------------------------------------------------------------------------------------------
   // INPUT  : 
   //          fits_vis_real, fits_vis_imag : visibilities (REAL and IMAG 2D arrays as FITS class) 
   //          fits_vis_u, fits_vis_v, fits_vis_w : UVW (real values baselines in units of wavelength - see TMS)
   //          delta_u, delta_v : size of the UV cell 
   //          frequency_mhz : frequency in MHz
   //
   // OUTPUT : 
   //          - uv_grid_real, uv_grid_imag : visibilities on UV grid (real and imag arrays)
   //          - uv_grid_counter : visibility counter and 
   //-----------------------------------------------------------------------------------------------------------------------------
   void gridding( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                  CBgFits& uv_grid_real, CBgFits& uv_grid_imag, CBgFits& uv_grid_counter, double delta_u, double delta_v, 
                  double frequency_mhz,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="" // weighting : U for uniform (others not implemented)
                );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // Executes imager:
   // INPUT  : 
   //          fits_vis_real, fits_vis_imag : visibilities (REAL and IMAG 2D arrays as FITS class) 
   //          fits_vis_u, fits_vis_v, fits_vis_w : UVW (real values baselines in units of wavelength - see TMS)     
   //-----------------------------------------------------------------------------------------------------------------------------
   bool run_imager( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, 
                    double frequency_mhz,
                    int    n_pixels,
                    double FOV_degrees,
                    double min_uv=-1000,        // minimum UV 
                    bool   do_gridding=true,    // excute gridding  (?)
                    bool   do_dirty_image=true, // form dirty image (?)
                    const char* weighting="",   // weighting : U for uniform (others not implemented)
                    const char* in_fits_file_uv_re="", // gridded visibilities can be provided externally
                    const char* in_fits_file_uv_im="",  // gridded visibilities can be provided externally
                    const char* szBaseOutFitsName=NULL
                  );

                  
   bool run_imager( float* data_real, 
                    float* data_imag,
                    int n_ant, 
                    int n_pol,
                    double frequency_mhz, 
                    int n_pixels,
                    double FOV_degrees,
                    double min_uv=-1000,      // minimum UV
                    bool do_gridding=true,    // excute gridding  (?)
                    bool do_dirty_image=true, // form dirty image (?)
                    const char* weighting="", // weighting : U for uniform (others not implemented)
                    const char* szBaseOutFitsName=NULL
                  );
};

#endif 