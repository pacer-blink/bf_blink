#include <bifrost/array.h>
#include <bifrost/common.h>
#include <bifrost/ring.h>
#include <assert.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <iostream>

#include "msok_imager.h"
#include "pacer_imager.h"

/*!
  \class bf_msimg
  \brief Plugger class for MsokImager
*/
class bf_msimg {
private:
    // Parameters passed during init() may be stored here, use _n_param 
    int _n_ant;
    int _n_pol;
    float _frequency_mhz;
    int _n_pixels;
    float _fov_deg;
    float _min_uv;
    bool _do_gridding;
    bool _do_dirty_image;
    char _weighting[255];
    char _base_out_fits_name[255];
    char _antenna_positions_file[255];
    CPacerImager _imager;
    
public:
    // Constructor
    bf_msimg() {

    }

    // Destructor
    ~bf_msimg() {

    }
    
    // Initialize your plugin here
    void init(const char* antenna_positions_file, int n_ant, int n_pol, float frequency_mhz, int n_pixels, float fov_deg, float min_uv, 
              bool do_gridding, bool do_dirty_image, const char* weighting, const char* base_out_fits_name) {
    
      _n_ant = n_ant;
      _n_pol = n_pol;
      _frequency_mhz = frequency_mhz;
      _n_pixels =  n_pixels;
      _fov_deg =  fov_deg; 
      _min_uv = min_uv;
      _do_gridding = do_gridding;
      _do_dirty_image = do_dirty_image;
      strcpy(_weighting, weighting);
      strcpy(_base_out_fits_name, base_out_fits_name);
      strcpy(_antenna_positions_file, antenna_positions_file);
      _imager.m_ImagerParameters.m_bConstantUVW = 1;
      _imager.m_ImagerParameters.SetGlobalParameters( _antenna_positions_file, 1);
    }

    // Do any zeroing / memset stuff here
    void reset_state() {
        
    }

    // execute your plugin
    void exec(BFarray const* in_data_real, BFarray const* in_data_imag) {
    float* d_real = (float *)in_data_real->data;
    float* d_imag = (float *)in_data_imag->data;

    bool imager_ret = _imager.run_imager( d_real, d_imag,
                                          _n_ant, _n_pol, _frequency_mhz, _n_pixels, _fov_deg, _min_uv,
                                          _do_gridding, _do_dirty_image, 
                                          _weighting, _base_out_fits_name);
    }
};

// Used by bifrost python wrapper at instantiation
BFstatus MsokImagerCreate(bfplugin* plugin_ptr) {
    BF_ASSERT(plugin_ptr, BF_STATUS_INVALID_POINTER);
    BF_TRY_RETURN_ELSE(*plugin_ptr = new bf_msimg(),
                       *plugin_ptr = 0);
}

// Initialisation for plugin 
BFstatus MsokImagerInit(bfplugin plugin, const char* antenna_positions_file, int n_ant, 
                        int n_pol, float frequency_mhz, int n_pixels, float fov_deg, 
                        float min_uv, bool do_gridding, bool do_dirty_image, 
                        const char* weighting, const char* base_out_fits_name) {
    BF_ASSERT(plugin, BF_STATUS_INVALID_HANDLE);
    BF_TRY_RETURN(plugin->init(antenna_positions_file, n_ant, n_pol, frequency_mhz, n_pixels, 
                               fov_deg, min_uv, do_gridding, do_dirty_image, 
                               weighting, base_out_fits_name));
}

// Assign to CUDA stream
BFstatus MsokImagerSetStream(bfplugin plugin, void const* stream) {
        //BF_ASSERT(plan, BF_STATUS_INVALID_HANDLE);
       // BF_ASSERT(stream, BF_STATUS_INVALID_POINTER);
       //BF_TRY_RETURN(plan->set_stream(*(cudaStream_t*)stream));
}

// Reset state of any internal memory 
BFstatus MsokImagerResetState(bfplugin plugin) {
        BF_ASSERT(plugin, BF_STATUS_INVALID_HANDLE);
        BF_TRY_RETURN(plugin->reset_state());
}

// Main method to execute data processing tasks
BFstatus MsokImagerExecute(bfplugin plugin,
                     BFarray const* in_data_real,
                     BFarray const* in_data_imag) {

        BF_ASSERT(plugin, BF_STATUS_INVALID_HANDLE);
        BF_TRY_RETURN(plugin->exec(in_data_real, in_data_imag));
}

// Called by python wrapper at deletion time
BFstatus MsokImagerDestroy(bfplugin plugin) {
    BF_ASSERT(plugin, BF_STATUS_INVALID_HANDLE);
    delete plugin;
    return BF_STATUS_SUCCESS;
}

//int main() {
//    printf("Uncomment this if you wan to compile an executable for debugging.");
//    return 0;
// }
