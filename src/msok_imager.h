#include <bifrost/array.h>
#include <bifrost/common.h>
#include <bifrost/ring.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct bf_msimg* bfplugin;

BFstatus MsokImagerCreate(bfplugin* plugin_ptr);
BFstatus MsokImagerInit(bfplugin plugin, const char* antenna_positions_file, int n_ant, 
                        int n_pol, float frequency_mhz, int n_pixels, float fov_deg, 
                        float min_uv, bool do_gridding, bool do_dirty_image, 
                        const char* weighting, const char* base_out_fits_name);
BFstatus MsokImagerSetStream(bfplugin plugin,
                       void const* stream);
BFstatus MsokImagerResetState(bfplugin plugin);
BFstatus MsokImagerExecute(bfplugin plugin,
                     BFarray const* in_data_real,
                     BFarray const* in_data_imag);
BFstatus MsokImagerDestroy(bfplugin plugin);

#ifdef __cplusplus
} // extern "C"
#endif
