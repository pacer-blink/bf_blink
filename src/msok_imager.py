from bifrost.libbifrost import _check, _get, BifrostObject
from bifrost.ndarray import asarray
try:
    from . import msok_imager_generated as _gen
except:
    import msok_imager_generated as _gen
    
class MsokImager(BifrostObject):
    def __init__(self):
        BifrostObject.__init__(self, _gen.MsokImagerCreate, _gen.MsokImagerDestroy)

    #void init(const char* antenna_positions_file, int n_ant, int n_pol, float frequency_mhz, int n_pixels, float fov_deg, float min_uv, 
    #          bool do_gridding, bool do_dirty_image, const char* weighting, const char* base_out_fits_name
    def init(self, antenna_positions_file, n_ant, n_pol, frequency_mhz, n_pixels,
             fov_deg, min_uv, do_gridding, do_dirty_image, weighting, base_out_fits_name):
        _check(_gen.MsokImagerInit(self.obj, antenna_positions_file, n_ant, n_pol, 
                                   frequency_mhz, n_pixels, fov_deg, min_uv, 
                                   do_gridding, do_dirty_image, weighting, base_out_fits_name))

    def execute(self, d_real, d_imag):
        _check(_gen.MsokImagerExecute(self.obj, asarray(d_real).as_BFarray(),
                                asarray(d_imag).as_BFarray()))
        return 0

    def set_stream(self, stream_ptr_generic):
        pass
        #_check(_gen.BTccSetStream(self.obj, stream_ptr_generic))

    def reset_state(self):
        _check(_gen.MsokImagerResetState(self.obj))

if __name__ == "__main__":
    import numpy as np 
    import bifrost as bf
    import time
    import h5py

    # Configuration
    filepath        = "../testdata/visibility_test_data.hdf5"
    n_ant           = 256
    n_pol           = 2
    frequency_mhz   = 159.375
    n_pixels        = 180
    fov_deg         = 180.0
    min_uv          = -1000
    do_gridding     = True
    do_dirty_image  = True
    weighting               = "U"
    base_out_fits_name      = "tcc_test"
    antenna_positions_file  = "../testdata/antenna_locations.txt"

    # Load data
    fh = h5py.File(filepath)
    corr_matrix = fh['data_matrix'][0,:,:,0,0]
    d_real  = bf.ndarray(np.ascontiguousarray(np.real(corr_matrix)), dtype='f32')
    d_imag  = bf.ndarray(np.ascontiguousarray(np.imag(corr_matrix)), dtype='f32')
    print("D_REAL PYTHON 0,1", d_real.ravel()[0:2])


    # Create imager instance
    msok = MsokImager()

    # Initialize
    msok.init(antenna_positions_file, n_ant, n_pol, frequency_mhz, n_pixels,
             fov_deg, min_uv, do_gridding, do_dirty_image, weighting, base_out_fits_name)

    # Execute
    t0 = time.time()
    msok.execute(d_real, d_imag)
    t1 = time.time()
    print(f'Time taken: {t1-t0:.3f}s')
    
