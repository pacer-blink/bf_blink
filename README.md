## msok_imager

Bifrost realtime imager plugin for CHASM

Author: Marcin Sokolowski

### Compiling

To build your plugin:

0) Setup your build environment (edit/run `source setup_env.sh`).
2) compile with meson by running:

```
meson setup build # add --wipe to recompile
cd build
meson compile
```

### Using msok_imager

From root directory, you can run:

```python
from build.msok_imager import *
import numpy as np 
import bifrost as bf
import time
import h5py


# Configuration
filepath        = "testdata/output.hdf5"
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
antenna_positions_file  = "/home/aavs/aavs-calibration/config/eda2/antenna_locations.txt"

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
```

