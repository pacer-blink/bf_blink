from bifrost.libbifrost import _check, _get, BifrostObject
from bifrost.ndarray import asarray

import msok_imager_generated as _gen

class MsokImager(BifrostObject):
    def __init__(self):
        BifrostObject.__init__(self, _gen.MsokImagerCreate, 
                               _gen.MsokImagerDestroy)
    def init(self):
        _check(_gen.MsokImagerInit(self.obj))

    def execute(self, in_BFarray, out_BFarray):
        _check(_gen.MsokImagerExecute(self.obj, asarray(in_BFarray).as_BFarray(),
                                asarray(out_BFarray).as_BFarray()))
        return out_BFarray

    def set_stream(self, stream_ptr_generic):
        _check(_gen.MsokImagerSetStream(self.obj, stream_ptr_generic))

    def reset_state(self):
        _check(_gen.MsokImagerResetState(self.obj))