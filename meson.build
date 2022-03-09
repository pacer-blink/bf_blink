# Bifrost wrapper file
project('msok_imager', 'cpp', version: '0.0.1',
         default_options : ['cpp_std=c++11'])
plugin_name = 'msok_imager'
cc = meson.get_compiler('cpp')

#add_global_arguments('-D_UNIX', language : 'cpp')


# CUDA setup
#cuda = import('unstable-cuda')
# cuda_dep = dependency('cuda', version : '>=10', modules : ['cudart', 'nvrtc'])

# Dependency: FFTW3
fftw3_dep = dependency('fftw3')

# Dependency: CFITSIO 
cfitsio_dep = dependency('cfitsio')

# Dependency: libnova
libnova_dep = cc.find_library('libnova', dirs : '/usr/lib/x86_64-linux-gnu') 

# Dependency: dl (linking)
ldl_dep     = cc.find_library('dl') 

# Dependency: bifrost
bf_dep = dependency('bifrost')
bifrost_src_path = '/home/dancpr/install/bifrost/src'
mslib_src_path   = '/home/aavs/Software/msfitslib/src'

mslib_dep = cc.find_library('msfitslib', dirs : '/home/aavs/Software/msfitslib/build') 


include_dirs = include_directories(
    bifrost_src_path, 
    bifrost_src_path / 'bifrost', 
    mslib_src_path
    )

# CUDA flags
#nvflags = cuda.nvcc_arch_flags('11.1', 'Volta')
#message(cuda.min_driver_version('11.1'))

# Generate python wrapper from .h file
builddir = meson.project_build_root()

# Generate shared object from imager cpp code (non-bifrost)
pacer_imager_lib = library('pacer_imager', ['src/pacer_imager.cpp', 'src/antenna_positions.cpp'], 
              version : '0.0.1', soversion : '0',  cpp_args :'-D_UNIX',
              dependencies: [fftw3_dep, mslib_dep, libnova_dep, ldl_dep, cfitsio_dep], 
              include_directories: include_dirs)

# Generate shared object from imager cpp code (bifrost wrapper)
msok_lib = library('msok_imager', 'src/msok_imager.cpp', 
              version : '0.0.1', soversion : '0', 
              link_with: pacer_imager_lib,
              dependencies: [bf_dep, mslib_dep],  
              include_directories: include_dirs)

#msok_exe = executable('msok_imager', 'src/msok_imager.cpp', 
#              link_with: pacer_imager_lib,
#              dependencies: [bf_dep, mslib_dep], 
#              include_directories: include_dirs)

message('Generating ctypes wrapper from msok_imager.h')
#r = run_command('./generate_wrapper.py', 'src' / plugin_name, '-o', builddir, check: true)
r = run_command('cp', 'msok_imager.py', builddir)
message(r.stdout().strip())
r = run_command('cp', 'msok_imager_generated.py', builddir)
message(r.stdout().strip())
