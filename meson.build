# Bifrost plugin for msok_imager
project('msok_imager', 'cpp', version: '0.0.1')

# Basic setup
plugin_name = 'msok_imager'

comp_flags = ['-Wno-unused-variable', '-Wno-dangling-else', '-Wno-sign-compare', '-Wno-parentheses',
              '-Wno-non-virtual-dtor', '-Wno-format', '-Wno-reorder', '-Wno-unused-but-set-variable']
add_global_arguments(comp_flags, language : 'cpp')

# Dependency: bifrost
bf_dep = dependency('bifrost')
bifrost_src_path = '/home/dancpr/software/bifrost/src'

# Subprojects: msfitslib and imager
msfits     = subproject('msfitslib')
imager     = subproject('imager')
correlator = subproject('correlator')

msfitslib_dep = msfits.get_variable('msfitslib_dep')
imager_dep = imager.get_variable('imager_dep')

# List of directories to include when building plugin
include_dirs = include_directories(
    bifrost_src_path, 
    bifrost_src_path / 'bifrost', 
    #mslib_src_path
    )

# Generate shared object from imager cpp code (non-bifrost)
#pacer_imager_lib = library('pacer_imager', ['src/pacer_imager.cpp', 'src/antenna_positions.cpp'], 
#              version : '0.0.1', soversion : '0',  cpp_args :'-D_UNIX',
#              dependencies: [fftw3_dep, mslib_dep, libnova_dep, ldl_dep, cfitsio_dep], 
#              include_directories: include_dirs)

# Generate shared object from imager cpp code (bifrost wrapper)
msok_lib = library('msok_imager', 'src/msok_imager.cpp', 
              version : '0.0.1', soversion : '0', 
              dependencies: [imager_dep, msfitslib_dep],  
              include_directories: include_dirs)

# Generate an executable (useful for debugging segfaults)
#msok_exe = executable('msok_imager', 'src/msok_imager.cpp', 
#              link_with: pacer_imager_lib,
#              dependencies: [bf_dep, mslib_dep], 
#              include_directories: include_dirs)

#message('Copying python wrapper for msok_imager')
#r = run_command('cp', 'src/msok_imager.py', builddir)
#r = run_command('cp', 'src/msok_imager_generated.py', builddir)
#r = run_command('cp', 'src/__init__.py', builddir)
#message(r.stdout().strip())

