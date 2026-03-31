project('dg-swem', ['cpp', 'fortran'],
        default_options: ['buildtype=debugoptimized'],
        meson_version: '>=1.1.0')

add_project_arguments('-DRKSSP', language: 'fortran')
add_project_arguments('-DSLOPE5', language: 'fortran')

compiler_name = meson.get_compiler('fortran').get_id()
message('Using compiler', compiler_name)

if get_option('precision') == 'double'
  sz = 8
else
  sz = 4
endif 

if compiler_name == 'gfortran'
  add_project_arguments('-ffixed-line-length-132', language: 'fortran')
  if sz == 8
    add_project_arguments('-fdefault-real-8', 
                          '-fdefault-double-8',
                          language: 'fortran')
  endif
elif compiler_name == 'nvidia_hpc'
  add_project_arguments('-Mextend',
                        f'-r@sz@',
                        language: 'fortran')
  if get_option('gpu')
    add_project_arguments('-cuda',
                          '-acc',
                          '-gpu=mem:unified,lineinfo',
                          '-Minfo=accel',
                          language: 'fortran')
    add_project_link_arguments('-cuda',
                               '-acc',
                               '-gpu=mem:unified',
                               language: 'fortran')
  endif
elif compiler_name == 'intel-llvm'
  add_project_arguments(f'-r@sz@',
                        '-xHost',
                        '-132',
                        language: 'fortran')

else
  error('Compiler not supported.')
endif

    


subdir('src')
