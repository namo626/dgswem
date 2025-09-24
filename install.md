Installation
=======
To obtain the developmental version, clone the repo. Then set the compiler of choice:

    export FC=ifx
    export FC=nvfortran

If not set, the default is `gfortran`. 
Run

    cd work
    make all sz=4 

to build all executables (`dgswem_serial, dgswem, adcprep, adcpost`) with single precision (use `sz=8` for double precision).

## GPU support
The code can be run on NVIDIA GPUs through OpenACC. Heterogenous memory management (HMM) must be enabled as it uses unified memory. Change the `gpu` flag to 1.

    make all sz=4 gpu=1
