ifneq (,$(findstring ifx,$(FC))) # Intel
  PPFC          :=  ifx
  PFC           :=  mpiifx
  CC            :=  icx
  FFLAGS1       :=  -r8 $(INCDIRS) -O3 -xHost -132
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  FFLAGS4       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA -DRKSSP -DSLOPE5
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DRKSSP -DSLOPE5
  DPRE          :=  -DREAL8 -DLINUX -DRKSSP -DSLOPE5
  DPRE2         :=  -DREAL8 -DLINUX -DCMPI 
  CFLAGS        :=  -O3 -I. -Wno-implicit-function-declaration
  IMODS         :=  -module
  LIBS          :=  -L ../metis -lmetis
  MSGLIBS       :=

else ifneq (,$(findstring nvfortran,$(FC)))   # NVIDIA
  sz            := 8
  PPFC	        :=  nvfortran
  PFC	        :=  mpif90
  ifeq ($(gpu),1)
	FFLAGS1	:=   -Mextend -r$(sz)   -cuda -traceback -g -O0 -acc -gpu=mem:unified,lineinfo,stacklimit:nostacklimit -Minfo=accel
  else
    FFLAGS1	:=  -Mextend -Mbounds -r$(sz) -traceback -g -O3 -tp=native
  endif
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  FFLAGS4	:=  $(FFLAGS1)
  DA  	        :=  -DREAL$(sz) -DLINUX -DCSCA -DRKSSP -DSLOPE5
  DP  	        :=  -DREAL$(sz) -DLINUX -DCSCA -DCMPI -DRKSSP -DSLOPE5
  DPRE	        :=  -DREAL$(sz) -DLINUX -DRKSSP -DSLOPE5
  DPRE2         :=  -DREAL$(sz) -DLINUX -DCMPI
  IMODS 	:=  -module
  CC            :=  nvc
  CXX    := nvc++
  CXXFLAGS := -O3 -g
  CFLAGS        :=  -O3 -I. -DLINUX
  LIBS  	:=  -L ../metis  -lmetis
  MSGLIBS	:=

else  # GCC
  sz            := 8
  ifeq ($(sz),8)
    RFLAG = -fdefault-real-8 -fdefault-double-8
  endif
  PPFC	        :=  gfortran
  PFC	        :=  mpif90
  FFLAGS1	:= $(RFLAG) -g -O3 -march=native -ffixed-line-length-132 -std=legacy -fallow-argument-mismatch -lz
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  FFLAGS4	:=  $(FFLAGS1)
  DA  	        :=  -DREAL$(sz) -DLINUX -DCSCA -DRKSSP -DSLOPE5
  DP  	        :=  -DREAL$(sz) -DLINUX -DCSCA -DCMPI -DRKSSP -DSLOPE5
  DPRE	        :=  -DREAL$(sz) -DLINUX -DRKSSP -DSLOPE5
  DPRE2         :=  -DREAL$(sz) -DLINUX -DCMPI
  IMODS 	:=  -J
  CC            :=  gcc
  CXX    := g++
  CXXFLAGS := -O3 -g
  CFLAGS        :=  -O3 -I.  -DLINUX -Wno-incompatible-pointer-types
  LIBS  	:=  -L ../metis  -lmetis
  MSGLIBS	:=
endif
