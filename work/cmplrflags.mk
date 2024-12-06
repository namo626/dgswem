# sb46.50.02 These flags work on the UT Austin Lonstar cluster.
ifneq (,$(findstring ifx,$(FC)))
  PPFC          :=  ifx
  PFC           :=  mpiifx
  CC            :=  icx
  FFLAGS1       :=  -r8 $(INCDIRS) -O3 -xHost -msse3 -132  #-traceback -check all #-prof-gen -prof-dir /bevo2/michoski/v21/work -pg -prof-use
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  FFLAGS4       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA -DRKSSP -DSLOPE5 #-DDGOUT #-DOUT_TEC #-DARTDIF  #  -DWETDR # -DSED_LAY -DOUT_TEC #-DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL -DFILTER
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DRKSSP -DSLOPE5 #-DDGOUT #-DOUT_TEC #-DARTDIF #-DWETDR  # -DSED_LAY -DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL
  DPRE          :=  -DREAL8 -DLINUX -DRKSSP -DSLOPE5 #-DOUT_TEC #-DSWAN #-DARTDIF  # -DWETDR #-DSED_LAY -DSWAN #-DOUT_TEC #-DSWAN -DRKC -DTRACE -DSED_LAY -DCHEM -DP0 -DP_AD -DSLOPEALL
  DPRE2         :=  -DREAL8 -DLINUX -DCMPI 
  CFLAGS        :=  -O3 -xSSSE3 -I. -Wno-implicit-function-declaration
  IMODS         :=  -module
  LIBS          :=  -L ../metis -lmetis
  MSGLIBS       :=

else ifneq (,$(findstring nvfortran,$(FC)))   # NVIDIA
  sz            := 8
  PPFC	        :=  nvfortran
  PFC	        :=  mpif90
ifeq ($(gpu),1)
  FFLAGS1	:=  -r$(sz) -Mextend -Mlarge_arrays -cuda -traceback -g -O3 -acc -gpu=unified,lineinfo -Minfo=accel
else
  FFLAGS1	:=  -r$(sz) -Mextend -traceback -g -O3 -tp=native
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
  CFLAGS        :=  -O3 -I. -fastsse -DLINUX
  LIBS  	:=  -L ../metis  -lmetis
  MSGLIBS	:=

else
  sz            := 8
  ifeq ($(sz),8)
    RFLAG = -fdefault-real-8 -fdefault-double-8
  else
    RFLAG = -freal-8-real-4
  endif
  PPFC	        :=  gfortran
  PFC	        :=  mpif90
  FFLAGS1	:= $(RFLAG) -g -O3 -march=native -ffixed-line-length-132 -std=legacy -fallow-argument-mismatch
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
  CFLAGS        :=  -O3 -I.  -DLINUX
  LIBS  	:=  -L ../metis  -lmetis
  MSGLIBS	:=
endif
