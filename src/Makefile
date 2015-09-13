CXX				:= gcc
LINKER				:= gcc -fPIC
OBJDIR				:= obj
CFLAGS				:= -O3 -g $(INCLUDES)
LDFLAGS				:= -lrt -lm

include common.mk

# files
TENSOR_SOURCES	:= \
	functions.c \
	energy.c \
	IO/dcdio.c \
	IO/pdbio.c \
	IO/topio.c \
	IO/config_reader.c \
	IO/tensorio.c \
	Tensors/cauchy_stress.c \
	Tensors/deform_gradient.c \
	Tensors/sphere_stress.c \
	Tensors/strain_tensor.c \
	Tensors/stretch_tensor.c \
	tensors.c 

all: tensorsCalc

TENSOR_OBJS			:= $(call objects, $(TENSOR_SOURCES))
-include $(TENSOR_OBJS:.o=.d)

tensorsCalc: $(TENSOR_OBJS)
	$(LINKER) -o $@ $(TENSOR_OBJS) $(LDFLAGS)

clean:
	rm -f tensorsCalc
	rm -rf "$(OBJDIR)"
	
.PHONY: makedirs clean all
