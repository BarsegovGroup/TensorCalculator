SOURCES=tensors.cpp IO/config_reader.cpp IO/dcdio.cpp IO/pdbio.cpp IO/topio.cpp functions.cpp IO/tensorio.cpp energy.cpp Tensors/cauchy_stress.cpp Tensors/sphere_stress.cpp Tensors/strain_tensor.cpp Tensors/deform_gradient.cpp Tensors/stretch_tensor.cpp

OBJDIR=obj

OBJECTS=$(OBJDIR)/%.cpp.o: %.cpp

EXECUTABLE=TensorCalculator

dir_guard=@mkdir -p $(@D)

$(OBJDIR)/%.cpp.o: %.cpp
	$(dir_guard)
	g++ -Wall -o $@ -c $<

objects = $(patsubst %.cpp, $(OBJDIR)/%.cpp.o, $(1))

SOPGPU_OBJS := $(call objects, $(SOURCES))

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(SOPGPU_OBJS)
	g++ -Wall $(SOPGPU_OBJS) -o $@

clean:
	rm -r obj/ tensors
