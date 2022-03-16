#
# this file was created by a computer. trust it.
#

# compiler tools
XILINX_VITIS ?= /home/chomper/workspace/Xilinx/Vitis/2021.2
XILINX_XRT ?= /opt/xilinx/xrt
XILINX_VIVADO ?= /home/chomper/workspace/Xilinx/Vivado/2021.2
XILINX_HLS ?= $(XILINX_VITIS)/Vivado_HLS

VPP_LINKER ?= ${XILINX_VITIS}/bin/v++
RM = rm -f
RMDIR = rm -rf

VITIS_PLATFORM = xilinx_u55n_gen3x4_xdma_1_202110_1
VITIS_PLATFORM_PATH = $(VITIS_PLATFORM)

VPP_OPTS = --target hw

#
# Hw kernel files
#

BINARY_CONTAINERS += compute.xclbin

BUILD_SUBDIRS += compute.build
COMPUTE_OBJS += ../../compute_kernels/Hardware/build/compute.xo

ALL_MESSAGE_FILES = $(subst .xclbin,.mdb,$(BINARY_CONTAINERS))

#
# primary build targets
#

.PHONY: all clean
all: $(BINARY_CONTAINERS)


clean:
	-$(RM) $(BINARY_CONTAINERS) $(ALL_MESSAGE_FILES) 
	-$(RM) *.xclbin.sh *.xclbin.info *.xclbin.link_summary*
	-$(RMDIR) $(BUILD_SUBDIRS)
	-$(RMDIR) .Xil

ifeq ($(MAKECMDGOALS),incremental)
.PHONY: compute.xclbin
.PHONY: ../../compute_kernels/Hardware/build/compute.xo
endif

.PHONY: incremental
incremental: all


nothing:

#
# binary container: compute.xclbin
#

compute.xclbin: $(COMPUTE_OBJS) compute-link.cfg
	-@echo $(VPP_LINKER) $(VPP_OPTS) --link --config compute-link.cfg -o"$@" $(COMPUTE_OBJS) > compute.xclbin.sh
	$(VPP_LINKER) $(VPP_OPTS) --link --config compute-link.cfg -o"$@" $(COMPUTE_OBJS)

