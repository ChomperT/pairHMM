# PairHMM Accelerator by FPGA

## Source Tree
```
├── E_coli_K12
│   ├── benchmark           # benchmark CPU && FPGA
│   ├── compute.xclbin      # xclbin
│   ├── do_index            # index E_coli_K12
│   ├── gatk_cpu            # run gatk with cpu mode
│   ├── gatk_fpga           # run gatk with fpga mode
│   ├── input               # source file
│   │   ├── fasta
│   │   └── fastq
│   └── libcompute.so       # jni impl by xilopencl
├── prebuilt
│   └── gatk.jar            # gatk prebuit
├── README.md                # Readme
└── src
    ├── gatk                 # gatk patch && build 
    │   ├── gatk.patch
    ├── jni
    │   └── compute         # libcompute source
    └── vitis
        ├── compute_kernels     # vitis kernel
        ├── gene_system         
        └── gene_system_hw_link
```

## Features
- Configurable compute units
- Configurable parallel scale
- Configurable computational refinement
- Variable input and output length
- Multi-core parallel computing
- Definable sequence length
- Flexible configuration according to hardware resource size
- Supports different scenarios such as Vitis CL/HLS/SD Accel


## Support Hardware
- Zynq/MP（with generic uio）
  - Ultra96 6cu @ 500MHz
- 7 series / UltraScale / Plus / HBM (xdma) 
  - VCU1525 24cu @ 250 MHz
- Vitis (OpenCL)
  - C1100 6cu @ 300 MHz

## Get Started
- depend
  - Ubuntu Linx
  - Xilinx XRT for C1100
  - cmake bwa samtools openjdk-11-jdk
  
- run test prebuilt
    ```bash
    git clone https://github.com/ChomperT/pairHMM
    sudo apt-get install cmake bwa samtools
    cd pairHMM/E_coli_K12

    ./do_index  # bwa index

    ./benchmark # benchmark...

    ./gatk_cpu  # gatk CPU

    ./gatk_fpga # gatk FPGA
    ```
- Build JNI Library libcompute.so
    ```bash
    cd pariHMM && mkdir build
    cd build && cmake ../src/jni/compute
    make
    ```
- Build Gatk.jar
    ```bash
    cd pairHMM/src/gatk
    git clone https://github.com/broadinstitute/gatk
    cd gatk
    git apply ../gatk.patch
    ./gradlew
    cp build/libs/gatk.jar ../../prebuilt/
    ```
## Perforamnce
```bash
# benchmark C1100 cu: 6, group: 11
# cat /proc/cpuinfo  | grep "model name" | head -n 1
model name      : 12th Gen Intel(R) Core(TM) i7-12700K

# ./benchmark gatk.in
INFO: Reading compute.xclbin
Loading: 'compute.xclbin'
181217 pairs need compute!
..........................
109, -1.675846099853515625
109, -1.675838470458984375
109, -1.800930023193359375
109, -1.800930023193359375
109, -1.660717010498046875
109, -1.660717010498046875
fpga cost: 2653.1589ms
cpu  cost: 90585.0078ms
```
## Known issues
- FP32 calculation accuracy error (0.01%)
- Part of the data FP32 calculation overflow
- DSP48 over-occupancy causes wiring difficulties
- C1100 xdma driver does not work properly

## TODO
- Testing Versal Series
  - Hardware's DSP58 performance improvement is huge and requires VCK5000
- Improve Gatk data input efficiency
  - Improve Gatk data input capability and increase parallel computing to facilitate the use of multiple FPGA computing units
- Break OpenCL frequency, scale limits
  - C1100 300/500 Mhz limits the number and frequency of deployed cores
  - XDMA/QDMA implementation performance would be better if the driver is stable
  - Current area footprint is less