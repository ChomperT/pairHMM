# PairHMM Acceleter by FPGA

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
- 可配置计算单元
- 可配置并行规模
- 可配置计算精读
- 可变输入、输出长度
- 多核心并行计算
- 序列长度可定义
- 可根据硬件资源规模灵活配置
- 支持Vitis CL/HLS/SD Accel等不同场景

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
- FP32 计算精度误差（0.01%)
- 部分数据FP32计算溢出
- DSP48 占用过高导致布线困难
- C1100 xdma 驱动无法正常工作

## TODO
- 测试 Versal 系列
  - 硬件的DSP58性能提升巨大，需求VCK5000
- 提高Gatk 数据输入效率
  - 提高Gatk的数据输入能力，增加并行计算，以便于利用多个FPGA计算单元
- 突破OpenCL 频率、规模限制
  - C1100 300/500 Mhz 限制了部署核心的数量和频率
  - 如果驱动稳定，XDMA/QDMA实现性能会更好一些
  - 目前面积占用较少
