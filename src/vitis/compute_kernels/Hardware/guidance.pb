
?
Hardware$df5a4716-f716-4ef8-b4ab-23ee55d21b75Vitis IDE session Hardware*O"K/home/chomper/workspace/Vivado/Vitis/compute_kernels/Hardware/guidance.html2M"I/home/chomper/workspace/Vivado/Vitis/compute_kernels/Hardware/guidance.pb$*??
????	Interface"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?Multiple burst reads of variable length and bit width 64 in loop 'loop_hr'(!%1!) has been inferred on bundle 'gmem0'. These burst requests might be further partitioned into multiple requests during RTL generation, based on max_read_burst_length or max_write_burst_length settings. (!%2!)
?
?2?Multiple burst reads of variable length and bit width 64 in loop 'loop_hr'(%REF) has been inferred on bundle 'gmem0'. These burst requests might be further partitioned into multiple requests during RTL generation, based on max_read_burst_length or max_write_burst_length settings. (%REF)

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=222
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=222Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;ZcomputeZ	InterfaceZ Acceleratorh 
????	Interface"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?Multiple burst reads of variable length and bit width 64 in loop 'loop_rr'(!%1!) has been inferred on bundle 'gmem0'. These burst requests might be further partitioned into multiple requests during RTL generation, based on max_read_burst_length or max_write_burst_length settings. (!%2!)
?
?2?Multiple burst reads of variable length and bit width 64 in loop 'loop_rr'(%REF) has been inferred on bundle 'gmem0'. These burst requests might be further partitioned into multiple requests during RTL generation, based on max_read_burst_length or max_write_burst_length settings. (%REF)

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=253
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=253Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;Z AcceleratorZcomputeZ	Interfaceh 
????Latency"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
dCannot unroll loop 'VITIS_LOOP_269_2' (!%1!) in function 'compute' completely: variable loop bound.
?
f2dCannot unroll loop 'VITIS_LOOP_269_2' (%REF) in function 'compute' completely: variable loop bound.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=269Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;ZcomputeZLatencyZ Acceleratorh 
????Latency"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?Cannot flatten loop 'loop_calc' in function 'cu<float, float, unsigned char, ap_uint<10> >' the outer loop is not a perfect loop because either the parent loop or the sub loop has no computeable trip count.
?
?2?Cannot flatten loop 'loop_calc' in function 'cu<float, float, unsigned char, ap_uint<10> >' the outer loop is not a perfect loop because either the parent loop or the sub loop has no computeable trip count.
R?
!%0!?
??
5See here for more help on vitis_hls 200-960 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-960.htmlZ AcceleratorZcomputeZLatencyh 
????Latency"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?Cannot flatten loop 'loop_hap' (!%1!) in function 'compute' the outer loop is not a perfect loop because there is nontrivial logic before entering the inner loop.
?
?2?Cannot flatten loop 'loop_hap' (%REF) in function 'compute' the outer loop is not a perfect loop because there is nontrivial logic before entering the inner loop.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=442R?
!%0!?
??
5See here for more help on vitis_hls 200-960 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-960.htmlZ AcceleratorZcomputeZLatencyh 
????Latency"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
aCannot flatten loop 'loop_rr' (!%1!) in function 'compute' the outer loop is not a perfect loop.
?
c2aCannot flatten loop 'loop_rr' (%REF) in function 'compute' the outer loop is not a perfect loop.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=253R?
!%0!?
??
5See here for more help on vitis_hls 200-960 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-960.htmlZ AcceleratorZcomputeZLatencyh 
????Latency"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
UCannot flatten loop 'loop_read' (!%1!) in function 'compute' more than one sub loop.
?
W2UCannot flatten loop 'loop_read' (%REF) in function 'compute' more than one sub loop.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=438R?
!%0!?
??
5See here for more help on vitis_hls 200-960 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-960.htmlZ AcceleratorZcomputeZLatencyh 
????Latency"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
UCannot flatten loop 'loop_read' (!%1!) in function 'compute' more than one sub loop.
?
W2UCannot flatten loop 'loop_read' (%REF) in function 'compute' more than one sub loop.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=436R?
!%0!?
??
5See here for more help on vitis_hls 200-960 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-960.htmlZ AcceleratorZcomputeZLatencyh 
????Latency"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?Cannot flatten loop 'loop_main' (!%1!) in function 'compute' the outer loop is not a perfect loop because there is nontrivial logic before entering the inner loop.
?
?2?Cannot flatten loop 'loop_main' (%REF) in function 'compute' the outer loop is not a perfect loop because there is nontrivial logic before entering the inner loop.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=416R?
!%0!?
??
5See here for more help on vitis_hls 200-960 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-960.htmlZ AcceleratorZcomputeZLatencyh 
?	??	?Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'compute_Pipeline_loop_hr' (loop 'loop_hr'): Unable to schedule 'store' operation ('hp_6_addr_write_ln229', !%1!) of variable 'trunc_ln229', !%2! on array 'hp_6' due to limited memory ports (II = 1). Please consider using a memory core with more ports or partitioning the array 'hp_6'.
?
?2?The II Violation in module 'compute_Pipeline_loop_hr' (loop 'loop_hr'): Unable to schedule 'store' operation ('hp_6_addr_write_ln229', %REF) of variable 'trunc_ln229', %REF on array 'hp_6' due to limited memory ports (II = 1). Please consider using a memory core with more ports or partitioning the array 'hp_6'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=229
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=229R?
!%0!?
??
5See here for more help on vitis_hls 200-885 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-885.htmlZ AcceleratorZcomputeZ
Kernelh 
?	??	?Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'compute_Pipeline_loop_hr' (loop 'loop_hr'): Unable to schedule 'store' operation ('hp_6_addr_1_write_ln229', !%1!) of variable 'trunc_ln229_1', !%2! on array 'hp_6' due to limited memory ports (II = 2). Please consider using a memory core with more ports or partitioning the array 'hp_6'.
?
?2?The II Violation in module 'compute_Pipeline_loop_hr' (loop 'loop_hr'): Unable to schedule 'store' operation ('hp_6_addr_1_write_ln229', %REF) of variable 'trunc_ln229_1', %REF on array 'hp_6' due to limited memory ports (II = 2). Please consider using a memory core with more ports or partitioning the array 'hp_6'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=229
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=229R?
!%0!?
??
5See here for more help on vitis_hls 200-885 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-885.htmlZ AcceleratorZcomputeZ
Kernelh 
?	??	?Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'compute_Pipeline_loop_hr' (loop 'loop_hr'): Unable to schedule 'store' operation ('hp_6_addr_2_write_ln229', !%1!) of variable 'trunc_ln229_2', !%2! on array 'hp_6' due to limited memory ports (II = 3). Please consider using a memory core with more ports or partitioning the array 'hp_6'.
?
?2?The II Violation in module 'compute_Pipeline_loop_hr' (loop 'loop_hr'): Unable to schedule 'store' operation ('hp_6_addr_2_write_ln229', %REF) of variable 'trunc_ln229_2', %REF on array 'hp_6' due to limited memory ports (II = 3). Please consider using a memory core with more ports or partitioning the array 'hp_6'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=229
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=229R?
!%0!?
??
5See here for more help on vitis_hls 200-885 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-885.htmlZ AcceleratorZcomputeZ
Kernelh 
?	??	?Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'compute_Pipeline_loop_hr' (loop 'loop_hr'): Unable to schedule 'store' operation ('hp_6_addr_3_write_ln229', !%1!) of variable 'trunc_ln229_3', !%2! on array 'hp_6' due to limited memory ports (II = 4). Please consider using a memory core with more ports or partitioning the array 'hp_6'.
?
?2?The II Violation in module 'compute_Pipeline_loop_hr' (loop 'loop_hr'): Unable to schedule 'store' operation ('hp_6_addr_3_write_ln229', %REF) of variable 'trunc_ln229_3', %REF on array 'hp_6' due to limited memory ports (II = 4). Please consider using a memory core with more ports or partitioning the array 'hp_6'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=229
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=229R?
!%0!?
??
5See here for more help on vitis_hls 200-885 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-885.htmlZ AcceleratorZcomputeZ
Kernelh 
?	??	?Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'compute_Pipeline_loop_hr' (loop 'loop_hr'): Unable to schedule 'store' operation ('hp_6_addr_6_write_ln229', !%1!) of variable 'trunc_ln229_6', !%2! on array 'hp_6' due to limited memory ports (II = 7). Please consider using a memory core with more ports or partitioning the array 'hp_6'.
?
?2?The II Violation in module 'compute_Pipeline_loop_hr' (loop 'loop_hr'): Unable to schedule 'store' operation ('hp_6_addr_6_write_ln229', %REF) of variable 'trunc_ln229_6', %REF on array 'hp_6' due to limited memory ports (II = 7). Please consider using a memory core with more ports or partitioning the array 'hp_6'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=229
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=229R?
!%0!?
??
5See here for more help on vitis_hls 200-885 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-885.htmlZcomputeZ
KernelZ Acceleratorh 
????
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
MPipelining result : Target II = NA, Final II = 8, Depth = 10, loop 'loop_hr'
Q
O2MPipelining result : Target II = NA, Final II = 8, Depth = 10, loop 'loop_hr'
Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;ZcomputeZ
ThroughputZ Acceleratorh 
????
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
UPipelining result : Target II = NA, Final II = 1, Depth = 1, loop 'VITIS_LOOP_241_1'
Y
W2UPipelining result : Target II = NA, Final II = 1, Depth = 1, loop 'VITIS_LOOP_241_1'
Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;Z AcceleratorZcomputeZ
Throughputh 
????Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'compute_Pipeline_VITIS_LOOP_269_2' (loop 'VITIS_LOOP_269_2'): Unable to schedule 'store' operation ('rb_6_addr_7_write_ln274', !%1!) of variable 'trunc_ln13_read' on array 'rb_6' due to limited memory ports (II = 1). Please consider using a memory core with more ports or partitioning the array 'rb_6'.
?
?2?The II Violation in module 'compute_Pipeline_VITIS_LOOP_269_2' (loop 'VITIS_LOOP_269_2'): Unable to schedule 'store' operation ('rb_6_addr_7_write_ln274', %REF) of variable 'trunc_ln13_read' on array 'rb_6' due to limited memory ports (II = 1). Please consider using a memory core with more ports or partitioning the array 'rb_6'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=274R?
!%0!?
??
5See here for more help on vitis_hls 200-885 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-885.htmlZ AcceleratorZcomputeZ
Kernelh 
????Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'compute_Pipeline_VITIS_LOOP_269_2' (loop 'VITIS_LOOP_269_2'): Unable to schedule 'store' operation ('rb_6_addr_6_write_ln274', !%1!) of variable 'trunc_ln269_5_read' on array 'rb_6' due to limited memory ports (II = 2). Please consider using a memory core with more ports or partitioning the array 'rb_6'.
?
?2?The II Violation in module 'compute_Pipeline_VITIS_LOOP_269_2' (loop 'VITIS_LOOP_269_2'): Unable to schedule 'store' operation ('rb_6_addr_6_write_ln274', %REF) of variable 'trunc_ln269_5_read' on array 'rb_6' due to limited memory ports (II = 2). Please consider using a memory core with more ports or partitioning the array 'rb_6'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=274R?
!%0!?
??
5See here for more help on vitis_hls 200-885 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-885.htmlZ AcceleratorZcomputeZ
Kernelh 
????Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'compute_Pipeline_VITIS_LOOP_269_2' (loop 'VITIS_LOOP_269_2'): Unable to schedule 'store' operation ('rb_6_addr_5_write_ln274', !%1!) of variable 'trunc_ln269_read' on array 'rb_6' due to limited memory ports (II = 3). Please consider using a memory core with more ports or partitioning the array 'rb_6'.
?
?2?The II Violation in module 'compute_Pipeline_VITIS_LOOP_269_2' (loop 'VITIS_LOOP_269_2'): Unable to schedule 'store' operation ('rb_6_addr_5_write_ln274', %REF) of variable 'trunc_ln269_read' on array 'rb_6' due to limited memory ports (II = 3). Please consider using a memory core with more ports or partitioning the array 'rb_6'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=274R?
!%0!?
??
5See here for more help on vitis_hls 200-885 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-885.htmlZcomputeZ
KernelZ Acceleratorh 
????Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'compute_Pipeline_VITIS_LOOP_269_2' (loop 'VITIS_LOOP_269_2'): Unable to schedule 'store' operation ('rb_6_addr_4_write_ln274', !%1!) of variable 'trunc_ln269_14_read' on array 'rb_6' due to limited memory ports (II = 4). Please consider using a memory core with more ports or partitioning the array 'rb_6'.
?
?2?The II Violation in module 'compute_Pipeline_VITIS_LOOP_269_2' (loop 'VITIS_LOOP_269_2'): Unable to schedule 'store' operation ('rb_6_addr_4_write_ln274', %REF) of variable 'trunc_ln269_14_read' on array 'rb_6' due to limited memory ports (II = 4). Please consider using a memory core with more ports or partitioning the array 'rb_6'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=274R?
!%0!?
??
5See here for more help on vitis_hls 200-885 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-885.htmlZ AcceleratorZcomputeZ
Kernelh 
????Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'compute_Pipeline_VITIS_LOOP_269_2' (loop 'VITIS_LOOP_269_2'): Unable to schedule 'store' operation ('rb_6_addr_1_write_ln274', !%1!) of variable 'trunc_ln269_29_read' on array 'rb_6' due to limited memory ports (II = 7). Please consider using a memory core with more ports or partitioning the array 'rb_6'.
?
?2?The II Violation in module 'compute_Pipeline_VITIS_LOOP_269_2' (loop 'VITIS_LOOP_269_2'): Unable to schedule 'store' operation ('rb_6_addr_1_write_ln274', %REF) of variable 'trunc_ln269_29_read' on array 'rb_6' due to limited memory ports (II = 7). Please consider using a memory core with more ports or partitioning the array 'rb_6'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=274R?
!%0!?
??
5See here for more help on vitis_hls 200-885 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-885.htmlZ AcceleratorZcomputeZ
Kernelh 
????
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
UPipelining result : Target II = NA, Final II = 8, Depth = 8, loop 'VITIS_LOOP_269_2'
Y
W2UPipelining result : Target II = NA, Final II = 8, Depth = 8, loop 'VITIS_LOOP_269_2'
Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;Z AcceleratorZcomputeZ
Throughputh 
????
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
UPipelining result : Target II = NA, Final II = 1, Depth = 1, loop 'VITIS_LOOP_187_1'
Y
W2UPipelining result : Target II = NA, Final II = 1, Depth = 1, loop 'VITIS_LOOP_187_1'
Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;Z AcceleratorZcomputeZ
Throughputh 
???
?
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 1, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', !%1!) of variable 'y', !%2! on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', !%3!) on local variable 'yp_11_out'.
?
?2?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 1, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', %REF) of variable 'y', %REF on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', %REF) on local variable 'yp_11_out'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=133R?
!%0!?
??
5See here for more help on vitis_hls 200-880 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-880.htmlZcomputeZ
ThroughputZ Acceleratorh 
???
?
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 2, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', !%1!) of variable 'y', !%2! on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', !%3!) on local variable 'yp_11_out'.
?
?2?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 2, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', %REF) of variable 'y', %REF on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', %REF) on local variable 'yp_11_out'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=133R?
!%0!?
??
5See here for more help on vitis_hls 200-880 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-880.htmlZ AcceleratorZcomputeZ
Throughputh 
???
?
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 3, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', !%1!) of variable 'y', !%2! on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', !%3!) on local variable 'yp_11_out'.
?
?2?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 3, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', %REF) of variable 'y', %REF on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', %REF) on local variable 'yp_11_out'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=133R?
!%0!?
??
5See here for more help on vitis_hls 200-880 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-880.htmlZ AcceleratorZcomputeZ
Throughputh 
???
?
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 4, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', !%1!) of variable 'y', !%2! on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', !%3!) on local variable 'yp_11_out'.
?
?2?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 4, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', %REF) of variable 'y', %REF on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', %REF) on local variable 'yp_11_out'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=133R?
!%0!?
??
5See here for more help on vitis_hls 200-880 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-880.htmlZ AcceleratorZcomputeZ
Throughputh 
???
?
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 7, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', !%1!) of variable 'y', !%2! on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', !%3!) on local variable 'yp_11_out'.
?
?2?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 7, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', %REF) of variable 'y', %REF on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', %REF) on local variable 'yp_11_out'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=133R?
!%0!?
??
5See here for more help on vitis_hls 200-880 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-880.htmlZ AcceleratorZcomputeZ
Throughputh 
???
?
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 9, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', !%1!) of variable 'y', !%2! on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', !%3!) on local variable 'yp_11_out'.
?
?2?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 9, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', %REF) of variable 'y', %REF on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', %REF) on local variable 'yp_11_out'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=133R?
!%0!?
??
5See here for more help on vitis_hls 200-880 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-880.htmlZcomputeZ
ThroughputZ Acceleratorh 
????
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 10, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', !%1!) of variable 'y', !%2! on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', !%3!) on local variable 'yp_11_out'.
?
?2?The II Violation in module 'cu_float_float_unsigned_char_ap_uint_10_Pipeline_loop_calc_sub' (loop 'loop_calc_sub'): Unable to enforce a carried dependence constraint (II = 10, distance = 1, offset = 0) between 'store' operation ('yp_11_out_write_ln140', %REF) of variable 'y', %REF on local variable 'yp_11_out' and 'load' operation ('yp_11_out_load', %REF) on local variable 'yp_11_out'.

_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=140
_
]gatk.cpp2O"A/home/chomper/workspace/Vivado/Vitis/compute_kernels/src/gatk.cpp2line=133R?
!%0!?
??
5See here for more help on vitis_hls 200-880 guidance.Mwww.xilinx.com"9/cgi-bin/docs/rdoc?v=2021.2;t=hls+guidance;d=200-880.htmlZ AcceleratorZcomputeZ
Throughputh 
????
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
UPipelining result : Target II = NA, Final II = 11, Depth = 170, loop 'loop_calc_sub'
Y
W2UPipelining result : Target II = NA, Final II = 11, Depth = 170, loop 'loop_calc_sub'
Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;Z AcceleratorZcomputeZ
Throughputh 
????
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
VPipelining result : Target II = NA, Final II = 1, Depth = 39, loop 'VITIS_LOOP_203_3'
Z
X2VPipelining result : Target II = NA, Final II = 1, Depth = 39, loop 'VITIS_LOOP_203_3'
Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;ZcomputeZ
ThroughputZ Acceleratorh 
????
Throughput"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
NPipelining result : Target II = NA, Final II = 1, Depth = 71, loop 'loop_out'
R
P2NPipelining result : Target II = NA, Final II = 1, Depth = 71, loop 'loop_out'
Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;ZcomputeZ
ThroughputZ Acceleratorh 
????	Interface"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
xDesign has inferred MAXI bursts and missed bursts, see Vitis HLS GUI synthesis summary report for detailed information.
|
z2xDesign has inferred MAXI bursts and missed bursts, see Vitis HLS GUI synthesis summary report for detailed information.
Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;Z AcceleratorZcomputeZ	Interfaceh 
????Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;J?
F**** Loop Constraint Status: All loop constraints were NOT satisfied.
J
H2F**** Loop Constraint Status: All loop constraints were NOT satisfied.
Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;Z AcceleratorZcomputeZ
Kernelh 
????Kernel"j
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;JH
 **** Estimated Fmax: 376.96 MHz
$
"2 **** Estimated Fmax: 376.96 MHz
Rp
 l
jh
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;Z AcceleratorZcomputeZ
Kernelh B?
?
	Interface?
	InterfaceHLS Interface RelatedHLS"%s: Accelerator:
Kernel:	InterfaceB Jj
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;P?  ?? 
?
Kernel?
KernelHLS Kernel RelatedHLS"%s: Accelerator:
Kernel:
KernelB Jj
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;P?? ?  
?

Throughput?

ThroughputHLS Throughput RelatedHLS"%s: Accelerator:
Kernel:
ThroughputB Jj
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;P?  ?? 
?
Latency?
LatencyHLS Latency RelatedHLS"%s: Accelerator:
Kernel:LatencyB Jj
h
Vitis HLS User Guide (UG1399)Gwww.xilinx.com"3/cgi-bin/docs/rdoc?v=2021.2;d=ug1399-vitis-hls.pdf;P?  ?? 