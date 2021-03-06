
?
Hardware$04135a9f-3d9b-4bea-9cf2-c9e2c3945f68Vitis IDE session Hardware*S"O/home/chomper/workspace/Vivado/Vitis/gene_system_hw_link/Hardware/guidance.html2Q"M/home/chomper/workspace/Vivado/Vitis/gene_system_hw_link/Hardware/guidance.pb *?
????AUTO-FREQ-SCALING-08"i
g
setting\www.xilinx.com"H/cgi-bin/docs/rdoc?t=vitis+guidance;v=2021.2;d=AUTO-FREQ-SCALING-08.htmlJ?
?For clock !%0!, the auto scaled frequency 304.7 MHz exceeds the original specified frequency. The compiler will select the original specified frequency of 300.0 MHz.?
?
?
clk_kernel_00_unbuffered_net2?"?/home/chomper/workspace/Vivado/Vitis/gene_system_hw_link/Hardware/compute.build/link/vivado/vpl/prj/prj.runs/impl_1/dr_timing_summary.rpt
2304.7
2300.0R?
?The automatic frequency scaling feature allows user kernels to operate in hardware, even if at a lower frequency than intended. In this case the clock may in fact be able to run at a higher frequency than specified. You may want to consider !URI%1! the clock frequency higher for better performance. The '--kernel_frequency' option is one way to control the frequency specification.k
ig
setting\www.xilinx.com"H/cgi-bin/docs/rdoc?t=vitis+guidance;v=2021.2;d=AUTO-FREQ-SCALING-08.htmlZ AcceleratorZcomputeZPerformanceh 
????AUTO-FREQ-SCALING-08"i
g
setting\www.xilinx.com"H/cgi-bin/docs/rdoc?t=vitis+guidance;v=2021.2;d=AUTO-FREQ-SCALING-08.htmlJ?
?For clock !%0!, the auto scaled frequency 464.4 MHz exceeds the original specified frequency. The compiler will select the original specified frequency of 450.0 MHz.?
?
?
hbm_aclk2?"?/home/chomper/workspace/Vivado/Vitis/gene_system_hw_link/Hardware/compute.build/link/vivado/vpl/prj/prj.runs/impl_1/dr_timing_summary.rpt
2464.4
2450.0R?
?The automatic frequency scaling feature allows user kernels to operate in hardware, even if at a lower frequency than intended. In this case the clock may in fact be able to run at a higher frequency than specified. You may want to consider !URI%1! the clock frequency higher for better performance. The '--kernel_frequency' option is one way to control the frequency specification.k
ig
setting\www.xilinx.com"H/cgi-bin/docs/rdoc?t=vitis+guidance;v=2021.2;d=AUTO-FREQ-SCALING-08.htmlZ AcceleratorZcomputeZPerformanceh :?	
X?S?	SYSLINK-1 BA
0Created top level block diagram design dr.bd.tcl
2	dr.bd.tcl
????PLATFORM-CLOCK-DOMAINS-01?
?
automatic frequency scalingawww.xilinx.com"M/cgi-bin/docs/rdoc?t=vitis+guidance;v=2021.2;d=PLATFORM-CLOCK-DOMAINS-01.html*= or >B?
?The compiler selected the following frequencies for the runtime controllable kernel clock(s) and scalable system clock(s): 
Kernel: ulp_ucs/aclk_kernel_01 = 500.0 MHz 
Kernel: ulp_ucs/aclk_kernel_00 = 300.0 MHz 
System: hbm_aclk = 450.0 MHz 
Scalable clock ulp_ucs/aclk_kernel_01 (Id = 1) is used for rtl kernels. This design has 0 rtl kernel(s).
Scalable clock ulp_ucs/aclk_kernel_00 (Id = 0) is used for hls kernels. This design has 1 hls kernel(s).?
?2?
Kernel: ulp_ucs/aclk_kernel_01 = 500.0 MHz 
Kernel: ulp_ucs/aclk_kernel_00 = 300.0 MHz 
System: hbm_aclk = 450.0 MHz 
Scalable clock ulp_ucs/aclk_kernel_01 (Id = 1) is used for rtl kernels. This design has 0 rtl kernel(s).
Scalable clock ulp_ucs/aclk_kernel_00 (Id = 0) is used for hls kernels. This design has 1 hls kernel(s).Jm
kThe !URI%1! feature allows user kernels to operate in hardware, even if at a lower frequency than intended.R AcceleratorR
SystemRPerformanceB?
?
AUTO-FREQ-SCALING-08?
AUTO-FREQ-SCALING-082Auto frequency scaling - Higher frequency possiblesdx"?For clock %REF, the auto scaled frequency %s MHz exceeds the original specified frequency. The compiler will select the original specified frequency of %s MHz.: Accelerator:
xclbin:PerformanceB?The automatic frequency scaling feature allows user kernels to operate in hardware, even if at a lower frequency than intended. In this case the clock may in fact be able to run at a higher frequency than specified. You may want to consider !URI%1! the clock frequency higher for better performance. The '--kernel_frequency' option is one way to control the frequency specification.Ji
g
setting\www.xilinx.com"H/cgi-bin/docs/rdoc?t=vitis+guidance;v=2021.2;d=AUTO-FREQ-SCALING-08.htmlP?  ?? J?	
u
	SYSLINK-1h
	SYSLINK-1!system_link Top Level BD Creationsystem_link"+Created top level block diagram design %STR
?
PLATFORM-CLOCK-DOMAINS-01?
PLATFORM-CLOCK-DOMAINS-01CRuntime controllable clock domains - Achieved clock frequency (MHz)sdx")One or more clocks failed a timing check.:
System:Performance: AcceleratorBkThe !URI%1! feature allows user kernels to operate in hardware, even if at a lower frequency than intended.J?
?
automatic frequency scalingawww.xilinx.com"M/cgi-bin/docs/rdoc?t=vitis+guidance;v=2021.2;d=PLATFORM-CLOCK-DOMAINS-01.htmlP?= or >?}The compiler selected the following frequencies for the runtime controllable kernel clock(s) and scalable system clock(s): %s?? :	text/htmlBcharset=UTF-8J?<html> Kernel clocks (and system clocks for some platforms) are scalable; they can preserve functionality at the cost of performance by running at a lower frequency than requested. To be scalable, a clock must be driven by an MMCM where the control registers for the MMCM can be set by the runtime over AXI4-Lite. This item shows the final runtime controlled frequencies for the scalable clocks.</html>?  ? ? 