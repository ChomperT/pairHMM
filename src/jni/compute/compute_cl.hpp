#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "common.h"
#define OCL_CHECK(error, call)                                                                   \
	call;                                                                                        \
	if (error != CL_SUCCESS) {                                                                   \
		printf("%s:%d Error calling " #call ", error code is: %d\n", __FILE__, __LINE__, error); \
		exit(EXIT_FAILURE);                                                                      \
	}

#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1

#include <CL/cl2.hpp>

//Customized buffer allocation for 4K boundary alignment
template <typename T>
struct aligned_allocator
{
  using value_type = T;
  T* allocate(std::size_t num)
  {
	void* ptr = nullptr;
	if (posix_memalign(&ptr,4096,num*sizeof(T)))
	  throw std::bad_alloc();
	return reinterpret_cast<T*>(ptr);
  }
  void deallocate(T* p, std::size_t num)
  {
	free(p);
  }
};

static const int DATA_SIZE = 4096*1024;

static const std::string error_message =
	"Error: Result mismatch:\n"
	"i = %d CPU result = %d Device result = %d\n";


class Compute {
public:
	Compute(const char *file) {
		bool found_device = false;
		std::string xclbinFilename = file;
		cl::Platform::get(&platforms);
		for(size_t i = 0; (i < platforms.size() ) & (found_device == false) ;i++){
			cl::Platform platform = platforms[i];
			std::string platformName = platform.getInfo<CL_PLATFORM_NAME>();
			if ( platformName == "Xilinx"){
				devices.clear();
				platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices);
			if (devices.size()){
				device = devices[0];
				found_device = true;
				break;
			}
			}
		}
		if (found_device == false){
		   std::cout << "Error: Unable to find Target Device "
			   << device.getInfo<CL_DEVICE_NAME>() << std::endl;
		   return;
		}

		// Creating Context and Command Queue for selected device
		OCL_CHECK(err, context = cl::Context(device, NULL, NULL, NULL, &err));

		OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err));
		//OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err));

		std::cout << "INFO: Reading " << xclbinFilename << std::endl;
		FILE* fp;
		if ((fp = fopen(xclbinFilename.c_str(), "r")) == nullptr) {
			printf("ERROR: %s xclbin not available please build\n", xclbinFilename.c_str());
			exit(EXIT_FAILURE);
		}
		// Load xclbin
		std::cout << "Loading: '" << xclbinFilename << "'\n";
		std::ifstream bin_file(xclbinFilename, std::ifstream::binary);
		bin_file.seekg (0, bin_file.end);
		unsigned nb = bin_file.tellg();
		bin_file.seekg (0, bin_file.beg);
		char *buf = new char [nb];
		bin_file.read(buf, nb);

		// Creating Program from Binary File
		cl::Program::Binaries bins;
		bins.push_back({buf,nb});
		devices.resize(1);
		OCL_CHECK(err, program = cl::Program(context, devices, bins, NULL, &err));

		// This call will get the kernel object from program. A kernel is an
		// OpenCL function that is executed on the FPGA.
		OCL_CHECK(err, krnl_compute = cl::Kernel(program,"compute", &err));
		this->isize = isize = 0;
		this->osize = osize = 0;
		imem = 0;
		omem = 0;
		ibuf = nullptr;
		obuf = nullptr;
	}
	void alloc(int input_size, int output_size) {
		OCL_CHECK(err, imem = cl::Buffer(context, CL_MEM_READ_ONLY , input_size, 0, &err));
		OCL_CHECK(err, omem = cl::Buffer(context, CL_MEM_WRITE_ONLY, output_size, 0, &err));
		OCL_CHECK(err, ibuf = q.enqueueMapBuffer (imem, CL_TRUE , CL_MAP_WRITE, 0, input_size, NULL, NULL, &err));
		OCL_CHECK(err, obuf = q.enqueueMapBuffer (omem, CL_TRUE , CL_MAP_READ , 0, output_size, NULL, NULL, &err));
		this->isize = input_size;
		this->osize = output_size;
	}
	void release() {
		if(ibuf || obuf) {
			q.enqueueUnmapMemObject(imem, ibuf);
			q.enqueueUnmapMemObject(omem, obuf);
			q.finish();
		}
	}
	void compute(intype *input, precision* output, u32 max_size){

		krnl_compute.setArg(0, imem);
		krnl_compute.setArg(1, omem);
		krnl_compute.setArg(2, max_size*2);

		q.enqueueMigrateMemObjects({imem}, 0/* 0 means from host*/);
		q.enqueueTask(krnl_compute);
		q.enqueueMigrateMemObjects({omem}, CL_MIGRATE_MEM_OBJECT_HOST/* 0 means from host*/);
		q.finish();

	}
	~Compute() {

	}

public:
	void *ibuf, *obuf;
	int isize, osize;

private:
	std::vector<cl::Device> devices;
	cl::Device device;
	cl_int err;
	cl::Context context;
	cl::CommandQueue q;
	cl::Kernel krnl_compute;
	cl::Program program;
	cl::Buffer imem, omem;
	std::vector<cl::Platform> platforms;
};
