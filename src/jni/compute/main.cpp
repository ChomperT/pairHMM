#include <stdlib.h>
#include <fstream>
#include <iostream>

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

int main(int argc, char* argv[]) {

	//TARGET_DEVICE macro needs to be passed from gcc command line
	if(argc != 2) {
		std::cout << "Usage: " << argv[0] <<" <xclbin>" << std::endl;
		return EXIT_FAILURE;
	}

	std::string xclbinFilename = argv[1];

	// Compute the size of array in bytes
	size_t size_in_bytes = DATA_SIZE * sizeof(int);

	// Creates a vector of DATA_SIZE elements with an initial value of 10 and 32
	// using customized allocator for getting buffer alignment to 4k boundary

	std::vector<cl::Device> devices;
	cl::Device device;
	cl_int err;
	cl::Context context;
	cl::CommandQueue q;
	cl::Kernel krnl_vector_add;
	cl::Program program;
	std::vector<cl::Platform> platforms;
	bool found_device = false;

	//traversing all Platforms To find Xilinx Platform and targeted
	//Device in Xilinx Platform
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
	   return EXIT_FAILURE;
	}

	// Creating Context and Command Queue for selected device
	OCL_CHECK(err, context = cl::Context(device, NULL, NULL, NULL, &err));
	OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err));

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
	OCL_CHECK(err, krnl_vector_add = cl::Kernel(program,"krnl_vadd", &err));

	// These commands will allocate memory on the Device. The cl::Buffer objects can
	// be used to reference the memory locations on the device.
	OCL_CHECK(err, cl::Buffer buffer_a(context, CL_MEM_READ_ONLY, size_in_bytes, NULL, &err));
	OCL_CHECK(err, cl::Buffer buffer_b(context, CL_MEM_READ_ONLY, size_in_bytes, NULL, &err));
	OCL_CHECK(err, cl::Buffer buffer_result(context, CL_MEM_WRITE_ONLY, size_in_bytes, NULL, &err));

	//set the kernel Arguments
	int narg=0;
	OCL_CHECK(err, err = krnl_vector_add.setArg(narg++,buffer_a));
	OCL_CHECK(err, err = krnl_vector_add.setArg(narg++,buffer_b));
	OCL_CHECK(err, err = krnl_vector_add.setArg(narg++,buffer_result));
	OCL_CHECK(err, err = krnl_vector_add.setArg(narg++,DATA_SIZE));

	//We then need to map our OpenCL buffers to get the pointers
	int *ptr_a;
	int *ptr_b;
	int *ptr_result;
	OCL_CHECK(err, ptr_a = (int*)q.enqueueMapBuffer (buffer_a , CL_TRUE , CL_MAP_WRITE , 0, size_in_bytes, NULL, NULL, &err));
	OCL_CHECK(err, ptr_b = (int*)q.enqueueMapBuffer (buffer_b , CL_TRUE , CL_MAP_WRITE , 0, size_in_bytes, NULL, NULL, &err));
	OCL_CHECK(err, ptr_result = (int*)q.enqueueMapBuffer (buffer_result , CL_TRUE , CL_MAP_READ , 0, size_in_bytes, NULL, NULL, &err));

	// Data will be migrated to kernel space
	OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_a,buffer_b},0/* 0 means from host*/));

	//Launch the Kernel
	OCL_CHECK(err, err = q.enqueueTask(krnl_vector_add));

	// The result of the previous kernel execution will need to be retrieved in
	// order to view the results. This call will transfer the data from FPGA to
	// source_results vector
	OCL_CHECK(err, q.enqueueMigrateMemObjects({buffer_result},CL_MIGRATE_MEM_OBJECT_HOST));

	OCL_CHECK(err, q.finish());

	//Verify the result
	int match = 0;
	for (int i = 0; i < DATA_SIZE; i++) {
		int host_result = ptr_a[i] + ptr_b[i];
		if (ptr_result[i] != host_result) {
			printf(error_message.c_str(), i, host_result, ptr_result[i]);
			match = 1;
			break;
		}
	}

	OCL_CHECK(err, err = q.enqueueUnmapMemObject(buffer_a , ptr_a));
	OCL_CHECK(err, err = q.enqueueUnmapMemObject(buffer_b , ptr_b));
	OCL_CHECK(err, err = q.enqueueUnmapMemObject(buffer_result , ptr_result));
	OCL_CHECK(err, err = q.finish());

	std::cout << "TEST " << (match ? "FAILED" : "PASSED") << std::endl;
	return (match ? EXIT_FAILURE :  EXIT_SUCCESS);

}
