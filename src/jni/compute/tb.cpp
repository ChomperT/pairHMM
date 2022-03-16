#include "common.h"

#include <vector>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "compute_cl.hpp"

using namespace std;

struct packet {
	int hs, rs;
	std::vector<u8> rd;
	std::vector<u8> hd;
	std::vector<u16> hls;
};

extern "C" void compute(intype *, precision*, u32 max);

template<typename btype=intype>
void load(const char *path, std::vector<packet> &pkts) {

	auto ifs = std::ifstream(path);
	while (!ifs.eof()) {
		string rss, hss;
		int rs, hs;
		ifs >> rss;
		ifs >> hss;
		rs = atoi(rss.c_str());
		hs = atoi(hss.c_str());

		packet pkt;
		pkt.rs = rs;
		pkt.hs = hs;

#define PAD(a, n)\
	for(u8 q = 0; q < n; q++)\
	a.push_back(0);


#define PLOAD(a, x, j) \
	for(u8 q = 0; q < sizeof(btype); q++) {\
	a.push_back(x[j+q]);\
	}

#define PLOAD2(a, x, j) \
	for(u8 q = 0; q < sizeof(btype); q++) {\
	a.push_back(x[j+q]-33);\
	}

		for (auto i = 0; i < rs; i++) {
			string rb, rq, iq, dq, gq;
			ifs >> rb;
			ifs >> rq;
			ifs >> iq;
			ifs >> dq;
			ifs >> gq;

			btype len = rb.size();
			u8 *d = (u8 *) &len;

			for(u8 m = 0; m < sizeof(btype); m++)
				pkt.rd.push_back(d[m]);

			u32 pn = len % sizeof(btype);

			if(pn != 0) {
				u32 vn = sizeof(btype) - pn;
				PAD(rb, vn);
				PAD(rq, vn);
				PAD(iq, vn);
				PAD(dq, vn);
				PAD(gq, vn);
			}

			for (int j = 0; j < rb.size(); j += sizeof(btype)) {
				PLOAD(pkt.rd, rb, j);
				PLOAD2(pkt.rd, rq, j);
				PLOAD2(pkt.rd, iq, j);
				PLOAD2(pkt.rd, dq, j);
				PLOAD2(pkt.rd, gq, j);
			}
		}

		for (auto i = 0; i < hs; i++) {
			string hss;
			ifs >> hss;
			btype v = hss.size();
			u8 *d = (u8 *) &v;
			pkt.hls.push_back(v);

			u32 pn = v % sizeof(btype);
			if(pn != 0) {
				u32 vn = sizeof(btype) - pn;
				PAD(hss, vn);
			}
			for(u8 m = 0; m < sizeof(btype); m++)
				pkt.hd.push_back(d[m]);

			for (int j = 0; j < hss.size(); j++) {
				pkt.hd.push_back(hss[j]);
			}
		}
		if (rs <= 0 || hs <= 0)
			continue;
		pkts.push_back(pkt);
	}
}

#include <sys/time.h>
#if QT
extern "C" int main(int c, char **argv) {
#else
int main(int argc, char **argv) {
#endif
	std::vector<packet> pkts;
	const char *path = "gatk.in";
	if(argc >= 2) {
		path = argv[1];
	}
	load<intype>(path, pkts);

	int nstop = argc >= 3 ? atoi(argv[2]) : pkts.size();
	int nstart = 0;
	u32 mrsize = 0;
	u32 mhsize = 0;
	for(int i = 0; i < nstop && pkts.size(); i++) {
		auto pkt = pkts[i];
	//for(auto pkt : pkts) {
		mrsize += 16;
		mrsize += pkt.rd.size();
		mrsize += pkt.hd.size();
		mhsize += pkt.hs*pkt.rs;
	}
	mrsize += 16;
	intype *ibuf = nullptr;
	precision *obuf = nullptr;
	Compute clc("compute.xclbin");
	if(mrsize > clc.isize || mhsize > clc.osize) {
		clc.release();
		clc.alloc(mrsize, mhsize*sizeof(precision));
	}
	ibuf = (intype *)clc.ibuf;
	obuf = (precision *)clc.obuf;
	u32 rcur = 0;
	u32 tosize = 0;
	intype *ep = ibuf;
	int cnt = 0;
	constexpr int hdr_len = 16/sizeof(intype);
	for(int i = nstart; i < nstop && i < pkts.size(); i++) {
		auto &pkt = pkts[i];
		intype *rmem = ibuf + rcur; u32 rlen = pkt.rd.size()/sizeof(intype);
		u32 hlen = pkt.hd.size()/sizeof(intype);
		tosize += pkt.rs * pkt.hs;
		rcur += (hdr_len) + rlen + hlen;
		u32 *mh = (u32 *)rmem;

		mh[0] = pkt.rs;
		mh[1]= pkt.hs;
		mh[2] = rlen;
		mh[3] = hlen;

		memcpy(&rmem[hdr_len], pkt.rd.data(), rlen*sizeof(intype));
		memcpy(&rmem[hdr_len+rlen], pkt.hd.data(), hlen*sizeof(intype));
		ep = (ibuf + rcur);
	}
	ep[0] = 0;

	struct timeval s, e, cs, ce;
	printf("%u pairs need compute!\n", tosize);
	gettimeofday(&cs, 0);
	compute(ibuf, obuf, mrsize);
	gettimeofday(&ce, 0);

	gettimeofday(&s, 0);
	clc.compute(ibuf, obuf, mrsize);
	gettimeofday(&e, 0);


	int ci = 0;
	cnt = 0;
	for(auto j = nstart; j < nstop && j < pkts.size(); j++) {
		auto &pkt = pkts[j];
		for(int i = 0; i < pkt.rs*pkt.hs; i++) {
			u32 vn = ci++;
			precision t = obuf[vn];//-LOG10(pkt.hls[i%pkt.hs]);
	//		precision o = LOG10(obuf[vn])-LOG10(pkt.hls[i%pkt.hs]) - LOG10(MAX_VALUE);
			printf("%d, %.18lf\n", j, t);
		}
	}
	float use = (e.tv_usec - s.tv_usec)/1000.0 + (e.tv_sec - s.tv_sec)*1000.0;
	printf("fpga cost: %.4fms\n", use);

	use = (ce.tv_usec - cs.tv_usec)/1000.0 + (ce.tv_sec - cs.tv_sec)*1000.0;
	printf("cpu  cost: %.4fms\n", use);

	return 0;
}
