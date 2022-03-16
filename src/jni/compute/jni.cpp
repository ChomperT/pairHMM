#include <jni.h>
#include <float.h>
#include <stdlib.h>

#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <sys/time.h>

#include "common.h"
#include "compute_cl.hpp"

struct packet {
	int hs, rs;
	std::vector<u8> rd;
	std::vector<u8> hd;
	std::vector<u16> hls;
};
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
template<typename btype=intype>
static void load_base(packet &pkt, std::string &rb, std::string &rq, std::string &iq, std::string &dq, std::string &gq) {
	btype len = rb.size();
	u8 *d = (u8 *)&len;

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
		PLOAD(pkt.rd, rq, j);
		PLOAD(pkt.rd, iq, j);
		PLOAD(pkt.rd, dq, j);
		PLOAD(pkt.rd, gq, j);
	}
}

template<typename btype=intype>
static void load_hap(packet &pkt, std::string &hp) {
	btype v = hp.size();
	pkt.hls.push_back(v);

	u8* d = (u8 *) &v;

	u32 pn = v % sizeof(btype);
	if(pn != 0) {
		u32 vn = sizeof(btype) - pn;
		PAD(hp, vn);
	}
	for(u8 m = 0; m < sizeof(btype); m++)
		pkt.hd.push_back(d[m]);

	for (int j = 0; j < hp.size(); j++) {
		pkt.hd.push_back(hp[j]);
	}
}

static jfieldID _rb, _rq, _iq, _dq, _gq, _hp;
static Compute *compute;
extern "C" JNIEXPORT void JNICALL Java_chomper_gatk_cfgengine_CFGEngine_initNative(JNIEnv *env, jclass clazz, jclass RC, jclass HC) {
	(void)clazz;
	_rb = env->GetFieldID(RC, "readBases", "[B");
	_rq = env->GetFieldID(RC, "readQuals", "[B");
	_iq = env->GetFieldID(RC, "insertionGOP", "[B");
	_dq = env->GetFieldID(RC, "deletionGOP", "[B");
	_gq = env->GetFieldID(RC, "overallGCP", "[B");
	_hp = env->GetFieldID(HC, "haplotypeBases", "[B");
	compute = new Compute("compute.xclbin");
	printf("hello init!\n");
}


static int rln, hln;
extern "C" JNIEXPORT void JNICALL Java_chomper_gatk_cfgengine_CFGEngine_computeNative(JNIEnv *env, jclass clazz, jobjectArray rbs, jobjectArray hbs, jdoubleArray out) {
	(void)clazz;
	u32 hs, rs;
	struct timeval st, ed;
	gettimeofday(&st, 0);
	rs = env->GetArrayLength(rbs);
	hs = env->GetArrayLength(hbs);
	packet pkt;
	//	printf("hello compute!: %d, %d\n", rs, hs);
	for(int i = 0; i < rs; i++) {
		auto rd = env->GetObjectArrayElement(rbs, i);
		auto rb = (jbyteArray)env->GetObjectField(rd, _rb);
		auto rq = (jbyteArray)env->GetObjectField(rd, _rq);
		auto iq = (jbyteArray)env->GetObjectField(rd, _iq);
		auto dq = (jbyteArray)env->GetObjectField(rd, _dq);
		auto gq = (jbyteArray)env->GetObjectField(rd, _gq);
		int n = env->GetArrayLength(rb);
		auto rbp = (char *)env->GetByteArrayElements(rb, JNI_FALSE);
		auto rqp = (char *)env->GetByteArrayElements(rq, JNI_FALSE);
		auto iqp = (char *)env->GetByteArrayElements(iq, JNI_FALSE);
		auto dqp = (char *)env->GetByteArrayElements(dq, JNI_FALSE);
		auto gqp = (char *)env->GetByteArrayElements(gq, JNI_FALSE);

		std::string rbv(rbp, n);
		std::string rqv(rqp, n);
		std::string iqv(iqp, n);
		std::string dqv(dqp, n);
		std::string gqv(gqp, n);
		if(n > rln)
			rln = n;
		//		std::cout << rbv <<" " << rqv << " " << dqv << std::endl;
		load_base(pkt, rbv, rqv, iqv, dqv, gqv);
		env->ReleaseByteArrayElements(rb, (jbyte *)rbp, 0);
		env->ReleaseByteArrayElements(rq, (jbyte *)rqp, 0);
		env->ReleaseByteArrayElements(iq, (jbyte *)iqp, 0);
		env->ReleaseByteArrayElements(dq, (jbyte *)dqp, 0);
		env->ReleaseByteArrayElements(gq, (jbyte *)gqp, 0);

		env->DeleteLocalRef(rb);
		env->DeleteLocalRef(rq);
		env->DeleteLocalRef(iq);
		env->DeleteLocalRef(dq);
		env->DeleteLocalRef(gq);
		env->DeleteLocalRef(rd);
	}

	for(int i = 0; i < hs; i++) {
		auto hd = env->GetObjectArrayElement(hbs, i);
		auto hb = (jbyteArray)env->GetObjectField(hd, _hp);
		int n = env->GetArrayLength(hb);
		if(n > hln)
			hln = n;
		auto hbp = (char *)env->GetByteArrayElements(hb, JNI_FALSE);
		std::string hbv(hbp, n);
		load_hap(pkt, hbv);
		env->ReleaseByteArrayElements(hb, (jbyte *)hbp, 0);
		env->DeleteLocalRef(hb);
		env->DeleteLocalRef(hd);
	}

	pkt.hs = hs;
	pkt.rs = rs;
	u32 mhsize = hs * rs * sizeof(precision);
	u32 mrsize = 16 + pkt.rd.size() + pkt.hd.size() + 16;
	if(mrsize > compute->isize || mhsize > compute->osize) {
		compute->release();
		compute->alloc(mrsize, mhsize*sizeof(precision));
	}
	intype *ibuf = (intype *)compute->ibuf;
	precision *obuf = (precision *)compute->obuf;
	u32 rcur = 0;
	u32 tosize = 0;
	intype *ep = ibuf;
	constexpr int hdr_len = 16/sizeof(intype);
	intype *rmem = &ibuf[rcur];
	u32 rlen = pkt.rd.size()/sizeof(intype);
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

	ep = &rmem[hdr_len+rlen+hlen];
	ep[0] = 0;
	ep[1] = 0;
	double *os = env->GetDoubleArrayElements(out, JNI_FALSE);
	//compute(ibuf, (precision*)os, mrsize, rs*hs);
	compute->compute(ibuf, obuf, mrsize);
	/*
	for(int i = 0; i < rs*hs; i++) {
		//precision ot = LOG10(obuf[i]) - LOG10(pkt.hls[i%pkt.hs])-lmax;
		os[i] = obuf[i];//ot;
	}
	*/
	env->SetDoubleArrayRegion(out, 0, hs*rs, os);
	env->ReleaseDoubleArrayElements(out, os, 0);
	//FREE(ibuf);
	gettimeofday(&ed, 0);
	float use = (ed.tv_sec - st.tv_sec)*1000.0 + (ed.tv_usec - st.tv_usec)/1000.0;
	printf("compute time: %.2fms rs: %u, hs: %u\n", use, rs, hs);
	//FREE(obuf);
}

extern "C" JNIEXPORT void JNICALL Java_chomper_gatk_cfgengine_CFGEngine_doneNative(JNIEnv *env, jclass clazz) {
	(void)env;
	(void)clazz;
	delete compute;
}
