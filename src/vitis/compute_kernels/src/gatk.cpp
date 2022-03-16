#include "common.h"
constexpr int num_core = CORE;
constexpr int max_read = MAX_SEQ_LEN/sizeof(intype);
constexpr u8 group = GROUP;
constexpr int maxseq = MAX_SEQ_LEN/group;

#ifndef __SYNTHESIS__
void dump(volatile void *addr, u32 size) {
		volatile u32 *a = (u32 *)addr;
		for(u32 i = 0; i < size; i++) {
				if(i && !(i % 4))
						putchar('\n');
				printf("%.8X ", a[i]);
		}
		puts("\n");
}
#endif

template<typename T, typename T2=u8>
T equal(T2 a) {
	return POW((T)10.0, a*((T)-0.1));
}

template<typename T, typename T2=u8>
T lequal(T2 a) {
	return (T)1.0 + POW((T)10.0, a*((T)0.1));
}

template<typename T, typename T2=u8>
static inline T mqual(T2 a, T2 b) {
#pragma HLS inline
	T2 va = a;
	T2 vb = b;
	T2 diff = 0;
	//T2 vmin = b;
	if (va > vb) {
		diff = va-vb;
		vb = va;
	} else {
		diff = vb-va;
	  //  vmin = a;
	}
	auto tabs_lequal = &tabs_equal[240];
	auto tequal = &tabs_equal[0];
	return tequal[vb] * tabs_lequal[diff];
}

template<typename T>
static T ABS_SUB(T a, T b) {
	return a > b ? a-b : b-a;
}

template<typename T, typename T2=T, typename TI = u8, typename len_t=u10>
T2 cu(const TI *hap, const TI *rbs, const TI *rqs, const TI *iq, const TI *dq, const TI *gq, T2 mv) {
	#pragma HLS EXPRESSION_BALANCE off
	auto tabs_lequal = &tabs_equal[240];
	auto tabs_lequal2 = &tabs_lequal[80];
	auto tequal_d1 = &tabs_equal[80];
	auto tequal_d3 = &tabs_equal[160];
	auto tequal = tabs_equal;

//#pragma HLS BIND_STORAGE variable=tequal type=ram_1wnr impl=uram
#pragma HLS inline off
	T2 res = 0;
	u16 hsize, rows;
	TI *hr = (TI *)&hsize;
	TI *rr = (TI *)&rows;
	hr[0] = hap[0];
	hr[1] = hap[1];
	rr[0] = rbs[0];
	rr[1] = rbs[1];

	T2 dxp[MAX_SEQ_LEN*3];
//#pragma HLS BIND_STORAGE variable=dxp type=ram_1wnr impl=uram
	//T2 dyp[MAX_SEQ_LEN];
	//T2 dzp[MAX_SEQ_LEN];
	T2 *dyp = &dxp[MAX_SEQ_LEN];
	T2 *dzp = &dxp[MAX_SEQ_LEN*2];

	auto rrb = &rbs[1];
	auto rrq = &rqs[1];
	auto riq = &iq[1];
	auto rdq = &dq[1];
	auto rgq = &gq[1];
	auto rhp = &hap[1];
	if(hsize == 0 || rows == 0 || rows > MAX_SEQ_LEN || hsize > MAX_SEQ_LEN)
		return 0;

	dzp[0] = 0;
	dxp[0] = 0;
	dyp[0] = mv;//MAX_VALUE;

	loop_calc:for(len_t j = 0; j < hsize;) {
#pragma HLS LOOP_TRIPCOUNT min=2 max=maxseq
		j++;
		T2 tnn = 0;
		T2 xp[group], yp[group], zp[group];
#pragma HLS BIND_OP variable=tnn op=hadd impl=fulldsp
		loop_calc_sub:for(len_t i = 1, i0 = 0; i <= rows; i++, i0++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=maxseq
#pragma HLS PIPELINE
			T2  mm, mi, md, im, ii, dd;
			TI rb = rrb[i];
			TI rq = rrq[i];
			TI iqi = riq[i];
			TI dqi = rdq[i];
			TI gqi = rgq[i];

			if(iqi == dqi)
				mm = tabs_lequal2[iqi];
			else {
				mm = (T2)1.0 - tequal[MIN(iqi, dqi)] * tabs_lequal[ABS_SUB(iqi, dqi)];
			}
			mi = tequal[iqi];
			md = tequal[dqi];
			ii = tequal[gqi];
			im = tequal_d1[gqi];

			dd = ii;
			T2 xm = j > 1 ? dxp[i] : 0;
			T2 zm = j > 1 ? dzp[i] : 0;
			T2 d1 = tequal_d1[rq];
			T2 d3 = tequal_d3[rq];
			T2 xpp, ypp, zpp;

			xpp = dxp[i0];
			ypp = dyp[i0];
			zpp = dzp[i0];

			tnn = 0; //
//			T2 xp[group], yp[group], zp[group];
			for(u8 e = 0; e < group; e++) {
				if(i == 1) {
					xp[e] = yp[e] = 0;
					zp[e] = mv;;
				}
				T2 prior = (rb == rhp[j+e]) || (rb == 'N') || (rhp[j+e] == 'N') ? d1 : d3;
				T2 t = prior * im;
				T2 x = prior * xpp * mm + t * ypp + t * zpp;
				T2 y = xp[e] * mi + yp[e] * ii;
				T2 z = xm * md + zm * ii;

				xpp = xp[e];
				ypp = yp[e];
				zpp = zp[e];
				xp[e] = x;
				yp[e] = y;
				zp[e] = z;
				xm = x;
				zm = z;
				if((j+e) <= hsize) {
					tnn += (x+y);
				}
			}
			tnn += res;

			dxp[i0] = xpp;
			dyp[i0] = ypp;
			dzp[i0] = zpp;
		}
		res = tnn;
		j += group - 1;
	}
	return res;
}

template<typename T=u8, typename len_t=u10>
void cublock(const u8 rb[CORE][MAX_SEQ_LEN],
		   const u8 rq[CORE][MAX_SEQ_LEN],
		   const u8 iq[CORE][MAX_SEQ_LEN],
		   const u8 dq[CORE][MAX_SEQ_LEN],
		   const u8 gq[CORE][MAX_SEQ_LEN],
		   const u8 hap[CORE][MAX_SEQ_LEN],
		   precision ot[CORE], u8 ve) {
#pragma HLS inline off
	typedef f32 alu_type;
	alu_type oo[CORE];
	alu_type hlf[CORE];
	u8 hls[CORE];

#pragma HLS ARRAY_PARTITION variable=oo dim=1 complete
#pragma HLS ARRAY_PARTITION variable=hls dim=1 complete
#pragma HLS ARRAY_PARTITION variable=hlf dim=1 complete

	constexpr alu_type mv = std::numeric_limits<alu_type>::max()/(alu_type)16.0;

	for(u8 j = 0; j < ve; j++) {
#pragma HLS LOOP_TRIPCOUNT max=num_core
#pragma HLS pipeline
		u16 hl;
		u8 *hl8 = (u8 *)&hl;
		hl8[0] = hap[j][0];
		hl8[1] = hap[j][1];
		//oo[j] = (alu_type)mv/(alu_type)hl;
		hls[j] = hl;
	}

	for(u8 i = 0; i < CORE; i++) {
#pragma HLS unroll
		hlf[i] = hls[i];
		oo[i] = cu<precision, alu_type>(hap[i], rb[i], rq[i], iq[i], dq[i], gq[i], MAX_VALUE);
	}
	for(u8 j = 0; j < ve; j++) {
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT max=num_core
		//ot[j] = LOG10(oo[j])-LOG10(mv);//-LOG10(hs[j]);
		ot[j] = LOG10(oo[j]/hlf[j]) - LMAX;
	}
}

template<typename rtype=intype, typename T=u8, typename len_t=u10>
u32 hload(rtype *d, T *dout) {
#pragma HLS inline
//#pragma HLS EXPRESSION_BALANCE off
	u32 n = 1;
	rtype v = d[0];
	T *vb = &dout[2];
	u16 len = v;
	u8 *l8 = (u8 *)&len;
	dout[0] = l8[0];
	dout[1] = l8[1];
	loop_hr:for(len_t i = 0; i < len; i+=sizeof(rtype)) {
#pragma HLS LOOP_TRIPCOUNT min=2 max=max_read
#pragma HLS pipeline
		rtype pn = d[n++];
		u8 *p = (u8 *)&pn;
		for(len_c_t e = 0; e < sizeof(rtype); e++) {
#pragma HLS unroll
			vb[i+e] = p[e];
		}
	}
	return n;
}
template<typename T=u8, typename bw_t=intype, typename len_t=u10>
u32 mrload(const bw_t *d, T rb[CORE][MAX_SEQ_LEN], T rq[CORE][MAX_SEQ_LEN], T iq[CORE][MAX_SEQ_LEN], T dq[CORE][MAX_SEQ_LEN], T gq[CORE][MAX_SEQ_LEN], u8 s, u8 e) {
#pragma HLS inline
//#pragma HLS EXPRESSION_BALANCE off
	u32 n = 1;
	u16 rl = d[0];
	T *l8 = (T *)&rl;
	for(len_c_t i = s; i < CORE; i++) {
#pragma HLS LOOP_TRIPCOUNT max=num_core
#pragma HLS pipeline
		if(i < e) {
			rb[i][0] = l8[0];
			rb[i][1] = l8[1];
		} else {
			rb[i][0] = 0;
			rb[i][1] = 0;
		}
	}
	bw_t *dt = (bw_t *)d;
	loop_rr:for(len_t i = 0; i < rl; i+=sizeof(bw_t)) {
#pragma HLS LOOP_TRIPCOUNT min=2 max=max_read
#pragma HLS pipeline
		bw_t r0, r1, r2, r3, r4;
		r0 = dt[n++];
		r1 = dt[n++];
		r2 = dt[n++];
		r3 = dt[n++];
		r4 = dt[n++];
		u8 *rb0, *rb1, *rb2, *rb3, *rb4;
		rb0 = (u8 *)&r0;
		rb1 = (u8 *)&r1;
		rb2 = (u8 *)&r2;
		rb3 = (u8 *)&r3;
		rb4 = (u8 *)&r4;

		for(len_c_t j = s;j < e; j++) {
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT max=num_core
			for(u8 m = 0; m < sizeof(bw_t); m++) {
#pragma HLS unroll
				rb[j][i+2+m] = rb0[m];
				rq[j][i+2+m] = rb1[m];
				iq[j][i+2+m] = rb2[m];
				dq[j][i+2+m] = rb3[m];
				gq[j][i+2+m] = rb4[m];
			}
		}
		//n += 5;
	}
	return n;
}

template<typename T=u8, typename bw_t=intype, typename len_t=u10>
u32 mrload2(const bw_t *d, T rb[CORE][MAX_SEQ_LEN*6], len_c_t s, len_c_t e) {
#pragma HLS inline
#pragma HLS EXPRESSION_BALANCE off
	u32 n = 1;
	u16 rl = d[0];
	T *l8 = (T *)&rl;
	for(len_c_t i = s; i < CORE; i++) {
#pragma HLS LOOP_TRIPCOUNT max=num_core
#pragma HLS pipeline
		if(i < e) {
			rb[i][0] = l8[0];
			rb[i][1] = l8[1];
		} else {
			rb[i][0] = 0;
			rb[i][1] = 0;
		}
	}
	bw_t *dt = (bw_t *)d;
	loop_rr:for(len_t i = 0; i < rl; i+=sizeof(bw_t)) {
#pragma HLS LOOP_TRIPCOUNT min=2 max=max_read
		bw_t r0, r1, r2, r3, r4;
		r0 = dt[n++];
		r1 = dt[n++];
		r2 = dt[n++];
		r3 = dt[n++];
		r4 = dt[n++];
		u8 *rb0, *rb1, *rb2, *rb3, *rb4;
		rb0 = (u8 *)&r0;
		rb1 = (u8 *)&r1;
		rb2 = (u8 *)&r2;
		rb3 = (u8 *)&r3;
		rb4 = (u8 *)&r4;

		loop_rfill:for(len_c_t j = s;j < e; j++) {
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT max=num_core
			T *rq = &rb[j][MAX_SEQ_LEN];
			T *iq = &rb[j][MAX_SEQ_LEN*2];
			T *dq = &rb[j][MAX_SEQ_LEN*3];
			T *gq = &rb[j][MAX_SEQ_LEN*4];
			for(u8 m = 0; m < sizeof(bw_t); m++) {
#pragma HLS unroll
				rb[j][i+2+m] = rb0[m];
				rq[i+2+m] = rb1[m];
				iq[i+2+m] = rb2[m];
				dq[i+2+m] = rb3[m];
				gq[i+2+m] = rb4[m];
			}
		}
		//n += 5;
	}
	return n;
}
#define swap64(a) (\
				((a & 0x00000000000000ff) << 56) |\
				((a & 0x000000000000ff00) << 40) |\
				((a & 0x0000000000ff0000) << 24) |\
				((a & 0x00000000ff000000) << 8 ) |\
				((a & 0x000000ff00000000) >> 8 ) |\
				((a & 0x0000ff0000000000) >> 24) |\
				((a & 0x00ff000000000000) >> 40) |\
				((a & 0xff00000000000000) >> 56))

#define swap32(a) (\
				((a & 0x000000ff) << 24) |\
				((a & 0x0000ff00) << 8) |\
				((a & 0x00ff0000) >> 8) |\
				((a & 0xff000000) >>24))

void test(u64 *a, u64 *b, u32 n) {
	for(u32 i = 0; i < n; i++) {
		a[i] = swap64(b[i]);
	}
}
constexpr u64 mem_size=(2*1024*1024*1024UL)-1UL;
extern "C" void compute(const intype *rmem, precision *omem, u32 max) {
//#pragma HLS INTERFACE m_axi bundle=gmem1 port=omem offset=slave
//#pragma HLS INTERFACE m_axi bundle=gmem0 port=rmem offset=slave

//#pragma HLS INTERFACE s_axilite port=max
//#pragma HLS INTERFACE s_axilite port=return

/*
	u8 rb[CORE][MAX_SEQ_LEN*6];
#pragma HLS BIND_STORAGE variable=rb type=ram_1wnr impl=uram
*/
	precision ot[CORE];
	u10 hs[CORE];
	u10 rs[CORE];
	u8 rb[CORE][MAX_SEQ_LEN];
	u8 rq[CORE][MAX_SEQ_LEN];
	u8 iq[CORE][MAX_SEQ_LEN];
	u8 dq[CORE][MAX_SEQ_LEN];
	u8 gq[CORE][MAX_SEQ_LEN];
	u8 hp[CORE][MAX_SEQ_LEN];


#pragma HLS BIND_STORAGE variable=rb type=ram_1wnr impl=bram
#pragma HLS BIND_STORAGE variable=rq type=ram_1wnr impl=bram
#pragma HLS BIND_STORAGE variable=iq type=ram_1wnr impl=bram
#pragma HLS BIND_STORAGE variable=dq type=ram_1wnr impl=bram
#pragma HLS BIND_STORAGE variable=gq type=ram_1wnr impl=bram
#pragma HLS BIND_STORAGE variable=hp type=ram_1wnr impl=bram


#pragma HLS ARRAY_PARTITION variable=rb dim=1 factor=num_core block
#pragma HLS ARRAY_PARTITION variable=rq dim=1 factor=num_core block
#pragma HLS ARRAY_PARTITION variable=iq dim=1 factor=num_core block
#pragma HLS ARRAY_PARTITION variable=dq dim=1 factor=num_core block
#pragma HLS ARRAY_PARTITION variable=gq dim=1 factor=num_core block
#pragma HLS ARRAY_PARTITION variable=hp dim=1 factor=num_core block
#pragma HLS ARRAY_PARTITION variable=hs dim=1 complete
#pragma HLS ARRAY_PARTITION variable=rs dim=1 complete
//#pragma HLS ARRAY_PARTITION variable=ot dim=1 complete
	for(len_c_t i = 0; i < CORE; i++) {
#pragma HLS unroll
		rs[i] = hs[i] = 0;
		ot[i] = 0;
	}
	u32 hep = 0;
	len_c_t ve = 0;
	u32 hc = 0;
	u32 iv = 0;
	constexpr u8 hdr_size = 16/sizeof(intype);
	loop_main: for(u32 rp = 0; rp < max;){
#pragma HLS LOOP_TRIPCOUNT max=1
		u32 hsize, rsize;
		u32 rlen, hlen, ln;
		intype ids[hdr_size];
		for(u8 u = 0; u < hdr_size; u++)
#pragma HLS unroll
			ids[u] = rmem[rp++];

		u32 *vd = (u32 *)ids;
		rsize = vd[0];
		hsize = vd[1];
		rlen = vd[2];
		hlen = vd[3];
		hep = 0;
		hc = 0;

		auto *imem = &rmem[rp];
		auto *hmem = &rmem[rp+rlen];
		rp += hlen+rlen;
		ln = rmem[rp];

		if(hsize == 0 || rsize == 0 || rp >= max) {
			return;
		}
		loop_read: for(u32 i = 0, ip=0; i < rsize;) {
#pragma HLS LOOP_TRIPCOUNT min=2 max=80
			len_c_t vs = ve;
			u32 rip = 0;
			loop_hap: for(u32 e = hc; e < hsize && ve < CORE; e++) {
#pragma HLS LOOP_TRIPCOUNT min=2 max=80
				hep += hload(&hmem[hep], hp[ve]);
				//auto hpp = &rb[ve][MAX_SEQ_LEN*5];
				//hep += hload(&hmem[hep], hpp, &hs[ve]);
				ve++;
				hc++;
			}
			if(vs != ve) {
				rip = mrload(&imem[ip], rb, rq, iq, dq, gq, vs, ve);
			}
			if(hc >= hsize) {
				hc = 0;
				hep = 0;
				i++;
				ip += rip;
			}
			if(!ln && ve < CORE && i < rsize) {
				continue;
			}
			cublock(rb, rq, iq, dq, gq, hp, ot, ve);
			//cublock2(rb, ot, hs, rs);
			f32 *ob = &omem[iv];
			loop_out:for(u8 j = 0; j < ve; j++) {
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT max=num_core
				ob[j] = ot[j];//LOG10(ot[j]/hs[ve]/MAX_VALUE);
				//printf("% 3u, % 3u, %.18f\n", rs[j], hs[j], ot[j]);
			}
			iv += ve;
			ve = 0;
		}
	}
}
