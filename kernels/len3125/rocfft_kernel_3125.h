// butterfly radix-5 constants
#define C5QA 0.30901699437494742410229341718282
#define C5QB 0.95105651629515357211643933337938
#define C5QC 0.50000000000000000000000000000000
#define C5QD 0.58778525229247312916870595463907
#define C5QE 0.80901699437494742410229341718282

__device__ void FwdRad5B1(float2* R0, float2* R1, float2* R2, float2* R3, float2* R4)
{

    float  TR0, TI0, TR1, TI1, TR2, TI2, TR3, TI3, TR4, TI4;

    TR0 = (*R0).x + (*R1).x + (*R2).x + (*R3).x + (*R4).x;
    TR1 = ((*R0).x - C5QC * ((*R2).x + (*R3).x)) + C5QB * ((*R1).y - (*R4).y)
          + C5QD * ((*R2).y - (*R3).y) + C5QA * (((*R1).x - (*R2).x) + ((*R4).x - (*R3).x));
    TR4 = ((*R0).x - C5QC * ((*R2).x + (*R3).x)) - C5QB * ((*R1).y - (*R4).y)
          - C5QD * ((*R2).y - (*R3).y) + C5QA * (((*R1).x - (*R2).x) + ((*R4).x - (*R3).x));
    TR2 = ((*R0).x - C5QC * ((*R1).x + (*R4).x)) - C5QB * ((*R2).y - (*R3).y)
          + C5QD * ((*R1).y - (*R4).y) + C5QA * (((*R2).x - (*R1).x) + ((*R3).x - (*R4).x));
    TR3 = ((*R0).x - C5QC * ((*R1).x + (*R4).x)) + C5QB * ((*R2).y - (*R3).y)
          - C5QD * ((*R1).y - (*R4).y) + C5QA * (((*R2).x - (*R1).x) + ((*R3).x - (*R4).x));

    TI0 = (*R0).y + (*R1).y + (*R2).y + (*R3).y + (*R4).y;
    TI1 = ((*R0).y - C5QC * ((*R2).y + (*R3).y)) - C5QB * ((*R1).x - (*R4).x)
          - C5QD * ((*R2).x - (*R3).x) + C5QA * (((*R1).y - (*R2).y) + ((*R4).y - (*R3).y));
    TI4 = ((*R0).y - C5QC * ((*R2).y + (*R3).y)) + C5QB * ((*R1).x - (*R4).x)
          + C5QD * ((*R2).x - (*R3).x) + C5QA * (((*R1).y - (*R2).y) + ((*R4).y - (*R3).y));
    TI2 = ((*R0).y - C5QC * ((*R1).y + (*R4).y)) + C5QB * ((*R2).x - (*R3).x)
          - C5QD * ((*R1).x - (*R4).x) + C5QA * (((*R2).y - (*R1).y) + ((*R3).y - (*R4).y));
    TI3 = ((*R0).y - C5QC * ((*R1).y + (*R4).y)) - C5QB * ((*R2).x - (*R3).x)
          + C5QD * ((*R1).x - (*R4).x) + C5QA * (((*R2).y - (*R1).y) + ((*R3).y - (*R4).y));

    ((*R0).x) = TR0;
    ((*R0).y) = TI0;
    ((*R1).x) = TR1;
    ((*R1).y) = TI1;
    ((*R2).x) = TR2;
    ((*R2).y) = TI2;
    ((*R3).x) = TR3;
    ((*R3).y) = TI3;
    ((*R4).x) = TR4;
    ((*R4).y) = TI4;
}


__device__ void
FwdPass0_len3125(const float2 *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, float2 *bufIn, float  *bufOutRe, float  *bufOutIm, float2 *R0, float2 *R1, float2 *R2, float2 *R3, float2 *R4, float2 *R5, float2 *R6, float2 *R7, float2 *R8, float2 *R9, float2 *R10, float2 *R11, float2 *R12, float2 *R13, float2 *R14, float2 *R15, float2 *R16, float2 *R17, float2 *R18, float2 *R19, float2 *R20, float2 *R21, float2 *R22, float2 *R23, float2 *R24)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*5 + 0 + 0 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*5 + 1 + 0 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*5 + 2 + 0 )*stride_in];
	(*R15) = bufIn[inOffset + ( 0 + me*5 + 3 + 0 )*stride_in];
	(*R20) = bufIn[inOffset + ( 0 + me*5 + 4 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*5 + 0 + 625 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*5 + 1 + 625 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*5 + 2 + 625 )*stride_in];
	(*R16) = bufIn[inOffset + ( 0 + me*5 + 3 + 625 )*stride_in];
	(*R21) = bufIn[inOffset + ( 0 + me*5 + 4 + 625 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*5 + 0 + 1250 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*5 + 1 + 1250 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*5 + 2 + 1250 )*stride_in];
	(*R17) = bufIn[inOffset + ( 0 + me*5 + 3 + 1250 )*stride_in];
	(*R22) = bufIn[inOffset + ( 0 + me*5 + 4 + 1250 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*5 + 0 + 1875 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*5 + 1 + 1875 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*5 + 2 + 1875 )*stride_in];
	(*R18) = bufIn[inOffset + ( 0 + me*5 + 3 + 1875 )*stride_in];
	(*R23) = bufIn[inOffset + ( 0 + me*5 + 4 + 1875 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*5 + 0 + 2500 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*5 + 1 + 2500 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*5 + 2 + 2500 )*stride_in];
	(*R19) = bufIn[inOffset + ( 0 + me*5 + 3 + 2500 )*stride_in];
	(*R24) = bufIn[inOffset + ( 0 + me*5 + 4 + 2500 )*stride_in];
	}



	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 1 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 2 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 3 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 4 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 1 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 2 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 3 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 4 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 1 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 2 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 3 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 4 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 1 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 2 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 3 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 4 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 1 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 2 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 3 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 4 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 1 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 2 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 3 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 4 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

__device__ void
FwdPass1_len3125(const float2 *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, float  *bufInRe, float  *bufInIm, float  *bufOutRe, float  *bufOutIm, float2 *R0, float2 *R1, float2 *R2, float2 *R3, float2 *R4, float2 *R5, float2 *R6, float2 *R7, float2 *R8, float2 *R9, float2 *R10, float2 *R11, float2 *R12, float2 *R13, float2 *R14, float2 *R15, float2 *R16, float2 *R17, float2 *R18, float2 *R19, float2 *R20, float2 *R21, float2 *R22, float2 *R23, float2 *R24)
{




	{
		float2 W = twiddles[4 + 4*((5*me + 0)%5) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 0)%5) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 0)%5) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 0)%5) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R4).x; ry = (*R4).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R4).x = TR;
		(*R4).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 1)%5) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 1)%5) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 1)%5) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 1)%5) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 2)%5) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R11).x; ry = (*R11).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R11).x = TR;
		(*R11).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 2)%5) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 2)%5) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 2)%5) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 3)%5) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R16).x; ry = (*R16).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R16).x = TR;
		(*R16).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 3)%5) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 3)%5) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R18).x; ry = (*R18).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R18).x = TR;
		(*R18).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 3)%5) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 4)%5) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 4)%5) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 4)%5) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		float2 W = twiddles[4 + 4*((5*me + 4)%5) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 5 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 10 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 15 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 20 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 5 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 10 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 15 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 20 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 5 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 10 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 15 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 20 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 5 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 10 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 15 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 20 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 5 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 10 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 15 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 20 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 5 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 10 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 15 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 20 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 5 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 10 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 15 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 20 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 5 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 10 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 15 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 20 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 5 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 10 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 15 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 20 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 5 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 10 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 15 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 20 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

__device__ void
FwdPass2_len3125(const float2 *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, float  *bufInRe, float  *bufInIm, float  *bufOutRe, float  *bufOutIm, float2 *R0, float2 *R1, float2 *R2, float2 *R3, float2 *R4, float2 *R5, float2 *R6, float2 *R7, float2 *R8, float2 *R9, float2 *R10, float2 *R11, float2 *R12, float2 *R13, float2 *R14, float2 *R15, float2 *R16, float2 *R17, float2 *R18, float2 *R19, float2 *R20, float2 *R21, float2 *R22, float2 *R23, float2 *R24)
{




	{
		float2 W = twiddles[24 + 4*((5*me + 0)%25) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 0)%25) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 0)%25) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 0)%25) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R4).x; ry = (*R4).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R4).x = TR;
		(*R4).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 1)%25) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 1)%25) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 1)%25) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 1)%25) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 2)%25) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R11).x; ry = (*R11).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R11).x = TR;
		(*R11).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 2)%25) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 2)%25) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 2)%25) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 3)%25) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R16).x; ry = (*R16).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R16).x = TR;
		(*R16).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 3)%25) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 3)%25) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R18).x; ry = (*R18).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R18).x = TR;
		(*R18).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 3)%25) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 4)%25) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 4)%25) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 4)%25) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		float2 W = twiddles[24 + 4*((5*me + 4)%25) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 25 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 50 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 75 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 100 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 25 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 50 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 75 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 100 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 25 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 50 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 75 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 100 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 25 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 50 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 75 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 100 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 25 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 50 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 75 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 100 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 25 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 50 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 75 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 100 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 25 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 50 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 75 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 100 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 25 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 50 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 75 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 100 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 25 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 50 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 75 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 100 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 25 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 50 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 75 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 100 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

__device__ void
FwdPass3_len3125(const float2 *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, float  *bufInRe, float  *bufInIm, float  *bufOutRe, float  *bufOutIm, float2 *R0, float2 *R1, float2 *R2, float2 *R3, float2 *R4, float2 *R5, float2 *R6, float2 *R7, float2 *R8, float2 *R9, float2 *R10, float2 *R11, float2 *R12, float2 *R13, float2 *R14, float2 *R15, float2 *R16, float2 *R17, float2 *R18, float2 *R19, float2 *R20, float2 *R21, float2 *R22, float2 *R23, float2 *R24)
{




	{
		float2 W = twiddles[124 + 4*((5*me + 0)%125) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 0)%125) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 0)%125) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 0)%125) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R4).x; ry = (*R4).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R4).x = TR;
		(*R4).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 1)%125) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 1)%125) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 1)%125) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 1)%125) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 2)%125) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R11).x; ry = (*R11).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R11).x = TR;
		(*R11).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 2)%125) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 2)%125) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 2)%125) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 3)%125) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R16).x; ry = (*R16).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R16).x = TR;
		(*R16).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 3)%125) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 3)%125) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R18).x; ry = (*R18).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R18).x = TR;
		(*R18).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 3)%125) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 4)%125) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 4)%125) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 4)%125) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		float2 W = twiddles[124 + 4*((5*me + 4)%125) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 125 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 250 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 375 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 500 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 125 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 250 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 375 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 500 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 125 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 250 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 375 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 500 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 125 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 250 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 375 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 500 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 125 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 250 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 375 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 500 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 125 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 250 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 375 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 500 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 125 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 250 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 375 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 500 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 125 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 250 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 375 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 500 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 125 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 250 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 375 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 500 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 125 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 250 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 375 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 500 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

__device__ void
FwdPass4_len3125(const float2 *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, float  *bufInRe, float  *bufInIm, float2 *bufOut, float2 *R0, float2 *R1, float2 *R2, float2 *R3, float2 *R4, float2 *R5, float2 *R6, float2 *R7, float2 *R8, float2 *R9, float2 *R10, float2 *R11, float2 *R12, float2 *R13, float2 *R14, float2 *R15, float2 *R16, float2 *R17, float2 *R18, float2 *R19, float2 *R20, float2 *R21, float2 *R22, float2 *R23, float2 *R24)
{




	{
		float2 W = twiddles[624 + 4*((5*me + 0)%625) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 0)%625) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 0)%625) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 0)%625) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R4).x; ry = (*R4).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R4).x = TR;
		(*R4).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 1)%625) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 1)%625) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 1)%625) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 1)%625) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 2)%625) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R11).x; ry = (*R11).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R11).x = TR;
		(*R11).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 2)%625) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 2)%625) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 2)%625) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 3)%625) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R16).x; ry = (*R16).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R16).x = TR;
		(*R16).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 3)%625) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 3)%625) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R18).x; ry = (*R18).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R18).x = TR;
		(*R18).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 3)%625) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 4)%625) + 0];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 4)%625) + 1];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 4)%625) + 2];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		float2 W = twiddles[624 + 4*((5*me + 4)%625) + 3];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOut[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R5);
	bufOut[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R10);
	bufOut[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R15);
	bufOut[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R20);
	bufOut[outOffset + ( 5*me + 0 + 625 )*stride_out] = (*R1);
	bufOut[outOffset + ( 5*me + 1 + 625 )*stride_out] = (*R6);
	bufOut[outOffset + ( 5*me + 2 + 625 )*stride_out] = (*R11);
	bufOut[outOffset + ( 5*me + 3 + 625 )*stride_out] = (*R16);
	bufOut[outOffset + ( 5*me + 4 + 625 )*stride_out] = (*R21);
	bufOut[outOffset + ( 5*me + 0 + 1250 )*stride_out] = (*R2);
	bufOut[outOffset + ( 5*me + 1 + 1250 )*stride_out] = (*R7);
	bufOut[outOffset + ( 5*me + 2 + 1250 )*stride_out] = (*R12);
	bufOut[outOffset + ( 5*me + 3 + 1250 )*stride_out] = (*R17);
	bufOut[outOffset + ( 5*me + 4 + 1250 )*stride_out] = (*R22);
	bufOut[outOffset + ( 5*me + 0 + 1875 )*stride_out] = (*R3);
	bufOut[outOffset + ( 5*me + 1 + 1875 )*stride_out] = (*R8);
	bufOut[outOffset + ( 5*me + 2 + 1875 )*stride_out] = (*R13);
	bufOut[outOffset + ( 5*me + 3 + 1875 )*stride_out] = (*R18);
	bufOut[outOffset + ( 5*me + 4 + 1875 )*stride_out] = (*R23);
	bufOut[outOffset + ( 5*me + 0 + 2500 )*stride_out] = (*R4);
	bufOut[outOffset + ( 5*me + 1 + 2500 )*stride_out] = (*R9);
	bufOut[outOffset + ( 5*me + 2 + 2500 )*stride_out] = (*R14);
	bufOut[outOffset + ( 5*me + 3 + 2500 )*stride_out] = (*R19);
	bufOut[outOffset + ( 5*me + 4 + 2500 )*stride_out] = (*R24);
	}

}

__device__ void 
fwd_len3125_device(const float2 *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, float2 *lwbIn, float2 *lwbOut, float  *lds)
{
	float2 R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24;
	FwdPass0_len3125(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass1_len3125(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass2_len3125(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass3_len3125(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass4_len3125(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
}

__global__ void 
fft_fwd_ip_len3125( const float2 * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, float2 * __restrict__ gb)
{

	__shared__ float  lds[3125];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	float2 *lwb;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

	size_t counter_mod = batch;
	if(dim == 1){
		ioOffset += counter_mod*stride_in[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		ioOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		ioOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			ioOffset += (counter_mod / currentLength)*stride_in[i];
			counter_mod = counter_mod % currentLength;
		}
		ioOffset+= counter_mod * stride_in[1];
	}
	lwb = gb + ioOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	fwd_len3125_device(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, lwb, lwb, lds);
}

