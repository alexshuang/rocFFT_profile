// butterfly radix-8 constants
#define C8Q 0.70710678118654752440084436210485

// butterfly radix-16 constants
#define C16A 0.923879532511286738
#define C16B 0.382683432365089837

__device__ void FwdRad16B1(float2* R0,
                           float2* R8,
                           float2* R4,
                           float2* R12,
                           float2* R2,
                           float2* R10,
                           float2* R6,
                           float2* R14,
                           float2* R1,
                           float2* R9,
                           float2* R5,
                           float2* R13,
                           float2* R3,
                           float2* R11,
                           float2* R7,
                           float2* R15)
{

    float2 res;

    (*R1)  = (*R0) - (*R1);
    (*R0)  = 2.0 * (*R0) - (*R1);
    (*R3)  = (*R2) - (*R3);
    (*R2)  = 2.0 * (*R2) - (*R3);
    (*R5)  = (*R4) - (*R5);
    (*R4)  = 2.0 * (*R4) - (*R5);
    (*R7)  = (*R6) - (*R7);
    (*R6)  = 2.0 * (*R6) - (*R7);
    (*R9)  = (*R8) - (*R9);
    (*R8)  = 2.0 * (*R8) - (*R9);
    (*R11) = (*R10) - (*R11);
    (*R10) = 2.0 * (*R10) - (*R11);
    (*R13) = (*R12) - (*R13);
    (*R12) = 2.0 * (*R12) - (*R13);
    (*R15) = (*R14) - (*R15);
    (*R14) = 2.0 * (*R14) - (*R15);

    (*R2)  = (*R0) - (*R2);
    (*R0)  = 2.0 * (*R0) - (*R2);
    (*R3)  = (*R1) + float2 (-(*R3).y, (*R3).x);
    (*R1)  = 2.0 * (*R1) - (*R3);
    (*R6)  = (*R4) - (*R6);
    (*R4)  = 2.0 * (*R4) - (*R6);
    (*R7)  = (*R5) + float2 (-(*R7).y, (*R7).x);
    (*R5)  = 2.0 * (*R5) - (*R7);
    (*R10) = (*R8) - (*R10);
    (*R8)  = 2.0 * (*R8) - (*R10);
    (*R11) = (*R9) + float2 (-(*R11).y, (*R11).x);
    (*R9)  = 2.0 * (*R9) - (*R11);
    (*R14) = (*R12) - (*R14);
    (*R12) = 2.0 * (*R12) - (*R14);
    (*R15) = (*R13) + float2 (-(*R15).y, (*R15).x);
    (*R13) = 2.0 * (*R13) - (*R15);

    (*R4)  = (*R0) - (*R4);
    (*R0)  = 2.0 * (*R0) - (*R4);
    (*R5)  = ((*R1) - C8Q * (*R5)) - C8Q * float2 ((*R5).y, -(*R5).x);
    (*R1)  = 2.0 * (*R1) - (*R5);
    (*R6)  = (*R2) + float2 (-(*R6).y, (*R6).x);
    (*R2)  = 2.0 * (*R2) - (*R6);
    (*R7)  = ((*R3) + C8Q * (*R7)) - C8Q * float2 ((*R7).y, -(*R7).x);
    (*R3)  = 2.0 * (*R3) - (*R7);
    (*R12) = (*R8) - (*R12);
    (*R8)  = 2.0 * (*R8) - (*R12);
    (*R13) = ((*R9) - C8Q * (*R13)) - C8Q * float2 ((*R13).y, -(*R13).x);
    (*R9)  = 2.0 * (*R9) - (*R13);
    (*R14) = (*R10) + float2 (-(*R14).y, (*R14).x);
    (*R10) = 2.0 * (*R10) - (*R14);
    (*R15) = ((*R11) + C8Q * (*R15)) - C8Q * float2 ((*R15).y, -(*R15).x);
    (*R11) = 2.0 * (*R11) - (*R15);

    (*R8) = (*R0) - (*R8);
    (*R0) = 2.0 * (*R0) - (*R8);
    (*R9) = ((*R1) - C16A * (*R9)) - C16B * float2 ((*R9).y, -(*R9).x);
    res   = (*R8);
    (*R1) = 2.0 * (*R1) - (*R9);

    (*R10) = ((*R2) - C8Q * (*R10)) - C8Q * float2 ((*R10).y, -(*R10).x);
    (*R2)  = 2.0 * (*R2) - (*R10);
    (*R11) = ((*R3) - C16B * (*R11)) - C16A * float2 ((*R11).y, -(*R11).x);
    (*R3)  = 2.0 * (*R3) - (*R11);

    (*R12) = (*R4) + float2 (-(*R12).y, (*R12).x);
    (*R4)  = 2.0 * (*R4) - (*R12);
    (*R13) = ((*R5) + C16B * (*R13)) - C16A * float2 ((*R13).y, -(*R13).x);
    (*R5)  = 2.0 * (*R5) - (*R13);

    (*R14) = ((*R6) + C8Q * (*R14)) - C8Q * float2 ((*R14).y, -(*R14).x);
    (*R6)  = 2.0 * (*R6) - (*R14);
    (*R15) = ((*R7) + C16A * (*R15)) - C16B * float2 ((*R15).y, -(*R15).x);
    (*R7)  = 2.0 * (*R7) - (*R15);

    res    = (*R1);
    (*R1)  = (*R8);
    (*R8)  = res;
    res    = (*R2);
    (*R2)  = (*R4);
    (*R4)  = res;
    res    = (*R3);
    (*R3)  = (*R12);
    (*R12) = res;
    res    = (*R5);
    (*R5)  = (*R10);
    (*R10) = res;
    res    = (*R7);
    (*R7)  = (*R14);
    (*R14) = res;
    res    = (*R11);
    (*R11) = (*R13);
    (*R13) = res;
}

__device__ void
FwdPass2_len4096(const float2 *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, float  *bufInRe, float  *bufInIm, float2 *bufOut, float2 *R0, float2 *R1, float2 *R2, float2 *R3, float2 *R4, float2 *R5, float2 *R6, float2 *R7, float2 *R8, float2 *R9, float2 *R10, float2 *R11, float2 *R12, float2 *R13, float2 *R14, float2 *R15)
{




	{
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 0];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 1];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 2];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 3];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 4];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R5).x; ry = (*R5).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R5).x = TR;
		(*R5).y = TI;
	}

	{
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 5];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 6];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 7];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 8];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 9];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R10).x; ry = (*R10).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R10).x = TR;
		(*R10).y = TI;
	}

	{
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 10];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 11];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 12];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 13];
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
		float2 W = twiddles[255 + 15*((1*me + 0)%256) + 14];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R15).x; ry = (*R15).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R15).x = TR;
		(*R15).y = TI;
	}

	FwdRad16B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15);


	if(rw)
	{
	bufOut[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 1*me + 0 + 256 )*stride_out] = (*R1);
	bufOut[outOffset + ( 1*me + 0 + 512 )*stride_out] = (*R2);
	bufOut[outOffset + ( 1*me + 0 + 768 )*stride_out] = (*R3);
	bufOut[outOffset + ( 1*me + 0 + 1024 )*stride_out] = (*R4);
	bufOut[outOffset + ( 1*me + 0 + 1280 )*stride_out] = (*R5);
	bufOut[outOffset + ( 1*me + 0 + 1536 )*stride_out] = (*R6);
	bufOut[outOffset + ( 1*me + 0 + 1792 )*stride_out] = (*R7);
	bufOut[outOffset + ( 1*me + 0 + 2048 )*stride_out] = (*R8);
	bufOut[outOffset + ( 1*me + 0 + 2304 )*stride_out] = (*R9);
	bufOut[outOffset + ( 1*me + 0 + 2560 )*stride_out] = (*R10);
	bufOut[outOffset + ( 1*me + 0 + 2816 )*stride_out] = (*R11);
	bufOut[outOffset + ( 1*me + 0 + 3072 )*stride_out] = (*R12);
	bufOut[outOffset + ( 1*me + 0 + 3328 )*stride_out] = (*R13);
	bufOut[outOffset + ( 1*me + 0 + 3584 )*stride_out] = (*R14);
	bufOut[outOffset + ( 1*me + 0 + 3840 )*stride_out] = (*R15);
	}

}

__device__ void
FwdPass1_len4096(const float2 *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, float  *bufInRe, float  *bufInIm, float  *bufOutRe, float  *bufOutIm, float2 *R0, float2 *R1, float2 *R2, float2 *R3, float2 *R4, float2 *R5, float2 *R6, float2 *R7, float2 *R8, float2 *R9, float2 *R10, float2 *R11, float2 *R12, float2 *R13, float2 *R14, float2 *R15)
{




	{
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 0];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 1];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 2];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 3];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 4];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R5).x; ry = (*R5).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R5).x = TR;
		(*R5).y = TI;
	}

	{
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 5];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 6];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 7];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 8];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 9];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R10).x; ry = (*R10).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R10).x = TR;
		(*R10).y = TI;
	}

	{
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 10];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 11];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 12];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 13];
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
		float2 W = twiddles[15 + 15*((1*me + 0)%16) + 14];
		float  TR, TI;
		float wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R15).x; ry = (*R15).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R15).x = TR;
		(*R15).y = TI;
	}

	FwdRad16B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15);


	if(rw)
	{
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 16 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 32 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 48 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 64 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 80 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 96 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 112 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 128 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 144 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 160 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 176 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 192 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 208 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 224 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 240 ) ] = (*R15).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 256 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 512 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 768 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 1024 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 1280 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 1536 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 1792 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 2048 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 2304 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 2560 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 2816 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 3072 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 3328 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 3584 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 3840 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 16 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 32 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 48 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 64 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 80 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 96 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 112 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 128 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 144 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 160 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 176 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 192 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 208 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 224 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((1*me + 0)/16)*256 + (1*me + 0)%16 + 240 ) ] = (*R15).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 256 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 512 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 768 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 1024 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 1280 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 1536 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 1792 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 2048 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 2304 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 2560 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 2816 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 3072 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 3328 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 3584 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 3840 ) ];
	}


	__syncthreads();

}

__device__ void
FwdPass0_len4096(const float2 *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, float2 *bufIn, float  *bufOutRe, float  *bufOutIm, float2 *R0, float2 *R1, float2 *R2, float2 *R3, float2 *R4, float2 *R5, float2 *R6, float2 *R7, float2 *R8, float2 *R9, float2 *R10, float2 *R11, float2 *R12, float2 *R13, float2 *R14, float2 *R15)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 256 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 512 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 768 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*1 + 0 + 1024 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*1 + 0 + 1280 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*1 + 0 + 1536 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*1 + 0 + 1792 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*1 + 0 + 2048 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*1 + 0 + 2304 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*1 + 0 + 2560 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*1 + 0 + 2816 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*1 + 0 + 3072 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*1 + 0 + 3328 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*1 + 0 + 3584 )*stride_in];
	(*R15) = bufIn[inOffset + ( 0 + me*1 + 0 + 3840 )*stride_in];
	}



	FwdRad16B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15);


	if(rw)
	{
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 10 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 11 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 12 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 13 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 14 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 15 ) ] = (*R15).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 256 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 512 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 768 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 1024 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 1280 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 1536 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 1792 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 2048 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 2304 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 2560 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 2816 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 3072 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 3328 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 3584 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 3840 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 10 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 11 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 12 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 13 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 14 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*16 + (1*me + 0)%1 + 15 ) ] = (*R15).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 256 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 512 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 768 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 1024 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 1280 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 1536 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 1792 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 2048 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 2304 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 2560 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 2816 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 3072 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 3328 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 3584 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 3840 ) ];
	}


	__syncthreads();

}

__device__ void 
fwd_len4096_device(const float2 *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, float2 *lwbIn, float2 *lwbOut, float  *lds)
{
	float2 R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15;
	FwdPass0_len4096(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15);
	FwdPass1_len4096(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15);
	FwdPass2_len4096(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15);
}


__global__ void 
fft_fwd_ip_len4096( const float2 * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, float2 * __restrict__ gb)
{

	__shared__ float  lds[4096];
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
	fwd_len4096_device(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, lwb, lwb, lds);
}

