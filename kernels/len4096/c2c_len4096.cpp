#include <cassert>
#include <complex>
#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <hip/hip_runtime_api.h>

#include "hip/hip_runtime.h"
#include "rocfft.h"

#include "rocfft_kernel_4096.h"

#define KERN_ARGS_ARRAY_WIDTH 16

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define MAX_WORK_GROUP_SIZE 1024
#define TWO_PI (-6.283185307179586476925286766559)
#define TWO_PI_FP32 (6.283185307179586476925286766559F)
#define TWO_PI_FP64 (6.283185307179586476925286766559)

/* radix table: tell the FFT algorithms for size <= 4096 ; required by twiddle, passes, and kernel*/
struct SpecRecord
{
    size_t length;
    size_t workGroupSize;
    size_t numTransforms;
    size_t numPasses;
    size_t radices[12]; // Setting upper limit of number of passes to 12
};
inline const std::vector<SpecRecord>& GetRecord()
{
    static const std::vector<SpecRecord> specRecord = 
	{
        //  Length, WorkGroupSize (thread block size), NumTransforms , NumPasses,
        //  Radices
        //  vector<size_t> radices; NUmPasses = radices.size();
        //  Tuned for single precision on OpenCL stack; double precsion use the
        //  same table as single
        {4096, 256,  1, 3, 16, 16, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // pow2
        {2048, 256,  1, 4,  8,  8,  8, 4, 0, 0, 0, 0, 0, 0, 0, 0},
        {1024, 128,  1, 4,  8,  8,  4, 4, 0, 0, 0, 0, 0, 0, 0, 0},
        { 512,  64,  1, 3,  8,  8,  8, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        { 256,  64,  1, 4,  4,  4,  4, 4, 0, 0, 0, 0, 0, 0, 0, 0},
        { 128,  64,  4, 3,  8,  4,  4, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {  64,  64,  4, 3,  4,  4,  4, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {  32,  64, 16, 2,  8,  4,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {  16,  64, 16, 2,  4,  4,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {   8,  64, 32, 2,  4,  2,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {   4,  64, 32, 2,  2,  2,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {   2,  64, 64, 1,  2,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    };

    return specRecord;
}
inline void DetermineSizes(const size_t& length, size_t& workGroupSize, size_t& numTrans)
{
    assert(MAX_WORK_GROUP_SIZE >= 64);

    if(length == 1) // special case
    {
        workGroupSize = 64;
        numTrans      = 64;
        return;
    }

    size_t baseRadix[]   = {13, 11, 7, 5, 3, 2}; // list only supported primes
    size_t baseRadixSize = sizeof(baseRadix) / sizeof(baseRadix[0]);

    size_t                   l = length;
    std::map<size_t, size_t> primeFactorsExpanded;
    for(size_t r = 0; r < baseRadixSize; r++)
    {
        size_t rad = baseRadix[r];
        size_t e   = 1;
        while(!(l % rad))
        {
            l /= rad;
            e *= rad;
        }

        primeFactorsExpanded[rad] = e;
    }

    assert(l == 1); // Makes sure the number is composed of only supported primes

    if(primeFactorsExpanded[2] == length) // Length is pure power of 2
    {
        if(length >= 1024)
        {
            workGroupSize = (MAX_WORK_GROUP_SIZE >= 256) ? 256 : MAX_WORK_GROUP_SIZE;
            numTrans      = 1;
        }
        else if(length == 512)
        {
            workGroupSize = 64;
            numTrans      = 1;
        }
        else if(length >= 16)
        {
            workGroupSize = 64;
            numTrans      = 256 / length;
        }
        else
        {
            workGroupSize = 64;
            numTrans      = 128 / length;
        }
    }
    else if(primeFactorsExpanded[3] == length) // Length is pure power of 3
    {
        workGroupSize = (MAX_WORK_GROUP_SIZE >= 256) ? 243 : 27;
        numTrans      = length >= 3 * workGroupSize ? 1 : (3 * workGroupSize) / length;
    }
    else if(primeFactorsExpanded[5] == length) // Length is pure power of 5
    {
        workGroupSize = (MAX_WORK_GROUP_SIZE >= 128) ? 125 : 25;
        numTrans      = length >= 5 * workGroupSize ? 1 : (5 * workGroupSize) / length;
    }
    else if(primeFactorsExpanded[7] == length) // Length is pure power of 7
    {
        workGroupSize = 49;
        numTrans      = length >= 7 * workGroupSize ? 1 : (7 * workGroupSize) / length;
    }
    else if(primeFactorsExpanded[11] == length) // Length is pure power of 11
    {
        workGroupSize = 121;
        numTrans      = length >= 11 * workGroupSize ? 1 : (11 * workGroupSize) / length;
    }
    else if(primeFactorsExpanded[13] == length) // Length is pure power of 13
    {
        workGroupSize = 169;
        numTrans      = length >= 13 * workGroupSize ? 1 : (13 * workGroupSize) / length;
    }
    else
    {
        size_t leastNumPerWI    = 1; // least number of elements in one work item
        size_t maxWorkGroupSize = MAX_WORK_GROUP_SIZE; // maximum work group size desired

        if(primeFactorsExpanded[2] * primeFactorsExpanded[3] == length)
        {
            if(length % 12 == 0)
            {
                leastNumPerWI    = 12;
                maxWorkGroupSize = 128;
            }
            else
            {
                leastNumPerWI    = 6;
                maxWorkGroupSize = 256;
            }
        }
        else if(primeFactorsExpanded[2] * primeFactorsExpanded[5] == length)
        {
            if(length % 20 == 0)
            {
                leastNumPerWI    = 20;
                maxWorkGroupSize = 64;
            }
            else
            {
                leastNumPerWI    = 10;
                maxWorkGroupSize = 128;
            }
        }
        else if(primeFactorsExpanded[2] * primeFactorsExpanded[7] == length)
        {
            leastNumPerWI    = 14;
            maxWorkGroupSize = 64;
        }
        else if(primeFactorsExpanded[3] * primeFactorsExpanded[5] == length)
        {
            leastNumPerWI    = 15;
            maxWorkGroupSize = 128;
        }
        else if(primeFactorsExpanded[3] * primeFactorsExpanded[7] == length)
        {
            leastNumPerWI    = 21;
            maxWorkGroupSize = 128;
        }
        else if(primeFactorsExpanded[5] * primeFactorsExpanded[7] == length)
        {
            leastNumPerWI    = 35;
            maxWorkGroupSize = 64;
        }
        else if(primeFactorsExpanded[2] * primeFactorsExpanded[3] * primeFactorsExpanded[5] == length)
        {
            leastNumPerWI    = 30;
            maxWorkGroupSize = 64;
        }
        else if(primeFactorsExpanded[2] * primeFactorsExpanded[3] * primeFactorsExpanded[7] == length)
        {
            leastNumPerWI    = 42;
            maxWorkGroupSize = 60;
        }
        else if(primeFactorsExpanded[2] * primeFactorsExpanded[5] * primeFactorsExpanded[7] == length)
        {
            leastNumPerWI    = 70;
            maxWorkGroupSize = 36;
        }
        else if(primeFactorsExpanded[3] * primeFactorsExpanded[5] * primeFactorsExpanded[7] == length)
        {
            leastNumPerWI    = 105;
            maxWorkGroupSize = 24;
        }
        else if(primeFactorsExpanded[2] * primeFactorsExpanded[11] == length)
        {
            leastNumPerWI    = 22;
            maxWorkGroupSize = 128;
        }
        else if(primeFactorsExpanded[2] * primeFactorsExpanded[13] == length)
        {
            leastNumPerWI    = 26;
            maxWorkGroupSize = 128;
        }
        else
        {
            leastNumPerWI    = 210;
            maxWorkGroupSize = 12;
        }

        if(maxWorkGroupSize > MAX_WORK_GROUP_SIZE)
            maxWorkGroupSize = MAX_WORK_GROUP_SIZE;
        assert(leastNumPerWI > 0 && length % leastNumPerWI == 0);

        for(size_t lnpi = leastNumPerWI; lnpi <= length; lnpi += leastNumPerWI)
        {
            if(length % lnpi != 0)
                continue;

            if(length / lnpi <= MAX_WORK_GROUP_SIZE)
            {
                leastNumPerWI = lnpi;
                break;
            }
        }

        numTrans      = maxWorkGroupSize / (length / leastNumPerWI);
        numTrans      = numTrans < 1 ? 1 : numTrans;
        workGroupSize = numTrans * (length / leastNumPerWI);
    }

    assert(workGroupSize <= MAX_WORK_GROUP_SIZE);
}
std::vector<size_t> GetRadices(size_t length)
{
    std::vector<size_t> radices;

    // get number of items in this table
    std::vector<SpecRecord> specRecord  = GetRecord();
    size_t tableLength = specRecord.size();

    printf("tableLength=%zu\n", tableLength);
    for(int i = 0; i < tableLength; i++)
    {
        if(length == specRecord[i].length)
        { // if find the matched size

            size_t numPasses = specRecord[i].numPasses;
            printf("numPasses=%zu, table item %d \n", numPasses, i);
            for(int j = 0; j < numPasses; j++)
            {
                radices.push_back((specRecord[i].radices)[j]);
            }
            break;
        }
    }

    // if not in the table, then generate the radice order with the algorithm.
    if(radices.size() == 0)
    {
        size_t R = length;

        // Possible radices
        size_t cRad[]   = {13, 11, 10, 8, 7, 6, 5, 4, 3, 2, 1}; // Must be in descending order
        size_t cRadSize = (sizeof(cRad) / sizeof(cRad[0]));

        size_t workGroupSize;
        size_t numTrans;
        // need to know workGroupSize and numTrans
        DetermineSizes(length, workGroupSize, numTrans);
        size_t cnPerWI = (numTrans * length) / workGroupSize;

        // Generate the radix and pass objects
        while(true)
        {
            size_t rad;

            // Picks the radices in descending order (biggest radix first) performanceCC
            // purpose
            for(size_t r = 0; r < cRadSize; r++)
            {

                rad = cRad[r];
                if((rad > cnPerWI) || (cnPerWI % rad))
                    continue;

                if(!(R % rad)) // if not a multiple of rad, then exit
                    break;
            }

            assert((cnPerWI % rad) == 0);

            R /= rad;
            radices.push_back(rad);

            assert(R >= 1);
            if(R == 1)
                break;

        } // end while
    } // end if(radices == empty)

	printf("radices size = %zu:\n",radices.size());
	for(std::vector<size_t>::const_iterator i = radices.begin(); i != radices.end(); i++)
		printf("\t%lu\n",*i);
	
    return radices;
}
class TwiddleTable
{
    size_t N;
    float2 * wc;

public:
    TwiddleTable(size_t length) : N(length)    {  wc = new float2[N];    }

    ~TwiddleTable()    { delete[] wc;    }

    float2 * GenerateTwiddleTable(const std::vector<size_t>& radices)
    {
		size_t radix;
        size_t L  = 1;
        size_t nt = 0;
		size_t idxcnt = 0;
        for(std::vector<size_t>::const_iterator i = radices.begin(); i != radices.end(); i++)
        {
            radix = *i;
            L *= radix;
			
			printf("radix = %zu, L = %zu\n", radix, L);

            // Twiddle factors
            for(size_t k = 0; k < (L / radix); k++)
            {
                double theta = TWO_PI * (k) / (L);
				printf("	[k/L] = [%zu/%zu], theta = %.3f\n", k, L, theta);
				
				{
					size_t j = 0;
					double c = cos((j)*theta);
					double s = sin ((j)*theta);
					printf("		j = %zu, alph = %.3f, c = %.3f, s = %.3f\n", j, j*theta, c,s);
				}
					
                for(size_t j = 1; j < radix; j++)
                {
                    double c = cos((j)*theta);
                    double s = sin ((j)*theta);
					printf("[%zu]		j = %zu, alph = %.3f, c = %.3f, s = %.3f\n", idxcnt, j, j*theta, c,s);
					idxcnt++;

                    wc[nt].x = c;
                    wc[nt].y = s;					
                    nt++;
                }
            }
        }

		printf("\n");
        return wc;
    }
};

int main(int argc, char* argv[])
{		
	// The problem size
	const std::vector<size_t> length{4096};
	const std::vector<size_t> inStride{1};
	const std::vector<size_t> outStride{1};
	const size_t batch_count = 1;
	const size_t dim = length.size();

    size_t kargs[3 * KERN_ARGS_ARRAY_WIDTH] = {0};
	size_t iDist, oDist;
	size_t *dev_kargs = NULL;

	size_t sz = batch_count;
	int i = 0;
    for (; i < dim; i++) {
        sz *= length[i];
        kargs[i + 0 * KERN_ARGS_ARRAY_WIDTH] = length[i];
        kargs[i + 1 * KERN_ARGS_ARRAY_WIDTH] = inStride[i];
        kargs[i + 2 * KERN_ARGS_ARRAY_WIDTH] = outStride[i];
	}
	iDist = oDist = sz;
    kargs[i + 1 * KERN_ARGS_ARRAY_WIDTH] = iDist;
    kargs[i + 2 * KERN_ARGS_ARRAY_WIDTH] = oDist;
	hipMalloc(&dev_kargs, sizeof(kargs) * sizeof(kargs[0]));
    hipMemcpy(dev_kargs, kargs, sizeof(kargs) * sizeof(kargs[0]), hipMemcpyHostToDevice);

	std::vector<std::complex<float>> bufIn(sz);
	std::vector<std::complex<float>> bufOut(sz);

	for(size_t i = 0; i < sz; ++i)
	{
		bufIn[i] = std::complex<float>(1.0f*i, -0.1f*i);
		bufOut[i] = std::complex<float>(0, 0);
	}

	// Create HIP device objects:
	float2 *dev_bufIn = NULL;
	float2 *dev_bufOut = NULL;
	hipMalloc(&dev_bufIn, bufIn.size() * sizeof(float2));
	hipMalloc(&dev_bufOut, bufOut.size() * sizeof(float2));
	hipMemcpy(dev_bufIn, bufIn.data(), bufIn.size() * sizeof(float2), hipMemcpyHostToDevice);

	printf("twiddles_create_pr \n");
	std::vector<size_t> radices;
	radices = GetRadices(length[0]);

	float2 * twtc;
	float2 * dev_tw = NULL;
	TwiddleTable twTable(length[0]);
	twtc = twTable.GenerateTwiddleTable(radices);
	hipMalloc(&dev_tw, length[0] * sizeof(float2));
	hipMemcpy(dev_tw, twtc, length[0] * sizeof(float2), hipMemcpyHostToDevice);	

	hipLaunchKernelGGL(fft_fwd_ip_len4096, dim3(1), dim3(256), 0, 0,
						dev_tw, dim, dev_kargs, dev_kargs + 1 * KERN_ARGS_ARRAY_WIDTH, batch_count, dev_bufIn);

	hipFree(dev_bufIn);
	hipFree(dev_bufOut);

	std::cout << "done\n";
}

