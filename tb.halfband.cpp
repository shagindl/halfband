#include "pch.h"
#include "halfband.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace tb_halfband {

TEST_CLASS(test_halfband) {
    static const int BUF_SIZE = 2048;
    const int MAX_FREQS = 32;

    void test_sine(double* buf, int n, double F, double Fs) {
        for (int i = 0; i < n; i++)       {
            auto t = (1 / Fs) * i;
            buf[i] = sin(2*M_PI*F * t);
        }
    }
public:
    TEST_METHOD(test_halfband_decimator) {
        double buf[BUF_SIZE], gain;
        Halfband<2, 8> cascade;

        for (int k = 1; k < MAX_FREQS; k++) {
            double filt_sum = 0, filt, orig_sum = 0, freq = 1000, Ks;

            Ks = 2 * (MAX_FREQS / k);
            test_sine(buf, BUF_SIZE, freq, Ks * freq);

            for (int j = 0; j < BUF_SIZE; j += 4) {
                filt = cascade.process(&buf[j]);
                // skip transient
                if (j > BUF_SIZE / 3) {
                    filt_sum += pow(filt, 2);
                    orig_sum += pow(buf[j], 2);
                }
            }

            gain = 20 * log10(sqrt(filt_sum) / sqrt(orig_sum));
            Assert::AreEqual(gain, 0, 0.001);
            gain = 0;
        }
    }
};

const int test_halfband::BUF_SIZE;

}/* namespace tb_halfband */

#if 0
void test_halfband(void)
{
    int orders[6] = {12,10,8,6,4,2};
    int n_orders = 6;
    int order;
    int steep, i, j, k;
    halfband *half;
    double buf[BUF_SIZE];
        
    
    for(steep=0;steep<2;steep++)
    {
        printf("Steep: %d\n", steep);
        for(i=0;i<n_orders;i++)
        {   
            order = orders[i];
            printf("\tOrder: %d\n", order, steep);
            half = create_halfband(order, steep);
            
            for(k=1;k<MAX_FREQS;k++)
            {
                
                double sum, filt, sine, orig_sum, gain, freq;
                freq = k/((double)MAX_FREQS);
                test_sine(buf, BUF_SIZE, freq);
                
                sum = orig_sum = 0;
                for(j=0;j<BUF_SIZE;j++)                           
                {
                    filt = process_halfband(half, buf[j]);
                    // skip transient
                    if(j>BUF_SIZE/3)
                    {
                        sum += filt * filt;
                        orig_sum += buf[j] * buf[j];
                    }                                    
                }                            
                
                gain = sqrt(sum) / sqrt(orig_sum);                                
                printf("\t\tFreq: %.2f Gain: %6.1fdB\n", freq, 20*log(gain)/log(10));                
            }                        
            destroy_halfband(half);            
        }
        
    }
}
#endif
