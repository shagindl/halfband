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
    TEST_METHOD(test_process_receive_array) {
        double buf[BUF_SIZE], gain;
        Halfband<2, 8> cascade;

        for (int k = 1; k < MAX_FREQS; k++) {
            double filt_sum = 0, filt, orig_sum = 0, freq = 1000, Ks;

            Ks = 3 * (MAX_FREQS / k);
            test_sine(buf, BUF_SIZE, freq, Ks * freq);

            for (int j = 0; j < BUF_SIZE; j += 4) {
                filt = cascade.process(&buf[j]);
                // skip transient
                if (j > BUF_SIZE / 3) {
                    filt_sum += pow(filt, 2);
                    for (int i = 0; i < 4; i++) {
                        orig_sum += pow(buf[j + i], 2);
                    }
                }
            }
            gain = 20 * log10(sqrt(filt_sum) / sqrt(orig_sum / 4));

            std::wstringstream str;
            str << L"\nKs = " << Ks << L"; orig_sum/4 = " << orig_sum/4 << L"; filt_sum = " << filt_sum << L"\n";

            if (Ks > 16) {
                Assert::AreEqual(0, gain, 0.05, str.str().c_str(), LINE_INFO());
            } else if (Ks >= 8) {
                Assert::AreEqual(0, gain, 0.06, str.str().c_str(), LINE_INFO());
            } else if (Ks > 4) {
                Assert::AreEqual(-110, gain, 20, str.str().c_str(), LINE_INFO());
            } else {
                Assert::IsTrue(-100 > gain, str.str().c_str(), LINE_INFO());
            }
        }
    }
};

const int test_halfband::BUF_SIZE;

}/* namespace tb_halfband */
