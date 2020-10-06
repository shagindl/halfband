/* Allpass cascade half band filtering. From http://www.musicdsp.org/showone.php?id=39, poretd to plain c */
#include <memory>
#include <stdlib.h>
#include <string.h>
#include "library/Delegate.h"

template<int n, int order, typename T = double, bool steep = false> 
class Halfband {
    typedef struct allpass_t {
        T x0, x1, x2;
        T y0, y1, y2;
        T a;
    public:
        allpass_t(const T& _a) {
            a = _a;
            x0 = x1 = x2 = 0;
            y0 = y1 = y2 = 0;
        }
        T process(const T& input) {
            x2=x1;
            x1=x0;
            x0=input;

            y2=y1;    
            y1=y0;
    
            y0 = x2+((input-y2)*a);

            return y0;
        }
    } allpass;

    class allpass_cascade {
        std::unique_ptr<allpass> filters[6];
        const int NF;
    public:
        allpass_cascade(const T *coefficients, int nf) : NF(nf) {
            for(int i = 0; i < NF; i++)
                filters[i] = std::make_unique<allpass>(coefficients[i]);
        }
        T process(const T& input) {
            T output = input;

            for (int i = 0; i < NF; i++)
                output = filters[i]->process(output);

            return output;    
        }
    };

    typedef struct halfband_t {
        std::unique_ptr<allpass_cascade> a;
        std::unique_ptr<allpass_cascade> b;
        T oldout;
        // --
        halfband_t* next;

        halfband_t();
        T process(const T& input){
            T output = input;

            output = (a->process(input) + oldout) * 0.5;
            oldout = b->process(input);

            return output;
        }
    }halfband;

    struct halfband_cascade {
        uint32_t i = 0;
        halfband_t hb;
        halfband_cascade *next;
        Delegate<5> Output;
    
        void process(const T& input) {
            
            T output = hb.process(input);
            
            if ( !(++i % 2) ) {
                Output(output);
                if(next) next->process(output);
            }
        }
    };

    /* cascade of half band filters, for 2^n times oversampling */
    halfband_cascade halfs[n];
    T buf[1 << n];
public:
    Halfband() {
        halfs[n - 1].next = nullptr;
        for (int i = 1; i < n; i++) {
            halfs[i - 1].next = &halfs[i];
        }
        for (auto& data : buf) data = 0.;
    }
    inline Delegate<5>& Subscribe(const int i) {
        assert(i < n);

        return halfs[i].Output;
    }
    // should receive an array of 2**n inputs to process
    T process(const T* input);
    void process(const T& input) {
        halfs[0].process(input);
    }
    //double process_half_cascade_v2(const double* input);
    //double process_half_cascade(halfband* half, const double input);
};

template<int n, int order, typename T, bool steep>
Halfband<n, order, T, steep>::halfband::halfband_t()
{
    if (steep)
    {
        if (order == 12)    //rejection=104dB, transition band=0.01
        {
            constexpr T a_coefficients[6] =
            { 0.036681502163648017
            ,0.2746317593794541
            ,0.56109896978791948
            ,0.769741833862266
            ,0.8922608180038789
            ,0.962094548378084
            };

            constexpr T b_coefficients[6] =
            { 0.13654762463195771
            ,0.42313861743656667
            ,0.6775400499741616
            ,0.839889624849638
            ,0.9315419599631839
            ,0.9878163707328971
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 6);
            b = std::make_unique<allpass_cascade>(b_coefficients, 6);
        }
        else if (order == 10)    //rejection=86dB, transition band=0.01
        {
            constexpr T a_coefficients[5] =
            { 0.051457617441190984
            ,0.35978656070567017
            ,0.6725475931034693
            ,0.8590884928249939
            ,0.9540209867860787
            };

            constexpr T b_coefficients[5] =
            { 0.18621906251989334
            ,0.529951372847964
            ,0.7810257527489514
            ,0.9141815687605308
            ,0.985475023014907
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 5);
            b = std::make_unique<allpass_cascade>(b_coefficients, 5);

        }
        else if (order == 8)    //rejection=69dB, transition band=0.01
        {
            constexpr T a_coefficients[4] =
            { 0.07711507983241622
            ,0.4820706250610472
            ,0.7968204713315797
            ,0.9412514277740471
            };

            constexpr T b_coefficients[4] =
            { 0.2659685265210946
            ,0.6651041532634957
            ,0.8841015085506159
            ,0.9820054141886075
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 4);
            b = std::make_unique<allpass_cascade>(b_coefficients, 4);

        }
        else if (order == 6)    //rejection=51dB, transition band=0.01
        {
            constexpr T a_coefficients[3] =
            { 0.1271414136264853
            ,0.6528245886369117
            ,0.9176942834328115
            };

            constexpr T b_coefficients[3] =
            { 0.40056789819445626
            ,0.8204163891923343
            ,0.9763114515836773
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 3);
            b = std::make_unique<allpass_cascade>(b_coefficients, 3);
        }
        else if (order == 4)    //rejection=53dB,transition band=0.05
        {
            constexpr T a_coefficients[2] =
            { 0.12073211751675449
            ,0.6632020224193995
            };

            constexpr T b_coefficients[2] =
            { 0.3903621872345006
            ,0.890786832653497
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 2);
            b = std::make_unique<allpass_cascade>(b_coefficients, 2);
        }

        else    //order=2, rejection=36dB, transition band=0.1
        {
            constexpr T a_coefficients[1] = { 0.23647102099689224 };
            constexpr T b_coefficients[1] = { 0.7145421497126001 };

            a = std::make_unique<allpass_cascade>(a_coefficients, 1);
            b = std::make_unique<allpass_cascade>(b_coefficients, 1);
        }
    }
    else    //softer slopes, more attenuation and less stopband ripple
    {
        if (order == 12)    //rejection=150dB, transition band=0.05
        {
            constexpr T a_coefficients[6] =
            { 0.01677466677723562
            ,0.13902148819717805
            ,0.3325011117394731
            ,0.53766105314488
            ,0.7214184024215805
            ,0.8821858402078155
            };

            constexpr T b_coefficients[6] =
            { 0.06501319274445962
            ,0.23094129990840923
            ,0.4364942348420355
            ,0.6329609551399348
            ,0.80378086794111226
            ,0.9599687404800694
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 6);
            b = std::make_unique<allpass_cascade>(b_coefficients, 6);
        }
        else if (order == 10)    //rejection=133dB, transition band=0.05
        {
            constexpr T a_coefficients[5] =
            { 0.02366831419883467
            ,0.18989476227180174
            ,0.43157318062118555
            ,0.6632020224193995
            ,0.860015542499582
            };

            constexpr T b_coefficients[5] =
            { 0.09056555904993387
            ,0.3078575723749043
            ,0.5516782402507934
            ,0.7652146863779808
            ,0.95247728378667541
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 5);
            b = std::make_unique<allpass_cascade>(b_coefficients, 5);
        }
        else if (order == 8)    //rejection=106dB, transition band=0.05
        {
            constexpr T a_coefficients[4] =
            { 0.03583278843106211
            ,0.2720401433964576
            ,0.5720571972357003
            ,0.827124761997324
            };

            constexpr T b_coefficients[4] =
            { 0.1340901419430669
            ,0.4243248712718685
            ,0.7062921421386394
            ,0.9415030941737551
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 4);
            b = std::make_unique<allpass_cascade>(b_coefficients, 4);
        }
        else if (order == 6)    //rejection=80dB, transition band=0.05
        {
            constexpr T a_coefficients[3] =
            { 0.06029739095712437
            ,0.4125907203610563
            ,0.7727156537429234
            };

            constexpr T b_coefficients[3] =
            { 0.21597144456092948
            ,0.6043586264658363
            ,0.9238861386532906
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 3);
            b = std::make_unique<allpass_cascade>(b_coefficients, 3);
        }
        else if (order == 4)    //rejection=70dB,transition band=0.1
        {
            constexpr T a_coefficients[2] =
            { 0.07986642623635751
            ,0.5453536510711322
            };

            constexpr T b_coefficients[2] =
            { 0.28382934487410993
            ,0.8344118914807379
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 2);
            b = std::make_unique<allpass_cascade>(b_coefficients, 2);
        }

        else    //order=2, rejection=36dB, transition band=0.1
        {
            constexpr T a_coefficients[1] = { 0.23647102099689224 };
            constexpr T b_coefficients[1] = { 0.7145421497126001 };

            a = std::make_unique<allpass_cascade>(a_coefficients, 1);
            b = std::make_unique<allpass_cascade>(b_coefficients, 1);
        }
    }
    oldout = 0;
}

#if 0
// should receive an array of 2**n inputs to process
// will return a single value following a cascade of filter->decimate steps
template<int n, int order, typename T, bool steep>
double Halfband<n, order, T, steep>::process(const T* input) {
    int i, j, k;
    T b;
    int nbuf;
    nbuf = 1 << n;

    // cascade 1
    k = 0;
    for (j = 0; j < nbuf; j += 2) {
        b = halfs[0].process(input[j]);
        b = halfs[0].process(input[j + 1]);
        // decimate and store
        buf[k++] = b;
    }
    // next buf size
    nbuf >>= 1;

    // subsequent cascades
    for (i = 1; i < n; i++) {
        k = 0;
        for (j = 0; j < nbuf; j += 2) {
            // note that we process the buffer in-place, writing
            // the decimated sample at the start of the buffer, always
            // behind the position read
            b = halfs[i].process(buf[j]);
            b = halfs[i].process(buf[j + 1]);
            buf[k++] = b;
        }
        // next buf size
        nbuf >>= 1;
    }

    // final fully decimated sample
    return buf[0];
}
#elif 1
    template<int n, int order, typename T, bool steep>
    T Halfband<n, order, T, steep>::process(const T* input) {
        T b;
        int nbuf = 1 << n;
        const T* input_buf = input;

        for (int i = 0; i < n; i++) {
            int k = 0;
            for (int j = 0; j < nbuf; j += 2) {
                // note that we process the buffer in-place, writing
                // the decimated sample at the start of the buffer, always
                // behind the position read
                b = halfs[i].hb.process(input_buf[j]);
                b = halfs[i].hb.process(input_buf[j + 1]);
                buf[k++] = b;
            }
            input_buf = buf;
            // next buf size
            nbuf >>= 1;
        }

        return input_buf[0];
    }
#endif
