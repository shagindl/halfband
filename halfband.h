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

            output = (a->process(input) + oldout) * 0.5f;
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
public:
    Halfband() {
        halfs[n - 1].next = nullptr;
        for (int i = 1; i < n; i++) {
            halfs[i - 1].next = &halfs[i];
        }
    }
    inline Delegate<5>& Subscribe(const int i) {
        assert(i < n);

        return halfs[i].Output;
    }
    // should receive an array of 2**n inputs to process
    __ramfunc void process(const T& input) {
        halfs[0].process(input);
    }
};

template<int n, int order, typename T, bool steep>
Halfband<n, order, T, steep>::halfband::halfband_t()
{
    if (steep)
    {
        if (order == 12)    //rejection=104dB, transition band=0.01
        {
            constexpr T a_coefficients[6] =
            { 0.036681502163648017f
            ,0.2746317593794541f
            ,0.56109896978791948f
            ,0.769741833862266f
            ,0.8922608180038789f
            ,0.962094548378084f
            };

            constexpr T b_coefficients[6] =
            { 0.13654762463195771f
            ,0.42313861743656667f
            ,0.6775400499741616f
            ,0.839889624849638f
            ,0.9315419599631839f
            ,0.9878163707328971f
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 6);
            b = std::make_unique<allpass_cascade>(b_coefficients, 6);
        }
        else if (order == 10)    //rejection=86dB, transition band=0.01
        {
            constexpr T a_coefficients[5] =
            { 0.051457617441190984f
            ,0.35978656070567017f
            ,0.6725475931034693f
            ,0.8590884928249939f
            ,0.9540209867860787f
            };

            constexpr T b_coefficients[5] =
            { 0.18621906251989334f
            ,0.529951372847964f
            ,0.7810257527489514f
            ,0.9141815687605308f
            ,0.985475023014907f
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 5);
            b = std::make_unique<allpass_cascade>(b_coefficients, 5);

        }
        else if (order == 8)    //rejection=69dB, transition band=0.01
        {
            constexpr T a_coefficients[4] =
            { 0.07711507983241622f
            ,0.4820706250610472f
            ,0.7968204713315797f
            ,0.9412514277740471f
            };

            constexpr T b_coefficients[4] =
            { 0.2659685265210946f
            ,0.6651041532634957f
            ,0.8841015085506159f
            ,0.9820054141886075f
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 4);
            b = std::make_unique<allpass_cascade>(b_coefficients, 4);

        }
        else if (order == 6)    //rejection=51dB, transition band=0.01
        {
            constexpr T a_coefficients[3] =
            { 0.1271414136264853f
            ,0.6528245886369117f
            ,0.9176942834328115f
            };

            constexpr T b_coefficients[3] =
            { 0.40056789819445626f
            ,0.8204163891923343f
            ,0.9763114515836773f
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 3);
            b = std::make_unique<allpass_cascade>(b_coefficients, 3);
        }
        else if (order == 4)    //rejection=53dB,transition band=0.05
        {
            constexpr T a_coefficients[2] =
            { 0.12073211751675449f
            ,0.6632020224193995f
            };

            constexpr T b_coefficients[2] =
            { 0.3903621872345006f
            ,0.890786832653497f
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 2);
            b = std::make_unique<allpass_cascade>(b_coefficients, 2);
        }

        else    //order=2, rejection=36dB, transition band=0.1
        {
            constexpr T a_coefficients[1] = { 0.23647102099689224f };
            constexpr T b_coefficients[1] = { 0.7145421497126001f };

            a = std::make_unique<allpass_cascade>(a_coefficients, 1);
            b = std::make_unique<allpass_cascade>(b_coefficients, 1);
        }
    }
    else    //softer slopes, more attenuation and less stopband ripple
    {
        if (order == 12)    //rejection=150dB, transition band=0.05
        {
            constexpr T a_coefficients[6] =
            { 0.01677466677723562f
            ,0.13902148819717805f
            ,0.3325011117394731f
            ,0.53766105314488f
            ,0.7214184024215805f
            ,0.8821858402078155f
            };

            constexpr T b_coefficients[6] =
            { 0.06501319274445962f
            ,0.23094129990840923f
            ,0.4364942348420355f
            ,0.6329609551399348f
            ,0.80378086794111226f
            ,0.9599687404800694f
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 6);
            b = std::make_unique<allpass_cascade>(b_coefficients, 6);
        }
        else if (order == 10)    //rejection=133dB, transition band=0.05
        {
            constexpr T a_coefficients[5] =
            { 0.02366831419883467f
            ,0.18989476227180174f
            ,0.43157318062118555f
            ,0.6632020224193995f
            ,0.860015542499582f
            };

            constexpr T b_coefficients[5] =
            { 0.09056555904993387f
            ,0.3078575723749043f
            ,0.5516782402507934f
            ,0.7652146863779808f
            ,0.95247728378667541f
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 5);
            b = std::make_unique<allpass_cascade>(b_coefficients, 5);
        }
        else if (order == 8)    //rejection=106dB, transition band=0.05
        {
            constexpr T a_coefficients[4] =
            { 0.03583278843106211f
            ,0.2720401433964576f
            ,0.5720571972357003f
            ,0.827124761997324f
            };

            constexpr T b_coefficients[4] =
            { 0.1340901419430669f
            ,0.4243248712718685f
            ,0.7062921421386394f
            ,0.9415030941737551f
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 4);
            b = std::make_unique<allpass_cascade>(b_coefficients, 4);
        }
        else if (order == 6)    //rejection=80dB, transition band=0.05
        {
            constexpr T a_coefficients[3] =
            { 0.06029739095712437f
            ,0.4125907203610563f
            ,0.7727156537429234f
            };

            constexpr T b_coefficients[3] =
            { 0.21597144456092948f
            ,0.6043586264658363f
            ,0.9238861386532906f
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 3);
            b = std::make_unique<allpass_cascade>(b_coefficients, 3);
        }
        else if (order == 4)    //rejection=70dB,transition band=0.1
        {
            constexpr T a_coefficients[2] =
            { 0.07986642623635751f
            ,0.5453536510711322f
            };

            constexpr T b_coefficients[2] =
            { 0.28382934487410993f
            ,0.8344118914807379f
            };

            a = std::make_unique<allpass_cascade>(a_coefficients, 2);
            b = std::make_unique<allpass_cascade>(b_coefficients, 2);
        }

        else    //order=2, rejection=36dB, transition band=0.1
        {
            constexpr T a_coefficients[1] = { 0.23647102099689224f };
            constexpr T b_coefficients[1] = { 0.7145421497126001f };

            a = std::make_unique<allpass_cascade>(a_coefficients, 1);
            b = std::make_unique<allpass_cascade>(b_coefficients, 1);
        }
    }
    oldout = 0;
}
