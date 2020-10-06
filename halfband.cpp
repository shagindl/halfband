/* Allpass cascade half band filtering. 
From http://www.musicdsp.org/showone.php?id=39
poretd to plain C, and with a simple cascade
function to allow n times oversampling for n a power of 2
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "halfband.h"

#if 0
template<int n, int order, typename T, bool steep>
T Halfband<n, order, T, steep>::process(const T& input) {
    int i, j, k;
    T b;
    int nbuf = 1 << n;

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
#endif
#if 0
// should receive an array of 2**n inputs to process
// will return a single value following a cascade of filter->decimate steps
template<int n, int order, typename T>
double Halfband<n, order, T>::process_half_cascade_v2(const double *input){
    int i, j, k;
    double b;
    int nbuf = 1<<cascade->n;
    const double* buf = input;
    
    for(i = 0; i < cascade->n; i++){
       k = 0;
       for(j = 0; j < nbuf; j += 2){
           // note that we process the buffer in-place, writing
           // the decimated sample at the start of the buffer, always
           // behind the position read
           b = process_halfband(halfs[i], buf[j]);
           b = process_halfband(cascade->halfs[i], buf[j+1]);
           cascade->buf[k++] = b;
       }
       buf = cascade->buf;
       // next buf size
       nbuf >>= 1;        
    }
    
    // final fully decimated sample
    return buf[0];
}
#endif
#if 0
template<int n, int order, typename T, bool steep>
T Halfband<n, order, T, steep>::process(const T& input){
    T b;
    int nbuf = 1 << n;
    int &i = cascade->i, &j = half->j, &k = cascade->k;
    int &mode = cascade->mode;


    int nbuf = 1 << n;
    const T* input_buf = input;

    for (int i = 0; i < n; i++) {
        int k = 0;
        for (int j = 0; j < nbuf; j += 2) {
            // note that we process the buffer in-place, writing
            // the decimated sample at the start of the buffer, always
            // behind the position read
            halfs[i].process(input_buf[j]);
            buf[k++] = halfs[i].process(input_buf[j + 1]);
        }
        input_buf = buf;
        // next buf size
        nbuf >>= 1;
    }

    return input_buf[0];



    return b;

    switch( mode ) {
    case 0:
        mode++;
        halfs[i].process(input);
        break;
    case 1:
        buf[k++] = halfs[i].process(input);
        j += 2;
        if( j < nbuf ) mode = 0;
        else {
            j = 0;
            nbuf >>= 1;
        }
        break;

        if( i < cascade->n ) {
            k = 0;
            if( j < nbuf ) {
                // note that we process the buffer in-place, writing
                // the decimated sample at the start of the buffer, always
                // behind the position read
                b = process_halfband(cascade->halfs[i], buf[j]);
                b = process_halfband(cascade->halfs[i], buf[j + 1]);
                cascade->buf[k++] = b;
            }
            buf = cascade->buf;
            // next buf size
            nbuf >>= 1;
        }
    }
    // final fully decimated sample
    return cascade->buf[0];
}
double process_half_cascade(half_cascade* cascade, const  double input) {
    return 0.;
}
#endif