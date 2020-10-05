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
// should receive an array of 2**n inputs to process
// will return a single value following a cascade of filter->decimate steps
template<int n, int order, typename T>
double Halfband<n, order, T>::process_half_cascade(halfband *half, const  double input){
    double b;
    int nbuf = 1<<n;
    int &i = cascade->i, &j = half->j, &k = cascade->k;
    int &mode = cascade->mode;

    b = process_halfband(half, input), j++;
    if( j % 2 == 0 ) {
        b = process_half_cascade(cascade, b);
    }

    return b;

    switch( mode ) {
    case 0:
        mode++;
        process_halfband(cascade->halfs[i], input);
        break;
    case 1:
        cascade->buf[k++] = process_halfband(cascade->halfs[i], input);
        j += 2;
        if( j < nbuf ) mode = 0;
        else {
            j = 0;

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