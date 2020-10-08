/* Allpass cascade half band filtering. 
From http://www.musicdsp.org/showone.php?id=39
poretd to plain C, and with a simple cascade
function to allow n times oversampling for n a power of 2
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "halfband.h"

// 
