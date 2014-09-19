#include "RandT3.h"
#include "TRandom3.h"

TRandom3 *randGen = new TRandom3(0);

double RandT3::randomOne(){
    return randGen->Uniform();
}
