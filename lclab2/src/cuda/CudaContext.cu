#include <stdio.h>


// Don't include anything to do with lclab2... just create the test func
namespace LC { namespace Cuda {

    extern "C" void test() {
        printf("Hello from cuda compiled function 'test'!\n");
    }

}}