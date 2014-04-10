#include <stdio.h>

#include "clar.h"

#include "test.h"

unsigned long global_test_counter;

int main(int argc, char *argv[])
{
    global_test_counter = 0u;
    /* Run the test suite */
    int result = clar_test(argc, argv);
    if (!result) {
        printf("Ran %lu test%s and %s passed!\n", global_test_counter, global_test_counter == 1 ? "" : "", global_test_counter == 1 ? "it" : "all");
    } else {
        printf("Aw snap! Some tests failed...\n");
    }
    return result;
}
