#include <iostream>

#include "doctest.h"

TEST_CASE("Hello world") {
    INFO("Hello there!");
    CHECK(1 == 2);
}
