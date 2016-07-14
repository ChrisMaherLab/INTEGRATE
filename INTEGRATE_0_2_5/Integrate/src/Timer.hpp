#pragma once

#include <ctime>

class Timer {
public:

    Timer() {
        reset();
    }

    void reset() {
        _start = clock();
    }

    double elapsed() const {
        return (clock() - _start) / double(CLOCKS_PER_SEC);
    }

private:
    clock_t _start;
};
