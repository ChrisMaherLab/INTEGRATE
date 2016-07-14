#pragma once

#include "SuffixArray.h"

#include <stdexcept>
#include <stdint.h>
#include <string>

template<typename T>
class BWIndex {
public:
    BWIndex(std::string const& sequence, SuffixArray<T> const& sa)
        : _bwt(sequence.size(), 0)
        , _sa(sa)
    {
        _bwt[0] = *sequence.rbegin();
        size_t i;
        for (i = 0; sa[i] != 0; ++i) {
            _bwt[i + 1] = sequence[sa[i] - 1];
        }

        for(++i; i < sequence.size(); ++i) {
            _bwt[i] = sequence[sa[i] - 1];
        }
    }

    std::string const& bwt() const {
        return _bwt;
    }

private:
    std::string _bwt;
    SuffixArray<T> const& _sa;
};
