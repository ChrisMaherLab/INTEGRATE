#pragma once

#include <divsufsort.h>
#include <divsufsort64.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// FIXME: pay attention to byte order in IO
namespace {
    uint64_t endian_swap(uint64_t x) {
        return x;
    }

    uint32_t endian_swap(uint32_t x) {
        return x;
    }

    void divsufsort_wrapper(std::string const& seq, uint64_t* data) {
        sauchar_t const* seqptr = reinterpret_cast<sauchar_t const*>(seq.data());
        saidx64_t* arr = reinterpret_cast<saidx64_t*>(data);
        if (divsufsort64(seqptr, arr, seq.size()) != 0) {
            throw std::runtime_error("Failed to build suffix array!");
        }
    }

    void divsufsort_wrapper(std::string const& seq, uint32_t* data) {
        sauchar_t const* seqptr = reinterpret_cast<sauchar_t const*>(seq.data());
        saidx_t* arr = reinterpret_cast<saidx_t*>(data);
        if (divsufsort(seqptr, arr, seq.size()) != 0) {
            throw std::runtime_error("Failed to build suffix array!");
        }
    }
}

template<typename T>
class SuffixArray {
public:
    typedef T value_type;

    SuffixArray() {}

    explicit SuffixArray(std::string const& sequence)
        : _sa(sequence.size(), 0)
        , _decimation(1)
    {
        divsufsort_wrapper(sequence, _sa.data());
    }

    void toStream(std::ostream& out) const {
        value_type sz = size();
        out.write(reinterpret_cast<char*>(&sz), sizeof(sz));
        out.write(reinterpret_cast<char const*>(_sa.data()), sizeof(value_type) * size());
    }

    static SuffixArray fromStream(std::istream& in) {
        value_type sz(0);
        in.read(reinterpret_cast<char*>(&sz), sizeof(sz));

        SuffixArray rv;
        rv._sa.resize(sz, 0);
        in.read(reinterpret_cast<char*>(rv._sa.data()), sizeof(value_type) * rv.size());

        return rv;
    }

    value_type operator[](size_t idx) const {
        return _sa[idx];
    }

    size_t size() const {
        return _sa.size();
    }

    bool empty() const {
        return _sa.empty();
    }

    size_t capacity() const {
        return _sa.capacity();
    }

    bool operator==(SuffixArray const& rhs) const {
        return _sa == rhs._sa;
    }

    void decimate(size_t factor) {
        if (_decimation != 1) {
            throw std::runtime_error("Suffix array is already decimated!");
        }

        if (factor <= 1) {
            return; // you're so silly!
        }

        _decimation = factor;

        size_t o = 1;
        for (size_t i = factor; i < _sa.size(); i += factor) {
            _sa[o++] = _sa[i];
        }
        _sa.resize(_sa.size() / factor + 1);
        std::vector<value_type>(_sa).swap(_sa);
    }

    friend std::ostream& operator<<(std::ostream& s, SuffixArray const& sa) {
        for (size_t i = 0; i < sa.size(); ++i) {
            s << i << ": " << sa[i] << "\n";
        }
        return s;
    }

    size_t decimation() const {
        return _decimation;
    }

    value_type const* data() const {
        return _sa.data();
    }

private:
    std::vector<value_type> _sa;
    size_t _decimation;
};
