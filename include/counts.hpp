#ifndef CHARON_COUNTS_H
#define CHARON_COUNTS_H

#pragma once

#include <iostream>
#include <vector>
#include <cassert>

// Extremly simple Counts class - represents a lower diagonal square matrix including the diagonal,
// expecting [row,col]==[col,row] so can store it once not twice
template<class T>
class Counts {
public:
    Counts() {};

    Counts(size_t size) {
        assert(size > 0);
        mRows = size;
        mData.resize(mRows * (mRows + 1) / 2);
        mData.shrink_to_fit();
    }

    void set_size(size_t size) {
        assert(size > 0);
        mRows = size;
        mData.resize(mRows * (mRows + 1) / 2);
        mData.shrink_to_fit();
    }

    T &operator()(size_t row, size_t col) {
        assert(row < mRows);
        assert(col < mRows);
        if (row < col) {
            col, row = row, col;
        }
        return mData[row * (row + 1) / 2 + col];
    }

    const T &operator()(size_t row, size_t col) const {
        assert(row >= 0 && row < mRows);
        assert(col >= 0 && row < mRows);
        return mData[row * (row + 1) / 2 + col];
    }

    size_t rows() const noexcept {
        return mRows;
    }

protected:
    size_t mRows;
    std::vector<T> mData;
};

#endif // CHARON_COUNTS_H
