#pragma once

#include <string.h> 
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <inttypes.h>
#include <iostream>
#include <string>
#include <stdexcept>

class RingMatrixException : public std::runtime_error 
{
public:
    explicit RingMatrixException(const std::string& message)
        : std::runtime_error(message) {}
};

class RingMatrix
{
public:
    int64_t rows;
    int64_t cols;
    int64_t mod;
    int64_t* data;

private:
    inline size_t total_elements() const
    {
        return (size_t)rows * (size_t)cols;
    }

    int64_t modInv(int64_t a) const
    {
        int64_t m = mod;
        int64_t aa = a % m;
        if (aa < 0) aa += m;
        if (aa == 0)
        {
            throw RingMatrixException("RingMatrix::modInv: zero element has no inverse modulo " + std::to_string(mod));
        }

        int64_t r0 = m, r1 = aa;
        int64_t s0 = 0, s1 = 1;

        while (r1 != 0)
        {
            int64_t q = r0 / r1;
            int64_t r2 = r0 - q * r1; r0 = r1; r1 = r2;
            int64_t s2 = s0 - q * s1; s0 = s1; s1 = s2;
        }

        if (r0 != 1)
        {
            throw RingMatrixException("RingMatrix::modInv: element " + std::to_string(a) +
                " is not invertible modulo " + std::to_string(mod));
        }

        int64_t inv = s0 % m;
        if (inv < 0) inv += m;
        return inv;
    }

    void swap_rows(int64_t r1, int64_t r2)
    {
        if (r1 == r2) return;
        for (int64_t j = 0; j < cols; ++j)
        {
            size_t i1 = index_pos(r1, j);
            size_t i2 = index_pos(r2, j);
            int64_t tmp = data[i1];
            data[i1] = data[i2];
            data[i2] = tmp;
        }
    }

    void swap_cols(int64_t c1, int64_t c2)
    {
        if (c1 == c2) return;
        for (int64_t i = 0; i < rows; ++i)
        {
            size_t i1 = index_pos(i, c1);
            size_t i2 = index_pos(i, c2);
            int64_t tmp = data[i1];
            data[i1] = data[i2];
            data[i2] = tmp;
        }
    }

    static void mat_mul_mod(const int64_t* A, const int64_t* B, int64_t* OUT, int64_t n, int64_t mod)
    {
        uint64_t M = (uint64_t)mod;
        uint64_t max_prod = (M - 1ULL) * (M - 1ULL);
        uint64_t threshold = UINT64_MAX - max_prod;

        for (int64_t i = 0; i < n; ++i)
        {
            const int64_t* Arow = A + (size_t)i * (size_t)n;
            int64_t* OutRow = OUT + (size_t)i * (size_t)n;
            for (int64_t j = 0; j < n; ++j)
            {
                uint64_t acc = 0;
                for (int64_t k = 0; k < n; ++k)
                {
                    uint64_t aval = (uint64_t)Arow[k];
                    uint64_t bval = (uint64_t)B[(size_t)k * (size_t)n + (size_t)j];
                    uint64_t prod = aval * bval;
                    if (acc > threshold) acc %= M;
                    acc += prod;
                }
                OutRow[j] = (int64_t)(acc % M);
            }
        }
    }

public:
    inline int64_t normalize(long long x) const
    {
        long long r = x % mod;
        if (r < 0) r += mod;
        return (int64_t)r;
    }

    inline size_t index_pos(int64_t i, int64_t j) const
    {
        return (size_t)i * (size_t)cols + (size_t)j;
    }

    static int digits_int64(int64_t x)
    {
        long long t = x;
        int d = (t <= 0);
        if (t < 0) t = -t;
        while (t) { t /= 10; ++d; }
        return d;
    }

    RingMatrix(int64_t r, int64_t c, int64_t m)
        : rows(r), cols(c), mod(m), data(nullptr)
    {
        if (r <= 0 || c <= 0 || m <= 0)
        {
            throw RingMatrixException("RingMatrix constructor: dimensions and modulus must be positive (got rows=" +
                std::to_string(r) + ", cols=" + std::to_string(c) + ", mod=" + std::to_string(m) + ")");
        }
        if ((size_t)r > SIZE_MAX / (size_t)c)
        {
            throw RingMatrixException("RingMatrix constructor: dimensions too large (rows=" +
                std::to_string(r) + ", cols=" + std::to_string(c) + ")");
        }
        size_t n = (size_t)r * (size_t)c;
        data = new int64_t[n];
        memset(data, 0, sizeof(int64_t) * n);
    }

    RingMatrix(const RingMatrix& o)
        : rows(o.rows), cols(o.cols), mod(o.mod), data(nullptr)
    {
        if (rows <= 0 || cols <= 0)
        {
            throw RingMatrixException("RingMatrix copy constructor: invalid dimensions (rows=" +
                std::to_string(rows) + ", cols=" + std::to_string(cols) + ")");
        }
        if ((size_t)rows > SIZE_MAX / (size_t)cols)
        {
            throw RingMatrixException("RingMatrix copy constructor: dimensions too large (rows=" +
                std::to_string(rows) + ", cols=" + std::to_string(cols) + ")");
        }
        size_t n = (size_t)rows * (size_t)cols;
        data = new int64_t[n];
        memcpy(data, o.data, sizeof(int64_t) * n);
    }

    friend void swap(RingMatrix& a, RingMatrix& b) noexcept
    {
        int64_t tmp_i64;
        int64_t* tmp_p;

        tmp_i64 = a.rows; a.rows = b.rows; b.rows = tmp_i64;
        tmp_i64 = a.cols; a.cols = b.cols; b.cols = tmp_i64;
        tmp_i64 = a.mod;  a.mod = b.mod;  b.mod = tmp_i64;
        tmp_p = a.data; a.data = b.data; b.data = tmp_p;
    }

    RingMatrix& operator=(RingMatrix other)
    {
        swap(*this, other);
        return *this;
    }

    ~RingMatrix()
    {
        delete[] data;
    }

    inline int64_t& at(int64_t i, int64_t j)
    {
        if ((uint64_t)i >= (uint64_t)rows || (uint64_t)j >= (uint64_t)cols)
        {
            throw RingMatrixException("RingMatrix::at: index out of bounds (i=" +
                std::to_string(i) + ", j=" + std::to_string(j) +
                ", rows=" + std::to_string(rows) + ", cols=" + std::to_string(cols) + ")");
        }
        return data[index_pos(i, j)];
    }

    inline int64_t at(int64_t i, int64_t j) const
    {
        if ((uint64_t)i >= (uint64_t)rows || (uint64_t)j >= (uint64_t)cols)
        {
            throw RingMatrixException("RingMatrix::at const: index out of bounds (i=" +
                std::to_string(i) + ", j=" + std::to_string(j) +
                ", rows=" + std::to_string(rows) + ", cols=" + std::to_string(cols) + ")");
        }
        return data[index_pos(i, j)];
    }

    RingMatrix operator+(const RingMatrix& b) const
    {
        if (rows != b.rows || cols != b.cols || mod != b.mod)
        {
            throw RingMatrixException("RingMatrix::operator+: matrix dimensions or modulus mismatch (" +
                std::to_string(rows) + "x" + std::to_string(cols) + " mod " + std::to_string(mod) +
                " vs " + std::to_string(b.rows) + "x" + std::to_string(b.cols) + " mod " + std::to_string(b.mod) + ")");
        }
        RingMatrix r(rows, cols, mod);
        size_t n = (size_t)rows * (size_t)cols;
        for (size_t i = 0; i < n; ++i)
        {
            int64_t s = data[i] + b.data[i];
            if (s >= mod) s -= mod;
            r.data[i] = s;
        }
        return r;
    }

    RingMatrix operator-(const RingMatrix& b) const
    {
        if (rows != b.rows || cols != b.cols || mod != b.mod)
        {
            throw RingMatrixException("RingMatrix::operator-: matrix dimensions or modulus mismatch (" +
                std::to_string(rows) + "x" + std::to_string(cols) + " mod " + std::to_string(mod) +
                " vs " + std::to_string(b.rows) + "x" + std::to_string(b.cols) + " mod " + std::to_string(b.mod) + ")");
        }
        RingMatrix r(rows, cols, mod);
        size_t n = (size_t)rows * (size_t)cols;
        for (size_t i = 0; i < n; ++i)
        {
            int64_t d = data[i] - b.data[i];
            if (d < 0) d += mod;
            r.data[i] = d;
        }
        return r;
    }

    RingMatrix operator*(int64_t k) const
    {
        if (mod <= 0)
        {
            throw RingMatrixException("RingMatrix::operator*: non-positive modulus (" + std::to_string(mod) + ")");
        }
        RingMatrix r(rows, cols, mod);
        int64_t kk = normalize(k);
        size_t n = (size_t)rows * (size_t)cols;
        uint64_t M = (uint64_t)mod;
        for (size_t i = 0; i < n; ++i)
        {
            uint64_t aval = (uint64_t)(data[i] < 0 ? (uint64_t)(data[i] + (int64_t)mod) : (uint64_t)data[i]);
            uint64_t prod = aval * (uint64_t)kk;
            r.data[i] = (int64_t)(prod % M);
        }
        return r;
    }

    friend RingMatrix operator*(int64_t k, const RingMatrix& m)
    {
        return m * k;
    }

    RingMatrix operator*(const RingMatrix& b) const
    {
        if (cols != b.rows)
        {
            throw RingMatrixException("RingMatrix::operator*: matrix dimensions incompatible for multiplication (" +
                std::to_string(rows) + "x" + std::to_string(cols) + " * " +
                std::to_string(b.rows) + "x" + std::to_string(b.cols) + ")");
        }
        if (mod != b.mod)
        {
            throw RingMatrixException("RingMatrix::operator*: modulus mismatch (" +
                std::to_string(mod) + " vs " + std::to_string(b.mod) + ")");
        }
        RingMatrix r(rows, b.cols, mod);

        uint64_t M = (uint64_t)mod;
        uint64_t max_prod = (M - 1ULL) * (M - 1ULL);
        uint64_t threshold = UINT64_MAX - max_prod;

        for (int64_t i = 0; i < rows; ++i)
        {
            const int64_t* a_row = data + (size_t)i * (size_t)cols;
            for (int64_t j = 0; j < b.cols; ++j)
            {
                uint64_t acc = 0;
                size_t b_col_index = (size_t)j;
                for (int64_t k = 0; k < cols; ++k)
                {
                    uint64_t aval = (uint64_t)(a_row[k] < 0 ? (uint64_t)(a_row[k] + mod) : (uint64_t)a_row[k]);
                    uint64_t bval = (uint64_t)(b.data[(size_t)k * (size_t)b.cols + b_col_index] < 0 ? (uint64_t)(b.data[(size_t)k * (size_t)b.cols + b_col_index] + b.mod) : (uint64_t)b.data[(size_t)k * (size_t)b.cols + b_col_index]);
                    uint64_t prod = aval * bval;
                    if (acc > threshold) acc %= M;
                    acc += prod;
                }
                r.data[(size_t)i * (size_t)b.cols + (size_t)j] = (int64_t)(acc % M);
            }
        }
        return r;
    }

    RingMatrix transpose() const
    {
        RingMatrix r(cols, rows, mod);
        for (int64_t i = 0; i < rows; ++i)
            for (int64_t j = 0; j < cols; ++j)
                r.data[(size_t)j * (size_t)r.cols + (size_t)i] = data[index_pos(i, j)];
        return r;
    }

    RingMatrix submatrix(int64_t top, int64_t bottom, int64_t left, int64_t right) const
    {
        if (top < 0 || left < 0 || bottom >= rows || right >= cols || top > bottom || left > right)
        {
            throw RingMatrixException("RingMatrix::submatrix: invalid submatrix indices (top=" +
                std::to_string(top) + ", bottom=" + std::to_string(bottom) +
                ", left=" + std::to_string(left) + ", right=" + std::to_string(right) +
                ", rows=" + std::to_string(rows) + ", cols=" + std::to_string(cols) + ")");
        }
        RingMatrix r(bottom - top + 1, right - left + 1, mod);
        for (int64_t i = 0; i < r.rows; ++i)
            for (int64_t j = 0; j < r.cols; ++j)
                r.data[(size_t)i * (size_t)r.cols + (size_t)j] = data[index_pos(top + i, left + j)];
        return r;
    }

    char* toString() const
    {
        size_t n = (size_t)rows * (size_t)cols;
        size_t len = (size_t)digits_int64(mod) + (size_t)digits_int64(rows) + (size_t)digits_int64(cols) + 2;
        for (size_t i = 0; i < n; ++i)
            len += (size_t)digits_int64(data[i]) + 1;

        char* s = new char[len + 1];
        int written = snprintf(s, (size_t)len + 1, "%" PRId64 " %" PRId64 " %" PRId64, (int64_t)mod, (int64_t)rows, (int64_t)cols);
        if (written < 0)
        {
            delete[] s;
            throw RingMatrixException("RingMatrix::toString: snprintf failed");
        }
        size_t p = (size_t)written;

        for (size_t i = 0; i < n; ++i)
        {
            written = snprintf(s + p, ((size_t)len + 1) - p, " %" PRId64, (int64_t)data[i]);
            if (written < 0)
            {
                delete[] s;
                throw RingMatrixException("RingMatrix::toString: snprintf failed");
            }
            p += (size_t)written;
        }
        s[p] = '\0';
        return s;
    }

    static RingMatrix fromString(const char* s)
    {
        if (!s)
        {
            throw RingMatrixException("RingMatrix::fromString: null pointer");
        }
        const char* p = s;
        char* endptr = nullptr;
        errno = 0;
        long long m_long = strtoll(p, &endptr, 10);
        if (endptr == p || errno == ERANGE)
        {
            throw RingMatrixException("RingMatrix::fromString: invalid or out-of-range modulus");
        }
        p = endptr;
        errno = 0;
        long long r_long = strtoll(p, &endptr, 10);
        if (endptr == p || errno == ERANGE)
        {
            throw RingMatrixException("RingMatrix::fromString: invalid or out-of-range row count");
        }
        p = endptr;
        errno = 0;
        long long c_long = strtoll(p, &endptr, 10);
        if (endptr == p || errno == ERANGE)
        {
            throw RingMatrixException("RingMatrix::fromString: invalid or out-of-range column count");
        }
        p = endptr;

        if (m_long <= 0 || r_long <= 0 || c_long <= 0)
        {
            throw RingMatrixException("RingMatrix::fromString: non-positive modulus or dimensions (mod=" +
                std::to_string(m_long) + ", rows=" + std::to_string(r_long) +
                ", cols=" + std::to_string(c_long) + ")");
        }
        int64_t m = (int64_t)m_long;
        int64_t r = (int64_t)r_long;
        int64_t c = (int64_t)c_long;

        RingMatrix mat(r, c, m);
        p = endptr;
        while (*p && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')) ++p;

        size_t total = (size_t)r * (size_t)c;
        for (size_t i = 0; i < total; ++i)
        {
            errno = 0;
            long long v = strtoll(p, &endptr, 10);
            if (endptr == p || errno == ERANGE)
            {
                throw RingMatrixException("RingMatrix::fromString: invalid or out-of-range matrix element at position " +
                    std::to_string(i));
            }
            mat.data[i] = mat.normalize(v);
            p = endptr;
            while (*p && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')) ++p;
        }
        return mat;
    }

    int64_t* characteristicPolynomialDanilevsky() const
    {
        if (rows != cols)
        {
            throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: matrix must be square (got " +
                std::to_string(rows) + "x" + std::to_string(cols) + ")");
        }
        int64_t n = rows;
        if (n <= 0)
        {
            throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: invalid dimension (n=" +
                std::to_string(n) + ")");
        }

        RingMatrix B(*this);

        size_t nn = (size_t)n * (size_t)n;
        int64_t* matU = new int64_t[nn];
        int64_t* matV = new int64_t[nn];
        int64_t* temp = new int64_t[nn];

        auto cleanup = [&]() {
            delete[] matU;
            delete[] matV;
            delete[] temp;
            };

        try
        {
            for (int64_t i = n - 1; i >= 1; --i)
            {
                int64_t pivot = B.at(i, i - 1);
                if (pivot == 0)
                {
                    int64_t found = -1;
                    for (int64_t p = 0; p <= i - 2; ++p)
                    {
                        if (B.at(i, p) != 0) { found = p; break; }
                    }
                    if (found == -1)
                    {
                        throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: cannot find non-zero pivot in row " +
                            std::to_string(i));
                    }
                    B.swap_cols(found, i - 1);
                    B.swap_rows(found, i - 1);
                    pivot = B.at(i, i - 1);
                    if (pivot == 0)
                    {
                        throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: pivot is zero after row/column swap");
                    }
                }

                int64_t inv_pivot = modInv(pivot);

                for (size_t t = 0; t < nn; ++t) matV[t] = 0;
                for (int64_t r = 0; r < n; ++r)
                {
                    size_t idx = (size_t)r * (size_t)n + (size_t)r;
                    if (idx >= nn) {
                        cleanup();
                        throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: index out of bounds in matV initialization");
                    }
                    matV[idx] = 1 % mod;
                }

                size_t row_offset = (size_t)(i - 1) * (size_t)n;
                if (row_offset >= nn) {
                    cleanup();
                    throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: row_offset out of bounds for matV");
                }

                for (int64_t j = 0; j < n; ++j)
                {
                    size_t idx = row_offset + (size_t)j;
                    if (idx >= nn) {
                        cleanup();
                        throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: index out of bounds for matV row");
                    }

                    if (j == i - 1)
                        matV[idx] = normalize(inv_pivot);
                    else
                    {
                        long long val = -(long long)B.at(i, j) * (long long)inv_pivot;
                        matV[idx] = normalize(val);
                    }
                }

                for (size_t t = 0; t < nn; ++t) matU[t] = 0;
                for (int64_t r = 0; r < n; ++r)
                {
                    size_t idx = (size_t)r * (size_t)n + (size_t)r;
                    if (idx >= nn) {
                        cleanup();
                        throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: index out of bounds in matU initialization");
                    }
                    matU[idx] = 1 % mod;
                }

                row_offset = (size_t)(i - 1) * (size_t)n;
                if (row_offset >= nn) {
                    cleanup();
                    throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: row_offset out of bounds for matU");
                }

                for (int64_t j = 0; j < n; ++j)
                {
                    size_t idx = row_offset + (size_t)j;
                    if (idx >= nn) {
                        cleanup();
                        throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: index out of bounds for matU row");
                    }
                    matU[idx] = normalize(B.at(i, j));
                }

                mat_mul_mod(B.data, matV, temp, n, mod);
                mat_mul_mod(matU, temp, B.data, n, mod);

                if (B.at(i, i - 1) != 1 % mod)
                {
                    throw RingMatrixException("RingMatrix::characteristicPolynomialDanilevsky: invariant violated at step i=" +
                        std::to_string(i));
                }
            }

            int64_t* coeffs = new int64_t[(size_t)n + 1];
            coeffs[0] = normalize(1);
            for (int64_t k = 1; k <= n; ++k)
            {
                coeffs[(size_t)k] = normalize(-(long long)B.at(0, k - 1));
            }

            delete[] matU;
            delete[] matV;
            delete[] temp;

            return coeffs;
        }
        catch (...)
        {
            cleanup();
            throw;
        }
    }
};