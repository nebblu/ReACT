#ifndef ARRAY_H
#define ARRAY_H

#include <vector>

#include "Common.h"

/* A simple multi-dimensional (up to 3) array, based on Numpy's array type.
 * The multi-dimensional array is represented in memory as a flat array
 * inherited from the STL vector class.  For 1-D arrays it may be used as
 * a plain old std::vector<real>.  For 2-D and 3-D arrays, the storage order
 * is row major. */
class array : public std::vector<real> {
public:
    /* Construct an empty (0-dimensional) array. */
    array();

    /* Construct an array with arbitrary initial values. */
    explicit array(size_type n0);
    array(size_type n0, size_type n1);
    array(size_type n0, size_type n1, size_type n2);

    /* Construct a 1-D array of length n, with all elements set to r. */
//    array(size_type n, real r);

    /* Construct an array with elements initialized from the C-array v.  For
     * 2-D arrays, the data elements are presumed to be laid out in memory as
     *   a_{ij} = v[n1*i + j] ,
     * i.e. in row-major order.  For 3-D arrays the memory layout is
     *   a_{ijk} = v[n2*(n1*i + j) + k] .  */
    array(size_type n0, const real* v);
    array(size_type n0, size_type n1, const real* v);
    array(size_type n0, size_type n1, size_type n2, const real* v);

    /* Construct a 1-D array from an STL vector */
    array(const std::vector<real>& v);

    /* Copy constructor. */
    array(const array& v);

    /* Assignment.  This array is resized, and the elements of v are copied. */
    array& operator=(const array& v);

    /* Destructor. */
    ~array();

    /* Construct an array with values initialized to zero. */
    static array zeros(size_type n0);
    static array zeros(size_type n0, size_type n1);
    static array zeros(size_type n0, size_type n1, size_type n2);

    /* Construct an array with values initialized to one. */
    static array ones(size_type n0);
    static array ones(size_type n0, size_type n1);
    static array ones(size_type n0, size_type n1, size_type n2);

    /* Constuct a 1-D array of length n with linearly or logarithmically spaced
     * values between min and max.  Note that if min is greater than max, the
     * values will be arranged in descending order.  Also note that for
     * logspace(), both min and max must be strictly positive. */
    static array linspace(real min, real max, size_type n = 50);
    static array logspace(real min, real max, size_type n = 50);

    /* Resize this array. */
    void resize(size_type n0);
    void resize(size_type n0, size_type n1);
    void resize(size_type n0, size_type n1, size_type n2);

    /* Transpose the current array over two axes. */
//    void transpose(int axis0 = 0, int axis1 = 1);

    void getshape(size_type* n0, size_type* n1, size_type* n2) const;

    real min() const;
    real max() const;
    array abs() const;

    /* Standard flat accessors. */
    const real& operator[](int i) const { return std::vector<real>::operator[](i); }
    real& operator[](int i) { return std::vector<real>::operator[](i); }

    /* Accessors for 1-D arrays */
    const real& operator()(int i0) const { return (*this)[i0]; }
    real& operator()(int i0) { return (*this)[i0]; }

    /* Accessors for 2-D arrays */
    const real& operator()(int i0, int i1) const { return (*this)[i0*n[1] + i1]; }
    real& operator()(int i0, int i1) { return (*this)[i0*n[1] + i1]; }

    /* Accessors for 3-D arrays */
    const real& operator()(int i0, int i1, int i2) const { return (*this)[(i0*n[1] + i1)*n[2] + i2]; }
    real& operator()(int i0, int i1, int i2) { return (*this)[(i0*n[1] + i1)*n[2] + i2]; }

    /* Automatic cast to C-array.  Note that the returned pointer is only valid
     * as long as this array is not resized. */
    operator real*() { return data(); }
    operator const real*() const { return data(); }

    /* Arithmetic operations.  Note that these methods do not check the shape
     * of the array v, only that its total size equals the total size of this
     * array. */
    array& operator+=(const array& v);  // element-wise addition
    array& operator-=(const array& v);  // element-wise subtraction
    array& operator*=(const array& v);  // element-wise multiplication
    array& operator/=(const array& v);  // element-wise division
    array& operator+=(real s);        // add a constant from each element
    array& operator-=(real s);        // subtract a constant from each element
    array& operator*=(real s);        // multiplication by a constant
    array& operator/=(real s);        // division by a constant

protected:
    size_type n[3];     // shape of array

    enum Spacing {
        LinearSpacing,
        LogarithmicSpacing
    };
    array(size_type n0, real min, real max, Spacing mode);

    array(size_type n0, size_type n1, size_type n2, real r);
};

array operator-(const array& v);    // unary negation
array operator+(const array& u, const array& v);
array operator-(const array& u, const array& v);
array operator*(const array& u, const array& v);
array operator/(const array& u, const array& v);
array operator+(const array& u, real s);
array operator+(real s, const array& u);
array operator-(const array& u, real s);
array operator*(const array& u, real s);
array operator*(real s, const array& u);
array operator/(const array& u, real s);

real min(const array& v);
real max(const array& v);
array abs(const array& v);

#endif // ARRAY_H
