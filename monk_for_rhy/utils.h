//! \file utils.h
//! Utilities
#ifndef _UTILS_H
#define _UTILS_H
#include <iostream>
#include <iterator>
#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

//! Print all elements in a std::vector
template <typename T>
std::ostream & operator << (std::ostream & out, const std::vector<T> & vec) {
  out << "[";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(out, ", "));
  out << "]";
  return out;
};

//! Print all elements in a std::array
template <typename T, size_t n>
std::ostream & operator << (std::ostream & out, const std::array<T, n> & vec) {
  out << "[";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(out, ", "));
  out << "]";
  return out;
}

//! Multiply std::arrays by a scalar
template <typename T, size_t n>
void mul(const double x, const std::array<T, n> & arr1, std::array<T, n> & arr2) {
    std::transform(arr1.cbegin(), arr1.cend(), arr2.begin(), std::bind2nd(std::multiplies<double>(), x));
}

//! Multiply one std::array by another
template <typename T, size_t n>
void mul(const std::array<T, n> & arr1, const std::array<T, n> & arr2, std::array<T, n> & arr3) {
    std::transform(arr1.cbegin(), arr1.cend(), arr2.cbegin(), arr3.begin(), std::multiplies<double>());
}

//! Array addition
template <typename T, size_t n>
void add(const std::array<T, n> & arr1, const std::array<T, n> & arr2, std::array<T, n> & arr3) {
    std::transform(arr1.cbegin(), arr1.cend(), arr2.cbegin(), arr3.begin(), std::plus<double>());
}

//! Array dot product
template <typename T, size_t n>
double arrdotprod(const std::array<T, n> & arr1, const std::array<T, n> & arr2) {
    double prod = 0.;
    for (size_t i = 0; i < n; ++i)
        prod += arr1[i] * arr2[i];
    return prod;
}

//! Dot product with arrays passed by pointer to fit all kinds of container
template <typename T>
double arrdotprod(size_t n, const T * arr1, const T * arr2) {
    double prod = 0.;
    for (size_t i = 0; i < n; ++i)
        prod += *(arr1 + i) * *(arr2 + i);
    return prod;
}

//! Cross product. Here we limit ourself in cross product in 3 dimension.
template <typename T>
void arrcrossprod(const std::array<T, 3> & u, const std::array<T, 3> & v, std::array<T, 3> & prod) {
    prod[0] = u[1] * v[2] - u[2] * v[1];
    prod[1] = u[2] * v[0] - u[0] * v[2];
    prod[2] = u[0] * v[1] - u[1] * v[0];
}

//! Overloaded `arrcrossprod()` with vectors being passed by pointer
template <typename T>
void arrcrossprod(const T * pu, const T * pv, T * pp) {
    *(pp) = *(pu + 1) * *(pv + 2) - *(pu + 2) * *(pv + 1);
    *(pp + 1) = *(pu + 2) * *(pv) - *(pu) * *(pv + 2);
    *(pp + 2) = *(pu) * *(pv + 1) - *(pu + 1) * *(pv);
}
    
//! Write std::vector to binary file
template <typename T>
void wdoublevec(std::vector<T> & vec, const std::string & path) {
    std::ofstream ofile(path, std::ios::out | std::ofstream::binary);
    vec.push_back(0.);
    std::copy(reinterpret_cast<const char*>(&vec.front()), reinterpret_cast<const char*>(&vec.back()),
    	std::ostreambuf_iterator<char>(ofile));
    vec.pop_back();
    ofile.flush();
    ofile.close();
}

//! Read a whole binary file and save data in a std::vector
template <typename T>
void rdoublevec(std::vector<T> & vec, const std::string & path) {
    std::ifstream ifile(path, std::ios::in | std::ios::binary | std::ios::ate);
    
    ifile.seekg(0, std::ios::end);
    auto fileSize = ifile.tellg();
    ifile.seekg(0, std::ios::beg);
    auto buffer = new char[fileSize];
    ifile.read(buffer, fileSize);
    ifile.close();
    
    T * double_values = (T*)buffer;
    vec =  std::vector<T>(double_values, double_values + (fileSize / sizeof(T)));
    delete[] buffer;
}

/*!
\brief Read a segment of a binary file and save data in a std::vector. This is necessary while reading a large file
@param vec the vector to contain the data
@param path path to binary file
@param rstart the offset, in number of elements (not in bytes)
@param ndata number of elements to be read
*/
template <typename T>
void rdoublevec(std::vector<T> & vec, const std::string & path, size_t rstart, size_t ndata) {
    size_t ndata_char = (sizeof(T)) / (sizeof(char)) * ndata;
    size_t pos_start = (sizeof(T)) / (sizeof(char)) * rstart;
    size_t pos_end = pos_start + ndata_char;
    std::ifstream ifile(path, std::ios::in | std::ios::binary | std::ios::ate);
    
    auto fileSize = ifile.tellg();
    if (pos_start > fileSize) {
        std::cerr << "Starting point out of range!" << std::endl;
        return;
    }
    
    if (pos_end > fileSize) {
        ndata_char = (size_t)(fileSize) - pos_start;
    }
    
    ifile.seekg(pos_start, std::ios::beg);
    auto buffer = new char[ndata_char];
    ifile.read(buffer, ndata_char);
    ifile.close();
    T * double_values = (T*)buffer;
    vec =  std::vector<T>(double_values, double_values + (ndata_char / sizeof(T)));
    delete[] buffer;
}

/*!
\brief Print progress on screen. Use this function **ONLY** when working in an interactive shell (other than submitting a job).
@param progress fraction of work having been done.
*/
inline void showprogress(double progress) {
    int barWidth = 70, pos;
    std::cout << "[";
    pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

/*!
\brief Find out file size in byte, with `std::string` argument.
*/
inline std::ifstream::pos_type get_filesize(const std::string filename) {
    std::ifstream ifile(filename, std::ios::in | std::ios::binary | std::ios::ate);
    auto fileSize = ifile.tellg();
    ifile.close();
    return fileSize;
}

/*!
\brief Find out file size in byte, with `const char *` argument.
*/
inline std::ifstream::pos_type get_filesize(const char * filename) {
    std::ifstream ifile(filename, std::ios::in | std::ios::binary | std::ios::ate);
    auto fileSize = ifile.tellg();
    ifile.close();
    return fileSize;
}

/*!
\brief Polynomial evaluation.
*/
template <size_t n>
double polyval(const std::array<double, n> & coeffs, const double x) {
    double res = coeffs[0];
    for (size_t i = 1; i < n; ++i) {
        res = (res * x + coeffs[i]);
    }
    return res;
}

/*!
\brief Polynomial evaluation.
*/
template <size_t n>
double polyval1(const std::array<double, n> & coeffs, const double x) {
    double res = coeffs[n - 1];
    for (int j = n - 2; j >= 0; --j) {
        res = (res * x + coeffs[j]);
    }
    return res;
}

/*!
\brief Return sign of `x` in double
*/
inline double dsign(const double x) {
    double y;
    y = (x >= 0.) ? 1 : -1;
    if (x == 0.) y = 0.;
    return y;
}

inline std::vector<std::string> str_split(const std::string & s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}
#endif
