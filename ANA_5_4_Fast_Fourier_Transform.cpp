// NUMERICAL INTEGRATION
// NEWTON-COTES FORMULAE
// SIMPSON'S RULE

// Date: 28/11/2020
// Author: Fernando Fernández del Cerro

// Cooley-Tukey formulation of the Fast Fourier Transform

#include <complex>
#include <cmath>
#include <vector>

// Define a shorthand for complex numbers
typedef std::complex<double> Complex;

// Define a shorthand for a vector of complex numbers
typedef std::vector<Complex> Vector;

// Recursive function to compute the FFT using the Cooley-Tukey algorithm
void fft(Vector& v, bool inverse = false) {
    int n = v.size();

    // Base case: if there is only one element in the input, it is already a FFT
    if (n == 1) return;

    // Divide and conquer: split the input into even and odd indices
    // and compute the FFT of each half recursively
    Vector even(n/2), odd(n/2);
    for (int i = 0; i < n; i += 2) {
        even[i/2] = v[i];
        odd[i/2] = v[i+1];
    }
    fft(even, inverse);
    fft(odd, inverse);

    // Combine the results of the subproblems
    double sign = (inverse ? 1 : -1); // sign of the exponent for the inverse FFT
    Complex w_n = exp(Complex(0, sign * 2 * M_PI / n)); // nth root of unity
    Complex w(1, 0); // start with w^0 = 1
    for (int i = 0; i < n/2; i++) {
        v[i] = even[i] + w * odd[i]; // combine even and odd indices
        v[i + n/2] = even[i] - w * odd[i]; // combine even and odd indices
        w *= w_n; // multiply w by the nth root of unity
    }

    // Divide the result by n for the inverse FFT
    if (inverse) {
        for (int i = 0; i < n; i++) {
            v[i] /= n;
        }
    }
}

// Function to compute the FFT of a vector of complex numbers
Vector fast_fourier_transform(const Vector& a) {
    Vector v = a;
    fft(v);
    return v;
}

// Function to compute the inverse FFT of a vector of complex numbers
Vector inverse_fast_fourier_transform(const Vector& v) {
    Vector a = v;
    fft(a, true);
    return a;
}

