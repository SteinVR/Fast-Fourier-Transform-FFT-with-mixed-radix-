#pragma once
#include <complex>
#include <vector>
#include <cstdlib>
#include <ctime>

class FFT {
public:
    // Forward FFT
    static std::vector<std::complex<double>> forward(const std::vector<std::complex<double>>& x);

    // Inverse FFT
    static std::vector<std::complex<double>> inverse(const std::vector<std::complex<double>>& x);

    // Random complex numbers for input
    static std::vector<std::complex<double>> random_complex_numbers(int n) {
        std::vector<std::complex<double>> result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = std::complex<double>(std::rand() / (double)RAND_MAX, std::rand() / (double)RAND_MAX);
        }
        return result;
    }
private:
    // Classic DFT
    static std::vector<std::complex<double>> classic_dft(const std::vector<std::complex<double>>& x, int direction);

    // Radix-2 FFT
    static std::vector<std::complex<double>> radix_2_fft(const std::vector<std::complex<double>>& x, int direction);

    // Radix-3 FFT
    static std::vector<std::complex<double>> radix_3_fft(const std::vector<std::complex<double>>& x, int direction);

    // Radix-5 FFT
    static std::vector<std::complex<double>> radix_5_fft(const std::vector<std::complex<double>>& x, int direction);

    // Main FFT function
    static std::vector<std::complex<double>> fft_main(const std::vector<std::complex<double>>& x, int direction = -1);
};

