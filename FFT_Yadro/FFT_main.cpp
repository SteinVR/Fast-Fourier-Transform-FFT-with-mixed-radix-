#include <iostream>
#include <algorithm>
#include "FFT.h"

int main() {
    // Input
    //std::vector<std::complex<double>> x = { 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<std::complex<double>> x = FFT::random_complex_numbers(20);

    // Forward FFT
    std::vector<std::complex<double>> fft_x = FFT::forward(x);

    // Inverse FFT
    std::vector<std::complex<double>> ifft_fft_x = FFT::inverse(fft_x);

    // Output
    std::cout << "\n" << std::endl;
    std::cout << "Forward: [";
    for (int i = 0; i < fft_x.size(); i++) {
        std::cout << fft_x[i];
        if (i != fft_x.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Invert: [";
    for (int i = 0; i < ifft_fft_x.size(); i++) {
        std::cout << ifft_fft_x[i];
        if (i != ifft_fft_x.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    std::cout << "\n" << std::endl;


    std::cout << "Origin: [";
    for (int i = 0; i < x.size(); i++) {
        std::cout << x[i];
        if (i != x.size() - 1) {
            std::cout << ", ";
        }
    }
    // Direct Compare
    /*
    double epsilon = 1e-4;
    bool equal = false;

    for (int i = 0; i < x.size(); i++) {
        double real_diff = std::abs(x[i].real() - ifft_fft_x[i].real());
        double imag_diff = std::abs(x[i].imag() - ifft_fft_x[i].imag());
        if (real_diff < epsilon || imag_diff < epsilon) {
            equal = true;
            break;
        }
    }

    if (equal) {
        std::cout << "\n" << std::endl;
        std::cout << "True\n" << std::endl;
    }
    else {
        std::cout << "\n" << std::endl;
        std::cout << "False\n" << std::endl;
    }
*/
    // RMSE
    double sum_squared_error = 0.0;
    for (int i = 0; i < x.size(); i++) {
        double error_real = x[i].real() - ifft_fft_x[i].real();
        double error_imag = x[i].imag() - ifft_fft_x[i].imag();
        sum_squared_error += error_real * error_real + error_imag * error_imag;
    }
    double mean_squared_error = sum_squared_error / x.size();
    double root_mean_squared_error = std::sqrt(mean_squared_error);

    std::cout << "\n" << std::endl;
    std::cout << "Root Mean Squared Error: " << root_mean_squared_error << std::endl;


    return 0;
}
