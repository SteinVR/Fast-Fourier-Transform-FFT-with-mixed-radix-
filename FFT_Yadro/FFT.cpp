#include "FFT.h"
#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


std::vector<std::complex<double>> FFT::forward(const std::vector<std::complex<double>>& x) {
    return fft_main(x, -1);
}

std::vector<std::complex<double>> FFT::inverse(const std::vector<std::complex<double>>& x) {
    std::vector<std::complex<double>> result = fft_main(x, 1);
    for (auto& value : result) {
        value /= result.size();
    }
    return result;
}
// FFT Main 
std::vector<std::complex<double>> FFT::fft_main(const std::vector<std::complex<double>>& x, int direction) {
    int N = x.size();
    if (N <= 1) {
        return x;
    }
    if (N % 5 == 0) {
        std::cout << "radix_5_fft" << std::endl;
        return radix_5_fft(x, direction);
    }
    else if (N % 3 == 0) {
        std::cout << "radix_3_fft" << std::endl;
        return radix_3_fft(x, direction);
    }
    else if (N % 2 == 0) {
        std::cout << "radix_2_fft" << std::endl;
        return radix_2_fft(x, direction);
    }
    else {
        std::cout << "classic_dft" << std::endl;
        return classic_dft(x, direction);
    }
}

std::vector<std::complex<double>> FFT::classic_dft(const std::vector<std::complex<double>>& x, int direction) {
    int N = x.size();
    std::vector<std::complex<double>> X(N);
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            double theta = direction * 2.0 * M_PI * k * n / N;
            X[k] += x[n] * std::complex<double>(cos(theta), sin(theta));
        }
    }
    return X;
}

std::vector<std::complex<double>> FFT::radix_2_fft(const std::vector<std::complex<double>>& x, int direction) {
    int N = x.size();
    if (N <= 1) {
        return x;
    }

    // Split input into 2 parts
    std::vector<std::complex<double>> even_elements(N / 2);
    std::vector<std::complex<double>> odd_elements(N / 2);
    for (int i = 0; i < N / 2; i++) {
        even_elements[i] = x[2 * i];
        odd_elements[i] = x[2 * i + 1];
    }

    //FFT for the 2 parts
    std::vector<std::complex<double>> fft_group_0 = fft_main(even_elements, direction);
    std::vector<std::complex<double>> fft_group_1 = fft_main(odd_elements, direction);

    // Combine FFT for the 2 parts
    std::vector<std::complex<double>> combined(N);
    for (int k = 0; k < N / 2; k++) {
    //    std::complex<double> t = std::polar(1.0, direction * 2.0 * M_PI * k / N) * fft_group_1[k];
        double theta = direction * 2.0 * M_PI * k / N;
        std::complex<double> t = std::complex<double>(cos(theta), sin(theta)) * fft_group_1[k];
        combined[k] = fft_group_0[k] + t;
        combined[k + N / 2] = fft_group_0[k] - t;
    }

    return combined;    
}

std::vector<std::complex<double>> FFT::radix_3_fft(const std::vector<std::complex<double>>& x, int direction) {
    int N = x.size();
    if (N <= 1) {
        return x;
    }

    // Split input into 3 parts
    std::vector<std::complex<double>> group_0(N / 3);
    std::vector<std::complex<double>> group_1(N / 3);
    std::vector<std::complex<double>> group_2(N / 3);
    for (int i = 0; i < N / 3; i++) {
        group_0[i] = x[3 * i];
        group_1[i] = x[3 * i + 1];
        group_2[i] = x[3 * i + 2];
    }

    //FFT for the 3 parts
    std::vector<std::complex<double>> fft_group_0 = fft_main(group_0, direction);
    std::vector<std::complex<double>> fft_group_1 = fft_main(group_1, direction);
    std::vector<std::complex<double>> fft_group_2 = fft_main(group_2, direction);

    // Combine FFT for the 3 parts
    std::vector<std::complex<double>> combined(N);
    for (int k = 0; k < N / 3; k++) {
        std::complex<double> t1 = std::polar(1.0, direction * 2.0 * M_PI * k / N) * fft_group_1[k];
        std::complex<double> t2 = std::polar(1.0, direction * 2.0 * M_PI * 2 * k / N) * fft_group_2[k];
        combined[k] = fft_group_0[k] + t1 + t2;
        combined[k + N / 3] = fft_group_0[k] + std::polar(1.0, direction * 2.0 * M_PI * (k + N / 3) / N) * fft_group_1[k] + std::polar(1.0, direction * 2.0 * M_PI * 2 * (k + N / 3) / N) * fft_group_2[k];
        combined[k + 2 * N / 3] = fft_group_0[k] + std::polar(1.0, direction * 2.0 * M_PI * (k + 2 * N / 3) / N) * fft_group_1[k] + std::polar(1.0, direction * 2.0 * M_PI * 2 * (k + 2 * N / 3) / N) * fft_group_2[k];
    }

    return combined;
}

std::vector<std::complex<double>> FFT::radix_5_fft(const std::vector<std::complex<double>>& x, int direction) {
    int N = x.size();
    if (N <= 1) {
        return x;
    }

    // Split input into 5 parts
    std::vector<std::complex<double>> group_0(N / 5);
    std::vector<std::complex<double>> group_1(N / 5);
    std::vector<std::complex<double>> group_2(N / 5);
    std::vector<std::complex<double>> group_3(N / 5);
    std::vector<std::complex<double>> group_4(N / 5);
    for (int i = 0; i < N / 5; i++) {
        group_0[i] = x[5 * i];
        group_1[i] = x[5 * i + 1];
        group_2[i] = x[5 * i + 2];
        group_3[i] = x[5 * i + 3];
        group_4[i] = x[5 * i + 4];
    }

    //FFT for the 5 parts
    std::vector<std::complex<double>> fft_group_0 = fft_main(group_0, direction);
    std::vector<std::complex<double>> fft_group_1 = fft_main(group_1, direction);
    std::vector<std::complex<double>> fft_group_2 = fft_main(group_2, direction);
    std::vector<std::complex<double>> fft_group_3 = fft_main(group_3, direction);
    std::vector<std::complex<double>> fft_group_4 = fft_main(group_4, direction);

    // Combine FFT for the 5 parts
    std::vector<std::complex<double>> combined(N);
    for (int k = 0; k < N / 5; k++) {
        std::complex<double> t1 = std::polar(1.0, direction * 2.0 * M_PI * k / N) * fft_group_1[k];
        std::complex<double> t2 = std::polar(1.0, direction * 2.0 * M_PI * 2 * k / N) * fft_group_2[k];
        std::complex<double> t3 = std::polar(1.0, direction * 2.0 * M_PI * 3 * k / N) * fft_group_3[k];
        std::complex<double> t4 = std::polar(1.0, direction * 2.0 * M_PI * 4 * k / N) * fft_group_4[k];

        combined[k] = fft_group_0[k] + t1 + t2 + t3 + t4;
        combined[k + N / 5] = fft_group_0[k] + std::polar(1.0, direction * 2.0 * M_PI * (k + N / 5) / N) * fft_group_1[k] + std::polar(1.0, direction * 2.0 * M_PI * 2 * (k + N / 5) / N) * fft_group_2[k] + std::polar(1.0, direction * 2.0 * M_PI * 3 * (k + N / 5) / N) * fft_group_3[k] + std::polar(1.0, direction * 2.0 * M_PI * 4 * (k + N / 5) / N) * fft_group_4[k];
        combined[k + 2 * N / 5] = fft_group_0[k] + std::polar(1.0, direction * 2.0 * M_PI * (k + 2 * N / 5) / N) * fft_group_1[k] + std::polar(1.0, direction * 2.0 * M_PI * 2 * (k + 2 * N / 5) / N) * fft_group_2[k] + std::polar(1.0, direction * 2.0 * M_PI * 3 * (k + 2 * N / 5) / N) * fft_group_3[k] + std::polar(1.0, direction * 2.0 * M_PI * 4 * (k + 2 * N / 5) / N) * fft_group_4[k];
        combined[k + 3 * N / 5] = fft_group_0[k] + std::polar(1.0, direction * 2.0 * M_PI * (k + 3 * N / 5) / N) * fft_group_1[k] + std::polar(1.0, direction * 2.0 * M_PI * 2 * (k + 3 * N / 5) / N) * fft_group_2[k] + std::polar(1.0, direction * 2.0 * M_PI * 3 * (k + 3 * N / 5) / N) * fft_group_3[k] + std::polar(1.0, direction * 2.0 * M_PI * 4 * (k + 3 * N / 5) / N) * fft_group_4[k];
        combined[k + 4 * N / 5] = fft_group_0[k] + std::polar(1.0, direction * 2.0 * M_PI * (k + 4 * N / 5) / N) * fft_group_1[k] + std::polar(1.0, direction * 2.0 * M_PI * 2 * (k + 4 * N / 5) / N) * fft_group_2[k] + std::polar(1.0, direction * 2.0 * M_PI * 3 * (k + 4 * N / 5) / N) * fft_group_3[k] + std::polar(1.0, direction * 2.0 * M_PI * 4 * (k + 4 * N / 5) / N) * fft_group_4[k];
    }

    return combined;
}