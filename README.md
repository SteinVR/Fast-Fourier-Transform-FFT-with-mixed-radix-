# Fast Fourier Transform (FFT) with 2,3,5 mixed radix
The Discrete Fourier Transform (DFT) is the version of the Fourier Transform used with digital signals.

The DFT takes a sequence of values (our discrete signal) and expresses it in terms of its frequency components. The formula for the DFT is:

$$
X(k) = \sum_{n=0}^{N-1} x(n) e^{-i2\pi kn/N}
$$

where:
- \(X(k)\) is the kth frequency bin of the DFT,
- \(x(n)\) is the nth sample of the input signal,
- \(N\) is the total number of samples,
- \(i\) is the imaginary unit,
- \(e\) is the base of the natural logarithm.

The DFT is a complex-valued function: each frequency bin \(X(k)\) is a complex number that describes the amplitude and phase of the kth frequency component.

**Fast Fourier Transform (FFT):**

The DFT is very computationally intensive. If signal have \(N\) samples, the DFT requires \(O(N^2)\) operations.

The Fast Fourier Transform (FFT) is an algorithm that computes the DFT in \(O(N \log N)\) operations, which is much faster for large \(N\).
The key insight behind the FFT is the "divide and conquer" strategy. It breaks the DFT of size \(N\) down into smaller DFTs of sizes \(N/2\), then combines those results to get the final DFT.

The most common FFT algorithm is the Cooley-Tukey algorithm, which works when \(N\) is a power of 2.

How it works:
- **Step 1: Divide:** Split the input sequence into two sequences: one for the even-indexed elements and one for the odd-indexed elements. This gives you two sequences of length \(N/2\).

- **Step 2: Conquer:** Compute the DFT of the two \(N/2\) sequences recursively. This is where the "divide and conquer" strategy comes in. The base case for the recursion is a DFT of size 1, which is just the input value itself.

- **Step 3: Combine:** Combine the two \(N/2\) DFTs to form the final DFT of size \(N\). This involves multiplying the DFT of the odd-indexed sequence by a "twiddle factor" $e^{-i2\pi k/N}\$ and adding the result to the DFT of the even-indexed sequence.

However, the standard Cooley-Tukey FFT algorithm is specifically designed for sequences whose lengths are powers of 2. This is because the algorithm works by recursively dividing the sequence into halves until it reaches sequences of length 1. If the length of the sequence is not a power of 2, this division process doesn't work.

**Mixed-Radix FFT:**

The mixed-radix FFT is a generalization of the Cooley-Tukey FFT that can handle sequences of any length, not just powers of 2. The term "mixed-radix" refers to the fact that the algorithm can divide the sequence into parts of different sizes (or "radices").

For example, if the length of the sequence is a multiple of 2, 3, and 5, the mixed-radix FFT could divide the sequence into parts of size 2, 3, or 5, depending on what's most efficient. This makes the mixed-radix FFT much more flexible than the standard Cooley-Tukey FFT.

The mixed-radix FFT still uses the same basic "divide and conquer" strategy as the Cooley-Tukey FFT. It divides the sequence into smaller parts, computes the FFT of each part, and then combines the results to get the FFT of the whole sequence. But because it can divide the sequence in different ways, it can handle a wider range of sequence lengths.
