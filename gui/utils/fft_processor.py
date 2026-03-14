"""
FFT Processor for NV Experiment GUI

Compute Fast Fourier Transform of timeseries data for live plotting.
CPU-based implementation with GPU acceleration placeholder for later.

Features:
- Multiple window functions (Hanning, Hamming, Blackman, etc.)
- Peak detection in frequency domain
- Normalization options
- Future GPU support via CuPy
"""

import numpy as np
from scipy.fft import fft, fftfreq, rfft, rfftfreq
from scipy.signal import find_peaks
from typing import Tuple, Optional, List
from enum import Enum


class WindowFunction(Enum):
    """Available window functions for FFT."""
    RECTANGULAR = "rectangular"
    HANNING = "hanning"
    HAMMING = "hamming"
    BLACKMAN = "blackman"
    BARTLETT = "bartlett"
    KAISER = "kaiser"


class FFTProcessor:
    """
    FFT computation with windowing and peak detection.

    Supports both CPU (NumPy/SciPy) and GPU (CuPy) backends.
    """

    def __init__(self, use_gpu: bool = False):
        """
        Initialize FFT processor.

        Args:
            use_gpu: If True, attempt to use GPU (CuPy). Falls back to CPU if unavailable.
        """
        self.use_gpu = use_gpu
        self._gpu_available = False

        if use_gpu:
            try:
                import cupy as cp
                self._gpu_available = True
                self.cp = cp
                print("✔ GPU acceleration available (CuPy)")
            except ImportError:
                print("⚠️ CuPy not available, using CPU for FFT")
                self.use_gpu = False

    def compute_fft(
        self,
        data: np.ndarray,
        sampling_rate: float,
        window: WindowFunction = WindowFunction.HANNING,
        real_only: bool = True,
        normalize: bool = True
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute FFT of timeseries data.

        Args:
            data: Input timeseries data (1D array)
            sampling_rate: Sampling rate in Hz
            window: Window function to apply
            real_only: If True, return only positive frequencies (for real signals)
            normalize: If True, normalize magnitude by length

        Returns:
            (frequencies, magnitudes) - Both 1D arrays
        """
        if len(data) == 0:
            return np.array([]), np.array([])

        # Apply window function
        windowed_data = self._apply_window(data, window)

        # Compute FFT
        if real_only:
            # Real FFT (more efficient for real signals)
            fft_result = rfft(windowed_data)
            freqs = rfftfreq(len(data), 1.0 / sampling_rate)
        else:
            # Complex FFT
            fft_result = fft(windowed_data)
            freqs = fftfreq(len(data), 1.0 / sampling_rate)
            # Take positive frequencies only
            mask = freqs >= 0
            freqs = freqs[mask]
            fft_result = fft_result[mask]

        # Compute magnitude
        magnitude = np.abs(fft_result)

        # Normalize
        if normalize:
            magnitude = magnitude / len(data)

        return freqs, magnitude

    def compute_power_spectrum(
        self,
        data: np.ndarray,
        sampling_rate: float,
        window: WindowFunction = WindowFunction.HANNING,
        db_scale: bool = False
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute power spectrum (magnitude squared).

        Args:
            data: Input timeseries data
            sampling_rate: Sampling rate in Hz
            window: Window function
            db_scale: If True, return in dB scale (10*log10)

        Returns:
            (frequencies, power) - Both 1D arrays
        """
        freqs, magnitude = self.compute_fft(data, sampling_rate, window, real_only=True)

        # Power = magnitude squared
        power = magnitude ** 2

        # Convert to dB if requested
        if db_scale:
            # Avoid log(0) by adding small epsilon
            power = 10 * np.log10(power + 1e-12)

        return freqs, power

    def find_peaks(
        self,
        frequencies: np.ndarray,
        magnitudes: np.ndarray,
        num_peaks: int = 5,
        min_prominence: Optional[float] = None
    ) -> List[Tuple[float, float]]:
        """
        Find peaks in frequency spectrum.

        Args:
            frequencies: Frequency array
            magnitudes: Magnitude array
            num_peaks: Maximum number of peaks to return
            min_prominence: Minimum peak prominence (relative to neighbors)

        Returns:
            List of (frequency, magnitude) tuples, sorted by magnitude descending
        """
        if len(magnitudes) == 0:
            return []

        # Detect peaks
        if min_prominence is None:
            min_prominence = np.max(magnitudes) * 0.1  # 10% of max

        peak_indices, properties = find_peaks(
            magnitudes,
            prominence=min_prominence
        )

        if len(peak_indices) == 0:
            return []

        # Extract peak frequencies and magnitudes
        peak_freqs = frequencies[peak_indices]
        peak_mags = magnitudes[peak_indices]

        # Sort by magnitude descending
        sorted_indices = np.argsort(peak_mags)[::-1]

        # Return top N peaks
        peaks = [
            (peak_freqs[i], peak_mags[i])
            for i in sorted_indices[:num_peaks]
        ]

        return peaks

    def _apply_window(self, data: np.ndarray, window: WindowFunction) -> np.ndarray:
        """
        Apply window function to data.

        Args:
            data: Input data
            window: Window function type

        Returns:
            Windowed data
        """
        n = len(data)

        if window == WindowFunction.RECTANGULAR:
            return data

        elif window == WindowFunction.HANNING:
            window_array = np.hanning(n)

        elif window == WindowFunction.HAMMING:
            window_array = np.hamming(n)

        elif window == WindowFunction.BLACKMAN:
            window_array = np.blackman(n)

        elif window == WindowFunction.BARTLETT:
            window_array = np.bartlett(n)

        elif window == WindowFunction.KAISER:
            # Kaiser window with beta=8.6 (good balance)
            window_array = np.kaiser(n, beta=8.6)

        else:
            # Default to Hanning
            window_array = np.hanning(n)

        return data * window_array

    def get_frequency_resolution(self, num_samples: int, sampling_rate: float) -> float:
        """
        Calculate frequency resolution.

        Args:
            num_samples: Number of samples in timeseries
            sampling_rate: Sampling rate in Hz

        Returns:
            Frequency resolution in Hz
        """
        return sampling_rate / num_samples

    def get_nyquist_frequency(self, sampling_rate: float) -> float:
        """
        Calculate Nyquist frequency (maximum observable frequency).

        Args:
            sampling_rate: Sampling rate in Hz

        Returns:
            Nyquist frequency in Hz
        """
        return sampling_rate / 2.0


# Convenience functions for quick usage

def compute_fft_simple(
    data: np.ndarray,
    sampling_rate: float,
    window: str = "hanning"
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simple FFT computation with default settings.

    Args:
        data: Timeseries data
        sampling_rate: Sampling rate in Hz
        window: Window function name

    Returns:
        (frequencies, magnitudes)
    """
    processor = FFTProcessor(use_gpu=False)
    window_enum = WindowFunction[window.upper()]
    return processor.compute_fft(data, sampling_rate, window=window_enum)


def find_dominant_frequency(
    data: np.ndarray,
    sampling_rate: float,
    window: str = "hanning"
) -> Tuple[float, float]:
    """
    Find the dominant frequency in the signal.

    Args:
        data: Timeseries data
        sampling_rate: Sampling rate in Hz
        window: Window function name

    Returns:
        (dominant_frequency, magnitude)
    """
    processor = FFTProcessor(use_gpu=False)
    window_enum = WindowFunction[window.upper()]

    freqs, mags = processor.compute_fft(data, sampling_rate, window=window_enum)

    if len(mags) == 0:
        return 0.0, 0.0

    # Find peak (excluding DC component at index 0)
    if len(mags) > 1:
        peak_idx = np.argmax(mags[1:]) + 1
    else:
        peak_idx = 0

    return freqs[peak_idx], mags[peak_idx]


# Example usage
if __name__ == "__main__":
    # Generate test signal: 10 Hz + 50 Hz sine waves
    sampling_rate = 1000  # Hz
    duration = 1.0  # seconds
    t = np.linspace(0, duration, int(sampling_rate * duration))

    signal = np.sin(2 * np.pi * 10 * t) + 0.5 * np.sin(2 * np.pi * 50 * t)
    noise = 0.1 * np.random.randn(len(t))
    noisy_signal = signal + noise

    # Compute FFT
    processor = FFTProcessor()
    freqs, mags = processor.compute_fft(
        noisy_signal,
        sampling_rate,
        window=WindowFunction.HANNING
    )

    # Find peaks
    peaks = processor.find_peaks(freqs, mags, num_peaks=5)

    print("Detected peaks:")
    for freq, mag in peaks:
        print(f"  {freq:.2f} Hz: magnitude = {mag:.4f}")

    # Expected: peaks at ~10 Hz and ~50 Hz
