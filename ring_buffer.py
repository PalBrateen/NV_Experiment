# ring_buffer.py
"""
Fixed-size ring buffer backed by a numpy array.

O(1) append, fixed memory footprint, fast chronological readout.
Designed for real-time DAQ streaming where display reads the last
N seconds while acquisition runs indefinitely.
"""

import numpy as np


class RingBuffer:
    """Fixed-size ring buffer backed by a contiguous numpy array.

    Parameters
    ----------
    capacity : int
        Maximum number of samples stored.
    dtype : numpy dtype
        Data type (default float64).

    Properties
    ----------
    total_count : int
        Cumulative number of samples written (monotonically increasing).
    is_full : bool
        True once the buffer has wrapped at least once.
    """

    def __init__(self, capacity: int, dtype=np.float64):
        self._buf = np.zeros(capacity, dtype=dtype)
        self._cap = capacity
        self._write = 0       # next write index (wraps)
        self._count = 0       # total samples ever written

    # ── write ───────────────────────────────────────────────────────

    def append(self, data: np.ndarray):
        """Append a chunk of data, overwriting oldest if necessary."""
        n = len(data)
        if n == 0:
            return
        if n >= self._cap:
            # Chunk larger than buffer — keep only the tail
            self._buf[:] = data[-self._cap:]
            self._write = 0
            self._count += n
            return

        end = self._write + n
        if end <= self._cap:
            self._buf[self._write:end] = data
        else:
            first = self._cap - self._write
            self._buf[self._write:] = data[:first]
            self._buf[:n - first] = data[first:]
        self._write = end % self._cap
        self._count += n

    # ── read ────────────────────────────────────────────────────────

    def get_ordered(self) -> np.ndarray:
        """Return data in chronological order (oldest → newest).

        Returns a copy — safe to hold across append calls.
        Before the buffer is full, returns only the valid samples.
        """
        if self._count == 0:
            return np.empty(0, dtype=self._buf.dtype)
        if self._count < self._cap:
            return self._buf[:self._write].copy()
        return np.roll(self._buf, -self._write).copy()

    @property
    def total_count(self) -> int:
        """Cumulative samples written since creation."""
        return self._count

    @property
    def is_full(self) -> bool:
        return self._count >= self._cap

    @property
    def capacity(self) -> int:
        return self._cap

    def __len__(self) -> int:
        """Number of valid samples currently stored."""
        return min(self._count, self._cap)

    def reset(self):
        """Clear buffer without reallocating."""
        self._buf[:] = 0
        self._write = 0
        self._count = 0
