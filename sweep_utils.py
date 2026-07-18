# sweep_utils.py
"""
Multi-parameter sweep utility with function-pointer optimization.

KEY CONCEPT — Two sweep tiers:

  INNER (first added): Changes fastest.  Same-shape acquisition.
      sweep.add('frequency', freq_values, setter=sg.set_freq)

  OUTER (subsequently added): Changes slowest.  Each outer value is a
      separate "experiment block" that gets Nruns of averaging before
      moving to the next outer value.
      sweep.add('mw_power', power_values, setter=sg.set_amp_rf)

The Sweep does NOT own the Nruns loop — the experiment class does.
The correct acquisition order is:

    for outer_combo in sweep.outer_combos():
        outer setters called automatically
        for i_run in range(Nruns):
            for inner_pt in sweep.inner_iter():
                inner setter called automatically
                acquire(...)

This ensures averaging happens at a fixed outer configuration.

CONFIGURATION parameters (daq_Nsamples, t_AOM, sampling_rate) that change
acquisition shape do NOT go into Sweep.  They are handled as an explicit
outermost loop that tears down / rebuilds the AI task.  See
DiodeExperiment._config_sweep() in experiment_base.py.
"""

import numpy as np
from itertools import product
from typing import Callable, Optional
from dataclasses import dataclass


@dataclass
class SweepIndex:
    """Rich index yielded per sweep point."""
    inner_idx: int              # index into inner values array
    outer_combo_idx: int        # flattened index of outer combination
    outer_values: dict          # {name: value} for all outer params


class Sweep:
    """
    Multi-parameter sweep.  First added = innermost (fastest changing).
    """

    def __init__(self):
        self._inner: Optional[tuple[str, np.ndarray, Optional[Callable]]] = None
        self._outers: list[tuple[str, np.ndarray, Optional[Callable]]] = []

    def add(self, name: str, values, setter: Optional[Callable] = None):
        """Add a sweep parameter.  First call sets the inner axis."""
        arr = np.asarray(values)
        if self._inner is None:
            self._inner = (name, arr, setter)
        else:
            self._outers.append((name, arr, setter))

    # ── shape queries ───────────────────────────────────────────────────

    @property
    def inner_name(self) -> str:
        return self._inner[0] if self._inner else ''

    @property
    def inner_values(self) -> np.ndarray:
        return self._inner[1] if self._inner else np.array([])

    @property
    def inner_count(self) -> int:
        return len(self._inner[1]) if self._inner else 0

    @property
    def outer_names(self) -> list[str]:
        return [n for n, _, _ in self._outers]

    @property
    def outer_shape(self) -> tuple[int, ...]:
        return tuple(len(v) for _, v, _ in self._outers) or (1,)

    @property
    def outer_count(self) -> int:
        n = 1
        for _, v, _ in self._outers:
            n *= len(v)
        return n

    @property
    def total_points(self) -> int:
        return self.inner_count * self.outer_count

    def __len__(self) -> int:
        return self.total_points

    def names(self) -> list[str]:
        """All parameter names, inner first."""
        out = []
        if self._inner:
            out.append(self._inner[0])
        out.extend(n for n, _, _ in self._outers)
        return out

    def values_for(self, name: str) -> np.ndarray:
        if self._inner and self._inner[0] == name:
            return self._inner[1]
        for n, v, _ in self._outers:
            if n == name:
                return v
        raise KeyError(name)

    # ── outer iteration ─────────────────────────────────────────────────

    def outer_combos(self):
        """
        Yield (outer_combo_idx, outer_values_dict) for each outer combo.
        Calls outer setters when values change.

        Usage:
            for oc_idx, oc_vals in sweep.outer_combos():
                for i_run in range(Nruns):
                    for si in sweep.inner_iter(oc_idx, oc_vals):
                        acquire(...)
        """
        if not self._outers:
            yield 0, {}
            return

        axes = [range(len(v)) for _, v, _ in self._outers]
        prev: list[Optional[int]] = [None] * len(self._outers)
        flat = 0

        for combo in product(*axes):
            indices = list(combo)
            # Call setters for changed outers
            for k in range(len(self._outers)):
                name, vals, setter = self._outers[k]
                if indices[k] != prev[k]:
                    if setter is not None:
                        setter(vals[indices[k]])
            prev = list(indices)

            ov = {self._outers[k][0]: self._outers[k][1][indices[k]]
                  for k in range(len(self._outers))}
            yield flat, ov
            flat += 1

    def inner_iter(self, outer_combo_idx: int = 0,
                   outer_values: Optional[dict] = None):
        """
        Yield SweepIndex for each inner point.  Calls the inner setter.
        """
        if self._inner is None:
            return
        name, vals, setter = self._inner
        ov = outer_values or {}
        for i, v in enumerate(vals):
            if setter is not None:
                setter(v)
            yield SweepIndex(inner_idx=i,
                             outer_combo_idx=outer_combo_idx,
                             outer_values=ov)

    # ── flat iteration (backward compat, wraps outer+inner) ─────────────

    def __iter__(self):
        """Yield params_dict for every point (flat, no run loop)."""
        for oc_idx, oc_vals in self.outer_combos():
            for si in self.inner_iter(oc_idx, oc_vals):
                d = dict(oc_vals)
                d[self.inner_name] = self.inner_values[si.inner_idx]
                yield d

    # ── display ─────────────────────────────────────────────────────────

    def info(self) -> str:
        lines = [f"Sweep: {self.total_points} pts = "
                 f"{self.inner_count} inner × {self.outer_count} outer"]
        if self._inner:
            n, v, s = self._inner
            sn = s.__qualname__ if s else "None"
            lines.append(f"  inner: {n}  [{v[0]:.6g}..{v[-1]:.6g}]  "
                         f"{len(v)} pts  setter={sn}")
        for i, (n, v, s) in enumerate(self._outers):
            sn = s.__qualname__ if s else "None"
            tag = "outer" if i == len(self._outers) - 1 else f"mid-{i}"
            if len(v) == 1:
                lines.append(f"  {tag}: {n} = {v[0]:.6g}  setter={sn}")
            else:
                lines.append(f"  {tag}: {n}  [{v[0]:.6g}..{v[-1]:.6g}]  "
                             f"{len(v)} pts  setter={sn}")
        return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════════
#  PulseBlaster smart setter factory
# ═══════════════════════════════════════════════════════════════════════════

def make_pb_setter(pb, param_index: int, instr: str, sequence: str,
                   seq_args_template: list, pb_channels: list) -> Callable:
    """
    Closure: overwrite seq_args[param_index], reprogram + restart PB.

    Works for any PB time parameter: t_MW (Rabi), tau (echo),
    relax_time (T1), t_precession (Ramsey), etc.

    The param_index corresponds to the position in sequenceArgs:
      Rabi:  [t_MW, t_AOM, ro_delay, AOM_lag, MW_lag]  → index 0
      Echo:  [tau, t_AOM, ro_delay, AOM_lag, MW_lag, t_pi] → index 0
      T1:    [relax_time, t_AOM, AOM_lag, ro_delay, MW_lag, t_pi] → index 0
      ESR:   PB doesn't change → don't use this, use sg.set_freq directly
    """
    from PBcontrol import PulseBlaster as PB

    def setter(value):
        seq_args_template[param_index] = value
        _, the_list = PB.PB_program(
            instr, sequence, seq_args_template + [pb_channels])
        instruction_list = [the_list[i][0] for i in range(len(the_list))]
        pb.run_sequence_for_diode(instruction_list)

    return setter
