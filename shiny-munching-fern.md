# NV Experiment Custom GUI - Implementation Plan

## Overview
Create a PySide6-based GUI for NV center experiment control with:
- **Main Window**: Instrument status, single-point/continuous monitoring, IPython console
- **Separate Experiment Windows**: Per-experiment sweep configuration with live plots and fitting
- **Resource Management**: One acquisition at a time with instrument claim/release
- **Live Plotting**: Voltage/Intensity with FFT capability
- **Dark Theme**: Consistent with existing Verdi V5 GUI

---

## Architecture

### Threading Model
```
Main GUI Thread
├─ Main Window UI updates
├─ Experiment Window UI updates (one or more windows)
└─ Embedded IPython Console

Worker Threads (one active at a time):
├─ Acquisition Thread (instruments I/O)
├─ Data Saving Thread (async file I/O)
└─ Fitting Thread (scipy.optimize)

Plotting: Main thread via Qt signals from workers
```

**Key Principle**: Only ONE acquisition thread runs at a time. Experiment windows can be configured in parallel but must claim instruments before running.

### Resource Management System
Centralized `InstrumentManager` class:
- Tracks instrument availability (FREE, CLAIMED, IN_USE)
- Handles claim/release by experiment windows
- Notifies on conflicts
- Ensures safe shutdown

```python
class InstrumentManager:
    def claim_instruments(self, requester, instruments_list) -> bool
    def release_instruments(self, requester)
    def get_instrument_status(self, instrument_name) -> (status, owner)
```

---

## Project Structure

```
gui/
├── main.py                          # Entry point, launch main window
├── main_window.py                   # Main GUI window class
├── experiment_window.py             # Separate window per experiment type
├── resource_manager.py              # Instrument claim/release system
├── gui_experiment_controller.py     # QObject wrapper for experiment controllers
│
├── widgets/
│   ├── instrument_status_card.py    # Connection status widget
│   ├── ipython_console_widget.py    # Embedded IPython console (dark mode)
│   ├── live_plot_widget.py          # Real-time plot with FFT option
│   ├── parameter_editor.py          # Unit-aware parameter inputs
│   ├── sweep_builder_widget.py      # Visual sweep configuration
│   ├── quick_fit_panel.py           # On-the-go fitting controls
│   └── scan_selector_widget.py      # Select specific scans to plot
│
├── dialogs/
│   ├── instrument_conflict_dialog.py
│   └── save_config_dialog.py
│
├── styles/
│   ├── dark_theme.qss               # Main stylesheet
│   └── colors.py                    # Color constants
│
└── utils/
    ├── unit_conversions.py          # Hz/GHz, ns/μs/ms conversion
    ├── fft_processor.py             # FFT computation (GPU later)
    └── plotting_utils.py            # Common plotting functions
```

---

## Phase 1: Foundation (Week 1-2)

### 1.1 Project Setup
**Files to create:**
- `gui/main.py` - Application entry point
- `gui/styles/dark_theme.qss` - Based on `verdi_v5_gui.py:32-75`
- `gui/styles/colors.py` - Color constants

**Dependencies:**
```python
# Add to requirements.txt or environment
PySide6
qtconsole
ipython
pyqtgraph
scipy
numpy
```

### 1.2 Resource Management System
**File**: `gui/resource_manager.py`

**Features:**
- Singleton `InstrumentManager` class
- Track instruments: `SignalGenerator`, `PulseBlaster`, `AnalogInputTask`, `AnalogOutputTask`, `Camera`
- States: `FREE`, `CLAIMED`, `IN_USE`
- Thread-safe with `QMutex`
- Emit signals on state changes

**Integration Points:**
- Wraps existing instruments from:
  - `SGcontrol.py`
  - `PBcontrol.py`
  - `DAQcontrol.py`
  - `Camcontrol.py`

### 1.3 GUI Experiment Controller Wrapper
**File**: `gui/gui_experiment_controller.py`

**Purpose**: Thread-safe wrapper around `DiodeExperimentController` and `CameraExperimentController`

```python
class GuiExperimentController(QObject):
    # Signals
    progress_updated = Signal(int, int)  # current, total
    data_point_acquired = Signal(dict)   # {param_values, signal, metadata}
    error_occurred = Signal(str)
    experiment_finished = Signal()

    # Methods
    def configure(self, config_dict)
    def start_acquisition(self)
    def pause_acquisition(self)
    def stop_acquisition(self)
    def release_instruments(self)
```

**Integration Points:**
- Uses `experiment_controller.py:DiodeExperimentController` / `CameraExperimentController`
- Uses `parameter_system.py:ParameterSweep`
- Uses `save_manager.py` for data saving

---

## Phase 2: Main Window (Week 2-3)

### 2.1 Main Window Shell
**File**: `gui/main_window.py`

**Layout:**
```
┌─────────────────────────────────────────────────────────────┐
│ [File] [Tools] [Help]                          [● SG] [● PB]│ ← MenuBar + Status
├─────────────────────────────────────────────────────────────┤
│ ┌─ Left Panel (60%) ────────┐  ┌─ IPython Console (40%) ──┐│
│ │  ┌── Tab: Instruments ───┐│  │  In [1]: gui.counter     ││
│ │  │ [Instrument Cards]    ││  │  Out[1]: 0               ││
│ │  └───────────────────────┘│  │                          ││
│ │  ┌── Tab: Config ────────┐│  │  In [2]: sg.frequency   ││
│ │  │ [Parameters]          ││  │  Out[2]: 2.87e9         ││
│ │  └───────────────────────┘│  │                          ││
│ │  ┌── Tab: Monitor ───────┐│  │  Variables:             ││
│ │  │ [▶ Start] [⏹ Stop]   ││  │  - gui, sg, pb, daq     ││
│ │  │ Mode: ⦿ Single Point  ││  │  - camera, np, plt      ││
│ │  │       ○ Continuous    ││  │                          ││
│ │  └───────────────────────┘│  │                          ││
│ │  ┌── Tab: Data ──────────┐│  │                          ││
│ │  │ [Save Settings]       ││  │                          ││
│ │  └───────────────────────┘│  │                          ││
│ └───────────────────────────┘  └──────────────────────────┘│
├─────────────────────────────────────────────────────────────┤
│ ┌─ Live Voltage (50%) ───────┐  ┌─ Live Intensity (50%) ──┐│
│ │ Display: ⦿ Timeseries       │  │ Display: ⦿ Timeseries   ││
│ │          ○ FFT              │  │          ○ FFT          ││
│ │ [Voltage vs Time Plot]      │  │ [Intensity vs Time]     ││
│ │ Buffer: 1000 points         │  │ [or Camera Image]       ││
│ └─────────────────────────────┘  └──────────────────────────┘│
├─────────────────────────────────────────────────────────────┤
│ Status: Idle │ Last acquisition: 2.45 V │ 14:32:15          │ ← StatusBar
└─────────────────────────────────────────────────────────────┘
```

**Key Features:**
- **Dockable IPython Console** (right side, resizable, can minimize/maximize)
- **Tabbed Left Panel**: Instruments, Config, Monitor, Data
- **Bottom Live Plots**: Dockable, expandable, timeseries/FFT toggle
- **Status Bar**: Instrument connection indicators

**Reference Files:**
- Layout structure from `verdi_v5_gui.py`
- Dark theme from `qtconsole_example.py:312-331`

### 2.2 IPython Console Widget
**File**: `gui/widgets/ipython_console_widget.py`

**Implementation:**
- Use `qtconsole.rich_ipython_widget.RichIPythonWidget`
- Configure dark mode from `qtconsole_example.py:174-289`
- Expose namespace: `gui`, `sg`, `pb`, `daq`, `camera`, `current_experiment`, `data`, `np`, `plt`
- Thread-safe: Console runs in kernel thread, GUI interactions via signals

**Features:**
- Syntax highlighting (pygments 'native' style)
- Auto-complete
- History navigation
- Variable inspection
- Access to live data during acquisition

### 2.3 Instrument Status Tab
**File**: `gui/widgets/instrument_status_card.py`

**Design:**
```
┌─ Signal Generator (SRS SG384) ──────────────────┐
│ Status: ● Connected (COM3)                       │
│ Frequency: 2.870 GHz │ Power: 8.0 dBm           │
│ Output: ON           │ Modulation: Pulse        │
│ Owner: [Main Window] │ State: FREE              │
│ [Reconnect] [Simulate] [Test Connection]        │
└──────────────────────────────────────────────────┘

┌─ PulseBlaster (SpinCore) ───────────────────────┐
│ Status: ● Connected                              │
│ Sequence: esr_seq    │ t_AOM: 8000 μs           │
│ Last Program: 14:32  │ Status: Idle             │
│ Owner: [ESR Window #1] │ State: CLAIMED         │
│ [Reconnect] [Verify Program]                    │
└──────────────────────────────────────────────────┘

[Similar cards for DAQ, Camera]
```

**Features:**
- Color-coded status: Green (connected), Yellow (claimed), Red (error), Gray (disconnected)
- Show current owner (main window or experiment window)
- Polling timer (1 Hz) for status updates
- Simulation mode toggle
- Reconnect on failure

**Integration:**
- Queries `InstrumentManager` for ownership status
- Reads from instrument instances for parameter display

### 2.4 Configuration Tab
**File**: `gui/widgets/parameter_editor.py`

**Purpose**: Configure parameters for single-point/continuous monitoring

**Layout:**
```
Experiment Type: [ESR ▼]  Mode: ⦿ Diode  ○ Camera

┌─ Microwave Parameters ─────────────────────────┐
│ Frequency:  [2.870] [GHz ▼]  Range: 0-20 GHz   │
│ Power:      [8    ] [dBm ▼]  Range: -130-25 dBm│
└─────────────────────────────────────────────────┘

┌─ Timing Parameters (PulseBlaster) ─────────────┐
│ t_AOM:      [8000 ] [μs ▼]   Min: 100 ns       │
│ ro_delay:   [1000 ] [ns ▼]   Min: 10 ns        │
│ Sequence:   [esr_seq ▼]                        │
└─────────────────────────────────────────────────┘

┌─ Acquisition Settings ─────────────────────────┐
│ Samples/point: [8    ]                          │
│ Sampling rate: [1.0  ] [MS/s ▼]                │
└─────────────────────────────────────────────────┘

[Load Template ▼] [Apply to Instruments]
```

**Features:**
- Custom `UnitSpinBox`: Dropdown for unit selection (Hz/kHz/MHz/GHz, ns/μs/ms/s)
- Auto-convert between units
- Validation: Red border if out of valid range, tooltip shows limits
- Load from `.py` or `.yaml` configs using `config_loader.py`
- "Apply" button sends to instruments

**Integration:**
- Loads configs from `esr_config_new.py`, `rabi_config_new.py`, etc.
- Uses `config_loader.py:load_config_from_yaml()` or direct Python import
- Updates instrument parameters via `InstrumentManager`

### 2.5 Monitor Tab (Single-Point/Continuous Mode)
**File**: Main window, monitor tab section

**Layout:**
```
┌─ Monitoring Control ────────────────────────────┐
│  Mode: ⦿ Single Point  ○ Continuous             │
│                                                  │
│  [▶ START]  [⏹ STOP]  [Test Single Point]      │
│                                                  │
│  Status: Acquiring...                           │
│  Signal: 2.341 V  │  Samples: 8                 │
│  Rate: 10 Hz      │  Buffer: 1000 pts           │
│                                                  │
│  ☑ Auto-save data (every 100 points)            │
└──────────────────────────────────────────────────┘

┌─ Quick Field Control (if using coils) ─────────┐
│  Bx: [0.00] mT  │  By: [0.00] mT  │  Bz: [0.00]│
│  [Field ON]  [Field OFF]  [Load Preset ▼]      │
└──────────────────────────────────────────────────┘
```

**Features:**
- **Single Point**: Acquire one measurement, display result
- **Continuous**: Stream data to live plots at specified rate
- Large START/STOP buttons
- Display current signal value
- Field control sliders (if coils connected)

**Integration:**
- Uses `GuiExperimentController` in single-point mode
- Emits data points to live plot widgets via signals
- Saves to file using `save_manager.py`

### 2.6 Live Plot Widget (Voltage/Intensity with FFT)
**File**: `gui/widgets/live_plot_widget.py`

**Design:**
```
┌─ Live Voltage ──────────────────────────────────┐
│ Display: ⦿ Timeseries  ○ FFT                    │
│                                                  │
│ [Real-time Line Plot]                           │
│ X: Time (s)  │  Y: Voltage (V)                  │
│ Buffer: [1000] points  │  Update: [10] Hz       │
│                                                  │
│ Current: 2.45 V │ Mean: 2.43 V │ Std: 0.12 V    │
│ [Clear] [Export PNG] [Pause Updates]            │
└──────────────────────────────────────────────────┘

[Same for Intensity plot]
```

**Features:**
- **Timeseries Mode**:
  - Rolling buffer (configurable size)
  - Real-time updates via `pyqtgraph`
  - Auto-scale or manual axes

- **FFT Mode**:
  - Compute FFT of buffer data
  - X: Frequency (Hz), Y: Magnitude
  - Peak detection and labeling
  - GPU acceleration placeholder (use CPU for now, GPU later via `cupy`)

- **Statistics Display**: Current, mean, std, min, max
- **Controls**: Clear, pause, export

**Implementation:**
- Use `pyqtgraph.PlotWidget` for speed
- Circular buffer for data storage
- FFT computed in separate thread (or GPU later)
- Update at 10-30 Hz max (balance responsiveness vs CPU)

**Integration:**
- Connected to `GuiExperimentController.data_point_acquired` signal
- For camera: Show intensity from ROI or full image statistics

### 2.7 Data Management Tab
**File**: Main window, data tab section

**Layout:**
```
┌─ Save Settings ─────────────────────────────────┐
│ Directory: [d:\Data\NV\2025-12-17          ]    │
│            [Browse...]                          │
│                                                  │
│ Filename Prefix: [NV1_monitor_]                 │
│ ☑ Auto-append timestamp                         │
│                                                  │
│ Save Format:                                     │
│ ☑ Parameters (YAML)                             │
│ ☑ Raw data (NumPy .npy)                         │
│ ☐ Plot image (PNG)                              │
│                                                  │
│ Auto-save interval: [100] points                │
└──────────────────────────────────────────────────┘

┌─ Recent Files ──────────────────────────────────┐
│ 📄 2025-12-17_1432_NV1_monitor.npy              │
│    1000 points │ 2.45 V avg                     │
│    [Load] [Delete]                              │
│ ...                                              │
└──────────────────────────────────────────────────┘
```

**Features:**
- Directory selection with auto-create
- Filename with timestamp
- Format selection (YAML params, NPY data, PNG plots)
- Auto-save during continuous monitoring
- Quick access to recent files

**Integration:**
- Uses `save_manager.py` backend
- Stores git commit hash in metadata for reproducibility

---

## Phase 3: Experiment Windows (Week 3-5)

### 3.1 Experiment Window Architecture
**File**: `gui/experiment_window.py`

**Purpose**: Separate floating window for each experiment (ESR, Rabi, T1, T2, etc.)

**How it's launched:**
1. User configures sweep parameters in main window
2. Clicks **"Open Sweep Window"** button
3. New `ExperimentWindow` instance created
4. Window attempts to **claim instruments** from `InstrumentManager`
5. If successful, window opens; if conflict, show dialog
6. User can configure multiple windows but only **run one at a time**

**Layout:**
```
┌─ ESR Experiment #1 ────────────────────────────────────────┐
│ [▶ START] [⏸ PAUSE] [⏹ STOP] [⚠ ABORT] [🗑 Close Window] │
├────────────────────────────────────────────────────────────┤
│ ┌─ Sweep Configuration (30%) ─┐  ┌─ Live Plot (70%) ─────┐│
│ │ ≡ Frequency                  │  │                        ││
│ │   Start: 2.82 GHz            │  │  [Real-time ESR Plot] ││
│ │   Stop:  2.92 GHz            │  │                        ││
│ │   Points: 101                │  │  X: Frequency (GHz)   ││
│ │                              │  │  Y: Signal (V)        ││
│ │ ≡ Power                      │  │                        ││
│ │   Values: [5, 8, 10] dBm    │  │  Contrast: 12.3%      ││
│ │                              │  │                        ││
│ │ Total: 303 points            │  │  ┌─ Quick Fit ───────┐││
│ │ Est. time: ~5m 3s            │  │  │ Model: Lorentzian ││
│ │                              │  │  │ Center: 2.8702 GHz││
│ │ [Add Parameter ▼]            │  │  │ Width: 8.2 MHz    ││
│ │                              │  │  │ R²: 0.9876        ││
│ │ ┌─ Scan Selector ──────────┐│  │  │ [Fit] [Copy]      ││
│ │ │ ☑ Run 1 (5 dBm)          ││  │  └───────────────────┘││
│ │ │ ☑ Run 2 (8 dBm)          ││  │                        ││
│ │ │ ☑ Run 3 (10 dBm)         ││  │  ┌─ Selected Scans ──┐││
│ │ │ ☐ Run 4 (pending...)     ││  │  │ ☑ Run 1 (5 dBm)   ││
│ │ │                          ││  │  │ ☑ Run 2 (8 dBm)   ││
│ │ │ [Plot Selected]          ││  │  │ ☐ Run 3 (10 dBm)  ││
│ │ └──────────────────────────┘│  │  │ [Overlay] [Stack] ││
│ │                              │  │  └───────────────────┘││
│ │ Progress: 47/303 (15.5%)    │  │                        ││
│ │ [████░░░░░] 4m 10s left     │  │  [Export PNG]         ││
│ │                              │  │  [Save Data]          ││
│ └──────────────────────────────┘  └────────────────────────┘│
├────────────────────────────────────────────────────────────┤
│ Status: Running │ Current: 2.8647 GHz, 8 dBm │ 14:32:45    │
└────────────────────────────────────────────────────────────┘
```

**Key Features:**
- **Sweep Builder** (left panel): Drag-reorder parameters, visual sweep config
- **Live Plot** (right panel): Real-time data with auto-fit overlay
- **Scan Selector**: After acquisition, select which runs/scans to plot
- **Quick Fit Panel**: On-the-go fitting with result display
- **Progress Tracking**: Bar, percentage, time remaining
- **Independent Configuration**: Copies parameters from main window at creation

**Window Lifecycle:**
1. **Created**: User clicks "Open Sweep Window" in main window
2. **Configured**: Inherits parameters, user can modify
3. **Claim Instruments**: Before START, claims SG, PB, DAQ/Camera
4. **Running**: Acquisition thread active, live updates
5. **Paused/Stopped**: User control, instruments still claimed
6. **Finished**: Auto-release instruments OR user keeps window open
7. **Closed**: Release all resources

**Multiple Windows:**
- User can open multiple experiment windows (e.g., ESR Window #1, Rabi Window #1)
- Only one can **run** at a time (enforced by `InstrumentManager`)
- Others are "configured" and ready to run in sequence
- User manually switches between windows to run experiments

### 3.2 Sweep Builder Widget
**File**: `gui/widgets/sweep_builder_widget.py`

**Design:**
```
┌─ Sweep Parameters ──────────────────────────────┐
│ Drag to reorder (top = outer loop)              │
│                                                  │
│ ┌────────────────────────────────────────────┐  │
│ │ ≡ Frequency                          [×]   │  │
│ │   ⦿ Linear   ○ Log   ○ Manual Array       │  │
│ │   Start: [2.82] [GHz▼]  Stop: [2.92] [GHz]│  │
│ │   Points: [101]  ➜  Step: 1.0 MHz         │  │
│ │   Preview: [2.82, 2.821, ..., 2.92] GHz   │  │
│ └────────────────────────────────────────────┘  │
│                                                  │
│ ┌────────────────────────────────────────────┐  │
│ │ ≡ Power                              [×]   │  │
│ │   ⦿ Manual Array   ○ Linear                │  │
│ │   Values: [5, 8, 10] dBm  (3 points)      │  │
│ └────────────────────────────────────────────┘  │
│                                                  │
│ [+ Add Parameter ▼]                              │
│                                                  │
│ Total measurements: 303 (101 × 3)               │
│ Estimated time: ~5 min 3 sec                    │
│ Memory: ~24 KB                                  │
└──────────────────────────────────────────────────┘
```

**Features:**
- **Drag-and-drop**: Reorder parameters to change loop nesting
- **Collapsible panels**: Each parameter in expandable card
- **Multiple input modes**: Linear sweep, log sweep, manual array
- **Smart parameter list**: Dropdown shows parameters from connected instruments
- **Real-time preview**: Show array values, total points, estimated time
- **Validation**: Check if values in valid ranges, highlight errors

**Implementation:**
- Use `QListWidget` with custom item delegates for drag-drop
- Custom widget class `SweepParameterCard` for each parameter
- Connect to `parameter_system.py:ParameterSweep` backend
- Estimate time based on sequence timing from PulseBlaster

**Integration:**
- Reads available parameters from `InstrumentManager`
- Generates config dict compatible with `experiment_controller.py`

### 3.3 Live Plot with Fitting
**File**: `gui/widgets/quick_fit_panel.py` + plot area in experiment window

**Plot Area:**
- Use `pyqtgraph` for real-time updates
- 1D plot for single parameter sweep
- 2D heatmap for dual parameter sweep
- Show current acquisition point as marker
- Auto-scale or manual axes
- Colormap selection for 2D

**Fit Panel:**
```
┌─ Quick Fit ─────────────────────────────────────┐
│ Model: [Lorentzian ▼]  [Fit Now]  ☑ Auto-fit   │
│                                                  │
│ ┌─ Fit Results ────────────────────────────────┐│
│ │ Center: 2.8702 ± 0.0003 GHz                 ││
│ │ Linewidth: 8.2 ± 0.5 MHz (FWHM)             ││
│ │ Contrast: 12.4 ± 0.8 %                      ││
│ │ R²: 0.9876                                   ││
│ │ Status: ✓ Converged                         ││
│ └──────────────────────────────────────────────┘│
│ [Copy Results] [Save Fit] [Show Residuals]     │
└──────────────────────────────────────────────────┘
```

**Available Fit Models:**
- **Lorentzian**: ESR, single resonance
- **Double Lorentzian**: Two resonances
- **Gaussian**: Inhomogeneous broadening
- **Exponential Decay**: T1, T2
- **Damped Sinusoid**: Rabi oscillations
- **Custom**: User-defined function (via IPython console)

**Features:**
- Auto-fit: Runs fitting every N points
- Manual fit: User clicks "Fit Now"
- Display fit curve overlay on plot
- Show fit parameters with uncertainties (from covariance matrix)
- R² goodness of fit
- Copy results to clipboard (formatted for lab notebook)
- Save fit parameters to file

**Implementation:**
- Use `scipy.optimize.curve_fit` for fitting
- Run in separate `QThread` to avoid blocking GUI
- Initial guess: Smart defaults based on data (e.g., peak finding for center)
- Handle fitting failures gracefully

### 3.4 Scan Selector Widget
**File**: `gui/widgets/scan_selector_widget.py`

**Purpose**: After multi-run or multi-parameter sweep, select specific scans to plot

**Design:**
```
┌─ Scan Selector ─────────────────────────────────┐
│ Acquired Scans:                                  │
│ ☑ Run 1: Power=5 dBm  (303 pts) ✓ Complete      │
│ ☑ Run 2: Power=8 dBm  (303 pts) ✓ Complete      │
│ ☑ Run 3: Power=10 dBm (303 pts) ✓ Complete      │
│ ☐ Run 4: Power=12 dBm (150 pts) ⚠ In Progress   │
│                                                  │
│ Display Mode:                                    │
│ ⦿ Overlay (same plot)  ○ Stack (separate axes)  │
│                                                  │
│ [Plot Selected] [Clear Selection] [Select All]  │
└──────────────────────────────────────────────────┘
```

**Features:**
- Checkbox list of all acquired scans
- Show metadata: Parameter values, number of points, status
- **Overlay mode**: Plot selected scans on same axes (different colors)
- **Stack mode**: Separate subplots for each scan
- Auto-update as new scans complete
- Export selected scans to file

**Use Case:**
User runs ESR at 3 different powers. After acquisition:
1. Checkboxes show all 3 scans
2. User selects scans 1 and 3
3. Clicks "Plot Selected"
4. Plot updates to show only those two scans overlaid
5. User can fit each individually or together

**Integration:**
- Reads scan data from `GuiExperimentController`
- Updates plot widget with selected data
- Connects to fit panel for per-scan fitting

### 3.5 Save Configuration in Experiment Window
**File**: Experiment window, save section

**Features:**
- Save button in experiment window
- Save current plot as PNG
- Save all scan data as NPY
- Save parameters as YAML
- Include fit results in metadata
- Auto-save after each run completion

**Directory Structure:**
```
d:\Data\NV\2025-12-17\
├── 1432_ESR_NV1\
│   ├── params.yaml          # All parameters
│   ├── data_run1.npy        # Raw data per run
│   ├── data_run2.npy
│   ├── data_run3.npy
│   ├── fit_results.yaml     # Fit parameters
│   ├── plot_run1.png        # Individual plots
│   ├── plot_overlay.png     # Overlay plot
│   └── metadata.json        # Git hash, timestamp, etc.
```

**Integration:**
- Uses `save_manager.py` backend
- Preserves experiment window state for later reload

---

## Phase 4: Advanced Features (Week 5-6)

### 4.1 FFT Processing
**File**: `gui/utils/fft_processor.py`

**Features:**
- Compute FFT of timeseries data (voltage or intensity)
- Window functions: Hanning, Hamming, Blackman
- Peak detection in frequency domain
- Normalization options
- GPU acceleration placeholder (use `cupy` later)

**Implementation (CPU for now):**
```python
import numpy as np
from scipy.fft import fft, fftfreq

def compute_fft(data, sampling_rate, window='hanning'):
    # Apply window
    if window == 'hanning':
        windowed = data * np.hanning(len(data))
    # Compute FFT
    fft_result = fft(windowed)
    freqs = fftfreq(len(data), 1/sampling_rate)
    # Return positive frequencies only
    mask = freqs > 0
    return freqs[mask], np.abs(fft_result[mask])
```

**GPU Acceleration (Phase 2, later):**
- Use `cupy` for FFT on GPU
- Auto-detect GPU availability
- Fallback to CPU if no GPU

### 4.2 Resource Conflict Handling
**File**: `gui/dialogs/instrument_conflict_dialog.py`

**Scenario**: User tries to start experiment window #2 while #1 is running

**Dialog:**
```
┌─ Instrument Conflict ───────────────────────────┐
│  ⚠️  Cannot start ESR Window #2                 │
│                                                  │
│  The following instruments are in use:          │
│  • Signal Generator (owner: ESR Window #1)      │
│  • PulseBlaster (owner: ESR Window #1)          │
│                                                  │
│  Options:                                        │
│  [ ] Wait for Window #1 to finish               │
│  [ ] Force stop Window #1 (⚠️ data may be lost) │
│  [Cancel]                  [OK]                 │
└──────────────────────────────────────────────────┘
```

**Features:**
- Show which instruments are unavailable
- Show current owner
- Options: Wait, force stop (with warning), cancel
- Queue system (optional): Auto-start when instruments free

### 4.3 Camera Tab (Placeholder)
**File**: Main window, camera tab (placeholder for now)

**To be designed later based on user requirements**

Suggested features:
- Live camera view
- ROI selection tool
- Exposure time control
- Trigger mode selection
- Snapshot capture
- Integration with intensity live plot

---

## Phase 5: Testing & Polish (Week 6-7)

### 5.1 Simulation Mode
**Feature**: Test GUI without hardware

**Implementation:**
- Use existing simulation classes:
  - `SignalGenerator_sim` from `SGcontrol.py`
  - Create similar sim classes for PB, DAQ, Camera
- Generate realistic fake data with noise
- Simulate acquisition delays
- Toggle via checkbox in main window

### 5.2 Error Handling
**Strategy:**
- Try-except around all instrument operations
- Display errors in status bar + dialog
- Log to IPython console
- Auto-reconnect on recoverable errors
- Safe shutdown on fatal errors

### 5.3 Performance Optimization
**Targets:**
- GUI refresh rate: 30 FPS minimum
- Plot update rate: 10-30 Hz
- Acquisition overhead: <100 μs per point
- Memory: No leaks during long acquisitions

**Profiling:**
- Use `cProfile` to identify bottlenecks
- Monitor memory with `tracemalloc`
- Optimize plot updates (downsample if >10k points)

### 5.4 Documentation
**Deliverables:**
- User manual with screenshots
- Keyboard shortcuts list
- Tooltips on all controls
- Example workflows (ESR, Rabi, T1, T2)
- Troubleshooting guide
- Video tutorial (optional)

---

## Critical Files to Integrate

### Existing Codebase:
1. **`mainControl.py`** - Entry point, adapt for GUI launch
2. **`experiment_controller.py`** - Wrap in `GuiExperimentController`
3. **`parameter_system.py`** - Use `ParameterSweep` for sweep logic
4. **`config_loader.py`** - Load YAML/Python configs
5. **`save_manager.py`** - Data saving backend
6. **`SGcontrol.py`** - Signal generator control
7. **`PBcontrol.py`** - PulseBlaster control
8. **`DAQcontrol.py`** - DAQ tasks (analog in/out)
9. **`Camcontrol.py`** - Camera control

### Reference Implementations:
10. **`verdi_v5_gui.py`** - GUI structure, dark theme, status indicators
11. **`qtconsole_example.py`** - IPython console integration, dark mode

---

## Implementation Priority

### Must-Have (Priority 1):
1. Main window shell with tabs
2. Instrument status cards
3. IPython console integration (dark mode)
4. Resource manager (claim/release)
5. Single experiment window (ESR)
6. Live plot widget (timeseries mode)
7. Basic sweep builder
8. Start/Stop/Abort controls
9. Save configuration

### Important (Priority 2):
10. Multiple experiment windows
11. Scan selector widget
12. On-the-go fitting (Lorentzian)
13. FFT mode in live plots
14. Parameter editor with unit conversion
15. Monitor tab (single-point/continuous)
16. Quick fit panel

### Nice-to-Have (Priority 3):
17. 2D heatmap plotting
18. Advanced fitting models
19. GPU-accelerated FFT
20. Conflict resolution dialog
21. Data browser
22. Camera tab
23. Dark theme polish
24. Keyboard shortcuts

---

## Development Workflow

### Phase-by-Phase:
1. **Week 1**: Project setup, resource manager, dark theme
2. **Week 2**: Main window shell, IPython console, instrument status
3. **Week 3**: First experiment window (ESR), live plot, sweep builder
4. **Week 4**: Fitting panel, scan selector, multiple experiment windows
5. **Week 5**: FFT integration, monitor tab, parameter editor
6. **Week 6**: Camera placeholder, conflict handling, testing
7. **Week 7**: Polish, documentation, user testing

### Testing Strategy:
- **Unit tests**: Resource manager, FFT processor, unit conversions
- **Integration tests**: Instrument claim/release, config loading
- **GUI tests**: Manual testing with simulation mode
- **Hardware tests**: Test with actual instruments (SG, PB, DAQ, Camera)

---

## Key Design Decisions

### Threading Model:
✅ **Use QThread** for acquisition, saving, fitting
✅ **Main thread** for GUI and plotting (via signals)
❌ **No multiprocessing** (unnecessary overhead for I/O-bound tasks)

### Plotting Library:
✅ **pyqtgraph** for real-time live plots (fast)
✅ **matplotlib** for export/publication quality (optional)

### GUI Framework:
✅ **PySide6** (better licensing than PyQt5, Qt6 features)

### Resource Management:
✅ **Centralized InstrumentManager** (singleton)
✅ **Claim/Release pattern** (explicit ownership)
✅ **One acquisition at a time** (enforced)

### Configuration Format:
✅ **Dual support**: Python `.py` and YAML `.yaml`
✅ **Use existing config_loader.py**

---

## Success Criteria

### Functional:
- ✅ Can connect to all instruments (SG, PB, DAQ, Camera)
- ✅ Can configure and run sweeps (ESR, Rabi, T1, T2)
- ✅ Can monitor single-point/continuous acquisitions
- ✅ Can open multiple experiment windows (run sequentially)
- ✅ Can fit data on-the-go (Lorentzian, exponential, etc.)
- ✅ Can plot timeseries and FFT
- ✅ Can save data and reload parameters
- ✅ IPython console works for troubleshooting

### Performance:
- ✅ GUI never freezes during acquisition
- ✅ Plot updates at 10+ Hz
- ✅ Acquisition overhead <100 μs per point
- ✅ No memory leaks during long runs

### UX:
- ✅ Dark theme consistent throughout
- ✅ Intuitive workflow (minimal clicks)
- ✅ Clear error messages
- ✅ Keyboard shortcuts for common actions
- ✅ Tooltips and documentation

---

## Next Steps

1. **User approval** of this plan
2. **Create project structure** (`gui/` folder)
3. **Implement resource manager** (foundation)
4. **Build main window shell** (GUI framework)
5. **Integrate IPython console** (dark mode)
6. **Develop first experiment window** (ESR as template)
7. **Test with simulation mode**
8. **Iterate based on user feedback**

---

## Questions for User (if any remain)

None - all clarified! Ready to proceed with implementation upon approval.
