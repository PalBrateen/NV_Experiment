# NV Experiment GUI - Quick Start Guide

## 🚀 Launch the GUI

```bash
cd d:\Brateen\NV_Experiment
python gui/main.py
```

---

## 🎯 What You'll See

### Main Window Layout

```
┌──────────────────────────────────────────────────────────────────┐
│ [File] [Tools] [View] [Help]              [🟢SG] [⚫PB] [⚫DAQ] │
├────────────────────────────────┬─────────────────────────────────┤
│                                │                                 │
│  [Instruments] [Config]        │    IPython Console              │
│  [Monitor] [Data]              │    In [1]: gui.                 │
│                                │                                 │
│  Instrument Status Cards:      │    Variables:                   │
│  • Signal Generator            │    - gui, sg, np, plt          │
│  • PulseBlaster               │    - data                       │
│  • DAQ                         │                                 │
│  • Camera                      │                                 │
│                                │                                 │
├────────────────────────────────┴─────────────────────────────────┤
│  Live Voltage (Timeseries)   │   Live Intensity (Timeseries)   │
│  [Real-time plot]            │   [Real-time plot]              │
│  Current: 2.45 V             │   Current: 1023 counts          │
│  [Pause] [Clear] [Export]    │   [Pause] [Clear] [Export]      │
├─────────────────────────────────────────────────────────────────┤
│ Status: Idle │ Last: 2.45 V │ 14:32:15                         │
└─────────────────────────────────────────────────────────────────┘
```

---

## ✨ Features to Try

### 1. **Test the Live Plots**

**Steps**:
1. Go to menu: **Tools → Generate Test Data** (check it)
2. Watch the bottom plots update in real-time!
3. **Voltage plot** (left): 1 Hz sine wave with noise
4. **Intensity plot** (right): 0.5 Hz sine wave with noise

**Controls**:
- **Timeseries ⟷ FFT**: Toggle between time-domain and frequency-domain
- **Buffer**: Adjust number of points (100-10000)
- **Update rate**: Change refresh rate (1-60 Hz)
- **Pause**: Freeze the display
- **Clear**: Erase all data
- **Export PNG**: Save plot as image

**Try FFT Mode**:
1. Let data accumulate for ~10 seconds
2. Click **FFT** radio button
3. See the dominant frequency peaks!
   - Voltage: Should see peak at ~1 Hz
   - Intensity: Should see peak at ~0.5 Hz

---

### 2. **Check Instrument Status**

**Steps**:
1. Click the **Instruments** tab
2. See 4 status cards:
   - **Signal Generator (SG384)** - Shows frequency, power, modulation
   - **PulseBlaster** - Shows sequence, timing parameters
   - **DAQ (NI P6363)** - Shows sampling rate, voltage range
   - **Camera (Hamamatsu)** - Shows exposure, ROI, trigger mode

**Status Colors**:
- 🟢 **Green dot** = Connected and FREE
- 🟡 **Yellow dot** = CLAIMED or IN_USE by an experiment window
- 🔴 **Red dot** = ERROR state
- ⚫ **Gray dot** = DISCONNECTED

**Buttons**:
- **Reconnect**: Attempt to reconnect to instrument
- **Simulate**: Toggle simulation mode
- **Test**: Test connection

---

### 3. **Use the IPython Console** (RIGHT SIDE)

The IPython console is **fully functional** and ready for troubleshooting!

**Pre-loaded Variables**:
```python
gui         # Main window instance
sg          # Signal Generator (simulated)
np          # NumPy
plt         # Matplotlib pyplot
data        # Last acquired dataset (None for now)
```

**Try These Commands**:
```python
# Access GUI properties
gui.test_time          # Current test data time
gui.test_data_enabled  # Is test data running?

# Access plot stats
gui.voltage_plot.stats
# Returns: {'current': 2.45, 'mean': 2.43, 'std': 0.12, ...}

# Get plot data
times, values = gui.voltage_plot.get_current_data()
print(f"Captured {len(values)} points")

# Quick numpy analysis
import numpy as np
data = np.random.randn(100)
print(f"Mean: {np.mean(data):.3f}, Std: {np.std(data):.3f}")

# Access signal generator (simulated)
sg.name  # 'sg1'

# List all variables
%whos

# Command history
%history
```

**Auto-Complete**:
- Type `gui.` and press **Tab** to see all methods
- Type `np.` and press **Tab** to see NumPy functions

**Syntax Highlighting**:
- Python keywords in blue
- Strings in orange
- Numbers in green
- Comments in green

---

### 4. **Explore the Menu Bar**

#### **File Menu**
- **Load Configuration...** (placeholder)
- **Save Configuration...** (placeholder)
- **Exit** - Close the GUI

#### **Tools Menu**
- **Open Sweep Window...** (placeholder - coming in Phase 2!)
- **Simulation Mode** ☑ - Toggle simulation for instruments
- **Generate Test Data** ☑ - Toggle test data for plots

#### **View Menu**
- **Toggle IPython Console** - Show/hide console
- **Toggle Live Plots** - Show/hide bottom plots

#### **Help Menu**
- **About** (placeholder)
- **Keyboard Shortcuts** (placeholder)

---

### 5. **Monitor the Status Bar** (BOTTOM)

The status bar shows:
- **Left**: Current activity status
  - "Status: Idle"
  - "Status: Generating test data..."
  - (Later: "Status: Running ESR sweep...")

- **Right**: Connection indicators for each instrument
  - **🟢 SG** = Signal Generator connected
  - **⚫ PB** = PulseBlaster disconnected
  - **⚫ DAQ** = DAQ disconnected
  - **⚫ CAM** = Camera disconnected

**Hover** over indicators to see tooltips with details!

---

## 🔧 Troubleshooting

### GUI won't launch?

**Check Python environment**:
```bash
python --version  # Should be Python 3.8+
```

**Check dependencies**:
```bash
pip list | grep -E "PySide6|qtconsole|pyqtgraph|scipy|numpy"
```

**Missing packages?**
```bash
pip install PySide6 qtconsole ipython pyqtgraph scipy numpy
```

### IPython console not working?

**Error**: "qtconsole not found"
```bash
pip install qtconsole ipython
```

### Live plots not updating?

1. Make sure **Tools → Generate Test Data** is **checked**
2. Check status bar says "Status: Generating test data..."
3. Try clicking **Clear** button on plots
4. Restart the GUI

### Dark theme not applied?

The QSS stylesheet should auto-load from `gui/styles/dark_theme.qss`.

**If it's missing**:
- Check the file exists
- The `load_stylesheet()` function in `gui/main.py` handles this

---

## 📚 What's Working (Phase 1 Complete!)

### ✅ Fully Functional:
1. **Main window** with dark theme
2. **IPython console** with dark mode and syntax highlighting
3. **Live plot widgets** (timeseries + FFT modes)
4. **Instrument status cards** (4 instruments)
5. **Resource manager** (claim/release system)
6. **Status bar indicators** with live updates
7. **Test data generator** for demonstrations
8. **FFT processor** with peak detection

### ⏳ Coming Next (Phase 2):
1. **Parameter editor** - Configure experiment parameters
2. **Experiment windows** - Separate windows for ESR, Rabi, T1, T2
3. **Sweep builder** - Visual parameter sweep configuration
4. **On-the-go fitting** - Lorentzian, exponential, Rabi fits
5. **Scan selector** - Choose which runs to plot
6. **Monitor tab** - Single-point and continuous acquisition modes

---

## 🎓 Understanding the Architecture

### Threading Model
```
Main GUI Thread
├─ UI updates (instant response)
├─ IPython console (separate kernel thread)
└─ Timer callbacks (status updates, test data)

Worker Threads (future):
├─ Acquisition thread (instrument I/O)
├─ Data saving thread (file I/O)
└─ Fitting thread (scipy calculations)
```

### Resource Management
```
InstrumentManager (Singleton)
├─ signal_generator: DISCONNECTED → FREE → CLAIMED → IN_USE
├─ pulseblaster: DISCONNECTED
├─ analog_input: DISCONNECTED
└─ camera: DISCONNECTED

Only ONE acquisition can run at a time!
```

### Data Flow (Test Mode)
```
Timer (50ms) → _generate_test_data()
    ↓
Generate sine wave + noise
    ↓
voltage_plot.add_data_point(t, value)
    ↓
LivePlotWidget updates buffer
    ↓
Timer (100ms) → _update_plot()
    ↓
pyqtgraph renders (GPU-accelerated)
```

---

## 🎯 Next Steps for You

1. **Explore the GUI**: Click around, toggle buttons, try the console
2. **Test live plotting**: Enable test data and watch FFT mode
3. **Check instrument cards**: See how status updates work
4. **Use IPython console**: Try the example commands
5. **Provide feedback**: What do you like? What needs improvement?

---

## 📖 Additional Resources

- **Full Implementation Progress**: See `GUI_IMPLEMENTATION_PROGRESS.md`
- **Complete Plan**: See `C:\Users\Smart\.claude\plans\shiny-munching-fern.md`
- **Color Constants**: `gui/styles/colors.py`
- **Dark Theme Stylesheet**: `gui/styles/dark_theme.qss`

---

**Ready to build Phase 2?** Let me know what feature you'd like next:
- Parameter editor for configuration?
- First experiment window (ESR)?
- Monitor tab for single-point acquisition?

🚀 **Happy experimenting!**
