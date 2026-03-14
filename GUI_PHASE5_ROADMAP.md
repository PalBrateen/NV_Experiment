# NV Experiment GUI - Phase 5 Roadmap

## Status: Planning Phase 📋

**Date**: 2025-12-17
**Goal**: Complete multi-parameter sweep integration with 2D/3D visualization
**Estimated Effort**: 3-5 hours of development

---

## Current State (After Phase 4)

### ✅ What's Complete:
- **Phase 1**: Foundation (resource manager, IPython console, live plots, main window)
- **Phase 2**: Experiment windows with fitting and scan comparison
- **Phase 3**: Real hardware integration (dual-mode acquisition)
- **Phase 4**: Multi-sweep widget integrated into Config tab

### 🟡 What's Partially Done:
- Multi-sweep parameters stored in `current_parameters['multi_sweep_params']`
- ExperimentWindow still uses single `SweepBuilderWidget` (1D only)
- No 2D/3D visualization yet

---

## Phase 5 Goals

### Goal 1: ExperimentWindow Multi-Sweep Support
**Priority**: High
**Effort**: 2-3 hours

Enable experiment windows to use multi-parameter sweeps from the Config tab.

### Goal 2: 2D Heatmap Visualization
**Priority**: High
**Effort**: 1-2 hours

Add heatmap display for 2D sweep results.

### Goal 3: Data Export
**Priority**: Medium
**Effort**: 1 hour

Export scans and fit parameters to CSV/YAML.

---

## Phase 5A: Multi-Sweep Integration in ExperimentWindow

### Objective
Modify `ExperimentWindow` to detect and use multi-sweep parameters from Config tab.

### Implementation Steps

#### Step 1: Detect Multi-Sweep Configuration

**File**: [gui/windows/experiment_window.py](gui/windows/experiment_window.py)

**Location**: In `__init__()` method

**Change**:
```python
def __init__(self, window_id: str, parameters: Dict[str, Any], parent=None):
    # ... existing code ...

    # Check for multi-sweep parameters
    if 'multi_sweep_params' in parameters and parameters['multi_sweep_params']:
        self.multi_sweep_params = parameters['multi_sweep_params']
        self.use_multi_sweep = True
        self.sweep_dimensions = len(self.multi_sweep_params)
    else:
        self.multi_sweep_params = []
        self.use_multi_sweep = False
        self.sweep_dimensions = 1

    # Store all parameters
    self.set_parameters(parameters)

    # ... rest of existing code ...
```

**Impact**: Window now knows if it should use multi-sweep

---

#### Step 2: Replace or Extend Sweep Builder

**Option A: Conditional Display** (Recommended)
```python
def _create_controls_panel(self) -> QWidget:
    widget = super()._create_controls_panel()

    # Display sweep configuration based on mode
    sweep_group = QGroupBox("Sweep Configuration")
    sweep_layout = QVBoxLayout()

    if self.use_multi_sweep:
        # Show multi-sweep summary (read-only)
        self._create_multi_sweep_summary()
    else:
        # Show editable sweep builder (current behavior)
        self.sweep_builder = SweepBuilderWidget()
        self.sweep_builder.sweep_changed.connect(self._on_sweep_changed)
        sweep_layout.addWidget(self.sweep_builder)

    sweep_group.setLayout(sweep_layout)
    self.parameters_layout.addWidget(sweep_group)

    return widget
```

**Option B: Unified Widget** (Future enhancement)
- Replace `SweepBuilderWidget` with `MultiSweepBuilderWidget` everywhere
- Single sweep = multi-sweep with one parameter
- More consistent but requires more changes

**Recommendation**: Use Option A for Phase 5 (minimal changes, backward compatible)

---

#### Step 3: Multi-Sweep Summary Display

**New Method**:
```python
def _create_multi_sweep_summary(self):
    """Display read-only summary of multi-sweep configuration."""
    summary_layout = QVBoxLayout()

    # Header
    header = QLabel(f"{self.sweep_dimensions}D Sweep Configuration (from Config tab)")
    header.setStyleSheet("font-weight: bold;")
    summary_layout.addWidget(header)

    # List each parameter
    for i, param in enumerate(self.multi_sweep_params):
        level = "Inner Loop" if i == 0 else f"Level {i+1}"
        text = f"{param.name}: {param.start:.3e} → {param.stop:.3e} ({param.points} pts, {param.mode}) [{level}]"
        label = QLabel(text)
        if i == 0:
            label.setStyleSheet("color: #00ff00;")  # Green for inner loop
        summary_layout.addWidget(label)

    # Total points
    total_points = 1
    for param in self.multi_sweep_params:
        total_points *= param.points

    total_label = QLabel(f"Total Points: {total_points:,}")
    total_label.setStyleSheet("font-weight: bold; color: #00ffff;")
    summary_layout.addWidget(total_label)

    self.parameters_layout.addLayout(summary_layout)
```

**Result**: User sees what sweep is configured, but can't change it in window (must change in Config tab)

---

#### Step 4: Build Multi-Parameter params_dict

**Modified Method**: `_build_params_dict()`

**Change**:
```python
def _build_params_dict(self, sweep_array: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """Build params_dict for hardware controller."""

    if self.use_multi_sweep:
        return self._build_multi_params_dict()
    else:
        # Current single-parameter approach
        return self._build_single_params_dict(sweep_array)

def _build_single_params_dict(self, sweep_array: np.ndarray) -> Dict[str, Any]:
    """Current implementation (Phase 3)."""
    sequence_map = {'ESR': 'esr', 'Rabi': 'rabi', ...}

    params_dict = {
        'seq': {'sequence': sequence_map.get(self.exp_type, 'esr')},
        'scan': {'Nruns': self.parameters.get('num_averages', 1000)},
        'mw': {
            'freq': sweep_array if self.exp_type == 'ESR' else self.parameters.get('frequency', 2.87e9),
            'power': self.parameters.get('power', 8.0)
        },
        'daq': {...},
        'sweep': {self.config['sweep_parameter']: sweep_array}
    }
    return params_dict

def _build_multi_params_dict(self) -> Dict[str, Any]:
    """Build params_dict for multi-parameter sweep."""
    sequence_map = {'ESR': 'esr', 'Rabi': 'rabi', ...}

    # Base params
    params_dict = {
        'seq': {'sequence': sequence_map.get(self.exp_type, 'esr')},
        'scan': {'Nruns': self.parameters.get('num_averages', 1000)},
        'mw': {
            'freq': self.parameters.get('frequency', 2.87e9),
            'power': self.parameters.get('power', 8.0)
        },
        'daq': {...},
        'sweep': {}  # Multi-parameter
    }

    # Add each sweep parameter
    for param in self.multi_sweep_params:
        param_array = param.get_array()

        # Map to correct key in params_dict
        if param.name == 'frequency':
            params_dict['sweep']['frequency'] = param_array
        elif param.name == 'power':
            params_dict['sweep']['power'] = param_array
        elif param.name in ['tau', 'delay']:
            params_dict['sweep'][param.name] = param_array
        # Add more mappings as needed

    return params_dict
```

**Impact**: Hardware controller receives multi-parameter sweep configuration

---

#### Step 5: Handle Multi-Dimensional Data Storage

**Modified Method**: `_on_hardware_data_point()`

**Change**:
```python
def _on_hardware_data_point(self, point_idx: int, param_values: dict, signal: float, metadata: dict):
    """Handle data point from hardware acquisition."""

    if self.use_multi_sweep:
        # Store full parameter values for multi-dimensional data
        self.sweep_data_points.append({
            'index': point_idx,
            'params': param_values.copy(),  # e.g., {'frequency': 2.87e9, 'power': 10}
            'signal': signal,
            'metadata': metadata
        })

        # For 2D, also store in structured arrays
        if self.sweep_dimensions == 2:
            self._store_2d_data_point(point_idx, param_values, signal)
    else:
        # Current 1D approach
        sweep_param = self.config['sweep_parameter']
        param_value = param_values.get(sweep_param, 0.0)

        self.sweep_params.append(param_value)
        self.sweep_signals.append(signal)

    # Update plot (method depends on dimensions)
    self._update_live_plot(point_idx, param_values, signal)
```

**New Method**:
```python
def _store_2d_data_point(self, point_idx: int, param_values: dict, signal: float):
    """Store data point for 2D sweep."""
    inner_param = self.multi_sweep_params[0]
    outer_param = self.multi_sweep_params[1]

    # Calculate 2D indices
    inner_idx = point_idx % inner_param.points
    outer_idx = point_idx // inner_param.points

    # Initialize 2D array on first call
    if not hasattr(self, 'data_2d'):
        self.data_2d = np.full((outer_param.points, inner_param.points), np.nan)

    self.data_2d[outer_idx, inner_idx] = signal
```

**Impact**: Data properly structured for 2D/3D visualization

---

#### Step 6: Adaptive Plot Display

**Modified Method**: `_update_live_plot()`

**New Logic**:
```python
def _update_live_plot(self, point_idx: int, param_values: dict, signal: float):
    """Update plot based on sweep dimensions."""

    if self.sweep_dimensions == 1:
        # Current behavior: line plot
        inner_param = self.multi_sweep_params[0] if self.use_multi_sweep else self.sweep_builder.get_sweep_array()
        param_value = param_values.get(inner_param.name if self.use_multi_sweep else self.config['sweep_parameter'])
        display_x = param_value * self.config['x_scale']
        self.live_plot.add_data_point(display_x, signal)

    elif self.sweep_dimensions == 2:
        # Update heatmap (implemented in Phase 5B)
        self._update_heatmap()

    else:
        # 3D+: Show latest 1D slice
        self._update_nd_slice()
```

**Impact**: Plot adapts to data dimensionality

---

### Files to Modify for Phase 5A:
1. `gui/windows/experiment_window.py`: ~100 lines of changes/additions
2. No new files needed

### Testing for Phase 5A:
1. Config tab: Add 1D frequency sweep → Window uses it
2. Config tab: Add 2D frequency×power sweep → Window detects and configures
3. Simulated acquisition: Data stored correctly
4. Hardware acquisition: Multi-sweep passed to controller

---

## Phase 5B: 2D Heatmap Visualization

### Objective
Replace line plot with heatmap for 2D sweeps.

### Implementation Steps

#### Step 1: Add ImageView Widget

**File**: [gui/windows/experiment_window.py](gui/windows/experiment_window.py)

**Location**: In `_create_plot_panel()` method

**Change**:
```python
def _create_plot_panel(self) -> QWidget:
    """Override to add fit panel and support 2D display."""
    widget = QWidget()
    layout = QVBoxLayout()
    layout.setContentsMargins(0, 0, 0, 0)

    # Splitter for plot and fit panel
    splitter = QSplitter(Qt.Vertical)

    # Create appropriate plot widget based on dimensions
    if self.sweep_dimensions == 2:
        # 2D heatmap
        from pyqtgraph import ImageView
        self.image_view = ImageView()
        splitter.addWidget(self.image_view)
        self.live_plot = None  # Not used for 2D
    else:
        # 1D line plot (current behavior)
        self.live_plot = LivePlotWidget(...)
        splitter.addWidget(self.live_plot)
        self.image_view = None  # Not used for 1D

    # Fit panel (works for both 1D and 2D)
    self.fit_panel = QuickFitPanel()
    self.fit_panel.set_fit_type(self.config['fit_type'])
    splitter.addWidget(self.fit_panel)

    splitter.setSizes([700, 300])
    layout.addWidget(splitter)

    widget.setLayout(layout)
    return widget
```

**Impact**: Window shows heatmap instead of line plot for 2D data

---

#### Step 2: Update Heatmap During Acquisition

**New Method**:
```python
def _update_heatmap(self):
    """Update 2D heatmap display."""
    if not hasattr(self, 'data_2d'):
        return

    # Update image
    self.image_view.setImage(self.data_2d.T, autoRange=True, autoLevels=True)

    # Set axis labels and scales
    inner_param = self.multi_sweep_params[0]
    outer_param = self.multi_sweep_params[1]

    # Scale axes
    x_scale = [inner_param.start * self.config['x_scale'],
               inner_param.stop * self.config['x_scale']]
    y_scale = [outer_param.start, outer_param.stop]

    # Set colormap (optional)
    from pyqtgraph import colormap
    cmap = colormap.get('viridis')
    self.image_view.setColorMap(cmap)
```

**Impact**: Real-time heatmap updates during acquisition

---

#### Step 3: Axis Labeling for Heatmap

**Enhancement**:
```python
def _configure_heatmap_axes(self):
    """Configure heatmap axes with labels and scales."""
    inner_param = self.multi_sweep_params[0]
    outer_param = self.multi_sweep_params[1]

    # Get ImageItem
    img_item = self.image_view.getImageItem()

    # Set transform for proper scaling
    inner_scale = (inner_param.stop - inner_param.start) / inner_param.points
    outer_scale = (outer_param.stop - outer_param.start) / outer_param.points

    transform = QTransform()
    transform.scale(inner_scale, outer_scale)
    transform.translate(inner_param.start, outer_param.start)
    img_item.setTransform(transform)

    # Add axis labels (requires custom PlotItem)
    # This is complex in pyqtgraph ImageView, consider using PlotItem + ImageItem instead
```

**Note**: pyqtgraph's `ImageView` has limited axis customization. For full control, use `PlotItem` with `ImageItem`.

---

#### Step 4: Alternative: Use PlotItem + ImageItem

**Better Approach** (more control):
```python
def _create_plot_panel(self) -> QWidget:
    """Create plot panel with full control over axes."""
    widget = QWidget()
    layout = QVBoxLayout()

    splitter = QSplitter(Qt.Vertical)

    if self.sweep_dimensions == 2:
        # Create PlotWidget for full axis control
        import pyqtgraph as pg
        self.plot_widget_2d = pg.PlotWidget()
        self.plot_widget_2d.setLabel('bottom',
                                     self.multi_sweep_params[0].name.capitalize(),
                                     units='')
        self.plot_widget_2d.setLabel('left',
                                     self.multi_sweep_params[1].name.capitalize(),
                                     units='')

        # Add ImageItem
        self.image_item = pg.ImageItem()
        self.plot_widget_2d.addItem(self.image_item)

        # Add colorbar
        self.colorbar = pg.ColorBarItem()
        self.colorbar.setImageItem(self.image_item)
        self.plot_widget_2d.addItem(self.colorbar)

        splitter.addWidget(self.plot_widget_2d)
    else:
        # 1D plot (current)
        self.live_plot = LivePlotWidget(...)
        splitter.addWidget(self.live_plot)

    # Fit panel
    splitter.addWidget(self.fit_panel)

    layout.addWidget(splitter)
    widget.setLayout(layout)
    return widget
```

**Update Method**:
```python
def _update_heatmap(self):
    """Update heatmap using ImageItem."""
    if not hasattr(self, 'data_2d'):
        return

    # Update image
    self.image_item.setImage(self.data_2d.T, autoLevels=True)

    # Set rect for proper scaling
    inner_param = self.multi_sweep_params[0]
    outer_param = self.multi_sweep_params[1]

    rect = QRectF(
        inner_param.start,
        outer_param.start,
        inner_param.stop - inner_param.start,
        outer_param.stop - outer_param.start
    )
    self.image_item.setRect(rect)
```

**Impact**: Full control over axes, labels, colorbar

---

### Files to Modify for Phase 5B:
1. `gui/windows/experiment_window.py`: ~80 lines of changes

### Testing for Phase 5B:
1. Config tab: Configure 2D sweep (frequency × power)
2. Open experiment window: See heatmap widget instead of line plot
3. Run simulated acquisition: Heatmap updates in real-time
4. Check axes labels and scaling
5. Verify colorbar shows correct signal range

---

## Phase 5C: Data Export

### Objective
Export scan data and fit parameters to files.

### Implementation Steps

#### Step 1: Add Export Buttons

**File**: [gui/widgets/scan_selector.py](gui/widgets/scan_selector.py)

**Location**: In control buttons section

**Change**:
```python
# Existing buttons
self.select_all_btn = QPushButton("Select All")
self.deselect_all_btn = QPushButton("Deselect All")
self.change_color_btn = QPushButton("Change Color")

# NEW: Export button
self.export_btn = QPushButton("Export Selected")
self.export_btn.setProperty("class", "success")
self.export_btn.clicked.connect(self._export_selected_scans)
button_layout.addWidget(self.export_btn)
```

---

#### Step 2: Export Dialog

**New Method**:
```python
def _export_selected_scans(self):
    """Export selected scans to file."""
    from PySide6.QtWidgets import QFileDialog, QDialog, QCheckBox

    selected = self.get_selected_scans()
    if not selected:
        return

    # Choose file
    file_path, _ = QFileDialog.getSaveFileName(
        self,
        "Export Scans",
        "",
        "CSV Files (*.csv);;NumPy Files (*.npy);;YAML Files (*.yaml)"
    )

    if not file_path:
        return

    # Export based on file type
    if file_path.endswith('.csv'):
        self._export_csv(selected, file_path)
    elif file_path.endswith('.npy'):
        self._export_npy(selected, file_path)
    elif file_path.endswith('.yaml'):
        self._export_yaml(selected, file_path)
```

---

#### Step 3: Export Formats

**CSV Export** (for 1D data):
```python
def _export_csv(self, scans: List[Dict], file_path: str):
    """Export scans to CSV."""
    import csv

    with open(file_path, 'w', newline='') as f:
        writer = csv.writer(f)

        # Header
        writer.writerow(['scan_id', 'timestamp', 'x', 'y'])

        # Data
        for scan_info in scans:
            scan_id = scan_info['index']
            timestamp = scan_info['timestamp'].isoformat()
            x_data, y_data = scan_info['data']

            for x, y in zip(x_data, y_data):
                writer.writerow([scan_id, timestamp, x, y])
```

**NumPy Export** (for arrays):
```python
def _export_npy(self, scans: List[Dict], file_path: str):
    """Export scans to NumPy file."""
    import numpy as np

    # Create structured array
    data_dict = {}
    for scan_info in scans:
        scan_id = f"scan_{scan_info['index']}"
        x_data, y_data = scan_info['data']
        data_dict[scan_id] = {'x': x_data, 'y': y_data}

    np.save(file_path, data_dict)
```

**YAML Export** (for metadata + data):
```python
def _export_yaml(self, scans: List[Dict], file_path: str):
    """Export scans with metadata to YAML."""
    import yaml

    export_data = []
    for scan_info in scans:
        x_data, y_data = scan_info['data']

        scan_export = {
            'scan_id': scan_info['index'],
            'timestamp': scan_info['timestamp'].isoformat(),
            'metadata': scan_info['metadata'],
            'data': {
                'x': x_data.tolist(),
                'y': y_data.tolist()
            }
        }
        export_data.append(scan_export)

    with open(file_path, 'w') as f:
        yaml.dump(export_data, f, default_flow_style=False)
```

---

#### Step 4: Export Fit Parameters

**File**: [gui/widgets/quick_fit_panel.py](gui/widgets/quick_fit_panel.py)

**New Method**:
```python
def export_fit_parameters(self, file_path: str):
    """Export fit parameters to YAML."""
    import yaml

    export_data = {
        'fit_type': self.current_fit_type,
        'parameters': self.fit_params,
        'uncertainties': self.fit_uncertainties,
        'r_squared': float(self.r2_label.text()) if self.r2_label.text() != "N/A" else None,
        'timestamp': datetime.now().isoformat()
    }

    with open(file_path, 'w') as f:
        yaml.dump(export_data, f, default_flow_style=False)
```

**Add Button**:
```python
# In _create_ui()
self.export_fit_btn = QPushButton("Export Fit")
self.export_fit_btn.setEnabled(False)
self.export_fit_btn.clicked.connect(self._export_fit)
control_layout.addWidget(self.export_fit_btn)

def _export_fit(self):
    """Export fit parameters."""
    file_path, _ = QFileDialog.getSaveFileName(
        self, "Export Fit", "", "YAML Files (*.yaml)"
    )
    if file_path:
        self.export_fit_parameters(file_path)
```

---

### Files to Modify for Phase 5C:
1. `gui/widgets/scan_selector.py`: ~100 lines
2. `gui/widgets/quick_fit_panel.py`: ~50 lines

### Testing for Phase 5C:
1. Run multiple scans
2. Select scans in scan selector
3. Click "Export Selected"
4. Choose CSV: Verify data exported correctly
5. Choose NPY: Verify arrays saved
6. Choose YAML: Verify metadata included
7. Export fit parameters: Verify YAML contains all fit info

---

## Phase 5 Timeline

### Week 1:
- **Day 1**: Phase 5A - Multi-sweep detection and display (~2 hours)
- **Day 2**: Phase 5A - params_dict building and data storage (~1 hour)
- **Day 3**: Phase 5A - Testing with simulated data (~1 hour)

### Week 2:
- **Day 4**: Phase 5B - 2D heatmap widget setup (~1 hour)
- **Day 5**: Phase 5B - Real-time heatmap updates (~1 hour)
- **Day 6**: Phase 5B - Testing with 2D simulated sweeps (~1 hour)

### Week 3:
- **Day 7**: Phase 5C - Export functionality (~1 hour)
- **Day 8**: Phase 5C - Testing all export formats (~0.5 hour)
- **Day 9**: Integration testing and documentation (~0.5 hour)

**Total Estimated Time**: ~9 hours

---

## Success Criteria

### Phase 5A Complete When:
- [  ] ExperimentWindow detects multi-sweep parameters
- [  ] Multi-sweep summary displayed in window
- [  ] Multi-parameter params_dict built correctly
- [  ] Data stored with full parameter values
- [  ] 1D and multi-sweep modes coexist without conflicts

### Phase 5B Complete When:
- [  ] 2D heatmap displays for 2D sweeps
- [  ] Heatmap updates in real-time during acquisition
- [  ] Axes labeled with parameter names and units
- [  ] Colorbar shows signal scale
- [  ] User can see current acquisition progress on heatmap

### Phase 5C Complete When:
- [  ] Export button in scan selector
- [  ] CSV export works for 1D data
- [  ] NumPy export works for array data
- [  ] YAML export includes metadata
- [  ] Fit parameters exportable to YAML
- [  ] Export dialog intuitive and functional

---

## Dependencies

### Phase 5A:
- No new dependencies (uses existing code)

### Phase 5B:
- pyqtgraph (already used)
- Consider: matplotlib for export plots (optional)

### Phase 5C:
- csv (standard library)
- numpy (already used)
- yaml (already used)

---

## Risk Assessment

### Low Risk:
- Phase 5A: Extends existing architecture, no breaking changes
- Phase 5C: Independent feature, doesn't affect acquisition

### Medium Risk:
- Phase 5B: pyqtgraph ImageView vs PlotItem+ImageItem decision
- Colormap customization may require trial and error
- Real-time heatmap performance for large arrays (>100×100)

### Mitigation:
- Start with pyqtgraph ImageView (simpler)
- Upgrade to PlotItem+ImageItem if needed (more control)
- Test performance with representative data sizes
- Add option to disable real-time updates for large sweeps

---

## Alternative Approaches

### Option 1: Separate 2D Window
Instead of adapting ExperimentWindow:
- Create `Experiment2DWindow` class
- Specialized for 2D/3D data
- Cleaner separation but more code duplication

**Pros**: Simpler logic, no conditionals
**Cons**: More maintenance, less flexible
**Recommendation**: Not preferred, use adaptive window

### Option 2: Use Matplotlib Instead of pyqtgraph
For heatmaps:
- matplotlib's `imshow()` is familiar
- Better for static plots
- Slower updates

**Pros**: More plotting options, better for publications
**Cons**: Slower real-time performance
**Recommendation**: Use pyqtgraph for GUI, matplotlib for exports

### Option 3: Web-Based Visualization
Use Plotly/Dash:
- Interactive 3D plots
- Modern web UI
- Remote access

**Pros**: Best visualization, cross-platform
**Cons**: Complete rewrite, steeper learning curve
**Recommendation**: Future consideration, not Phase 5

---

## Post-Phase 5 Enhancements

### Phase 6 Ideas:
1. **3D Visualization**: Slice viewer for 3D+ sweeps
2. **Adaptive Sweeps**: Zoom into resonance automatically
3. **Real-Time Fitting for 2D**: Fit each row/column during acquisition
4. **Batch Experiments**: Queue multiple experiments
5. **Database Integration**: Store results in SQLite/PostgreSQL
6. **Remote Control**: API for external control
7. **Automated Reports**: PDF generation with plots and analysis

---

## Questions for User

Before starting Phase 5:

1. **Priority**: Which is most urgent?
   - A) Multi-sweep integration (enable 2D hardware sweeps)
   - B) 2D visualization (see 2D data)
   - C) Data export (save results)

2. **2D Sweep Use Case**: What 2D sweeps do you run most?
   - Frequency × Power?
   - Frequency × Magnetic field?
   - Tau × Frequency?

3. **Export Format**: Which format is most important?
   - CSV for Excel/Origin?
   - NumPy for Python analysis?
   - YAML for config + data?
   - HDF5 for large datasets?

4. **Hardware Testing**: Are instruments available for testing Phase 5A?
   - Need working SG, PB, DAQ
   - Test 2D sweep on real hardware

---

## Summary

**Phase 5 Goal**: Complete multi-parameter sweep support with visualization and export

**Three Sub-Phases**:
- **5A**: ExperimentWindow multi-sweep integration (~3 hours)
- **5B**: 2D heatmap visualization (~2 hours)
- **5C**: Data export functionality (~1 hour)

**Total Effort**: ~6 hours development + 3 hours testing = 9 hours

**Current Blockers**: None - ready to start

**Recommendation**: Implement in order 5A → 5B → 5C for logical progression

---

*Created*: 2025-12-17
*Author*: Claude (Sonnet 4.5)
*Project*: NV Experiment GUI
*Status*: Planning complete, ready for implementation
