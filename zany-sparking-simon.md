# Phase 5B Revision: Multi-Sweep Line Plot Display

## User Request Summary
**Change from original Phase 5B plan**:
- NO 2D heatmap visualization
- Use line plots ONLY for all sweep displays
- For 2D sweeps: Show selection box with list of all slices
- User selects which slices to display as line plots
- Heatmap deferred to post-processing (outside GUI)

## Example: 2D Frequency × Power Sweep

### Original Plan (REJECTED):
- Show heatmap (frequency vs power with color indicating signal)
- User sees entire 2D data at once

### New Plan (APPROVED):
- Store 11 frequency sweeps (one per power value: 0, 2, 4, ..., 20 dBm)
- Show selection box listing: "Power = 0 dBm", "Power = 2 dBm", ..., "Power = 20 dBm"
- User checks which power values to display
- Selected sweeps shown as overlaid line plots (frequency vs signal)
- Each sweep has different color (from scan selector color palette)

## Implementation Approach

### Phase 5B: Slice Selector for Multi-Dimensional Sweeps

**Goal**: Display 2D/3D sweep data as selectable 1D line plots

**Key Concept**:
- Inner loop parameter → X-axis of line plot (e.g., frequency)
- Outer loop parameter(s) → Selection list (e.g., power values)
- User selects which outer values to display
- Multiple selections = multiple overlaid line plots

---

## Step 1: Store Slices During Acquisition

**File**: `gui/windows/experiment_window.py`

**For 2D Sweep** (e.g., frequency × power):
```python
def __init__(self, window_id: str, parameters: Dict[str, Any], parent=None):
    # ... existing code ...

    # For multi-sweep: store slices
    if self.use_multi_sweep and self.sweep_dimensions == 2:
        self.sweep_slices = []  # List of slices
        self.current_slice_data = {
            'outer_value': None,
            'x_data': [],
            'y_data': []
        }

def _on_hardware_data_point(self, point_idx: int, param_values: dict, signal: float, metadata: dict):
    """Handle data point - organize into slices for 2D."""

    if self.use_multi_sweep and self.sweep_dimensions == 2:
        inner_param = self.multi_sweep_params[0]
        outer_param = self.multi_sweep_params[1]

        # Get parameter values
        inner_value = param_values[inner_param.name]
        outer_value = param_values[outer_param.name]

        # Check if starting new slice (outer parameter changed)
        if (self.current_slice_data['outer_value'] is None or
            self.current_slice_data['outer_value'] != outer_value):

            # Save previous slice if exists
            if self.current_slice_data['outer_value'] is not None:
                self.sweep_slices.append(self.current_slice_data.copy())

            # Start new slice
            self.current_slice_data = {
                'outer_param': outer_param.name,
                'outer_value': outer_value,
                'x_data': [],
                'y_data': []
            }

        # Add point to current slice
        self.current_slice_data['x_data'].append(inner_value)
        self.current_slice_data['y_data'].append(signal)

        # Update live plot with current slice
        display_x = inner_value * self.config['x_scale']
        self.live_plot.add_data_point(display_x, signal)

    else:
        # 1D sweep: existing behavior
        # ... current code ...

def _on_hardware_finished(self):
    """Handle completion - save final slice."""

    if self.use_multi_sweep and self.sweep_dimensions == 2:
        # Save final slice
        if self.current_slice_data['outer_value'] is not None:
            self.sweep_slices.append(self.current_slice_data.copy())

    super()._on_hardware_finished()
```

**Result**: After 2D acquisition, `self.sweep_slices` contains:
```python
[
    {'outer_param': 'power', 'outer_value': 0, 'x_data': [2.85e9, ...], 'y_data': [2.3, ...]},
    {'outer_param': 'power', 'outer_value': 2, 'x_data': [2.85e9, ...], 'y_data': [2.4, ...]},
    ...
    {'outer_param': 'power', 'outer_value': 20, 'x_data': [2.85e9, ...], 'y_data': [2.7, ...]},
]
```

---

## Step 2: Create Slice Selector Widget

**New File**: `gui/widgets/slice_selector.py` (~200 lines)

**Purpose**: List all slices with checkboxes, allow selection for display

```python
"""
Slice Selector Widget for Multi-Dimensional Sweep Display

For 2D/3D sweeps, allows user to select which slices to display as line plots.

Example for frequency × power sweep:
- Lists all power values: 0 dBm, 2 dBm, 4 dBm, ..., 20 dBm
- User checks which to display
- Selected slices shown as overlaid line plots
"""

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QLabel, QListWidget, QListWidgetItem, QPushButton
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor

class SliceSelectorWidget(QWidget):
    """
    Slice selector for multi-dimensional sweep display.

    Signals:
        selection_changed: Emitted when slice selection changes (selected_indices: List[int])
    """

    selection_changed = Signal(list)  # List of selected slice indices

    def __init__(self, parent=None):
        super().__init__(parent)

        self.slices = []  # List of slice metadata
        self.slice_colors = {}  # slice_index -> QColor

        self._create_ui()

    def _create_ui(self):
        """Create slice selector UI."""
        layout = QVBoxLayout()

        # Header
        header_layout = QHBoxLayout()
        self.header_label = QLabel("Available Slices")
        self.header_label.setStyleSheet("font-weight: bold;")
        header_layout.addWidget(self.header_label)

        self.count_label = QLabel("0 slices")
        header_layout.addWidget(self.count_label)
        header_layout.addStretch()
        layout.addLayout(header_layout)

        # Slice list (checkboxes)
        self.slice_list = QListWidget()
        self.slice_list.setMaximumHeight(200)
        self.slice_list.itemChanged.connect(self._on_selection_changed)
        layout.addWidget(self.slice_list)

        # Control buttons
        btn_layout = QHBoxLayout()

        self.select_all_btn = QPushButton("Select All")
        self.select_all_btn.clicked.connect(self._select_all)
        btn_layout.addWidget(self.select_all_btn)

        self.deselect_all_btn = QPushButton("Deselect All")
        self.deselect_all_btn.clicked.connect(self._deselect_all)
        btn_layout.addWidget(self.deselect_all_btn)

        layout.addLayout(btn_layout)

        self.setLayout(layout)

    def set_slices(self, slices: List[Dict[str, Any]]):
        """
        Set slices to display.

        Args:
            slices: List of slice dicts with 'outer_param', 'outer_value', 'x_data', 'y_data'
        """
        self.slices = slices
        self.slice_list.clear()
        self.slice_colors.clear()

        # Color palette (same as scan selector)
        colors = ["#00ff00", "#ff0000", "#0000ff", "#ffff00", "#ff00ff",
                  "#00ffff", "#ffa500", "#ff1493", "#00ff7f", "#ff4500"]

        for i, slice_data in enumerate(slices):
            # Create label
            outer_param = slice_data['outer_param']
            outer_value = slice_data['outer_value']

            # Format based on value magnitude
            if abs(outer_value) < 0.01 or abs(outer_value) > 1000:
                value_str = f"{outer_value:.3e}"
            else:
                value_str = f"{outer_value:.2f}"

            label = f"{outer_param.capitalize()} = {value_str}"

            # Create list item
            item = QListWidgetItem(label)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked)  # Auto-select all
            item.setData(Qt.UserRole, i)  # Store slice index

            # Assign color
            color = QColor(colors[i % len(colors)])
            self.slice_colors[i] = color
            item.setForeground(color)

            self.slice_list.addItem(item)

        self.count_label.setText(f"{len(slices)} slices")
        self.header_label.setText(f"Available Slices ({slices[0]['outer_param'].capitalize()})")

    def get_selected_indices(self) -> List[int]:
        """Get indices of selected slices."""
        selected = []
        for i in range(self.slice_list.count()):
            item = self.slice_list.item(i)
            if item.checkState() == Qt.Checked:
                selected.append(item.data(Qt.UserRole))
        return selected

    def _select_all(self):
        """Select all slices."""
        for i in range(self.slice_list.count()):
            self.slice_list.item(i).setCheckState(Qt.Checked)

    def _deselect_all(self):
        """Deselect all slices."""
        for i in range(self.slice_list.count()):
            self.slice_list.item(i).setCheckState(Qt.Unchecked)

    def _on_selection_changed(self, item: QListWidgetItem):
        """Handle selection change."""
        selected = self.get_selected_indices()
        self.selection_changed.emit(selected)

    def get_slice_color(self, slice_index: int) -> QColor:
        """Get color for a slice."""
        return self.slice_colors.get(slice_index, QColor("#ffffff"))
```

---

## Step 3: Integrate Slice Selector into ExperimentWindow

**File**: `gui/windows/experiment_window.py`

**Add to Controls Panel**:
```python
def _create_controls_panel(self) -> QWidget:
    widget = super()._create_controls_panel()

    # ... existing sweep builder code ...

    # For 2D+ sweeps: Add slice selector
    if self.use_multi_sweep and self.sweep_dimensions >= 2:
        slice_group = QGroupBox("Slice Selection")
        slice_layout = QVBoxLayout()

        self.slice_selector = SliceSelectorWidget()
        self.slice_selector.selection_changed.connect(self._on_slice_selection_changed)
        slice_layout.addWidget(self.slice_selector)

        slice_group.setLayout(slice_layout)
        self.parameters_layout.addWidget(slice_group)

    return widget
```

**Populate After Acquisition**:
```python
def _on_hardware_finished(self):
    """Handle completion - populate slice selector."""

    if self.use_multi_sweep and self.sweep_dimensions == 2:
        # Save final slice
        if self.current_slice_data['outer_value'] is not None:
            self.sweep_slices.append(self.current_slice_data.copy())

        # Populate slice selector
        self.slice_selector.set_slices(self.sweep_slices)

        # Initial display: show all slices
        self._display_selected_slices(list(range(len(self.sweep_slices))))

    # Call parent for auto-fit, etc.
    super()._on_hardware_finished()
```

---

## Step 4: Display Selected Slices as Line Plots

**File**: `gui/windows/experiment_window.py`

```python
def _on_slice_selection_changed(self, selected_indices: List[int]):
    """Handle slice selection change - update plot."""
    self._display_selected_slices(selected_indices)

def _display_selected_slices(self, selected_indices: List[int]):
    """Display selected slices as overlaid line plots."""

    # Clear plot
    self.live_plot.clear_data()

    if not selected_indices:
        return

    # Plot each selected slice
    for slice_idx in selected_indices:
        slice_data = self.sweep_slices[slice_idx]
        color = self.slice_selector.get_slice_color(slice_idx)

        x_data = np.array(slice_data['x_data'])
        y_data = np.array(slice_data['y_data'])

        # Apply display scaling
        x_display = x_data * self.config['x_scale']

        # Get label for legend
        outer_param = slice_data['outer_param']
        outer_value = slice_data['outer_value']
        label = f"{outer_param}={outer_value:.2f}"

        # Plot line with color
        self.live_plot.plot_widget.plot(
            x_display,
            y_data,
            pen={'color': color.name(), 'width': 2},
            symbol='o',
            symbolSize=4,
            symbolBrush=color,
            name=label
        )

    # Add legend if multiple slices
    if len(selected_indices) > 1:
        self.live_plot.plot_widget.addLegend()
```

---

## Step 5: Memory Management (Optional Enhancement)

**For Large Datasets** (e.g., 100 × 100 points = 10,000 values):

**Approach**: Store slices in memory by default, add option to save to disk if needed

```python
def __init__(self, window_id: str, parameters: Dict[str, Any], parent=None):
    # ... existing code ...

    # Memory management
    self.max_memory_slices = 100  # Keep up to 100 slices in memory
    self.use_disk_cache = False
    self.disk_cache_dir = None

def _on_hardware_finished(self):
    """Handle completion - check if disk cache needed."""

    if len(self.sweep_slices) > self.max_memory_slices:
        # Prompt user to save to disk
        reply = QMessageBox.question(
            self,
            "Large Dataset",
            f"Sweep has {len(self.sweep_slices)} slices. Save to disk for better performance?",
            QMessageBox.Yes | QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            self._save_slices_to_disk()
            self.use_disk_cache = True

    # Continue normal flow
    # ...

def _save_slices_to_disk(self):
    """Save slices to temporary directory."""
    import tempfile
    import pickle

    self.disk_cache_dir = tempfile.mkdtemp(prefix="nv_sweep_")

    for i, slice_data in enumerate(self.sweep_slices):
        file_path = os.path.join(self.disk_cache_dir, f"slice_{i}.pkl")
        with open(file_path, 'wb') as f:
            pickle.dump(slice_data, f)

    # Clear from memory, keep only metadata
    self.sweep_slices_metadata = [
        {'outer_param': s['outer_param'], 'outer_value': s['outer_value'], 'index': i}
        for i, s in enumerate(self.sweep_slices)
    ]
    self.sweep_slices = []  # Free memory

def _load_slice_from_disk(self, slice_idx: int) -> Dict:
    """Load slice from disk cache."""
    import pickle

    file_path = os.path.join(self.disk_cache_dir, f"slice_{slice_idx}.pkl")
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def _display_selected_slices(self, selected_indices: List[int]):
    """Display selected slices - load from disk if needed."""

    self.live_plot.clear_data()

    for slice_idx in selected_indices:
        # Load slice (from memory or disk)
        if self.use_disk_cache:
            slice_data = self._load_slice_from_disk(slice_idx)
        else:
            slice_data = self.sweep_slices[slice_idx]

        # Plot as before
        # ...
```

---

## Files to Modify/Create

### New File:
1. `gui/widgets/slice_selector.py` (~200 lines)

### Modified Files:
1. `gui/windows/experiment_window.py`:
   - Add `sweep_slices` storage (~30 lines)
   - Integrate `SliceSelectorWidget` (~20 lines)
   - Implement `_display_selected_slices()` (~40 lines)
   - Optional: disk cache for large datasets (~60 lines)

**Total**: ~350 lines across 2 files

---

## Advantages of This Approach

✅ **Familiar**: Users already understand line plots
✅ **Comparative**: Easy to compare different parameter values (overlay multiple lines)
✅ **Flexible**: Select any combination of slices to display
✅ **Simple**: No need to implement heatmap, colorbar, axis transformations
✅ **Compatible**: Works with existing LivePlotWidget infrastructure
✅ **Exportable**: Each slice can be exported separately or together
✅ **Post-processing**: User can create heatmaps in Origin/matplotlib/etc. with exported data

---

## Example User Workflow

### 2D ESR Power Sweep:

1. **Configure in Config Tab**:
   - Frequency: 2.85-2.89 GHz, 51 points
   - Power: 0-20 dBm, 11 points
   - Total: 561 points

2. **Open Window**: ESR window shows "2D sweep configuration"

3. **Run Acquisition**:
   - Live plot shows current frequency sweep
   - Each time power changes, plot clears and shows new sweep
   - Status: "Power = 4 dBm, acquiring..." (updates as outer loop progresses)

4. **After Completion**:
   - Slice Selector appears with 11 checkboxes:
     ```
     ✓ Power = 0.00 dBm    (green)
     ✓ Power = 2.00 dBm    (red)
     ✓ Power = 4.00 dBm    (blue)
     ...
     ✓ Power = 20.00 dBm   (orange)
     ```
   - Plot shows all 11 sweeps overlaid with different colors

5. **Compare Specific Powers**:
   - User unchecks all
   - User checks only "Power = 0 dBm" and "Power = 10 dBm"
   - Plot updates to show only those two sweeps
   - Easy to see how resonance changes with power

6. **Export**:
   - User can export each slice individually
   - Or export all slices together as multi-column CSV

---

## Testing Plan

### Test 1: 1D Sweep (Regression)
- Configure 1D frequency sweep
- Run acquisition
- Verify: No slice selector appears (1D mode)
- Verify: Standard line plot works as before

### Test 2: 2D Sweep Display
- Configure 2D frequency × power sweep (11 powers)
- Run acquisition
- Verify: Slice selector appears with 11 items
- Verify: All 11 sweeps displayed with different colors
- Verify: Legend shows power values

### Test 3: Slice Selection
- Deselect all slices → Plot clears
- Select single slice → Single line appears
- Select 3 slices → Three colored lines appear
- Verify: Colors match list items

### Test 4: Large Dataset (Optional)
- Configure 2D sweep with 100 outer loop values
- Verify: Prompt for disk cache appears
- Accept disk cache
- Verify: Slices load on-demand when selected
- Verify: No memory issues

### Test 5: Export
- Select 2 slices
- Export to CSV
- Verify: Both slices in file
- Verify: Can plot in Excel/Origin

---

## Success Criteria

Phase 5B complete when:
- [  ] SliceSelectorWidget created and functional
- [  ] 2D sweep data organized into slices during acquisition
- [  ] Slice selector populated after 2D acquisition completes
- [  ] Multiple slices display as overlaid line plots
- [  ] Color-coded slices with legend
- [  ] Select/deselect individual slices updates plot
- [  ] 1D sweeps still work as before (no regression)
- [  ] Memory management works for large datasets (optional)

---

## Timeline

**Phase 5B Revised**: ~2-3 hours
- Hour 1: Create SliceSelectorWidget (~200 lines)
- Hour 2: Integrate into ExperimentWindow, slice storage (~90 lines)
- Hour 3: Testing and refinement
- Optional: Memory management (+1 hour if needed for large datasets)

**Much simpler than original heatmap plan!**

---

## Summary

**Changed Approach**:
- ❌ 2D heatmap with colorbar
- ✅ Line plot display with slice selector

**Benefits**:
- Simpler implementation
- More familiar to users
- Better for comparison
- Easier to export
- No new dependencies

**User can still create heatmaps later**:
- Export 2D data to CSV/NPY
- Load in Origin/matplotlib/Jupyter
- Create publication-quality heatmaps with full control

---

*Plan Status*: Ready for review and approval
*Estimated Effort*: 2-3 hours implementation + 1 hour testing
*Risk Level*: Low (extends existing architecture)
