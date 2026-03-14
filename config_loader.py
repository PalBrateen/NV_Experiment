# config_loader.py
"""
Config loader that converts YAML to params_dict format.
Also generates sweep values from start/stop/step.
"""

import yaml
import numpy as np
from pathlib import Path
from typing import Dict, Any, List
from collections import OrderedDict


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load YAML config and convert to experiment params_dict.
    
    Args:
        config_path: Path to .yaml config file
    
    Returns:
        params_dict compatible with experiment system
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return process_config(config)


def process_config(config: Dict) -> Dict[str, Any]:
    """
    Process raw YAML config into params_dict format.
    
    - Converts units to SI base units
    - Generates sweep arrays from start/stop/step
    - Builds structured params_dict
    """
    params = {}
    
    # --------------------------------------------
    # Experiment metadata
    # --------------------------------------------
    params['experiment'] = {
        'name': config['experiment']['name'],
        'description': config['experiment'].get('description', ''),
        'sequence': config['experiment']['sequence'],
        'mode': config['experiment'].get('mode', 'diode')
    }
    
    # --------------------------------------------
    # Sequence parameters (for PulseBlaster)
    # --------------------------------------------
    pb_config = config.get('pulse_blaster', {})
    params['seq'] = {
        'sequence': config['experiment']['sequence'],
        't_AOM': pb_config.get('t_AOM_us', 100) * 1e-6,  # Convert to seconds
        'ro_delay': pb_config.get('ro_delay_ns', 0) * 1e-9,  # Convert to seconds
        'AOM_lag': pb_config.get('AOM_lag_ns', 0) * 1e-9,  # Convert to seconds
        'MW_lag': pb_config.get('MW_lag_ns', 0) * 1e-9,  # Convert to seconds
        'Nsamples': config.get('daq', {}).get('Nsamples', 1),
        'PBchannels': {},  # To be filled by connection config or user
    }

    # PB parameters
    params['pb'] = {
        'clk_cyc': 1e3 / pb_config.get('clk_MHz', 100),  # in ns
    }
    
    params['pb'] = {
        'channels': pb_config.get('channels', {})
    }

    # --------------------------------------------
    # Signal generator parameters
    # --------------------------------------------
    sg_config = config.get('signal_generator', {})
    params['mw'] = {
        'power': sg_config.get('power_dBm', 0),
        'freq': sg_config.get('frequency_Hz', 2.87e9),
        'modulation': sg_config.get('modulation', 'pulse')
    }
    
    # --------------------------------------------
    # DAQ parameters
    # --------------------------------------------
    daq_config = config.get('daq', {})
    params['daq'] = {
        'Nsamples': daq_config.get('Nsamples', 1),
        'sampling_rate': daq_config.get('sampling_rate_Hz', 2e6),
        'voltage_range': tuple(daq_config.get('voltage_range_V', [-10, 10]))
    }
    
    # --------------------------------------------
    # Camera parameters (if present)
    # --------------------------------------------
    if 'camera' in config:
        cam_config = config['camera']
        params['camera'] = {
            'exposure_ms': cam_config.get('exposure_ms', 100),
            'trigger_mode': cam_config.get('trigger_mode', 'level'),
            'roi': cam_config.get('roi', [0, 0, 512, 512])
        }
    
    # --------------------------------------------
    # Acquisition settings
    # --------------------------------------------
    acq_config = config.get('acquisition', {})
    params['scan'] = {
        'Nruns': acq_config.get('Nruns', 1),
        'contrast_mode': acq_config.get('contrast_mode', 'ratio_signal_over_reference'),
        'randomize': acq_config.get('randomize', False)
    }
    
    # --------------------------------------------
    # Sweep configuration (ORDER PRESERVED!)
    # --------------------------------------------
    params['sweep'] = OrderedDict()
    
    for sweep_item in config.get('sweep', []):
        param_name = sweep_item['parameter']
        
        if 'values' in sweep_item:
            # Explicit values provided
            values = np.array(sweep_item['values'])
        else:
            # Generate from start/stop/step
            start = float(sweep_item['start'])
            stop = float(sweep_item['stop'])
            step = float(sweep_item['step'])
            values = np.arange(start, stop + step/2, step)  # +step/2 to include endpoint
        
        params['sweep'][param_name] = values
    
    # Also store sweep metadata for reference
    params['sweep_config'] = config.get('sweep', [])
    
    # --------------------------------------------
    # Output settings
    # --------------------------------------------
    output_config = config.get('output', {})
    params['output'] = {
        'save_path': output_config.get('save_path', 'Saved_Data/'),
        'filename_prefix': output_config.get('filename_prefix', 'data'),
        'auto_save_interval': output_config.get('auto_save_interval', 0)
    }
    
    params['plot'] = output_config.get('plot', {
        'x_label': 'Parameter',
        'x_scale': 1.0,
        'y_label': 'Signal',
        'live_update': True
    })
    
    return params


def save_config(params_dict: Dict, output_path: str):
    """Save params_dict back to YAML format."""
    # Convert numpy arrays to lists for YAML serialization
    def convert_for_yaml(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_for_yaml(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [convert_for_yaml(item) for item in obj]
        else:
            return obj
    
    yaml_dict = convert_for_yaml(params_dict)
    
    with open(output_path, 'w') as f:
        yaml.dump(yaml_dict, f, default_flow_style=False, sort_keys=False)


def flatten_params(params_dict: Dict) -> Dict[str, Any]:
    """
    Flatten nested params_dict for saving to simple format.
    Preserves category as prefix: 'mw.power', 'seq.t_AOM', etc.
    """
    flat = {}
    
    def _flatten(d, prefix=''):
        for key, value in d.items():
            new_key = f"{prefix}.{key}" if prefix else key
            if isinstance(value, dict) and not isinstance(value, OrderedDict):
                _flatten(value, new_key)
            elif isinstance(value, np.ndarray):
                flat[new_key] = value.tolist()
            else:
                flat[new_key] = value
    
    _flatten(params_dict)
    return flat


# ============================================================
# Convenience functions for common configs
# ============================================================

def create_esr_config(
    freq_start: float = 2.82e9,
    freq_stop: float = 2.92e9,
    freq_step: float = 1e6,
    power_dBm: float = 8,
    Nsamples: int = 8,
    Nruns: int = 30,
    t_AOM_us: float = 8000
) -> Dict:
    """Create ESR config programmatically."""
    return {
        'experiment': {'name': 'ESR', 'sequence': 'esr_seq', 'mode': 'diode'},
        'signal_generator': {'power_dBm': power_dBm, 'modulation': 'pulse'},
        'pulse_blaster': {'t_AOM_us': t_AOM_us, 'clk_MHz': 100},
        'daq': {'Nsamples': Nsamples, 'sampling_rate_Hz': 2e6},
        'sweep': [{'parameter': 'frequency', 'instrument': 'signal_generator',
                   'start': freq_start, 'stop': freq_stop, 'step': freq_step, 'units': 'Hz'}],
        'acquisition': {'Nruns': Nruns, 'contrast_mode': 'ratio_signal_over_reference'},
        'output': {'filename_prefix': 'ESR'}
    }


def create_rabi_config(
    duration_start: float = 10e-9,
    duration_stop: float = 500e-9,
    duration_step: float = 2e-9,
    mw_freq: float = 2.87e9,
    power_dBm: float = 8,
    Nsamples: int = 2,
    Nruns: int = 1
) -> Dict:
    """Create Rabi config programmatically."""
    return {
        'experiment': {'name': 'Rabi', 'sequence': 'rabi_seq', 'mode': 'diode'},
        'signal_generator': {'power_dBm': power_dBm, 'frequency_Hz': mw_freq, 'modulation': 'pulse'},
        'pulse_blaster': {'t_AOM_us': 20, 'ro_delay_ns': 1800, 'AOM_lag_ns': 800, 'MW_lag_ns': 150, 'clk_MHz': 100},
        'daq': {'Nsamples': Nsamples, 'sampling_rate_Hz': 2e6},
        'sweep': [{'parameter': 'pulse_duration', 'instrument': 'pulse_blaster',
                   'start': duration_start, 'stop': duration_stop, 'step': duration_step, 'units': 's'}],
        'acquisition': {'Nruns': Nruns, 'contrast_mode': 'ratio_signal_over_reference'},
        'output': {'filename_prefix': 'Rabi'}
    }