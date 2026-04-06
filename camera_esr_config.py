# camera_esr_config.py — Camera ESR experiment config
"""
Widefield NV ESR using Hamamatsu ORCA with level-triggered acquisition.

PulseBlaster drives the camera trigger, laser, and MW channels.
Camera acquires 2 frames per scan point (signal + reference).
Frequency sweep via SG384.

Usage:
    from camera_esr_config import config
    exp, data = run(instruments, config)
"""

from experiment_config import (
    ExperimentConfig, ScanAxis, CameraConfig, FieldConfig,
    MicrowaveConfig, PulseBlasterConfig, TimeConfig,
    SequenceConfig, PlotConfig, SaveConfig, RuntimeFlags,
    ns, us, ms, Hz, kHz, MHz, GHz,
)
from connectionConfig import (PBclk, laser, samp_clk, start_trig,
                               MW, camera, bx, by, bz)

t_AOM = 70 * ms
channels = {
    'samp': samp_clk, 'mw': MW, 'laser': laser,
    'start': start_trig, 'camera': camera,
    'bx': bx, 'by': by, 'bz': bz,
}

config = ExperimentConfig(
    scans={
        'freq': ScanAxis(
            name='freq',
            start=2.82 * GHz,
            stop=2.92 * GHz,
            step=1 * MHz,
        ),
    },

    mw=MicrowaveConfig(power=-24, freq=2.87 * GHz),
    times=TimeConfig(t_AOM=t_AOM),
    pb=PulseBlasterConfig(
        clock_MHz=PBclk,
        channels=dict(sorted(channels.items(), key=lambda kv: kv[1])),
    ),
    seq=SequenceConfig(
        name='esr_seq',
        args_names=['t_AOM'],
        args_values=[t_AOM],
        Nsamples=2,             # 2 frames per scan point (signal + reference)
    ),

    camera=CameraConfig(
        instr_mode='cam_levelm',
        exposure_s=0.020,       # 20 ms
        roi=[1020, 1124, 236, 240],
        frames_per_cycle=2,
        focus_interval=1,       # focus every run
        min_focus_time_s=5.0,
    ),
    field_cfg=FieldConfig(
        align_field=[0, 0, 0],
        test_field=[0, 0, 0],
        t_align_dc_ms=200,
    ),

    plot=PlotConfig(x_units=GHz, x_label='Frequency (GHz)'),
    save_opts=SaveConfig(prefix='CamESR'),
    runtime=RuntimeFlags(Nruns=10, detector='camera'),
).use('mw', 'camera')
