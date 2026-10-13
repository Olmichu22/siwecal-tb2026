#
# Gaudi/k4FWCore steering: siwecaldecoded ROOT tree(s) -> pedestal or MIP
# calibration table, via the PedestalMipCalibrator component.
#
#   CALIB_INPUT_FILES=/path/run.root \
#   CALIB_MODE=Pedestal CALIB_GAIN=high \
#   CALIB_OUTPUT_PEDESTAL_FILE=/path/Pedestal_run_highgain.txt \
#       k4run gaudi_source/options/run_pedestal_mip.py
#
#   CALIB_INPUT_FILES=/path/run1.root,/path/run2.root \
#   CALIB_MODE=Mip CALIB_GAIN=high \
#   CALIB_OUTPUT_MIP_FILE=/path/MIP_pedestalsubmode1_run_highgain.txt \
#       k4run gaudi_source/options/run_pedestal_mip.py
#
# Mip mode computes its own on-the-fly per-SCA pedestal internally (the
# ported method's pedestal_mode=1 -- see PedestalMipCalibrator.cpp) and does
# NOT read a separate pedestal file.
#
# CALIB_DIAGNOSTICS_FILE (optional): cross-check ROOT file with 2D
# layer/chip vs. sca/channel maps + a fit-status code map + 1D value
# distributions (see PedestalMipCalibrator.cpp's DiagnosticsFile doc).
# Empty (default) = no plots written.
#
# Condor fill/merge/fit pipeline modes (see gaudi_jobs/calibration/condor/):
#   CALIB_MODE=Fill: builds the ped_high/ped_low/mip_high/mip_low histogram
#     grids from CALIB_INPUT_FILES (no fit) and writes them, zero-suppressed,
#     to CALIB_OUTPUT_HISTOGRAM_FILE.
#   CALIB_MODE=FitPedestal|FitMip: reads CALIB_INPUT_HISTOGRAM_FILE (a
#     hadd-merged Fill output) instead of CALIB_INPUT_FILES, and runs the
#     same fit/robustness/diagnostics code as Mode=Pedestal|Mip.
#   CALIB_MODE=FitTrigger: the trigger turn-on per channel/chip/slab from the
#     same merged Fill output, with the set's high-gain pedestal and MIP tables
#     (CALIB_INPUT_PEDESTAL_FILE, CALIB_INPUT_MIP_FILE), written to
#     CALIB_OUTPUT_TRIGGER_{CHIP,CHANNEL,SLAB}_FILE (see TriggerTurnOn.h):
#
#   CALIB_MODE=FitTrigger CALIB_INPUT_HISTOGRAM_FILE=merged_th210.root \
#   CALIB_INPUT_PEDESTAL_FILE=.../pedestals/th210/Pedestal_..._highgain.txt \
#   CALIB_INPUT_MIP_FILE=.../mips/th210/MIP_..._highgain.txt \
#   CALIB_OUTPUT_TRIGGER_CHIP_FILE=.../trigger/th210/Trigger_chip_highgain.txt \
#   CALIB_OUTPUT_TRIGGER_CHANNEL_FILE=.../trigger/th210/Trigger_channel_highgain.txt \
#       k4run gaudi_source/options/run_pedestal_mip.py
#
# All work happens in PedestalMipCalibrator::initialize(); EvtMax=1 is enough
# since execute() is a no-op (this is a global reduction over an entire run,
# not a per-event k4FWCore transform -- see PedestalMipCalibrator.cpp).
import os

import Gaudi.Configuration as _gaudi_cfg
from Configurables import EventDataSvc
from Configurables import PedestalMipCalibrator
from k4FWCore import ApplicationMgr

# Log verbosity: WARNING by default (Condor fill/fit jobs are silent unless
# something is wrong -- keeps per-job .out logs tiny across thousands of jobs);
# override with CALIB_OUTPUT_LEVEL=INFO/DEBUG for interactive debugging.
_output_level = getattr(_gaudi_cfg, os.environ.get("CALIB_OUTPUT_LEVEL", "WARNING"), _gaudi_cfg.WARNING)

mode = os.environ.get("CALIB_MODE", "Pedestal")
if mode not in ("Pedestal", "Mip", "Fill", "FitPedestal", "FitMip", "FitTrigger"):
    raise SystemExit(f"CALIB_MODE must be one of Pedestal|Mip|Fill|FitPedestal|FitMip|FitTrigger, got: {mode}")

input_files_raw = os.environ.get("CALIB_INPUT_FILES", "")
if mode in ("Pedestal", "Mip", "Fill") and not input_files_raw:
    raise SystemExit(f"CALIB_MODE={mode} requires CALIB_INPUT_FILES (comma-separated siwecaldecoded ROOT files)")
input_files = [f for f in input_files_raw.split(",") if f.strip()]

_default_slab_z_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                                    "mappings", "slab_z_positions.yml")
if not os.path.isfile(_default_slab_z_file):
    _default_slab_z_file = ""

calib = PedestalMipCalibrator(
    "PedestalMipCalibrator",
    InputFiles=input_files,
    TreeName=os.environ.get("CALIB_TREE", "siwecaldecoded"),
    Mode=mode,
    Gain=os.environ.get("CALIB_GAIN", "high"),
    OutputPedestalFile=os.environ.get("CALIB_OUTPUT_PEDESTAL_FILE", ""),
    OutputMipFile=os.environ.get("CALIB_OUTPUT_MIP_FILE", ""),
    MaxNhit=int(os.environ.get("CALIB_MAX_NHIT", "1")),
    NSlabsHit=int(os.environ.get("CALIB_NSLABS_HIT", "8")),
    SpillSelect=os.environ.get("CALIB_SPILL_SELECT", "all"),
    SpillOccupancyCut=int(os.environ.get("CALIB_SPILL_OCC_CUT", "200")),
    MinMipIntegral=float(os.environ.get("CALIB_MIN_MIP_INTEGRAL", "200")),
    ChipFallbackMinIntegral=float(os.environ.get("CALIB_CHIP_FALLBACK_MIN_INTEGRAL", "2000")),
    PedestalMaxMeanAdc=float(os.environ.get("CALIB_PEDESTAL_MAX_MEAN_ADC", "300")),
    MipFallbackMaxAdc=float(os.environ.get("CALIB_MIP_FALLBACK_MAX_ADC", "300")),
    # Fit quality. There is no chi2/ndf gate any more: chi2/ndf is not comparable
    # across statistics (it grew 1.2 -> 17 on the same spectrum as it went from
    # 639 to 74,555 entries, with the MPV unmoved), so the old ceilings rejected
    # every fit above ~600 entries and every chip/slab fit. These two are
    # scale-free. See PedestalMipCalibrator.cpp's fit-quality property block.
    MipHistMaxTol=float(os.environ.get("CALIB_MIP_HISTMAX_TOL", "0.20")),
    MipPeakRebin=int(os.environ.get("CALIB_MIP_PEAK_REBIN", "4")),
    MaxMipRelErr=float(os.environ.get("CALIB_MAX_MIP_REL_ERR", "0.10")),
    MipLowLim=float(os.environ.get("CALIB_MIP_LOW_LIM", "10")),
    MipHighLim=float(os.environ.get("CALIB_MIP_HIGH_LIM", "50")),
    MipLowLimLowGain=float(os.environ.get("CALIB_MIP_LOW_LIM_LOW_GAIN", "2")),
    MipHighLimLowGain=float(os.environ.get("CALIB_MIP_HIGH_LIM_LOW_GAIN", "20")),
    # Per-slab MIP windows.  The slab file's technology block marks the slabs
    # that ran at a discriminator of their own (slab 12, the chip-on-board, at
    # DAC 243); their truncated spectra fit at 31-35 ADC and need a wider window.
    # CALIB_SLAB_Z_FILE="" disables it; CALIB_MIP_WINDOW_OVERRIDES="12:10:90,..."
    # sets explicit windows.
    SlabZFile=os.environ.get("CALIB_SLAB_Z_FILE", _default_slab_z_file),
    MipWindowThresholdSlabHigh=float(os.environ.get("CALIB_MIP_WINDOW_THRESHOLD_SLAB_HIGH", "90")),
    MipWindowSlabOverrides=[e for e in os.environ.get("CALIB_MIP_WINDOW_OVERRIDES", "").split(",") if e.strip()],
    DiagnosticsFile=os.environ.get("CALIB_DIAGNOSTICS_FILE", ""),
    OutputHistogramFile=os.environ.get("CALIB_OUTPUT_HISTOGRAM_FILE", ""),
    InputHistogramFile=os.environ.get("CALIB_INPUT_HISTOGRAM_FILE", ""),
    InputPedestalFile=os.environ.get("CALIB_INPUT_PEDESTAL_FILE", ""),
    InputMipFile=os.environ.get("CALIB_INPUT_MIP_FILE", ""),
    OutputTriggerChipFile=os.environ.get("CALIB_OUTPUT_TRIGGER_CHIP_FILE", ""),
    OutputTriggerChannelFile=os.environ.get("CALIB_OUTPUT_TRIGGER_CHANNEL_FILE", ""),
    OutputTriggerSlabFile=os.environ.get("CALIB_OUTPUT_TRIGGER_SLAB_FILE", ""),
    TriggerMinPassChip=float(os.environ.get("CALIB_TRIGGER_MIN_PASS_CHIP", "300")),
    TriggerMinPassChannel=float(os.environ.get("CALIB_TRIGGER_MIN_PASS_CHANNEL", "30")),
    TriggerMinTurnOnChannel=float(os.environ.get("CALIB_TRIGGER_MIN_TURNON_CHANNEL", "30")),
)

ApplicationMgr(TopAlg=[calib],
               EvtSel="NONE",
               EvtMax=1,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=_output_level)
