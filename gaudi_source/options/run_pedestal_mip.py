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
if mode not in ("Pedestal", "Mip", "Fill", "FitPedestal", "FitMip"):
    raise SystemExit(f"CALIB_MODE must be one of Pedestal|Mip|Fill|FitPedestal|FitMip, got: {mode}")

input_files_raw = os.environ.get("CALIB_INPUT_FILES", "")
if mode in ("Pedestal", "Mip", "Fill") and not input_files_raw:
    raise SystemExit(f"CALIB_MODE={mode} requires CALIB_INPUT_FILES (comma-separated siwecaldecoded ROOT files)")
input_files = [f for f in input_files_raw.split(",") if f.strip()]

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
    MinMipIntegral=float(os.environ.get("CALIB_MIN_MIP_INTEGRAL", "200")),
    ChipFallbackMinIntegral=float(os.environ.get("CALIB_CHIP_FALLBACK_MIN_INTEGRAL", "2000")),
    PedestalMaxMeanAdc=float(os.environ.get("CALIB_PEDESTAL_MAX_MEAN_ADC", "300")),
    MipFallbackMaxAdc=float(os.environ.get("CALIB_MIP_FALLBACK_MAX_ADC", "300")),
    MaxMipChi2Ndf=float(os.environ.get("CALIB_MAX_MIP_CHI2NDF", "1.1")),
    MaxMipChi2NdfFallback=float(os.environ.get("CALIB_MAX_MIP_CHI2NDF_FALLBACK", "25")),
    MipHistMaxTol=float(os.environ.get("CALIB_MIP_HISTMAX_TOL", "0.25")),
    MipLowLim=float(os.environ.get("CALIB_MIP_LOW_LIM", "10")),
    MipHighLim=float(os.environ.get("CALIB_MIP_HIGH_LIM", "50")),
    MipLowLimLowGain=float(os.environ.get("CALIB_MIP_LOW_LIM_LOW_GAIN", "2")),
    MipHighLimLowGain=float(os.environ.get("CALIB_MIP_HIGH_LIM_LOW_GAIN", "20")),
    DiagnosticsFile=os.environ.get("CALIB_DIAGNOSTICS_FILE", ""),
    OutputHistogramFile=os.environ.get("CALIB_OUTPUT_HISTOGRAM_FILE", ""),
    InputHistogramFile=os.environ.get("CALIB_INPUT_HISTOGRAM_FILE", ""),
)

ApplicationMgr(TopAlg=[calib],
               EvtSel="NONE",
               EvtMax=1,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=_output_level)
