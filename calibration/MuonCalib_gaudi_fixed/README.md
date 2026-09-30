# MuonCalib_gaudi_fixed

Pedestal, MIP and LG→HG anchor tables re-derived on the chunks re-decoded with the SCA data blocks paired to their
own BCIDs (`Data/rundata_converted_fixed`, fix 48d903b). Same layout as `MuonCalib_gaudi`, so
`generate_reco_dag.py --calib-dir calibration/MuonCalib_gaudi_fixed` reconstructs with it.

| set   | input runs (muons, normal incidence)                        | produced with |
|-------|-------------------------------------------------------------|---------------|
| th210 | eudaq 152–155, 158–160                                      | `run_calibration_batch.py --no-raw2root` |
| th220 | 60–64, 85–89, 142, 143, eudaq 145, 146 (eudaq 144 undecodable) | condor DAG, `--decoded-as-found` |
| th230 | 4                                                           | condor DAG |

Anchors (`anchor/th<N>/gain_anchor_th<N>.txt`), all with the fixed th210 pedestals, as `Reconstructed_final` uses
them: th210 eudaq 166 k 0.0924 c 2.35 (legacy 0.0925, 2.24); th220 run 72 0.0964, 1.81; th230 run 12 0.0958, 1.62.
The 4 % between th210 and th220/th230 is a property of the data-taking period, not of the fix.

## Comparison with the legacy tables (`compare/`)

Produced by `diagnostics/compare_pedestal_tables.py`, `compare_mip_tables.py` and `compare_thresholds.py`.

| high gain                          | th210 legacy → fixed | th220 legacy → fixed | th230 legacy → fixed |
|------------------------------------|----------------------|----------------------|----------------------|
| pedestal width, median [ADC]       | 1.73 → 1.37          | 2.36 → 1.42          | 2.15 → 1.44          |
| RMS of per-SCA offsets [ADC]       | 1.25 → 2.95          | 1.76 → 3.12          | 1.58 → 2.90          |
| pedestal mean new − old, median    | −0.06 (RMS 4.0)      | —                    | 0.01 (RMS 3.8)       |
| MIP MPV, median [ADC]              | 18.39 → 18.69        | 21.68 → 22.69        | 32.17 → 32.98        |
| MIP Landau width, median [ADC]     | 2.10 → 1.94          | 3.17 → 2.95          | 6.06 → 6.11          |
| channels with a MIP fit            | 12 290 → 13 436      | 14 611 → 14 518      | 5 871 → 3 452        |

The per-slab pedestal mean does not move; the real per-SCA structure appears and the widths shrink, because each
SCA's histogram now holds only its own physical cell. The three sets agree with each other on the fixed tables
(pedestal difference per cell th220 − th210: median +0.13, 68 % within 0.49 ADC; legacy 1.46 ADC). The pedestal
"multi-peaks" (20 % of cells) were the same mixing: 0.2 % on the fixed data (run 88).

The expected pedestal agreement of `ordenes_opus.txt` (|Δped| < 0.5 ADC) holds per slab, not per cell: per cell the
legacy tables were wrong by ~4 ADC RMS, in the direction of mixing cells.

Tool defects found while producing them: the Fill accepted ROOT-recovered truncated group files on retry (fixed,
55e540d; 78 of 290 th220 groups had to be redone), and the DAG's merge check rejects valid merged files this large.
