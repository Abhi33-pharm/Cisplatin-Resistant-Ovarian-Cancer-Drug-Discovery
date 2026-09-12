# Cisplatin-Resistant-Ovarian-Cancer-Drug-Discovery
Finding the appropriate targets for cisplatin resistant ovarian cancer patients.

# Cisplatin Resistance Biomarker/Target Pipeline — 8 Stages

Run these in order, from the same working directory. Each stage reads what
the previous one saved to `pipeline_outputs/`, and writes new files there —
you never need to keep variables in memory between scripts.

```
pipeline/
├── stage1_data_acquisition.py       # download GSE47856, QC plots
├── stage2_deg_analysis.py           # DEGs, volcano plot, heatmap, boxplots
├── stage3_wgcna.py                  # WGCNA, soft-power/dendrogram/module-trait plots
├── stage4_cross_validation.py       # validate DEGs in independent GEO cohorts
├── stage5_tcga_survival.py          # TCGA-OV survival, KM curves, forest plot
├── stage6_intersection_ppi.py       # Venn + PPI network + hub ranking
├── stage7_ml_feature_selection.py   # LASSO + SVM-RFE + Boruta/RF, all diagnostic plots
├── stage8_risk_score_validation.py  # risk score, ROC, calibration
└── README.md                        # this file
```

## Before you start

1. `pip install -r requirements.txt` (see bottom of this file)
2. Download **Supplementary Table 1** from Miow et al. 2014 (Oncogene),
   save the cell_line + GI50 columns as `miow2014_supplementary_table1.csv`
   in your working directory.
3. Run each script **one at a time**, in order. Every script prints
   "NEXT: run stage*.py" at the end.

## Where everything lands

- Tables  : `"C:\Users\naemi\Documents\Cis-platin\Cisplatin-Resistant-Ovarian-Cancer-Drug-Discovery\Results\Tables"`
- All plots: `"C:\Users\naemi\Documents\Cis-platin\Cisplatin-Resistant-Ovarian-Cancer-Drug-Discovery\Results\Figures"`

## requirements.txt

```
geoparse
pandas
numpy
scipy
statsmodels
matplotlib
seaborn
scikit-learn
PyWGCNA
networkx
matplotlib-venn
boruta
lifelines
requests
mygene
```

