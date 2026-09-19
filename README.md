# Cisplatin-Resistant-Ovarian-Cancer-Drug-Discovery
Finding the appropriate targets for cisplatin-resistant ovarian cancer patients.

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
# Significance of each step

Stage 1 — Data acquisition & QC
Before trusting any biology, you need to know the data itself is clean. This stage builds the expression matrix and checks for outlier arrays or batch effects (via PCA/correlation heatmap) — skipping it means any downstream signal could just be a scanning artifact rather than real biology.

Stage 2 — DEG analysis
This answers the most basic question: which genes actually differ between resistant and sensitive cell lines? It's your first, broadest candidate list — necessary but not sufficient on its own, because a gene can be statistically different without being biologically central to resistance.

Stage 3 — WGCNA
DEGs alone treat every gene independently and miss coordinated biology. WGCNA groups genes that rise and fall together into co-expression modules, then asks which whole module correlates with resistance (GI50) — this captures pathway-level, systems-biology signal that a simple two-group DEG test can't see, and gives you a second, independent line of evidence.

Stage 4 — Cross-validation in independent GEO cohorts
A single dataset can produce false positives from chance or dataset-specific quirks. Checking whether your DEGs move in the same direction in other, independent experiments is the cheapest, fastest filter against those false positives — before you invest further analysis in a signal that doesn't replicate.

Stage 5 — TCGA-OV survival validation
Everything so far comes from cell lines — useful for mechanism, but a 'biomarker' has to mean something in actual patients. This stage tests whether your surviving genes are prognostic for real patient outcomes, which is what separates a cell-culture finding from a clinically meaningful one.

Stage 6 — Intersection & PPI network
By now you have genes that are differentially expressed, replicated, co-expressed as a module, AND prognostic — intersecting these lists keeps only genes that satisfy every criterion at once. The PPI network then ranks survivors by how central they are in the interaction network, since core/hub proteins are generally more promising drug targets than peripheral ones.

Stage 7 — ML feature selection (LASSO/SVM-RFE/Boruta)
Even after all that filtering, you may still have more candidates than is practical to validate experimentally. Three different machine-learning algorithms, each with different statistical assumptions, independently ask 'which of these genes best predicts resistance?' — genes chosen by all three are far less likely to be an artifact of any single method's bias.

Stage 8 — Risk-score model & ROC validation
This converts your final gene panel into something testable: a quantitative score you can evaluate (ROC/AUC, calibration) and eventually validate in an independent cohort. It's the step that turns 'a list of interesting genes' into 'a candidate diagnostic/prognostic biomarker' ready for wet-lab confirmation.


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

