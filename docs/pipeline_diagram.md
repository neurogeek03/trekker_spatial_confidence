# Trekker Spatial Confidence Pipeline

```mermaid
flowchart TD
    %% Inputs
    DATA["📁 Trekker output dirs\n../data/{SAMPLE}/"]
    H5AD["📄 Annotated h5ad\nPCT_test_QC_merged_filtered_…h5ad"]
    CONF["⚙️ config/pipeline.conf\nSamples: BC3 BC9 BC13 BC14 BC15 BC28"]

    %% Orchestrator
    CONF --> SHELL["run_pipeline.sh\n--stage all | scoring | annotations | validation | visualization"]

    %% Stage 1
    SHELL -->|"per sample (loop)"| S1["Stage 1 · Scoring\nrun_scoring.py"]
    DATA --> S1
    S1 -->|"spatial_confidence.io\nload_trekker_data"| S1A["Load coords\n& cell groups"]
    S1A -->|"scoring.score_all_cells\nDBSCAN cascade minPts=4,3,2"| S1B["Score all cells\n(rescue unpositioned)"]
    S1B -->|"scoring.compute_composite_score"| S1C["Composite score\n+ penalties"]
    S1C -->|"scoring.add_trekker_status"| S1D["Merge Trekker status\nconfident / ambiguous / unpositioned"]
    S1D --> CSV1["📄 results/confidence_scores_rescued_{SAMPLE}.csv"]
    S1D --> FIG1["🖼️ figures/rescue_unpositioned_{SAMPLE}.png"]

    %% Stage 2
    SHELL --> S2["Stage 2 · Annotations\nrun_annotations.py"]
    H5AD -->|"io.load_annotations"| S2
    CSV1 -->|"io.load_confidence_scores"| S2
    S2 -->|"annotations.merge_annotations"| CSV2["📄 results/confidence_with_annotations_{SAMPLE}.csv"]

    %% Stage 3
    SHELL --> S3["Stage 3 · Validation\nrun_validation.py"]
    CSV2 -->|"io.load_merged_data"| S3
    S3 -->|"validation.sweep_thresholds\nthresholds 0.0–0.7"| S3A["Spatial coherence metrics\nseparation ratio, median kNN"]
    S3A --> CSV3["📄 results/spatial_coherence_by_threshold.csv"]
    S3A --> FIG3A["🖼️ figures/validation_hippocampal_spatial_maps.png"]
    S3A --> FIG3B["🖼️ figures/validation_coherence_metrics.png"]

    %% Stage 4
    SHELL --> S4["Stage 4 · Visualization\nrun_visualization.py"]
    CSV2 -->|"io.load_merged_data"| S4
    S4 -->|"--type neurons"| FIG4A["🖼️ figures/all_neurons_spatial_combined.png"]
    S4 -->|"--type dcx\n(requires dcx_expression_{SAMPLE}.csv)"| FIG4B["🖼️ figures/dcx_dg_spatial_clean.png"]

    %% Styling
    classDef input  fill:#2d4a6e,color:#fff,stroke:#4a90d9
    classDef stage  fill:#3b3b5c,color:#fff,stroke:#8888cc
    classDef output fill:#1e4d2b,color:#fff,stroke:#4caf50
    classDef fig    fill:#4d2b1e,color:#fff,stroke:#ff7043
    classDef shell  fill:#4a3b00,color:#fff,stroke:#ffc107

    class DATA,H5AD,CONF input
    class SHELL shell
    class S1,S1A,S1B,S1C,S1D,S2,S3,S3A,S4 stage
    class CSV1,CSV2,CSV3 output
    class FIG1,FIG3A,FIG3B,FIG4A,FIG4B fig
```

## Stage summary

| Stage | Script | Input | Output |
|---|---|---|---|
| 1 · Scoring | `run_scoring.py` | Trekker output dir (per sample) | `confidence_scores_rescued_{SAMPLE}.csv`, diagnostic figure |
| 2 · Annotations | `run_annotations.py` | h5ad + stage-1 CSVs | `confidence_with_annotations_{SAMPLE}.csv` |
| 3 · Validation | `run_validation.py` | stage-2 CSVs | `spatial_coherence_by_threshold.csv`, validation figures |
| 4 · Visualization | `run_visualization.py` | stage-2 CSVs (+ optional DCX CSV) | spatial map figures |

### Key scoring details (Stage 1)
- DBSCAN minPts cascade `4 → 3 → 2` rescues unpositioned cells
- Composite score penalized for low signal fraction and ambiguous Trekker status
- Output columns defined in `spatial_confidence/config.py → OUTPUT_COLUMNS`