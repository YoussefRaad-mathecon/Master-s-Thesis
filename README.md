TEST CODE DO NOT CONSIDER

```
├── README.md                          
├── data/                              
│   ├── players.json
│   ├── events_England.json
│   ├── events_France.json
│   ├── events_Germany.json
│   ├── events_Italy.json
│   └── events_Spain.json
├── R/                                 # R scripts and functions
│   ├── Get shots.R                   # Function to extract shot data
│   └── Get shots extra.R             # Extended shot data extraction
├── analysis/                          # Analysis R Markdown files
│   ├── Load rawdata.Rmd          # Data loading and preprocessing
│   ├── Exploratory analysis.Rmd   # EDA and feature engineering
│   ├── Fit models.Rmd            # Main model fitting (GBM & GLM)
│   ├── Fit models big.Rmd        # Extended model with more features
│   ├── Fit right footed models.Rmd
│   ├── Fit left footed models.Rmd
│   ├── Evaluate models.Rmd        # Main model evaluation
│   ├── Evaluate models big.Rmd
│   ├── Evaluate right footed models.Rmd
│   ├── Evaluate left footed models.Rmd
│   ├── Calculations.Rmd           # xG calculations by counter-attack
│   ├── Hypothesis test.Rmd        # Statistical tests
│   └── Mixture models analysis.Rmd # Distribution fitting
```
