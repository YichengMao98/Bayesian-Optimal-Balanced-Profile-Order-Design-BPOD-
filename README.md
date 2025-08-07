## Script Naming Convention

All R script files follow the naming pattern:

ModelType_OrderEffect_Kappa.R

Where:

- `ModelType` ∈ {`MNL`, `PML`}  
  Represents the model used in estimation:  
  - `MNL`: Multinomial Logit  
  - `PML`: Panel Mixed Logit

- `OrderEffect` ∈ {`primacy`, `recency`, `central`}  
  Represents the type of profile order effect simulated.

- `Kappa` ∈ {`1`, `0.3`, `0.1`}  
  Represents the value of the parameter κ.
