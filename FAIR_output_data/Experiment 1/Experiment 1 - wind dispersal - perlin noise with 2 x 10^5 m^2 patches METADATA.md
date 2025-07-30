## Experiment description:

Initial tree cover is varied in the range \[0.1, 0.87], and the spatial domain size is set to 960 by 960 m. For each level of initial tree cover, we obtain the corresponding self ignition rate as follows. We observe tree cover slope for 20 model runs and average these to obtain an indicator of simulated ecosystem instability. If the absolute value of the average of tree cover slopes is smaller, we associate this with a higher degree of model instability. After all, this indicates that the fraction of model runs that had a positive tree cover slope versus those with a negative tree cover slope was closer to 50%. We minimize this absolute value using a gradient descent algorithm to obtain the corresponding fire ignition rate.



\#### File creation date: 29-07-2025



## Parameter settings

**Dispersal mode**: wind  
**Resolution**: 1 m  
**Spatial domain size**: 960 m  
**Initial tree distribution**:  perlin\_noise
**patch\_width**: 500m
**Resource grid resolution**: 32m  
**Time steps**: 100  
**Batch mode**: saddle\_search

