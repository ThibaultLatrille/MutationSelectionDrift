# Replicate experiments 

Renaming necessary for compatibility with v1.0 of SimuEvol (executable renaming)


```python
!cp ../SimuEvol/build/SimuProfile ../SimuEvol/build/SimuDiv
!cp ../SimuEvol/build/PolyProfile ../SimuEvol/build/SimuPoly
```

## Simulated dataset

Simulations where traits (effective population size, mutation rate, generation time) are fluctuating along each branch during the simulation (traits are changing every generation).


```python
!python3 simulated_experiment.py --experiment Gen5Ma150L5000 --nbr_cpu 4
```

Other simulations where traits (effective population size, mutation rate, generation time) are considered constant for a specific branch during the simulation (traits are changing for each branch).


```python
!sed -i 's#BRANCH_WISE_CORRELATION: false#BRANCH_WISE_CORRELATION: true#g' config.yaml
```


```python
!python3 simulated_experiment.py --experiment Gen5Ma150L5000BranchWise --nbr_cpu 4
```
