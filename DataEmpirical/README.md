# Replicate experiments 

Parameter of inference (Ncat, number of points, burn-in) are located in _config.yaml_

All the dataset (tree, alignments, ...) *must be* in the folder given in the parameter '--name'

The filename of the rootedtree (in newick format) *must be* provided throught the parameter '--tree'

The user defined prefix of the experiment can be provided trought the parameter '--prefix'

## Extract alignnments (Isopods & OrthoMam)

CDS files (in .ali format) are extracted in the subfolder 'singlegene_alignments'. 



```python
!sh ./extract_alignments.sh
```

## mammalian dataset (with replicates)

The list of the alignment filenames *must be* provided trought the parameter '--cds' 


```python
!python3 replicate_experiments.py --prefix Cat50 --name OrthoMam --tree rootedtree.lht.nhx --cds cds.highcoverage.list --sample 18 --replicate 4 --lht life_history_traits.tsv --nbr_cpu 4

```

## Isopods (with replicates)

The list of the alignment filenames *must be* provided trought the parameter '--cds' 


```python
!python3 replicate_experiments.py --prefix Ncat50 --name Isopods --tree rootedtree.nhx --cds cds.highcoverage.list --sample 12 --replicate 4 --nbr_cpu 4

```

## Primates (no replicate)


```python
!python3 empirical_experiment.py --prefix IWCat50LHT --name Primates --cds CDS.ali --tree rootedtree.nhx --lht life_history_traits.tsv --nbr_cpu 4

```


```python
!python3 empirical_experiment.py --prefix IWCat50LHTCalib --name Primates --cds CDS.ali --tree rootedtree.nhx --lht life_history_traits.tsv --calibs calibs.tsv --nbr_cpu 4

```

## Custom made experiment

With replicates:


```python
!python3 replicate_experiments.py --help
```

    usage: replicate_experiments.py [-h] [-p PREFIX] -n NAME --tree TREE --cds CDS
                                    [--sample SAMPLE] [--replicate REPLICATE]
                                    [--lht LHT] [--calibs CALIBS]
                                    [--intersection INTERSECTION] [-s SCREEN]
                                    [-b SBATCH] [-c NBR_CPU]
                                    [--random_state RANDOM_STATE]
    
    optional arguments:
      -h, --help            show this help message and exit
      -p PREFIX, --prefix PREFIX
      -n NAME, --name NAME
      --tree TREE
      --cds CDS
      --sample SAMPLE
      --replicate REPLICATE
      --lht LHT
      --calibs CALIBS
      --intersection INTERSECTION
      -s SCREEN, --screen SCREEN
      -b SBATCH, --sbatch SBATCH
      -c NBR_CPU, --nbr_cpu NBR_CPU
      --random_state RANDOM_STATE


Without replicates: 


```python
!python3 empirical_experiment.py --help
```

    usage: empirical_experiment.py [-h] [-p PREFIX] -n NAME --cds CDS --tree TREE
                                   [--lht LHT] [--calibs CALIBS] [-s SCREEN]
                                   [-b SBATCH] [-c NBR_CPU]
    
    optional arguments:
      -h, --help            show this help message and exit
      -p PREFIX, --prefix PREFIX
      -n NAME, --name NAME
      --cds CDS
      --tree TREE
      --lht LHT
      --calibs CALIBS
      -s SCREEN, --screen SCREEN
      -b SBATCH, --sbatch SBATCH
      -c NBR_CPU, --nbr_cpu NBR_CPU


# Analysis of the different replicates


```python
!cd ../scripts && python3 replicate_analysis.py
```