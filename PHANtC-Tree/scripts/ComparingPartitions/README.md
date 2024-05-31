# ComparingPartitions v2

This repository contains *comparing_partitions_v2.py*, an automated script to run the original [ComparingPartitions](https://github.com/jacarrico/ComparingPartitions) (Collection of python functions for Partition Comparison - See www.comparingpartitions.info for more information).


Currently, *comparing_partitions_v2.py* compares subsequential partitions of a single method to determine stability regions (i.e. partition threshold intervals in which cluster composition is similar). Cluster similarity between different partitions is assessed with the Adjusted Wallace Coefficient. Other coefficients such as Simpson's Index of Diversity or the Adjusted Rand are also calculated.


## Usage

```bash
-h, --help            show this help message and exit
  -i1 INPUT1, --input1 INPUT1
                        [MANDATORY] Input matrix 1 (table with partitions)
  -i2 INPUT2, --input2 INPUT2
                        Input matrix 2 (table with partitions)
  -o1 ORDER1, --order1 ORDER1
                        Partitions order in matrix 1 (0: min -> max; 1: max ->
                        min) [0]
  -o2 ORDER2, --order2 ORDER2
                        Partitions order in matrix 2 (0: min -> max; 1: max ->
                        min) [0]
  -a ANALYSIS, --analysis ANALYSIS
                        Type of analysis (options: stability or
                        between_methods) [between_methods]
  -t TAG, --tag TAG     [MANDATORY] Tag for output file name
  -n N_OBS, --n_obs N_OBS
                        Minimum number of sequencial observations to consider
                        an interval for method stability analysis [5]
  -thr THRESHOLD, --threshold THRESHOLD
                        Threshold of Adjusted Wallace score to consider an
                        observation for method stability analysis [0.99]
  -log LOG, --log LOG   Log file
```


## Examples

### 1) obtain stability regions 
E.g. - comparing partitions of a single method to determine stability regions (where stability regions correspond to partition threshold ranges where at least 5 subsequent partitions have an Adjusted Wallace > 0.99)
```bash
comparing_partitions_v2.py -i1 PARTITIONS.tsv -o1 0 -a stability -t OUTPUT_NAME -n 5 -thr 0.99 -log LOG
```
