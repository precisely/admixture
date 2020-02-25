# admixture

Tool for running admixture builder


## Setup

Python 3.6+ required

Binaries for plink and admixture need to be on your PATH. Plink is expected to be called as "plink1.9" and admixture is expected as "admixture"

For local installation, a virtual environment is highly recommended. You can create one and activate it with the following:
```
$ python3.7 -m venv ~/Envs/precisely
$ source ~/Envs/precisely/bin/activate
```

When finished, deactivate with
```
$ deactivate
```


## Installation

The go to the top level of the code (with the requirements.txt file), and run:
```
$ pip install -r requirements.txt
```

Running the following, should now produce a help output:
```
$ python ancestry.py admixture --help
```


## Running

A single, starter admixture run on a VCF can be performed with:
```
$ python ancestry.py admixture start CONFIG_FILE TEST_VCF OUTPUT
```

The config file here follows the format described by templates/admix.json and points to populations to use.

The test vcf is the sample you wish to find admixture for, and the filename prefix MUST match the sample name. This means that input VCFs here are single sample. For example: sample name "aneil" is in VCF named "aneil.vcf.gz"


To run the full hierarchical admixture pipeline, use this:
```
$ python ancestry.py admixture full -d -t 36 /path/to/input.vcf.gz /path/to/output.json
```

The -d flag will create JSON dumps of intermediate admixtures as the subpopulation tests are completed. The -t flag, for threads, indicates the amount of threads to use when computing admixtures.

An example of a real command used for a 1x Gencove imputed VCF:

```
$ python ancestry.py admixture full -d -t 36 /scratch/giab/aneil/1x/aneil/78a75379-db2b-4beb-ab0e-539ab0293492/aneil.vcf.gz /scratch/giab/aneil/1x/aneil/78a75379-db2b-4beb-ab0e-539ab0293492/aneil.out.json
```




## Tests

Run all tests with:
```
$ nosetests
```

To run tests with more output:
```
$ nosetests --verbose
$ nosetests -v
```

To prevent Nose from capturing stdout:
```
$ nosetests --nocapture
$ nosetests -s
```
