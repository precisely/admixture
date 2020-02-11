# admixture
Tool for running admixture builder


##Setup

Python 3.6+ required


Binaries for plink and admixture need to be on your PATH. Plink is expected to be called as "plink1.9" and admixture is expected as "admixture"

For local installation, a virtual environment is highly recommended. You can create one and activate it with the following:

`$ python3.7 -m venv ~/Envs/precisely` 

`$ source ~/Envs/precisely/bin/activate`

When finished, deactivate with 
`$ deactivate`


##Installation
The go to the top level of the code (with the setup.py file), and run:

`pip install .`

This should now make the command `$ ancestry` available to you via the command line


##Running

A single, starter admixture run on a VCF can be performed with:

`$ ancestry admixture start CONFIG_FILE TEST_VCF OUTPUT`

The config file here follows the format described by templates/admix.json and points to populations to use.

The test vcf is the sample you wish to find admixture for, and the filename prefix MUST match the sample name. This means that input VCFs here are single sample. For example: sample name "aneil" is in VCF named "aneil.vcf.gz"

TODO: admixture pipeline incorporating subpop tests
