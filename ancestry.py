import json
import logging
import click
import os
import re
from datetime import datetime
import ancestry.admixture

logger = logging.getLogger("ancestry")
handler = logging.StreamHandler()
formatter = logging.Formatter(
    '%(asctime)s %(name)-8s [%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

@click.group()
def admixture():
    pass


@admixture.command()
# JSON file with the params listed below
@click.argument('config_file', type=click.File('r'))
@click.argument('test_ped')  # MUST BE INDIVIDUAL VCF FILE
@click.option('-t', '--threads', default=1)
@click.argument('output', type=click.File('w'), required=True)
def single(config_file, test_ped, threads, output):
    """
    The threading option affects supervised admixture subprocess ONLY

    Must make sure binaries for plink and admixture can be found by admixture.py
    Run admixture given the files and the k
    The config file requires the following:
    {
        "plink": "pathtobinary",
        "admixture": "pathtobinary",
        "ref_ped": "path",
        "ref_bim": "path",
        "ref_fam": "path",
        "k": int
    }
    """
    params = json.load(config_file)
    admixture = ancestry.admixture.run_admix(params, test_ped, threads)
    results = ancestry.admixture.postprocess(
        admixture[0], admixture[1], admixture[2])
    json.dump(results, output, indent=2)


@admixture.command()
@click.option('-d', '--debug', is_flag=True,
              default=False, help='Enables debug mode.')
# pop file created for admixture (see steps of admix in admixture.py
@click.argument('pop')
@click.argument('q')  # Q file generated by admixture for the same sample set
# k number use to run admixture with (the same number used in the *.config
# file from admix command)
@click.argument('k')
# output file to dump JSON to
@click.argument('output', type=click.File('w'), default='-')
def postprocess(debug, pop, q, k, output):
    """
    Only the postprocessing step of the admixture pipeline. This requires a *.pop file and a Q file generated by
    ADMIXTURE (for the same sample set). You can generate this files by running "ancestry admix"
    """
    results = ancestry.admixture.postprocess(pop, q, k)
    logger.debug(results)


@admixture.command()
@click.argument('global_admix', type=click.File('r'))
@click.argument('test_ped')  # MUST BE INDIVIDUAL VCF FILE
@click.argument('output', type=click.File('w'), default='-')
def subpopadmix(global_admix, test_ped, output):
    """
    Determines what subpopulations need to be tested based on the admixture from global run (see admix command).
    The config file path for each subpopulation test are provided within this code.

    """

    global_json = json.load(global_admix)
    to_do = ancestry.admixture.subpoptest(global_json)
    total_json = {
        "global": global_json
    }

    for pop in to_do:
        prefix, k = ancestry.admixture.create_reference(pop, "data/GlobalMerge")
        params = {
            "ref_ped": prefix + ".bed",
            "ref_bim": prefix + ".bim",
            "ref_fam": prefix + ".fam",
            "k": k
        }
        preresults = ancestry.admixture.run_admix(params, test_ped)
        results = ancestry.admixture.postprocess(
            preresults[0], preresults[1], preresults[2])
        total_json[pop] = results

    json.dump(total_json, output, indent=2)


@admixture.command()
# json file produced by subpopadmix
@click.argument('full_json', type=click.File('r'))
@click.argument('output', type=click.File('w'), default='-')
def handleedges(debug, full_json, output):
    full_json = json.load(full_json)
    edged_json = ancestry.admixture.filters(full_json)
    json.dump(edged_json, output, indent=2)


@admixture.command()
@click.argument("poptest")
@click.argument("global_prefix")
def create_ref(poptest, global_prefix):
    """
    Use create_reference() to create a new bed/bim/fam group using the population defined by the poptest name.
    This name should match to POPULATIONS in populations.py
    """
    prefix, k = ancestry.admixture.create_reference(poptest, global_prefix)
    params = {
        "ref_ped": prefix + ".bed",
        "ref_bim": prefix + ".bim",
        "ref_fam": prefix + ".fam",
        "k": k
    }
    logger.debug(params)


@admixture.command()
@click.argument('test_ped')  # MUST BE INDIVIDUAL VCF FILE
@click.option('-t', '--threads', default=1)
@click.argument('output', type=click.File('w'), required=True)
def init_global(test_ped, threads, output):
    """
    Use create_reference() to create a new bed/bim/fam group using the population defined by the poptest name.
    This name should match to POPULATIONS in populations.py

    TEST SAMPLE NAME MUST BE THE SAME AS SAMPLE FILENAME
    """
    test_prefix = re.sub('\..+', '', test_ped)
    sample_name = re.sub('\/.+\/', '', test_prefix)
    prefix, k = ancestry.admixture.create_reference("global", "data/GlobalMerge", sample_name)
    fullpath = os.path.abspath(prefix + ".bed")
    path = re.sub("\.bed", "", fullpath)
    params = {
        "ref_ped": path + ".bed",
        "ref_bim": path + ".bim",
        "ref_fam": path + ".fam",
        "k": k
    }
    logger.debug(params)
    preresults = ancestry.admixture.run_admix(params, test_ped, threads)
    results = ancestry.admixture.postprocess(
        preresults[0], preresults[1], preresults[2])
    json.dump(results, output, indent=2)


@admixture.command()
@click.argument('test_ped')  # MUST BE INDIVIDUAL VCF FILE
@click.option('-t', '--threads', default=1)
@click.argument('output', type=click.File('w'), required=True)
def full(test_ped, threads, output):
    """
    The full hierarchical ancestry caller for a given sample VCF.

    TEST SAMPLE NAME MUST BE THE SAME AS SAMPLE FILENAME.

    The threading option affects supervised admixture subprocess ONLY
    """
    test_prefix = re.sub('\..+', '', test_ped)
    sample_name = re.sub('\/.+\/', '', test_prefix)
    prefix, k = ancestry.admixture.create_reference("global", "data/GlobalMerge", sample_name)
    fullpath = os.path.abspath(prefix + ".bed")
    path = re.sub("\.bed", "", fullpath)
    params = {
        "ref_ped": path + ".bed",
        "ref_bim": path + ".bim",
        "ref_fam": path + ".fam",
        "k": k
    }
    logger.debug(params)
    preresults = ancestry.admixture.run_admix(params, test_ped, threads)
    global_json = ancestry.admixture.postprocess(
        preresults[0], preresults[1], preresults[2])

    to_do = ancestry.admixture.subpoptest(global_json)
    total_json = {
        "global": global_json
    }

    for pop in to_do:
        prefix, k = ancestry.admixture.create_reference(pop, "data/GlobalMerge", sample_name)
        params = {
            "ref_ped": prefix + ".bed",
            "ref_bim": prefix + ".bim",
            "ref_fam": prefix + ".fam",
            "k": k
        }
        preresults = ancestry.admixture.run_admix(params, test_ped)
        results = ancestry.admixture.postprocess(
            preresults[0], preresults[1], preresults[2])
        total_json[pop] = results

    json.dump(total_json, output, indent=2)

@click.group()
def cli():
    pass


cli.add_command(admixture)

if __name__ == '__main__':
    cli()