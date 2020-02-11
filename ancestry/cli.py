import json
import logging
import click
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
@click.option('-d', '--debug', is_flag=True,
              default=False, help='Enables debug mode.')
# JSON file with the params listed below
@click.argument('config_file', type=click.File('r'))
@click.argument('test_ped')  # MUST BE INDIVIDUAL VCF FILE
@click.argument('output', type=click.File('w'), default='-')
def start(debug, config_file, test_ped, output):
    """
    #must make sure binaries for plink and admixture can be found by admixture.py
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
    admixture = ancestry.admixture.run_admix(params, test_ped)
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
    print(results)





@click.group()
def cli():
    pass


cli.add_command(admixture)
