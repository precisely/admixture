from nose.tools import *
import unittest
import shutil


import ancestry.admixture


# setup before running this entire file
# function name matters
def setup():
    print("setting up nose test")
    # for example:
    global path_admixture
    global path_plink19
    path_admixture = shutil.which("admixture")
    path_plink19 = shutil.which("plink1.9")


# teardown after running this entire file
# function name matters
def teardown():
    # as needed:
    print("tearing down nose test")

# setup for specific tests, for use with @with_setup annotation
# function name does not matter
def setup_individual():
    print("setting up individual test")

# teardown for specific tests, for potential use with @with_setup annotation
# function name does not matter
def teardown_individual():
    print("tearing down individual test")


@with_setup(setup_individual, teardown_individual)
def test_aneil():
    # maybe expected should be in a separate file
    expected = {
        "AA_Ref_AndamanIsland": "0.000010",
        "AA_Ref_Baloch": "0.000379",
        "AA_Ref_Bengali": "0.000010"
        # and whatever else
    }
    admixture = ancestry.admixture.run_admix(...)
    actual = ancestry.admixture.postprocess(
        admixture[0], admixture[1], admixture[2])
    assert expected == actual
