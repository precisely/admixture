import logging
import os
import re
import subprocess
from data.populations import POPULATIONS

from ancestry import get_data_filename
from ancestry.utils import cwd

log = logging.getLogger("ancestry")

# rules determining threshold sums and included populations for that sum
# for each population test. pop test names must match across functions and
# output
ruleset = {
    "africa": {
        "threshold": 0.02,
        "populations": ["AA_Ref_African"]
    },
    "america": {
        "threshold": 0.05,
        "populations": ["AA_Ref_Amerindian"]
    },
    "eastasia": {
        "threshold": 0.03,
        "populations": ["AA_Ref_EastAsian"]
    },
    "southasia": {
        "threshold": 0.25,
        "populations": ["AA_Ref_SouthAsian"]
    },
    "europe": {
        "threshold": 0.10,
        "populations": ["AA_Ref_NorthernEuropean", "AA_Ref_SouthernEuropean"]
    }

}


def admix_prep(params, test_ped):
    """

    This will use the given file prefix to generate the inputs for ADMIXTURE supervised

    - make ped file with incoming customer VCF
    - use plink to do the merge (it will handle triallelics and other errors)
    - if plink finds an error in the merge, remove the offending snp columns from BOTH the ref and test peds
    - try previous step until it works (possibly >50 removes)
    - add row for test file to the fam and pop files
    - run admixture

    """

    # check if we have the files available
    for key, value in params.items():
        if key == "k":
            continue
        elif os.path.exists(value) is False:
            log.debug("{} did not exist at {}".format(key, value))
            return False

    if os.path.exists(test_ped) is False:
        log.debug("{} did not exist".format(test_ped))
        return False

    return True


# performs the creation of pop file for admixture and identifies the test
# sample with _
def identify_test_sample(fam_file, sample_name):
    pop_name = re.sub('fam', 'pop', fam_file)
    o = open(pop_name, 'w')
    with open(fam_file, 'r') as f:
        for line in f:
            line = re.sub('^'+sample_name, '_'+sample_name, line)
            o.write(line)
    o.close()
    if os.path.exists(pop_name):
        return True
    else:
        return False


def run(*args, **kwargs):
    log.debug(args)
    args = ' '.join(args)
    return subprocess.run(args, shell=True, stdout=subprocess.PIPE)


def plink(*args, **kwargs):
    cmd = ['plink1.9'] + list(args)
    try:
        proc = run(*cmd, **kwargs)
    except subprocess.CalledProcessError:
        log.error(f"Call to plink failed: {proc.stderr}")
        raise
    return proc.stdout


def admixture(*args, **kwargs):
    cmd = ['admixture'] + list(args)
    try:
        proc = run(*cmd, **kwargs)
    except subprocess.CalledProcessError:
        log.error(f"Call to admixture failed: {proc.stderr}")
        raise
    return proc.stdout


def create_reference(poptest_name, global_prefix, sample):
    """

    Use dict in populations.py to slice the necessary populations out of the Global bed defined by global_prefix
    Input: name of pop test that matches populations.py, the prefix-only path to the bed file being sliced out of
    Output: path + prefix to new bed/bim/fam files

    """
    pops = POPULATIONS[poptest_name]["pops"]
    k = POPULATIONS[poptest_name]["k"]

    with open("keeper.txt", "w") as f:
        for pop in pops:
            print(pop, file=f)

    try:
        plink("--bfile", global_prefix, "--make-bed", "--keep-fam", "keeper.txt", "--out",
              global_prefix + "." + poptest_name + "." + sample)
    except subprocess.CalledProcessError:
        log.debug("Error in creating new reference.")
        raise

    output_prefix = global_prefix + "." + poptest_name + "." + sample
    return output_prefix, int(k)


# perform merge of test and reference samples and run ADMIXTURE
def run_admix(params, test_ped, threads):
    # assign params to vars
    #ref_ped = get_data_filename(params['ref_ped'])
    ref_ped = params['ref_ped']
    k = str(params['k'])

    test_prefix = re.sub('\..+', '', test_ped)
    ref_prefix = re.sub('\.bed', '', ref_ped)
    out_prefix = test_prefix + ".out"
    sample_name = re.sub('\/.+\/', '', test_prefix)

    with cwd(os.path.dirname(os.path.abspath(test_prefix))):
        plink("--bfile", ref_prefix, "--write-snplist")
        plink("--vcf", test_ped, "--extract", "plink.snplist", "--make-bed", "--out", test_prefix)

        # plink merge test with refs
        # first cycle will catch a whole lot of shit
        plink("--bfile", ref_prefix, "--bmerge", test_prefix, "--out", out_prefix)

        # use the count of the missnps file to determine next steps
        missnp_file = out_prefix + ".missnp"
        if os.path.exists(missnp_file):
            process_missnp(missnp_file, ref_prefix, test_prefix, out_prefix)

        # if no missnp file and BED file created on first pass
        elif os.path.exists(missnp_file) is False and os.path.exists(out_prefix + ".bed"):
            log.debug("Successful merge at {}".format(out_prefix + ".bed"))
        # if no missnp file or BED file created on first pass
        else:
            log.debug(
                "Could not find the missnp file at {} and could not find output file of {}".format(
                    missnp_file, out_prefix))
            raise RuntimeError(".missnp file could no be found.")

        # a stage of filtering to make sure there aren't any completely
        # ungenotyped loci
        plink("--bfile", out_prefix + ".NoMulti.merged", "--geno", "0.999", "--make-bed", "--out", out_prefix + ".NoMulti.merged.filtered")

        if os.path.exists(out_prefix + ".NoMulti.merged.filtered.bed"):
            out_prefix = out_prefix + ".NoMulti.merged.filtered"
        else:
            raise RuntimeError(
                "Could not locate the filtered BED at {}".format(
                    out_prefix + ".NoMulti.merged.filtered.bed"))

        # prepare the pop file for admixture by adding "_" to the test sample name
        # in the fam file may
        pop_create = identify_test_sample(out_prefix + ".fam", sample_name)
        if pop_create is False:
            raise RuntimeError(
                "Failed to create popfile using {}".format(
                    out_prefix + ".fam"))

        # begin admixture run
        admixture("--supervised", out_prefix + ".bed", k, "-j"+str(threads))

    # catch for successful admixture run
    if os.path.exists(out_prefix + "." + k + "." + "Q"):
        log.debug(
            "Admixture Q file completed successfully. Time for postprocessing.")
    else:
        raise RuntimeError(
            "Could not locate the admixture Q file at {}".format(
                out_prefix + "." + k + "." + "Q"))

    # return the Q file and the pop file to the calling variable
    Q = out_prefix + "." + k + "." + "Q"
    pop = out_prefix + ".pop"
    return [pop, Q, k]


def bim_check(name):
    """
    boolean
    See if input file exists, use this function in PLINK processing to alert file read errors
    """
    # use the bim file
    if os.path.exists(name + ".bim") is False:
        return False
    else:
        return True


def process_missnp(missnp_file, ref_prefix, test_prefix, out_prefix):
    missnp_max = 5000  # max number of allowed variants in missnp file
    with open(missnp_file, 'r') as f:
        lines = sum(1 for _ in f)

    # if the missnp file contains too many markers, try to flip strands for
    # them
    if lines > missnp_max:
        log.debug(
            "Strand errors detected on merge of {} and {} with out prefix {}. Attempting repair".format(
                ref_prefix, test_prefix, out_prefix))
        # flipping command

        if bim_check(ref_prefix) is False:
            raise RuntimeError("Could not locate {}, possible file read error".format(ref_prefix + ".bim"))

        try:
            plink("--bfile", ref_prefix, "--flip", out_prefix + ".missnp", "--make-bed", "--out",
                  ref_prefix + "Flipped")
        except subprocess.CalledProcessError:
            log.debug("Error in strand repair.")
            raise

        # provide fam file for new flip naming
        try:
            run("cp", ref_prefix + ".fam", ref_prefix + "Flipped.fam")
        except subprocess.CalledProcessError:
            log.debug("Error creating the fam file {}".format(ref_prefix + "Flipped.fam"))
            raise

        # attempt to merge again
        try:
            plink("--bfile", ref_prefix + "Flipped", "--bmerge", test_prefix, "--out", out_prefix)
        except subprocess.CalledProcessError:
            log.debug(
                "Error in running merge after flipping {}".format(
                    ref_prefix + "Flipped"))
            raise

        # noinspection PyAugmentAssignment
        ref_prefix = ref_prefix + "Flipped"
        if bim_check(ref_prefix) is False:
            raise RuntimeError("Could not locate {}, possible file read error".format(ref_prefix + ".bim"))

        # check line count of missnps after flipped run
        with open(missnp_file, 'r') as f:
            l2 = sum(1 for _ in f)

        # if variants in missnp below the threshold, simply remove them
        if l2 < missnp_max:

            try:
                plink("--bfile", ref_prefix, "--exclude", out_prefix + ".missnp", "--make-bed", "--out",
                      ref_prefix + "NoMulti")
                plink("--bfile", test_prefix, "--exclude", out_prefix + ".missnp", "--make-bed", "--out",
                      out_prefix + ".NoMulti")
            except subprocess.CalledProcessError:
                log.debug(
                    "Error attempting to remove missnp variants from {} and {}.".format(
                        ref_prefix, out_prefix + ".missnp"))
                raise
            ref_prefix += "NoMulti"
        else:
            log.debug(
                "Too many variants in missnp file after flipping and exclusion")
            raise RuntimeError(
                "Too many variants remaining after flipping and exclusion")
    # if markers below threshold count, try to remove them
    else:
        try:
            plink("--bfile", ref_prefix, "--exclude", out_prefix + ".missnp", "--make-bed", "--out",
                  ref_prefix + "NoMulti")
            plink("--bfile", test_prefix, "--exclude", out_prefix + ".missnp", "--make-bed", "--out",
                  out_prefix + ".NoMulti")
        except subprocess.CalledProcessError:
            log.debug("Error attempting to remove missnp variants from {} and {}.".format(ref_prefix,
                                                                                          out_prefix + ".missnp"))
            raise
        ref_prefix = ref_prefix + "NoMulti"

    if bim_check(ref_prefix) is False:
        raise RuntimeError("Could not locate {}, possible file read error".format(ref_prefix + ".bim"))

    # try to perform a merge after the above repairs. If this merge fails,
    # sink the whole try
    log.debug("Starting a final merge attempt")
    try:
        plink("--bfile", ref_prefix, "--bmerge", out_prefix + ".NoMulti",  "--out", out_prefix + ".NoMulti.merged")
    except subprocess.CalledProcessError:
        log.debug("Error in the final merge on {} and {}".format(test_prefix, ref_prefix))
        raise


# postprocess admixture results into scoring format
def postprocess(pop, q, k):
    """
    read in the two files, measure line counts, and output some JSON about the admixure results of the test sample
    """
    ref_min = 0.99

    if os.path.exists(pop) is False or os.path.exists(q) is False:
        raise RuntimeError(
            "Could not read the input files {} and {}".format(
                pop, q))

    with open(pop, 'r') as pf:
        pcount = sum(1 for _ in pf)

    with open(q, 'r') as qf:
        qcount = sum(1 for _ in qf)

    if pcount != qcount:
        raise RuntimeError(
            "The line counts of the pop and Q files did not match each other.")
    header = [None] * int(k)
    test = []
    with open(pop) as bf1:
        with open(q) as bf2:
            for pl, ql in zip(bf1, bf2):
                # the delimiters are doing some weird shit, hopefully they just
                # use spaces and stay that way
                pllist = pl.split(' ')
                name = pllist[0]
                ql = ql.strip()
                qllist = ql.split(' ')
                # some exceptions to skip to next line. If its a test sample
                # and if the population header position has already been
                # determined
                if re.match('^_', name):
                    test = qllist
                    continue
                if name in header:
                    continue
                # searching for the column number for the ref pops using high
                # percentage matches in the Q file columns
                colcount = 0
                for col in qllist:
                    if float(col) > ref_min:
                        header[colcount] = name
                        break
                    else:
                        colcount += 1

                # check to see if the name is now present in the header dict.
                # If it's not, there's an issue.
                if name in header:
                    continue
                else:
                    raise RuntimeError(
                        "Reference pop with no samples named {}".format(name))

    test_len = len(test)
    header_len = len(header)
    out = {}
    for i in range(header_len):
        if header[i] is None:
            raise RuntimeError(
                "Population of k column {} is unknown. Please examine.".format(i))
        out[header[i]] = test[i]

    return out


# determine the global tests which need to be run based on input from
# global admixture run
def subpoptest(global_admix):
    """
    input is the dict-converted JSON from the initial (global) admix output file. Returns list of subgroup tests to run.

    Define the tests and admixture thresholds that determine what subpops need to be run through admixture again

    Test ascertainment is currently implemented by the sum of the admixture in the relevant populations.
    If we see it is needed, we can change to more population specific thresholds, or even a machine learning approach
    if we can find a decent number of samples.
    """

    to_do = []

    for test, pops in ruleset.items():
        summer = 0.0

        for pop in pops['populations']:
            try:
                admix = float(global_admix[pop])
            except BaseException:
                raise RuntimeError(
                    'Unable to find admixture for {} in global_admix'.format(pop))
            summer += admix

        if summer >= pops["threshold"]:
            to_do.append(test)
        else:
            continue

    if len(to_do) > 0:
        log.debug("Returning the following tests: {}".format(to_do))
        return to_do
    elif len(to_do) == 0:
        log.debug("NO SUBPOPULATION TESTS MATCH. You may want to inspect this sample manually!")
        return to_do


def filters(full_json):
    """
    Define edge cases for ancestry tests, and re-process results to fit them
    Input JSON made by subpoptest, output json to dump
    """

    # organize results as needed for app. search oceania and south asian
    # differently. Oceania attempts to pull first from any oceanian test, and
    # then from asia if its not present.
    log.debug("Beginning filtering and edge case handling.")

    # readjust global percentages based on threshold ruleset before processing
    # subgroup weights
    log.debug("Adjusting global percentages based on thresholds before processing")
    for group, rules in ruleset.items():
        summer = 0.0
        for rule in rules["populations"]:
            try:
                admix = float(full_json["global"][rule])
            except BaseException:
                raise RuntimeError(
                    "Unable to find admixture for {} in global populations".format(rule))
            summer += admix
        if summer >= rules["threshold"]:
            # special rule to re-weight ashkenazi if above zero threshold and below an adjustment threshold
            if group == "ashkenazi" and float(full_json["global"]["AA_Ref_Ashkenazi"]) <= 0.35:
                full_json["global"]["AA_Ref_Ashkenazi"] = float(full_json["global"]["AA_Ref_Ashkenazi"]) * 0.67
            continue
        else:
            # zero in global only
            for p in rules["populations"]:
                full_json["global"][p] = 0.0
    # redistribute global percentages so they add up to 1.0 again
    log.debug("Redistributing global percentages")
    tosum = list(float(x) for x in full_json["global"].values())
    denom = sum(tosum)
    for pop in full_json["global"]:
        full_json["global"][pop] = float(full_json["global"][pop]) / denom

    # holds the global results from the input JSON we will use to adjust
    # subpop percentages
    global_from_json = full_json["global"]

    # determine where to place south asian and oceanian admix values,
    # depending on presence of asia test
    has_oceania = 0
    has_asia = 0
    has_europe = 0


    # now that results are structured in correct population groups, we need to first make each sub group sum to 1.0
    # because we are taking subsets of each test
    # and then normalize each population by the total of their population
    # container defined by the ruleset
    log.debug("Adjusting subgroups to fit container and normalizing by globals")
    for group, pops in out_json.items():
        gsum = 0.0
        gpops = gsum_ruleset[group]

        # exception: if there is an oceanian test, make sure you normalize
        # using the oceanion populations, not asian
        if group == "oceania" and has_oceania == 1:
            gsum = float(global_from_json["AA_Ref_Oceanian"])
            # out_json is filtered for only "in group" pops, so its fine to sum
            # here
            psum = sum(out_json[group].values())
            log.debug("Oceanian edge case value={} group={} pop={} Psum={} Gsum={} ".format(
                out_json[group]["AA_Ref_Oceanian"], group, "AA_Ref_Oceanian", psum, gsum))
            if psum * gsum != 0.0:
                out_json[group]["AA_Ref_Oceanian"] = (float(
                    out_json[group]["AA_Ref_Oceanian"]) * gsum) / psum
            else:
                out_json[group]["AA_Ref_Oceanian"] = 0.0
            continue

        # if there is no oceanian test, normalize oceania with asia
        if group == "oceania" and has_oceania == 0 and has_asia == 1:
            gpops = gsum_ruleset["asia"]
            for gpop in gpops:
                gsum += float(global_from_json[gpop])
            psum = sum(out_json["asia"].values())
            log.debug(
                "Oceanian with-asia edge case value={} group={} pop={} Psum={} Gsum={} ".format(
                    out_json[group]["AA_Ref_Oceanian"],
                    group,
                    "AA_Ref_Oceanian",
                    psum,
                    gsum))
            if psum * gsum != 0.0:
                out_json[group]["AA_Ref_Oceanian"] = (float(
                    out_json[group]["AA_Ref_Oceanian"]) * gsum) / psum
            else:
                out_json[group]["AA_Ref_Oceanian"] = 0.0
            continue

        # if no asia test, just use the global oceania pop
        if group == "oceania" and has_asia == 0 and has_oceania == 0:
            out_json[group]["AA_Ref_Oceanian"] = global_from_json["AA_Ref_Oceanian"]
            continue

        # if there is asian test, for south asia use that test to normalize
        if group == "south_asia" and has_asia == 1:
            gpops = gsum_ruleset["asia"]
            for gpop in gpops:
                gsum += float(global_from_json[gpop])
            psum = sum(out_json["asia"].values())
            log.debug(
                "South asian with-asia edge case value={} group={} pop={} Psum={} Gsum={} ".format(
                    out_json[group]["AA_Ref_South_Asian"],
                    group,
                    "AA_Ref_South_Asian",
                    psum,
                    gsum))
            if psum * gsum != 0.0:
                out_json[group]["AA_Ref_South_Asian"] = (float(
                    out_json[group]["AA_Ref_South_Asian"]) * gsum) / psum
            else:
                out_json[group]["AA_Ref_South_Asian"] = 0.0
            continue

        # if there is no asian test, just use global for south asians
        if group == "south_asia" and has_asia == 0:
            out_json[group]["AA_Ref_South_Asian"] = global_from_json["AA_Ref_South_Asian"]
            continue

        # if the group is americas, CSA = CSA + (0.5 (CSA+NA)), then NA = (0.5 (CSA+NA)) - CSA
        if group == "americas":
            amdenom = out_json[group]["AA_Ref_Central_South_American"] + out_json[group]["AA_Ref_North_American"]
            out_json[group]["AA_Ref_Central_South_American"] += (0.5 * amdenom)
            if out_json[group]["AA_Ref_Central_South_American"] >= amdenom:
                out_json[group]["AA_Ref_Central_South_American"] = amdenom
            out_json[group]["AA_Ref_North_American"] = amdenom - out_json[group]["AA_Ref_Central_South_American"]

        for gpop in gpops:
            gsum += float(global_from_json[gpop])
        psum = sum(out_json[group].values())

        for pop in pops:
            log.debug("value={} group={} pop={} Psum={} Gsum={} ".format(
                out_json[group][pop], group, pop, psum, gsum))
            if psum * gsum != 0.0:
                out_json[group][pop] = (out_json[group][pop] * gsum) / psum
            else:
                out_json[group][pop] = 0.0

    # delete duplicate south asian entry from asian continent
    del out_json["asia"]["AA_Ref_South_Asian"]

    # move subgroups to first level of JSON, delete super groups
    final_json = {}
    for val in out_json.values():
        for k, v in val.items():
            new = re.sub("AA_Ref_", "", k)
            final_json[new] = v
    # rewieght everything again after the edge cases and dedups are resolved
    denom = sum(final_json.values())
    for pop, val in final_json.items():
        final_json[pop] = (val / denom) * 100


    # fix the pop names by removing AA_Ref

    # make sure the final percentages round out to 100
    # if int() ever changes from rounding down behavior, this will stop working
    n = 100
    dec = {}
    rounded = {}
    for key, val in final_json.items():
        dec[key] = val - int(val)
        rounded[key] = int(val)

    if sum(rounded.values()) == n:
        return rounded
    else:
        s = [(k, dec[k]) for k in sorted(dec, key=dec.get, reverse=True)]
        while sum(rounded.values()) != n:
            log.debug("Adjusting rounded values")
            for key, val in s:
                if sum(rounded.values()) != n:
                    rounded[key] += 1
        return rounded