import sys
import re
import os

def identify_test_sample(pop_name, fam_file, sample_name):
    #pop_name = re.sub('fam', 'pop', fam_file)
    o = open(pop_name, 'w')
    with open(fam_file, 'r') as f:
        for line in f:
            line = re.sub('^'+sample_name, '_'+sample_name, line)
            o.write(line)
            print("Modified sample name at {}".format(line))
    o.close()
    if os.path.exists(pop_name):
        return True
    else:
        return False

print(identify_test_sample(sys.argv[1], sys.argv[2], sys.argv[3]))