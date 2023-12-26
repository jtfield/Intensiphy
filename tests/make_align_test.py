import os
import sys
sys.path.insert(0, os.path.abspath(".."))

from modules.fetch_and_align import make_align

def make_align_test():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    print(split_path_and_name)

    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    test_align = split_path_and_name[0] + '/combo.fas'

    make_align(False, test_output_dir, test_align)

    symlink_path = test_output_dir + '/intermediate_files/alignment.fas'
    
    assert os.path.islink(symlink_path)