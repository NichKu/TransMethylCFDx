__author__ = "Nicholas Kueng"
__copyright__ = "Copyright 2023, Nicholas Kueng"
__license__ = "MIT"


def get_interval_list_path(bed, target_dir):
    """
    This function creates the path and the name of the interval files for the target and bait region BED files
    """
    base = os.path.basename(bed)
    split_base = os.path.splitext(base)
    file_name = os.path.splitext(base)[0]
    output_file = target_dir + "/" + file_name + ".interval_list"
    return output_file
