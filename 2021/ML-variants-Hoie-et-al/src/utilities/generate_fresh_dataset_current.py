import sys
import prism_machine_learning as prismml

from argparse import ArgumentParser
parser = ArgumentParser(description="Select training, validation proteins")
parser.add_argument("-o", dest="OUTPUT",
                    default="data/preprocessed.pkl", help="Print error messages")
args = parser.parse_args()

outpath = str(args.OUTPUT)
print("Outpath:", outpath)

print("Generating norm dataset ...")
dfs_proc, dfs_names, dfs_raw = prismml.generate_load_preprocessed_datasets("data/preprocessed/prism*.txt",
                                            prevent_create_new = False, force_create_new = True, normalize_mave_only = True, ranknorm = True,  filetype = "prismdata",
                                            preprocessed_path=outpath)
