from ComBat_seq import *
import argparse

parser = argparse.ArgumentParser()
# ~~~~Module Required Arguments~~~~ #
parser.add_argument("-f", "--filename", type=str,
                    help="Name of file to be read")
parser.add_argument("-o", "--output_filename", type=str,
                    help="The basename used for file output")

args = parser.parse_args()

# check file exists
if (args.filename is not None):
    f = open(args.filename)  # open input file

    out_filename = args.output_filename  # open output file

    with open(out_filename, 'w') as out:
        ########################
        # write some code that calls combat seq code
        ########################

    f.close()

else:
    print('Not given a valid file.')
