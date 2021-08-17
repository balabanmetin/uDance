import argparse
import sys
from multiprocessing import cpu_count
from os.path import abspath, expanduser
from sys import stderr
from didactic.decompose import decompose
from didactic.refine import refine
from didactic.stitch import stitch


def options_config():
    # Input arguments parser
    parser = argparse.ArgumentParser(description='Massively scalable divide-and-conquer phylogenetic inference',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--version', action='store_true', help='print the current version')
    # parser.add_argument('--debug', action='store_true', help='Print the traceback when an exception is raised')
    subparsers = parser.add_subparsers(title='commands',
                                       description='decompose   Create local inference (RAxML) tasks\n'
                                                   'refine      Refine partitions via consensus (ASTRAL) \n'
                                                   'stitch      Stitch back locally refined trees',
                                       help='Run run_didactic.py {commands} [-h] for additional help',
                                       dest='{commands}')

    # To make sure that subcommand is required in python >= 3.3
    python_version = sys.version_info
    if (python_version[0] * 10 + python_version[1]) >= 33:
        subparsers.required = True

    # decompose command subparser
    parser_decompose = subparsers.add_parser('decompose',
                                       description='Create local refinement tasks')
    parser_decompose.add_argument("-j", "--jplace", dest="jplace_fp",
                            help="path to the jplace placement file", metavar="FILE")
    parser_decompose.add_argument("-f", "--threshold", dest="threshold", default="600",
                            help="maximum number of elements in each cluster")
    parser_decompose.add_argument("-o", "--output", dest="output_fp",
                            help="path for the output directory where files will be placed",
                            metavar="DIRECTORY")
    parser_decompose.add_argument("-s", "--alignment-dir", dest="alignment_dir_fp",
                            help="path for input directory which contains "
                                 "extended reference alignment files (FASTA), "
                                 "containing reference and query sequences.",
                            metavar="FILE")
    parser_decompose.add_argument("-p", "--protein", dest="protein_seqs", action='store_true', default=False,
                            help="input sequences are protein sequences")
    parser_decompose.add_argument("-T", "--threads", type=int, dest="num_thread", default=0,
                            help="number of cores used in placement. "
                                 "0 to use all cores in the running machine", metavar="NUMBER")
    parser_decompose.add_argument("-l", "--overlap", type=int, dest="overlap_length", default=50,
                            help="minimum alignment overlap length needed to use the subalignment"
                                 "in subtree refinement", metavar="NUMBER")
    parser_decompose.add_argument("-c", "--constrain-outgroups", dest="constrain_outgroups", action='store_true',
                            default=False,
                            help="enforce outgroup topology on gene tree estimation stage.")
    parser_decompose.add_argument("-n", "--numtasks", type=int, dest="num_tasks", metavar='NUMBER', default=1,
                            help="number of tasks where local refinement jobs will be split.")
    parser_decompose.add_argument("-C", "--occupancy", type=float, dest="occupancy_threshold", default=0.66,
                            help="minimum fraction of species needed to call a gene. "
                                 "Must be a value between 0 and 1. highly occupant.")
    parser_decompose.add_argument("-e", "--edge-threshold", type=float, dest="edge_threshold", default=0.02,
                            help="maximun edge length in a cluster.")
    parser_decompose.add_argument("-m", "--method", dest="method", choices=['raxml-ng', 'iqtree', 'raxml-8'],
                                  default=False, help="method for subtree inference.")

    parser_decompose.set_defaults(func=decompose)

    # refine command subparser
    parser_ref = subparsers.add_parser('refine',
                                       description='Refine partitions via consensus (ASTRAL)')
    parser_ref.add_argument("-o", "--output", dest="output_fp",
                            help="path for the output directory where files will be placed",
                            metavar="DIRECTORY")
    parser_ref.add_argument("-T", "--threads", type=int, dest="num_thread", default=0,
                            help="number of cores used in placement. "
                                 "0 to use all cores in the running machine", metavar="NUMBER")
    parser_ref.add_argument("-m", "--method", dest="method", choices=['raxml-ng', 'iqtree', 'raxml-8'],
                                  default=False, help="method for subtree inference.")

    parser_ref.set_defaults(func=refine)

    # stitch command subparser
    parser_sti = subparsers.add_parser('stitch',
                                       description='Stitch back locally refined trees')
    parser_sti.add_argument("-o", "--output", dest="output_fp",
                            help="path for the output directory where files will be placed",
                            metavar="DIRECTORY")
    parser_sti.add_argument("-T", "--threads", type=int, dest="num_thread", default=0,
                            help="number of cores used in placement. "
                                 "0 to use all cores in the running machine", metavar="NUMBER")

    parser_sti.set_defaults(func=stitch)

    options = parser.parse_args()

    if options.num_thread == 0:
        options.num_thread = cpu_count()
    options.output_fp = abspath(expanduser(options.output_fp))

    return options
