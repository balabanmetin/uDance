import argparse
import sys
from multiprocessing import cpu_count
from os.path import abspath, expanduser

from uDance.decompose import decompose
from uDance.mainlines import mainlines
from uDance.refine import refine
from uDance.stitch import stitch


def options_config():
    # Input arguments parser
    parser = argparse.ArgumentParser(description='Massively scalable divide-and-conquer phylogenetic inference',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--version', action='store_true', help='print the current version')
    # parser.add_argument('--debug', action='store_true', help='Print the traceback when an exception is raised')
    subparsers = parser.add_subparsers(title='commands',
                                       description='mainlines   Strategically choose backbone taxa\n'
                                                   'decompose   Create local inference (RAxML) tasks\n'
                                                   'refine      Refine partitions via consensus (ASTRAL) \n'
                                                   'stitch      Stitch back locally refined trees',
                                       help='Run run_udance.py {commands} [-h] for additional help',
                                       dest='{commands}')

    # To make sure that subcommand is required in python >= 3.3
    python_version = sys.version_info
    if (python_version[0] * 10 + python_version[1]) >= 33:
        subparsers.required = True

    # mainlines command subparser
    parser_mainlines = subparsers.add_parser('mainlines',
                                             description='Strategically choose backbone taxa')
    parser_mainlines.add_argument("-s", "--alignment-dir", dest="alignment_dir_fp",
                                  help="path for input directory which contains "
                                       "extended reference alignment files (FASTA), "
                                       "containing reference and query sequences.",
                                  metavar="DIRECTORY")
    parser_mainlines.add_argument("-n", "--number", type=int, dest="target_num", metavar='NUMBER', default=1000,
                                  help="number of taxa to be selected.")
    parser_mainlines.add_argument("-p", "--protein", dest="protein_seqs", action='store_true', default=False,
                                  help="input sequences are protein sequences")
    parser_mainlines.add_argument("-l", "--length", type=int, dest="concat_length", default=5000,
                                  help="number of sites in the concatenated alignment created. "
                                       "This value is proportional with running time and memory use.", metavar="NUMBER")
    parser_mainlines.add_argument("-g", "--gap-threshold", type=float, dest="gap_threshold", default=0.95,
                                  help="Alignment filtering threshold. "
                                       "Sites with a gappiness value larger than 1-gap_threshold will be removed.")
    parser_mainlines.set_defaults(func=mainlines)

    # decompose command subparser
    parser_decompose = subparsers.add_parser('decompose',
                                             description='Create local refinement tasks')
    parser_decompose.add_argument("-j", "--jplace", dest="jplace_fp",
                                  help="path to the jplace placement file", metavar="FILE")
    parser_decompose.add_argument("-t", "--threshold", dest="threshold", default="1000",
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
    parser_decompose.add_argument("-l", "--sublength", type=int, dest="subalignment_length", default=100,
                                  help="minimum alignment length needed to use the subalignment"
                                       , metavar="NUMBER")
    parser_decompose.add_argument("-f", "--fraglength", type=int, dest="fragment_length", default=75,
                                  help="minimum seqence(fragment) length needed to use in the subalignment"
                                       , metavar="NUMBER")
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
    parser_decompose.add_argument("-m", "--method", dest="method", choices=['raxml-ng', 'iqtree', 'raxml-8', 'copy'],
                                  default='raxml-8', help="method for subtree inference.")

    parser_decompose.set_defaults(func=decompose)

    # refine command subparser
    parser_ref = subparsers.add_parser('refine',
                                       description='Refine partitions via consensus (ASTRAL)')
    parser_ref.add_argument("-p", "--partition", dest="partition_dir",
                            help="path for the directory of the partition to be refined. ",
                            metavar="DIRECTORY")
    parser_ref.add_argument("-T", "--threads", type=int, dest="num_thread", default=0,
                            help="number of cores used in the refinement. "
                                 "0 to use all cores in the running machine", metavar="NUMBER")
    parser_ref.add_argument("-M", "--memory", type=int, dest="memory", default=1000,
                            help="memory used in the refinement (MB). "
                                 "0 to use all cores in the running machine", metavar="NUMBER")
    parser_ref.add_argument("-m", "--method", dest="method", choices=['raxml-ng', 'iqtree', 'raxml-8'],
                            default=False, help="method for subtree inference.")
    parser_ref.add_argument("-g", "--use-gpu", dest="use_gpu", action='store_true',
                            default=False,
                            help="disable gpu usage (currently defuct)")
    parser_ref.add_argument("-c", "--contract", type=float, dest="contract_threshold", default=0.9,
                            help="contract branches with support less than given threshold"
                                 "in the inferred gene trees", metavar="NUMBER")
    parser_ref.add_argument("-o", "--occupancy", type=int, dest="occupancy_threshold", default=1,
                                  help="gene occupancy threshold for inclusion in ASTRAL step."
                                       , metavar="NUMBER")
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

    if hasattr(options, "num_thread") and options.num_thread == 0:
        options.num_thread = cpu_count()
    if hasattr(options, "output_fp"):
        options.output_fp = abspath(expanduser(options.output_fp))
    if hasattr(options, "cluster_dir"):
        options.cluster_dir = abspath(expanduser(options.cluster_dir))

    return options
