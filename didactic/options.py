from optparse import OptionParser
from multiprocessing import cpu_count
from os.path import abspath, expanduser
from sys import stderr


def options_config():
    parser = OptionParser()
    parser.add_option("-j", "--jplace", dest="jplace_fp",
                      help="path to the jplace placement file", metavar="FILE")
    parser.add_option("-f", "--threshold", dest="threshold", default="600",
                      help="maximum number of elements in each cluster")
    parser.add_option("-o", "--output", dest="output_fp",
                      help="path for the output directory where files will be placed",
                      metavar="DIRECTORY")
    parser.add_option("-s", "--alignment-dir", dest="alignment_dir_fp",
                      help="path for input directory which contains "
                           "extended reference alignment files (FASTA), "
                           "containing reference and query sequences.",
                      metavar="FILE")
    parser.add_option("-p", "--protein", dest="protein_seqs", action='store_true', default=False,
                      help="input sequences are protein sequences")
    parser.add_option("-T", "--threads", dest="num_thread", default="0",
                      help="number of cores used in placement. "
                           "0 to use all cores in the running machine", metavar="NUMBER")
    parser.add_option("-l", "--overlap", dest="overlap_length", default="50",
                      help="minimum alignment overlap length needed to use the subalignment"
                           "in subtree refinement", metavar="NUMBER")
    parser.add_option("-c", "--constrain-outgroups", dest="constrain_outgroups", action='store_true', default=False,
                      help="enforce outgroup topology on gene tree estimation stage.")
    parser.add_option("-n", "--numtasks", dest="num_tasks", metavar='NUMBER', default="1",
                      help="number of tasks where local refinement jobs will be split.")

    options, args = parser.parse_args()

    options.num_thread = int(options.num_thread)
    options.overlap_length = int(options.overlap_length)
    options.num_tasks = int(options.num_tasks)
    if not options.num_thread:
        options.num_thread = cpu_count()
    options.output_fp = abspath(expanduser(options.output_fp))
    if options.num_tasks < 1:
        stderr.write("Invalid number of tasks. Number of tasks is set to the minimum value: 1.\n")
        options.num_tasks = 1
    return (options, args)