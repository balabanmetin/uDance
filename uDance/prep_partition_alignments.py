import sys
import multiprocessing as mp

from uDance.PoolAlignmentWorker import PoolAlignmentWorker
from uDance.fasta2dic import fasta2dic

from os import listdir
from os.path import isfile, join, splitext


def prep_partition_alignments(alndir, protein_flag, species_path_list, num_thread, overlap):
    only_files = [f for f in listdir(alndir) if isfile(join(alndir, f))]
    all_scripts = []
    for aln in only_files:
        aln_input_file = join(alndir, aln)
        basename = splitext(aln)[0]
        # try:
        fa_dict = fasta2dic(aln_input_file, protein_flag, False)
        alignment_worker = PoolAlignmentWorker()
        alignment_worker.set_class_attributes(overlap, fa_dict, basename)
        pool = mp.Pool(num_thread)
        scripts = pool.map(alignment_worker.worker, species_path_list)
        pool.close()
        pool.join()
        valid_scripts = [s for s in scripts if s is not None]
        all_scripts += valid_scripts
        # main_script.write("\n".join(valid_scripts))
        # main_script.write("\n")
        # except Exception as e:
        #     print(e)
        #     print("Alignment %s is not a valid fasta alignment" % aln_input_file, file=sys.stderr)
    return all_scripts
