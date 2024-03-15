from abc import ABC, abstractmethod
from os.path import join
import json


class StitchStrat(ABC):
    def __init__(self, suffix, blen):
        self.suffix = suffix
        if blen:
            self.blen = '.bl'
        else:
            self.blen = ''

    @abstractmethod
    def get_astral_treename(self, output_dir, cluster):
        pass

    def get_suffix(self):
        return self.suffix


class IncrementalStrat(StitchStrat):
    def __init__(self, blen):
        super(IncrementalStrat, self).__init__('incremental', blen)

    def get_astral_treename(self, output_dir, cluster):
        return join(output_dir, cluster, 'astral_output.%s.nwk%s' % (self.suffix, self.blen))


class UpdatesStrat(StitchStrat):
    def __init__(self, blen):
        super(UpdatesStrat, self).__init__('updates', blen)

    def get_astral_treename(self, output_dir, cluster):
        return join(output_dir, cluster, 'astral_output.%s.nwk%s' % (self.suffix, self.blen))


class MaxqsStrat(StitchStrat):
    def __init__(self, blen):
        super(MaxqsStrat, self).__init__('maxqs', blen)

    def get_astral_treename(self, output_dir, cluster):
        qsdict = {}
        for s in ['incremental', 'updates']:
            logfile = join(output_dir, cluster, 'astral.%s.log' % s)
            with open(logfile, 'r') as f:
                for lines in f.readlines():
                    if lines.startswith('Final quartet score'):
                        qsdict[s] = int(lines.strip().split(' ')[4])
                        break
        with open(join(output_dir, cluster, 'quartet_scores.json'), 'w') as f:
            json.dump(qsdict, f)
        if qsdict['incremental'] > qsdict['updates']:
            return join(output_dir, cluster, 'astral_output.%s.nwk%s' % ('incremental', self.blen))
        else:
            return join(output_dir, cluster, 'astral_output.%s.nwk%s' % ('updates', self.blen))


def strategy_dealer(blen):
    return [IncrementalStrat(blen), UpdatesStrat(blen), MaxqsStrat(blen)]
