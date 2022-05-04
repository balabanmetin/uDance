from abc import ABC, abstractmethod
from os.path import join
import json


class StitchStrat(ABC):

    def __init__(self, suffix):
        self.suffix = suffix

    @abstractmethod
    def get_astral_treename(self, output_dir, cluster):
        pass

    def get_suffix(self):
        return self.suffix


class IncrementalStrat(StitchStrat):

    def __init__(self):
        super(IncrementalStrat, self).__init__("incremental")

    def get_astral_treename(self, output_dir, cluster):
        return join(output_dir, cluster, "astral_output.%s.nwk" % self.suffix)


class UpdatesStrat(StitchStrat):

    def __init__(self):
        super(UpdatesStrat, self).__init__("updates")

    def get_astral_treename(self, output_dir, cluster):
        return join(output_dir, cluster, "astral_output.%s.nwk" % self.suffix)


class MaxqsStrat(StitchStrat):

    def __init__(self):
        super(MaxqsStrat, self).__init__("maxqs")

    def get_astral_treename(self, output_dir, cluster):
        qsdict = {}
        for s in ["incremental", "updates"]:
            logfile = join(output_dir, cluster, "astral.%s.log" % s)
            with open(logfile, "r") as f:
                for lines in f.readlines():
                    if lines.startswith("Final quartet score"):
                        qsdict[s] = int(lines.strip().split(" ")[4])
                        break
        with open(join(output_dir, cluster, "quartet_scores.json"), "w") as f:
            json.dump(qsdict, f)
        if qsdict["incremental"] > qsdict["updates"]:
            return join(output_dir, cluster, "astral_output.%s.nwk" % "incremental")
        else:
            return join(output_dir, cluster, "astral_output.%s.nwk" % "updates")


def strategy_dealer():
    return [IncrementalStrat(), UpdatesStrat(), MaxqsStrat()]
