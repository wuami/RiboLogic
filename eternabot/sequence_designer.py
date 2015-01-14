import ensemble_utils
import ensemble_design

class SequenceDesigner:

    def __init__(self, id, secstruct, constraints):
        self.id = id
        self.secstruct = secstruct
        if constraints == None:
            self.constraints = 'N'*len(secstruct)
        else:
            self.constraints = constraints
        self.score_cutoff = 90
        if len(secstruct) <= 50:
            self.score_cutoff = 70
        elif len(secstruct) <= 80:
            self.score_cutoff = 80
        strategy_names = ['example_gc60', 'penguian_clean_dotplot', 'berex_simplified_berex_test']
        self.ensemble = ensemble_utils.Ensemble("conventional", strategy_names, None)
        self.solution = []

    def optimize_sequence(self):
        res = ensemble_design.inverse_fold_whole(self.secstruct, self.constraints, self.ensemble.score, self.score_cutoff, "conventional")
        self.solution = res['end']
        return

    def get_solution(self):
        return self.solution

    def check_current_secstructs(self):
        return True
