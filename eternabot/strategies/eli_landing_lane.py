from eterna_utils import *
import inv_utils
import strategy_template

class Strategy(strategy_template.Strategy):
    def __init__(self):
        
        strategy_template.Strategy.__init__(self)
        
        self.title_ = "Landing Lane"
        self.author_ = "Eli Fisker"
        self.url_ = "NA"
        #self.default_params_ = [2,1,1,2,1,1,0.5, 0.55, 100]
        #self.code_length_ = 40
        #self.publishable_ = True
        #self.denormalized_ = True
        #self.comprehensive_ = True
        #self.martin_weight_ = 4.361964
        #self.satisfying_point_ = 80

    def score(self, designs):            
        score = 0

        lane1 = [60,76]
        MS2 = [97,116]
        oligo1 = [0,22]
        oligo2 = [26,48]

        pairmaps = []
        for design in designs:
            pairmaps.append(get_pairmap_from_secstruct(design['secstruct']))

        for i in range(lane1[0],lane1[1]):
            if pairmaps[0][i] in range(MS2[0],MS2[1]):
                score += 1
            if pairmaps[1][i] in range(oligo2[0],oligo2[1]):
                score += 1
            if pairmaps[2][i] in range(oligo1[0],oligo1[1]):
                score += 1
            if pairmaps[3][i] in range(oligo1[0],oligo2[1]):
                score += 1
        for i in range(MS2[0],MS2[1]):
            if pairmaps[3][i] in range(oligo1[0],oligo2[1]):
                score += 1
        
        return float(score) / ((lane1[1]-lane1[0])*4 + MS2[1]-MS2[0])
         
