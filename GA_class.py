import math
import random
import copy


def cumsum(arg_list):
    for i in range(len(arg_list)):
        if i == 0:
            arg_list[i] = arg_list[0]
        else:
            arg_list[i] = arg_list[i] + arg_list[i - 1]


def decode(bin_list):
    t = 0
    for index in range(len(bin_list)):
        t += bin_list[index] * math.pow(2, index)
    return t * 5 / 1023


def my_fit(x):
    #return 10 * math.sin(5 * x) + 7 * math.cos(4 * x)
    return 5-(x-1)*(x-5)

def get_pfitvalue(fun_value):
    temp = 0
    fixvalue = []
    number = len(fun_value)
    for i in range(number):
        if fun_value[i] > 0:
            temp = fun_value[i]
        else:
            temp = 0
        fixvalue.append(temp)
    total = sum(fixvalue)
    return [fixvalue[i] / total for i in range(number)]


def exchange(indiv_1, indiv_2):
    """
    :type indiv_1: individual
    :type indiv_2: individual
    """
    chromlength = len(indiv_1.gene)
    exchange_point = random.randint(0, chromlength - 1)
    gene1 = indiv_1.gene[:]
    gene2 = indiv_2.gene[:]
    indiv_1.gene = gene1[0:exchange_point] + gene2[exchange_point:]
    indiv_2.gene = gene2[0:exchange_point] + gene1[exchange_point:]


class individual:
    def __init__(self, arg_chromlength: int, arg_pMutate: float):
        """
        :param arg_chromlength: 基因片段的长度
        :param arg_pMutate: 基因突变的概率
        """
        self.chromlength = arg_chromlength
        self.gene = [random.randint(0, 1) for i in range(self.chromlength)]
        # self.gene=[0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        self.pMutate = arg_pMutate

    def mutate(self):
        if random.random() < self.pMutate:
            point = random.randint(0, self.chromlength - 1)
            if self.gene[point] == 0:
                self.gene[point] = 1
            else:
                self.gene[point] = 0


class Population:
    def __init__(self, arg_number, arg_pCross=0.6, arg_chromlength=10, arg_pMutate=0.001):
        """
        :param arg_number: 初始群体数目
        :param arg_pCross: 交换率
        :param arg_chromlength: 基因片段的长度
        :param arg_pMutate: 基因突变的概率
        """

        self.number = arg_number
        self.pCross = arg_pCross
        self.chromlength = arg_chromlength

        self.fun_value = []
        self.individual_list = []
        for i in range(self.number):
            self.individual_list.append(individual(self.chromlength, arg_pMutate))

        self.bestgene = []
        self.bestfunvalue = 0

    def cross(self):
        for index in range(self.number - 1):
            if random.random() < self.pCross:
                exchange(self.individual_list[index], self.individual_list[index + 1])

    def mutate(self):
        for i in range(self.number):
            self.individual_list[i].mutate()

    def get_fun_value(self):
        self.fun_value = [my_fit(decode(self.individual_list[i].gene[:])) for i in range(self.number)]
        return self.fun_value[:]

    def select(self):
        p_fitvalue = get_pfitvalue(self.get_fun_value())
        cumsum(p_fitvalue)

        p_select = [random.random() for i in range(self.number)]
        p_select.sort()

        select_index = 0
        new_indiv_index = 0
        new_individual_list = []

        while select_index < self.number:
            if p_select[select_index] < p_fitvalue[new_indiv_index]:
                new_individual_list.append(copy.deepcopy(self.individual_list[new_indiv_index]))
                select_index += 1
            else:
                new_indiv_index += 1
        self.individual_list = new_individual_list[:]

    def best(self):
        self.fun_value = [my_fit(decode(self.individual_list[i].gene[:])) for i in range(self.number)]
        self.bestfunvalue = max(self.fun_value)
        self.best_index = self.fun_value.index(self.bestfunvalue)
        self.bestgene = self.individual_list[self.best_index].gene[:]

if __name__ == "__main__":

    pop = Population(50)
    bestgene = []
    best_value = []
    for i in range(100):
        pop.select()

        pop.best()
        bestgene.append(pop.bestgene[:])
        best_value.append(pop.bestfunvalue)
        pop.cross()
        pop.mutate()

    best_index = best_value.index(max(best_value))
    best = bestgene[best_index]

    print(decode(best))
