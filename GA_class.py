import math
import random

class individual():
    def __init__(self,arg_gene,arg_chromlength,arg_pMutate):
        """
        :param arg_gene: 基因编码
        :param arg_chromlength: 基因片段的长度
        :param arg_pMutate: 基因突变的概率
        """
        self.gene=arg_gene
        self.chromlength = arg_chromlength
        self.pMutate = arg_pMutate

    def b2d(self,bin):
        t = 0
        for i in range(len(bin)):
            t += bin[i] * math.pow(2, i)
        return t

    def objfun(self,x_list):
        x=self.b2d(x_list)*10/1023
        return 100 - (x - 2) * (x - 7)

    def fun_value(self):
        return self.objfun(self.gene)

    @property
    def get_gene(self):
        return self.gene

    @property
    def set_gene(self,arg_gene):
        self.gene=arg_gene

    def mutate(self):
        if random.random()<self.pMutate:
            point=random.randint(0,self.chromlength-1)
            if self.gene[point]==0:
                self.gene[point] =1
            else:
                self.gene[point]=0



class Population():
    def __init__(self,arg_number,arg_pCross=0.6,arg_chromlength=10,arg_pMutate=0.001):
        """
        :param arg_number: 初始群体数目
        :param arg_pCross: 交换率
        :param arg_chromlength: 基因片段的长度
        :param arg_pMutate: 基因突变的概率
        """
        self.number=arg_number
        self.pCross=arg_pCross
        self.chromlength = arg_chromlength
        self.pMutate = arg_pMutate
        self.fun_value=[]

        self.init_gene=[random.randint(0, 1) for i in range(self.chromlength)]
        self.init_gene
        self.individual_list=[individual(self.init_gene,self.chromlength,self.pMutate) for i in range(self.number)]


    def __exchange__(self, indiv_1,indiv_2):
        exchange_point=random.randint(0, self.chromlength)
        gene1 = indiv_1.get_gene
        gene2 = indiv_2.get_gene
        indiv_1.gene=gene1[0:exchange_point] + gene2[exchange_point:self.chromlength]
        indiv_2.gene=gene2[0:exchange_point] + gene1[exchange_point:self.chromlength]

    def cross(self):
        for i in range(self.number):
            if random.random()<self.pCross:
                self.__exchange__(self.individual_list[i],self.individual_list[random.randint(0, self.number-1)])

    def mutate(self):
        for i in range(self.number):
            self.individual_list[i].mutate()

    def __get_pfitvalue__(self):
        self.fun_value=[self.individual_list[i].fun_value() for i in range(self.number)]
        fixvalue=[self.fun_value[i] for i in range(self.number)]
        total=sum(fixvalue)
        return [fixvalue[i]/total for i in range(self.number)]

    def __cumsum__(self,arg_fitvalue):
        for i in range(len(arg_fitvalue)):
            if i == 0:
                arg_fitvalue[i] = arg_fitvalue[0]
            else:
                arg_fitvalue[i] = arg_fitvalue[i] + arg_fitvalue[i - 1]

    def select(self):
        p_fitvalue=self.__get_pfitvalue__()
        self.__cumsum__(p_fitvalue)

        p_select=[random.random() for i in range(self.number)]
        p_select.sort()

        select_index=0
        new_indiv_index=0
        new_individual_list=[]

        while(select_index<self.number):
            if p_select[select_index]<p_fitvalue[new_indiv_index]:
                new_individual_list.append(self.individual_list[new_indiv_index])
                select_index+=1
            else:
                new_indiv_index+=1
        self.individual_list=new_individual_list

    @property
    def get_best(self):
        bestfunvalue=self.fun_value[0]
        bestone=self.individual_list[0]

        for i in range(self.number):
            if bestfunvalue<=self.fun_value[i]:
                bestfunvalue=self.fun_value[i]
                bestone=self.individual_list[i]

        return bestfunvalue

if __name__=="__main__":

    pop=Population(500)

    for i in range(1000):
        pop.select()
        pop.cross()
        pop.mutate()
    bestfunvalue=pop.get_best
    print(bestfunvalue)