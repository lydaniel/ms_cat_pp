#!/usr/bin/python
# import libraries
import sys
sys.path.append('../oed')

import oed
import pmf
import church

import copy
import itertools
import sys
import pp
import datetime
import numpy
import scipy
import scipy.stats
import scipy.misc
import math
import os
import random

####
def random_comb(lst, num, rand_size) :
    #check limit
    #print lst
    if (scipy.misc.comb(len(lst), num) < rand_size) :
        return list(itertools.combinations(lst, num))
    else :
        flst = []
        while (len(flst) < rand_size) :
            #tlst = [random.randint(0, len(lst)) for i in ]
            tlst = []
            while (len(tlst) < num) :
                rr = random.randint(0, len(lst)-1)
                if rr not in tlst :
                    tlst.append(rr)
            tlst.sort()
            if tlst not in flst :
                flst.append(tuple(tlst))

        flst.sort()
                
        return flst

####

part_list_limit = 2**15
#part_list_limit = 950

#a  = random_comb(range(46), 2, part_list_limit)
#print a
#print len(a)

#sys.exit(0)

models = pmf.PMF(2)
models[0] = pmf.P("independent_cue.church", 1.0)
models[1] = pmf.P("context_theory.church", 1.0)
models.normalize()

query = "raw"
church_exec_path = "/opt/webchurch/church"
#church_exec_path = "/home/danly/webchurch/church"

t_outputs = [pmf.PMF() for i in range(len(models))]

data_dir = "data/gen_mturk_r7"
kl_dir = "data/kl_mturk_r7"

for filename in os.listdir(data_dir) :
    data_file = open(data_dir + "/" + filename, "r")
    out_file = open(kl_dir + "/" + filename, "w")

    for (i, line) in enumerate(data_file) :
        if (i == 0) :
            partcpnt_data = [{} for i in range(int(line))]
        elif (i == 1) :
            inputs = line
        else :
            stimulus = line.split("|")[0]
            temp_line = line.split("|")[1].strip("\n").split(",")
            for tl in range(len(temp_line)) :
                partcpnt_data[tl][stimulus] = temp_line[tl]

    inputs_list = []
    tinputs_list = inputs.strip("\n").split("|")
    for til in range(len(tinputs_list)) :
        tinputs_list[til] = tinputs_list[til].split(";")
        for ttil in tinputs_list[til] :
            inputs_list.append(ttil)

    for m in range(len(models)) :
        church.exec_model(models[m].x, inputs, t_outputs[m], query, church_exec_path)

    output_pmf = pmf.PMF(2)
    output_pmf[0] = pmf.P("A", 0.0)
    output_pmf[1] = pmf.P("B", 0.0)

    b_outputs = [0]*16
    for s in range(len(b_outputs)) :
        b_outputs[s] = [0]*len(models)
        for m in range(len(b_outputs[s])) :
            b_outputs[s][m] = copy.deepcopy(output_pmf)

    for s in range(len(b_outputs)) :
        for m in range(len(b_outputs[s])) :
            b_outputs[s][m][0].p = t_outputs[m][s].p
            b_outputs[s][m][1].p = 1 - t_outputs[m][s].p

    #print inputs_list
    #print b_outputs

    for np in range(1,len(partcpnt_data)+1) :
        #part_list = list(itertools.combinations(range(len(partcpnt_data)), np))
        #part_list = itertools.combinations(range(len(partcpnt_data)), np)
        part_list = random_comb(range(len(partcpnt_data)), np, part_list_limit)
        #print part_list
        #print np, len(part_list)
        #sys.stdout.flush()
        #for a in len(part_list) :
        #    print a
        #    sys.stdout.flush()
        #sys.exit(0)
        #if (len(part_list) < part_list_limit) :
        #    rpart_list = list(part_list)
        #else :
        #    rpart_list = list(random.sample(part_list, part_list_limit))
        #print len(part_list), len(rpart_list)
        #sys.stdout.flush()
        #print np, len(part_list)
        outputs = [[0 for i in range(len(models))] for j in range(len(part_list))]
        log_outputs = [[0 for i in range(len(models))] for j in range(len(part_list))]
        for pl in range(len(part_list)) :
            #print part_list[pl], outputs[pl]
            for m in range(len(models)) :
                for pll in part_list[pl]:
                    #print "  ", "pll: ", pll, "m: ", m

                    for il in inputs_list :
                        index = inputs_list.index(il)
                        #print il, index

                        if (partcpnt_data[pll][il]  == 'A') :
                            log_outputs[pl][m] += math.log(b_outputs[index][m][0].p)
                            #print "    ", b_outputs[index][m][0].p, math.log(b_outputs[index][m][0].p), log_outputs[pl][m]
                        else :
                            log_outputs[pl][m] += math.log(b_outputs[index][m][1].p)
                            #print "    ", b_outputs[index][m][1].p, math.log(b_outputs[index][m][1].p), log_outputs[pl][m]

        for pl in range(len(part_list)) :
            for m in range(len(models)) :
                outputs[pl][m] = math.exp(log_outputs[pl][m])

        #print log_outputs
        #print outputs

        foutputs = [[0 for i in range(len(models))] for j in range(len(part_list))]
        for pl in range(len(part_list)) :
            soutputs = sum(outputs[pl])
            for m in range(len(models)) :
                foutputs[pl][m] = outputs[pl][m]/soutputs
        #print foutputs

        kl = [0 for j in range(len(part_list))]
        for pl in range(len(part_list)) :
            for m in range(len(models)) :
                kl[pl] += foutputs[pl][m]*math.log(foutputs[pl][m]/models[m].p)
        #print kl
        #sys.stdout.flush()

        kl_mean = sum(kl)/len(kl)

        #print kl_mean

        if (len(kl) != 1) :
            kl_var = 0
            for l in range(len(kl)) :
                kl_var += (kl[l] - kl_mean)**2
            kl_var = kl_var/(len(kl)-1)

            kl_se = math.sqrt(kl_var/len(kl))

            #print np, kl_mean, kl_se
            out_file.write(str(np) + ", " + str(kl_mean) + ", " + str(kl_se) + "\n")
        else :
            #print np, kl_mean, "undef"
            out_file.write(str(np) + ", " + str(kl_mean) + ", " + "undef" + "\n")

        out_file.flush()

    data_file.close()
    out_file.close()



