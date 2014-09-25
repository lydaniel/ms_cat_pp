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

####

class Result:
        def __init__(self, inp, outp):
                self.inp = inp 
                self.outp = outp
        def __repr__(self):
                return repr((self.inp, self.outp))

####

def check_cats(cat_a, cat_b) :
    for ca in range(len(cat_a[0])) :
        a_count = 0
        b_count = 0
        for i in range(len(cat_a)) :
            if (cat_a[i][ca] == '1') :
                a_count += 1
            else :
                b_count += 1

        if (a_count <= b_count) :
            #print "failed a", ca, a_count, b_count
            return False

    for cb in range(len(cat_b[0])) :
        a_count = 0
        b_count = 0
        for i in range(len(cat_b)) :
            if (cat_b[i][cb] == '0') :
                b_count += 1
            else :
                a_count += 1

        if (b_count <= a_count) :
            #print "failed b", cb, a_count, b_count
            return False

    return True

####

in_sequence_raw = ["".join(seq) for seq in itertools.product("01", repeat=4)]
in_sequence = [0]*len(in_sequence_raw)
for i in range(len(in_sequence_raw)) :
    in_sequence[i] = str(in_sequence_raw[i]).replace("", ",")[1:-1]

num_sequence = range(len(in_sequence))

inputs_list = []
joint_cat_size_list = []
#for a in range(3,6):
#    for b in range(3,6):
a = 5
b = 4
a_comb = itertools.combinations(num_sequence, a)
for ac in a_comb :
    num_sequence_b = copy.deepcopy(num_sequence)
    for aac in ac :
        num_sequence_b.remove(aac)
    b_comb = itertools.combinations(num_sequence_b, b)
    for bc in b_comb :
        #print ac, bc
        a_cat_raw = []
        a_cat = []
        for aac in ac :
            a_cat_raw.append(in_sequence_raw[aac])
            a_cat.append(in_sequence[aac])
        b_cat_raw = []
        b_cat = []
        for bbc in bc :
            b_cat_raw.append(in_sequence_raw[bbc])
            b_cat.append(in_sequence[bbc])
            
        if (check_cats(a_cat_raw, b_cat_raw)) :
            a_string = ""
            t_joint_cat_size = 0
            for acs in a_cat :
                a_string += acs + ";"
                t_joint_cat_size += 1
            a_string = a_string[:-1]
            
            b_string = ""
            for bcs in b_cat :
                b_string += bcs + ";"
                t_joint_cat_size += 1
            b_string = b_string[:-1]

            inputs_list.append(a_string + "|" + b_string)
            joint_cat_size_list.append(t_joint_cat_size)
            #print a_cat, b_cat

#for (i,ilk) in enumerate(inputs_list) :
#    print ilk, joint_cat_size_list[i]

print "runs = ", len(inputs_list)

samples = 250

print "samples = ", samples
print ""

models = pmf.PMF(2)
models[0] = pmf.P("independent_cue.church", 1.0)
models[1] = pmf.P("context_theory.church", 1.0)
models.normalize()

query = "rejection-query"
church_exec_path = "/opt/webchurch/church"

results = []

print "################################"

for il in range(len(inputs_list)) :
#for il in range(20) :
    joint_cat_size = joint_cat_size_list[il]
    bit_sequence_raw = ["".join(seq) for seq in itertools.product("01", repeat=joint_cat_size)]
    bit_sequence = [0]*len(bit_sequence_raw)
    for i in range(len(bit_sequence_raw)) :
        bit_sequence[i] = str(bit_sequence_raw[i]).replace("", " ")[1:-1]

    t_output_pmf = pmf.PMF(len(bit_sequence))
    for i in range(len(bit_sequence_raw)) :
        t_output_pmf[i] = pmf.P(bit_sequence[i], 1.0)
    t_output_pmf.normalize() 

    t_outputs = [0]*len(models)
    for m in range(len(models)) :
        t_outputs[m] = copy.deepcopy(t_output_pmf)

    #inputs = "1,1,1,0;1,0,1,0;0,1,1,1|1,1,0,0;0,0,0,0"
    #inputs = "1,1,1,0;1,0,1,0;1,0,1,1;1,1,0,1;0,1,1,1|1,1,0,0;0,1,1,0;0,0,0,1;0,0,0,0"
    inputs = inputs_list[i]

    print "inputs = ", inputs

    for m in range(len(models)) :
        church.exec_model(models[m].x, inputs, t_outputs[m], query, \
                          church_exec_path, church.Type._string)

    for m in range(len(models)) :
        for o in range(len(t_outputs[m])) :
            t_outputs[m][o].p = int(t_outputs[m][o].p*samples)

    print "t_outputs = ", t_outputs

    output_pmf = pmf.PMF(2)
    output_pmf[0] = pmf.P("A", 0.0)
    output_pmf[1] = pmf.P("B", 0.0)

    outputs = [0]*joint_cat_size
    for s in range(len(outputs)) :
        outputs[s] = [0]*len(models)
        for m in range(len(outputs[s])) :
            outputs[s][m] = copy.deepcopy(output_pmf)

    for s in range(len(outputs)) :
        for m in range(len(outputs[s])) :
            for o in range(len(t_outputs[m])) :
                if (t_outputs[m][o].p != 0) :
                    if (t_outputs[m][o].x.replace(" ","")[s] == '0') :
                        outputs[s][m][0].p += t_outputs[m][o].p
                    else :
                        outputs[s][m][1].p += t_outputs[m][o].p

    for s in range(len(outputs)) :
        for m in range(len(outputs[s])) :
            outputs[s][m].normalize()

    print "outputs  = ", outputs

    expected_kl = 0
    for s in range(len(outputs)) :
        t_expected_kl = oed.get_expected_kl(models, outputs[s])
        expected_kl += t_expected_kl
        print s, ", t_ekl =", t_expected_kl, ", total ekl =", expected_kl

    temp_result = Result(inputs, expected_kl)
    print "result =", temp_result

    results.append(temp_result)

    sys.stdout.flush()

    print "################################"

sorted_results = sorted(results, key=lambda result: result.outp, reverse=True)
for s in sorted_results :
    print s

