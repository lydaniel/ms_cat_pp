#!/usr/bin/python
# import libraries
import sys
sys.path.append('../oed')

import oed
import pmf
import church
import lib_result

import copy
import itertools
import sys
import pp
import datetime
import numpy
import scipy
import scipy.stats
import math

####

def check_cats(cat_a, cat_b) :
    #print cat_a, cat_b
    for ca in range(len(cat_a[0])) :
        a_mode = 0
        b_mode = 0
        right_count = 0
        wrong_count = 0

        for i in range(len(cat_a)) :
            if (cat_a[i][ca] == '1') :
                a_mode += 1
                right_count += 1
            else :
                a_mode -= 1
                wrong_count += 1

        if (a_mode < 0) :
            #print "    a mode failed:", a_mode, ca
            return False

        for i in range(len(cat_b)) :
            if (cat_b[i][ca] == '0') :
                b_mode += 1
                right_count += 1
            else :
                b_mode -= 1
                wrong_count += 1

        if (b_mode < 0) :
            #print "    b mode failed:", b_mode, ca
            return False

        if ((a_mode == 0) or (b_mode == 0)) :
            #print "    count:", right_count, wrong_count, ca
            if (right_count < wrong_count) :
                #print "    count failed:", right_count - wrong_count, ca
                return False

    return True

def list_share_count(ac, bc, at, bt) :
    diff = 0
    for a in ac :
        if a not in at :
            diff += 1
    for b in bc :
        if b not in bt :
            diff += 1
    return diff

def mask(v1, v2) :
    x = [1]*len(v1)
    for i in range(len(v1)) :
        if (v1[i] == v2[i]) :
            x[i] = 0
    return x

def does_intersect(a1, a2, b1, b2) :
    ma = mask(a1, a2)
    mb = mask(b1, b2)
    #print "        ci:", ma, mb, sum(ma), sum(mb)
    if (sum(ma) == sum(mb)) :
        for i in range(len(ma)) :
            if (ma[i] != mb[i]) :
                return False

        #print "        di: true"
        for i in range(len(ma)) :
            if (ma[i] == 0) :
                if (a1[i] != b1[i]) :
                    return False
        return True
    else :
        return False

def make_key(a_cat, b_cat) :
    string = ""
    for acs in a_cat :
        string += acs + ","
    string = string[:-1]
    string += ";"

    for bcs in b_cat :
        string += bcs + ","
    string = string[:-1]

    return string

def add_to_hash(perm_hash, a_cat, b_cat) :
    a_swap_list = list(itertools.permutations(range(len(a_cat)), len(a_cat)))
    b_swap_list = list(itertools.permutations(range(len(b_cat)), len(b_cat)))
    row_swap_list = list(itertools.product(a_swap_list, b_swap_list))

    col_swap_list = list(itertools.permutations(range(len(a_cat[0])), len(a_cat[0])))

    for (asl, bsl) in row_swap_list :
        tt_a_cat = copy.deepcopy(a_cat)
        tt_b_cat = copy.deepcopy(b_cat)

        for i in range(len(asl)) :
            tt_a_cat[i] = a_cat[asl[i]]
        for i in range(len(bsl)) :
            tt_b_cat[i] = b_cat[bsl[i]]

        for csl in col_swap_list :
            t_a_cat = copy.deepcopy(tt_a_cat)
            t_b_cat = copy.deepcopy(tt_b_cat)
            for i in range(len(t_a_cat)) :
                t_a_cat_i = [int(j) for j in t_a_cat[i]]
                tt_a_cat_i = [int(j) for j in tt_a_cat[i]]
                for j in range(len(csl)) :
                    t_a_cat_i[j] = tt_a_cat_i[csl[j]]
                t_a_cat[i] = ''.join(str(j) for j in t_a_cat_i)

            for i in range(len(t_b_cat)) :
                t_b_cat_i = [int(j) for j in t_b_cat[i]]
                tt_b_cat_i = [int(j) for j in tt_b_cat[i]]
                for j in range(len(csl)) :
                    t_b_cat_i[j] = tt_b_cat_i[csl[j]]
                t_b_cat[i] = ''.join(str(j) for j in t_b_cat_i)
            
            t_a_cat.sort()
            t_b_cat.sort()

            perm_hash[make_key(t_a_cat, t_b_cat)] = True

    return True

####

def pp_execute_church(models, cat_size, inputs, query, church_exec_path, il, run) :
    t_outputs = [pmf.PMF() for i in range(len(models))]
    output_string = "" 
    output_string += "inputs = " + str(inputs) + "\n"
    #output_string += "check = " + str(check_cats_array) + "\n"

    for m in range(len(models)) :
        church.exec_model(models[m].x, inputs, t_outputs[m], query, \
                          church_exec_path)

    inputs_list = []
    tinputs_list = inputs.strip("\n").split("|")
    for til in range(len(tinputs_list)) :
        tinputs_list[til] = tinputs_list[til].split(";")
        for ttil in tinputs_list[til] :
            inputs_list.append(ttil)

    output_pmf = pmf.PMF(2)
    output_pmf[0] = pmf.P("A", 0.0)
    output_pmf[1] = pmf.P("B", 0.0)

    #### b_outputs
    b_outputs = [0]*(cat_size[0]+cat_size[1]+cat_size[2])
    for s in range(len(b_outputs)) :
        b_outputs[s] = [0]*len(models)
        for m in range(len(b_outputs[s])) :
            b_outputs[s][m] = copy.deepcopy(output_pmf)

    for s in range(len(b_outputs)) :
        for m in range(len(b_outputs[s])) :
            b_outputs[s][m][0].p = t_outputs[m][s].p
            b_outputs[s][m][1].p = 1 - t_outputs[m][s].p

    #output_string += "\n============\n"
    output_string += "b_outputs = " + str(b_outputs) + "\n"
    output_string += "\n============\n"

    #### recall
    bit_sequence = ["".join(seq) for seq in itertools.product("01", repeat=(cat_size[0]+cat_size[1]))]
    outputs_recall = [pmf.PMF(len(bit_sequence)) for i in range(len(models))]
    for m in range(len(outputs_recall)) :
        for o in range(len(outputs_recall[m])) :
            outputs_recall[m][o].x = bit_sequence[o]
            temp_b_output_prod = 0
            for b in range(len(outputs_recall[m][o].x)) :
                if (outputs_recall[m][o].x[b] == '0') :
                    temp_b_output_prod += math.log(b_outputs[b][m][0].p)
                else :
                    temp_b_output_prod += math.log(b_outputs[b][m][1].p)
            outputs_recall[m][o].p = math.exp(temp_b_output_prod)
    #output_string += "\n============\n"
    #output_string += "outputs_recall = " + str(outputs_recall) + "\n"

    expected_kl_recall = oed.get_expected_kl(models, outputs_recall)

    temp_result_recall = lib_result.Result(str(inputs) + " log" + str(il) + "_r" + str(run), expected_kl_recall)
    output_string += "result_recall = " + str(temp_result_recall) + "\n"

    #### transfer
    bit_sequence = ["".join(seq) for seq in itertools.product("01", repeat=(cat_size[2]))]
    outputs_transfer = [pmf.PMF(len(bit_sequence)) for i in range(len(models))]
    for m in range(len(outputs_transfer)) :
        for o in range(len(outputs_transfer[m])) :
            outputs_transfer[m][o].x = bit_sequence[o]
            temp_b_output_prod = 0
            for b in range(len(outputs_transfer[m][o].x)) :
                if (outputs_transfer[m][o].x[b] == '0') :
                    temp_b_output_prod += math.log(b_outputs[b+cat_size[0]+cat_size[1]][m][0].p)
                else :
                    temp_b_output_prod += math.log(b_outputs[b+cat_size[0]+cat_size[1]][m][1].p)
            outputs_transfer[m][o].p = math.exp(temp_b_output_prod)
    #output_string += "\n============\n"
    #output_string += "outputs_transfer = " + str(outputs_transfer) + "\n"

    expected_kl_transfer = oed.get_expected_kl(models, outputs_transfer)

    temp_result_transfer = lib_result.Result(str(inputs) + " log" + str(il) + "_r" + str(run), expected_kl_transfer)
    output_string += "result_transfer = " + str(temp_result_transfer) + "\n"

    #### total
    bit_sequence = ["".join(seq) for seq in itertools.product("01", repeat=(cat_size[0]+cat_size[1]+cat_size[2]))]
    outputs_total = [pmf.PMF(len(bit_sequence)) for i in range(len(models))]
    for m in range(len(outputs_total)) :
        for o in range(len(outputs_total[m])) :
            outputs_total[m][o].x = bit_sequence[o]
            temp_b_output_prod = 0
            for b in range(len(outputs_total[m][o].x)) :
                if (outputs_total[m][o].x[b] == '0') :
                    temp_b_output_prod += math.log(b_outputs[b][m][0].p)
                else :
                    temp_b_output_prod += math.log(b_outputs[b][m][1].p)
            outputs_total[m][o].p = math.exp(temp_b_output_prod)
    #output_string += "\n============\n"
    #output_string += "outputs_total = " + str(outputs_total) + "\n"

    expected_kl_total = oed.get_expected_kl(models, outputs_total)

    temp_result_total = lib_result.Result(str(inputs) + " log" + str(il) + "_r" + str(run), expected_kl_total)
    output_string += "result_total = " + str(temp_result_total) + "\n"

    return(temp_result_recall, temp_result_transfer, temp_result_total, output_string)

####

in_sequence_raw = ["".join(seq) for seq in itertools.product("01", repeat=4)]
in_sequence = [0]*len(in_sequence_raw)
for i in range(len(in_sequence_raw)) :
    in_sequence[i] = str(in_sequence_raw[i]).replace("", ",")[1:-1]

num_sequence = range(len(in_sequence))

#inputs_list_raw = []
inputs_perm_hash = {}
inputs_list = []
cat_size_list = []
#for a in range(3,6):
#    for b in range(3,6):
a_size = 5
b_size = 4
t_size = len(num_sequence) - a_size - b_size
a_comb = itertools.combinations(num_sequence, a_size)
for ac in a_comb :
    num_sequence_b = copy.deepcopy(num_sequence)
    for aac in ac :
        num_sequence_b.remove(aac)

    b_comb = itertools.combinations(num_sequence_b, b_size)

    for bc in b_comb :
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
            
        check_cats_valid = check_cats(a_cat_raw, b_cat_raw)
        if (check_cats_valid) :
            #print a_cat_raw, b_cat_raw
            valid_linsep = True
            for (a1, a2) in itertools.combinations(a_cat_raw, 2):
                ma = mask(a1, a2)
                #print "a:", a1, a2, ma
                if (sum(ma) != 1) :
                    for (b1, b2) in itertools.combinations(b_cat_raw, 2) :
                        #print "    b:", b1, b2
                        if does_intersect(a1, a2, b1, b2) :
                            valid_linsep = False
                            break

                #print "    ", valid_linsep
                if (valid_linsep == False) :
                    break
                #print ""

            if (valid_linsep) :
                #print make_key(a_cat_raw, b_cat_raw)
                #print inputs_perm_hash

                if (make_key(a_cat_raw, b_cat_raw) not in inputs_perm_hash):
                    tc = copy.deepcopy(num_sequence)
                    for aac in ac :
                        tc.remove(aac)
                    for bbc in bc :
                        tc.remove(bbc)
                    #print ac, bc, tc

                    t_cat_raw = []
                    t_cat = []
                    for ttc in tc :
                        t_cat_raw.append(in_sequence_raw[ttc])
                        t_cat.append(in_sequence[ttc])

                    a_string = ""
                    for acs in a_cat :
                        a_string += acs + ";"
                    a_string = a_string[:-1]

                    b_string = ""
                    for bcs in b_cat :
                        b_string += bcs + ";"
                    b_string = b_string[:-1]

                    t_string = ""
                    for tcs in t_cat :
                        t_string += tcs + ";"
                    t_string = t_string[:-1]

                    #inputs_list_raw.append([a_cat_raw, b_cat_raw])
                    #print "committing: ", [a_cat_raw, b_cat_raw]

                    add_to_hash(inputs_perm_hash, a_cat_raw, b_cat_raw)

                    #print inputs_perm_hash
                    #print len(inputs_perm_hash)

                    inputs_list.append(a_string + "|" + b_string + "|" + t_string)
                    cat_size_list.append([a_size, b_size, t_size])
                    #sys.exit(0)
                #else :
                #    print "hit :", make_key(a_cat_raw, b_cat_raw)
                    
##                check_cats_list.append(check_cats_array)

#a_target = ['0,1,1,1','1,0,1,0','1,0,1,1','1,1,0,1','1,1,1,0']
#b_target = ['0,0,0,0','0,0,0,1','0,1,1,0','1,1,0,0']

#inputs_list = [ \
#    '0,0,0,1;0,0,1,1;1,1,0,0;1,1,1,0;1,1,1,1|0,1,0,0;0,1,1,0;1,0,0,0;1,0,1,0|0,0,0,0;0,0,1,0;0,1,0,1;0,1,1,1;1,0,0,1;1,0,1,1;1,1,0,1', \
#    '0,0,1,1;0,1,0,0;1,0,1,1;1,1,0,0;1,1,1,0|0,0,0,0;0,0,1,0;1,0,0,0;1,0,1,0|0,0,0,1;0,1,0,1;0,1,1,0;0,1,1,1;1,0,0,1;1,1,0,1;1,1,1,1', \
#    '0,0,1,1;0,1,0,0;0,1,0,1;1,0,1,0;1,1,0,0|0,0,0,0;0,0,0,1;0,0,1,0;1,0,0,0|0,1,1,0;0,1,1,1;1,0,0,1;1,0,1,1;1,1,0,1;1,1,1,0;1,1,1,1', \
#    '0,1,1,1;1,0,1,0;1,0,1,1;1,1,0,1;1,1,1,0|0,0,0,0;0,0,0,1;0,1,1,0;1,1,0,0|0,0,1,0;0,0,1,1;0,1,0,0;0,1,0,1;1,0,0,0;1,0,0,1;1,1,1,1'] # MS
#
#cat_size_list = [[5, 4, 7], \
#                 [5, 4, 7], \
#                 [5, 4, 7], \
#                 [5, 4, 7]]

inputs_list.append('1,1,1,0;1,0,1,0;1,0,1,1;1,1,0,1;0,1,1,1|1,1,0,0;0,1,1,0;0,0,0,1;0,0,0,0|0,0,1,0;0,0,1,1;0,1,0,0;0,1,0,1;1,0,0,0;1,0,0,1;1,1,1,1')
cat_size_list.append([5, 4, 7])

inputs_list_f = inputs_list
cat_size_list_f = cat_size_list

############################################

t_init = datetime.datetime.now()
print "time start: ", str(t_init)
sys.stdout.flush()

samples = 5000

f_log = open("results/mh_s" + str(samples) + "_l1.txt", "w")
f_log_s_recall         = open("results/mh_s" + str(samples) + "_l1_sorted_recall.txt", "w")
f_log_stats_recall     = open("results/mh_s" + str(samples) + "_l1_stats_recall.txt", "w")
f_log_stats_s_recall   = open("results/mh_s" + str(samples) + "_l1_stats_recall_sorted.txt", "w")
f_log_s_transfer       = open("results/mh_s" + str(samples) + "_l1_sorted_transfer.txt", "w")
f_log_stats_transfer   = open("results/mh_s" + str(samples) + "_l1_stats_transfer.txt", "w")
f_log_stats_s_transfer = open("results/mh_s" + str(samples) + "_l1_stats_transfer_sorted.txt", "w")
f_log_s_total          = open("results/mh_s" + str(samples) + "_l1_sorted_total.txt", "w")
f_log_stats_total      = open("results/mh_s" + str(samples) + "_l1_stats_total.txt", "w")
f_log_stats_s_total    = open("results/mh_s" + str(samples) + "_l1_stats_total_sorted.txt", "w")

f_log.write("samples = " + str(samples) + "\n")
f_log.write("number of experiments = " + str(len(inputs_list_f)) + "\n\n")
f_log.flush()

models = pmf.PMF(2)
models[0] = pmf.P("independent_cue.church", 1.0)
models[1] = pmf.P("context_theory.church", 1.0)
models.normalize()

query = "raw"
#TODO
#church_exec_path = "/home/danly/webchurch/church"
church_exec_path = "/opt/webchurch/church"

results_recall = []
stat_results_recall = []
results_transfer = []
stat_results_transfer = []
results_total = []
stat_results_total = []

ppservers = ("128.0.0.1",)
#TODO
#ncpus = 24
ncpus = 4
job_server = pp.Server(ncpus, ppservers=ppservers, secret="cocolab")

#TODO
#input_len = len(inputs_list_f)
input_len = 4
run_len = 4
job_results = [0]*input_len
for il in range(input_len) :
    job_results[il] = [0]*run_len

#for il in range(input_len) :
#    for r in range(run_len) :
#        job_results[il][r] = pp_execute_church(models, cat_size_list_f[il], inputs_list_f[il], query, church_exec_path, check_cats_list[il], il, r)
#a = pp_execute_church(models, cat_size_list_f[0], inputs_list_f[0], query, church_exec_path, check_cats_list[il], il, r)
for il in range(input_len) :
    for r in range(run_len) :
        job_results[il][r] = job_server.submit(pp_execute_church, \
                                               (models, cat_size_list_f[il], inputs_list_f[il], query, church_exec_path, il, r), \
                                               (), ("oed", "pmf", "church", "copy", "itertools", "lib_result", "math"))

for il in range(input_len) :
    temp_ekl_recall   = [0]*run_len
    temp_ekl_transfer = [0]*run_len
    temp_ekl_total    = [0]*run_len
    for r in range(run_len) :
        f_log_tmp = open("results/log/log" + str(il) + "_r" + str(r) + "_mh_s" + str(samples) + "_l1.txt", "w")
        (temp_result_recall, temp_result_transfer, temp_result_total, output_string) = job_results[il][r]()

        results_recall.append(temp_result_recall)
        temp_ekl_recall[r] = temp_result_recall.outp
        results_transfer.append(temp_result_transfer)
        temp_ekl_transfer[r] = temp_result_transfer.outp
        results_total.append(temp_result_total)
        temp_ekl_total[r] = temp_result_total.outp

        f_log_tmp.write(output_string + "\n")
        f_log_tmp.flush()
        f_log_tmp.close()

        f_log.write(str(temp_result_recall) + " r " + "\n")
        f_log.write(str(temp_result_recall) + " t " + "\n")
        f_log.write(str(temp_result_recall) + " a " + "\n")
        f_log.flush()

    temp_ekl_mean_recall = numpy.mean(temp_ekl_recall)
    temp_ekl_se_recall   = scipy.stats.sem(temp_ekl_recall)
    temp_ekl_var_recall  = numpy.var(temp_ekl_recall)
    temp_stat_result_recall = lib_result.Stat_result(str(inputs_list_f[il]) + " log" + str(il), temp_ekl_mean_recall, temp_ekl_se_recall, temp_ekl_var_recall)
    stat_results_recall.append(temp_stat_result_recall)
    f_log_stats_recall.write(str(temp_stat_result_recall) + "\n")
    f_log_stats_recall.flush()

    temp_ekl_mean_transfer = numpy.mean(temp_ekl_transfer)
    temp_ekl_se_transfer   = scipy.stats.sem(temp_ekl_transfer)
    temp_ekl_var_transfer  = numpy.var(temp_ekl_transfer)
    temp_stat_result_transfer = lib_result.Stat_result(str(inputs_list_f[il]) + " log" + str(il), temp_ekl_mean_transfer, temp_ekl_se_transfer, temp_ekl_var_transfer)
    stat_results_transfer.append(temp_stat_result_transfer)
    f_log_stats_transfer.write(str(temp_stat_result_transfer) + "\n")
    f_log_stats_transfer.flush()

    temp_ekl_mean_total = numpy.mean(temp_ekl_total)
    temp_ekl_se_total   = scipy.stats.sem(temp_ekl_total)
    temp_ekl_var_total  = numpy.var(temp_ekl_total)
    temp_stat_result_total = lib_result.Stat_result(str(inputs_list_f[il]) + " log" + str(il), temp_ekl_mean_total, temp_ekl_se_total, temp_ekl_var_total)
    stat_results_total.append(temp_stat_result_total)
    f_log_stats_total.write(str(temp_stat_result_total) + "\n")
    f_log_stats_total.flush()

sorted_results_recall = sorted(results_recall, key=lambda result: result.outp, reverse=True)
for s in sorted_results_recall :
    f_log_s_recall.write(str(s) + "\n")
sorted_stat_results_recall = sorted(stat_results_recall, key=lambda result: result.outp_mean, reverse=True)
for sr in sorted_stat_results_recall :
    f_log_stats_s_recall.write(str(sr) + "\n")

sorted_results_transfer = sorted(results_transfer, key=lambda result: result.outp, reverse=True)
for s in sorted_results_transfer :
    f_log_s_transfer.write(str(s) + "\n")
sorted_stat_results_transfer = sorted(stat_results_transfer, key=lambda result: result.outp_mean, reverse=True)
for sr in sorted_stat_results_transfer :
    f_log_stats_s_transfer.write(str(sr) + "\n")

sorted_results_total = sorted(results_total, key=lambda result: result.outp, reverse=True)
for s in sorted_results_total :
    f_log_s_total.write(str(s) + "\n")
sorted_stat_results_total = sorted(stat_results_total, key=lambda result: result.outp_mean, reverse=True)
for sr in sorted_stat_results_total :
    f_log_stats_s_total.write(str(sr) + "\n")

t_finish = datetime.datetime.now()

f_log.close()
f_log_s_recall.close()
f_log_stats_recall.close()
f_log_stats_s_recall.close()
f_log_s_transfer.close()
f_log_stats_transfer.close()
f_log_stats_s_transfer.close()
f_log_s_transfer.close()
f_log_stats_transfer.close()
f_log_stats_s_transfer.close()

print "time finish: ", str(t_finish)
print "time elapsed: ", (t_finish - t_init).total_seconds()
sys.stdout.flush()


