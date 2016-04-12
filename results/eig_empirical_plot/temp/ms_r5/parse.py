#!/usr/bin/python

import csv

subj_file = open("2014-10-31-subjects_2.csv", 'r')
tran_file = open("2014-10-31-transfer_2.csv", 'r')

out_file = open("ms.txt", 'w')
stats_file = open("ms_stats.txt", 'w')

subj_csv = csv.DictReader(subj_file, delimiter=',')
subjid_success = []
for line in subj_csv :
    if (line["learn.success"] == "1") :
        subjid_success.append(line["subject"])

#print subjid_success
#print len(subjid_success)
out_file.write(str(len(subjid_success)) + "\n")

tran_csv = csv.DictReader(tran_file, delimiter=',')
tran_data = {}
for line in tran_csv :
    #print line["timeStart"]
    #print line["subject"]
    if line["subject"] in subjid_success :
        if line["subject"] not in tran_data :
            tran_data[line["subject"]] = {}
        tran_data[line["subject"]][line["stim"]] = line["category"]
    #    print line["trial.num"]
    #tran_csv[line[""]]


#inputs = '0,0,0,1;0,0,1,1;1,1,0,0;1,1,1,0;1,1,1,1|0,1,0,0;0,1,1,0;1,0,0,0;1,0,1,0|0,0,0,0;0,0,1,0;0,1,0,1;0,1,1,1;1,0,0,1;1,0,1,1;1,1,0,1'
inputs = '1,1,1,0;1,0,1,0;1,0,1,1;1,1,0,1;0,1,1,1|1,1,0,0;0,1,1,0;0,0,0,1;0,0,0,0|0,0,1,0;0,0,1,1;0,1,0,0;0,1,0,1;1,0,0,0;1,0,0,1;1,1,1,1'

inputs_list = []
tinputs_list = inputs.strip("\n").split("|")
for til in range(len(tinputs_list)) :
    tinputs_list[til] = tinputs_list[til].split(";")
    for ttil in tinputs_list[til] :
        inputs_list.append(ttil)

#print inputs
out_file.write(inputs + "\n")
#print inputs_list

#print tran_data
for il in inputs_list :
    out_str = il + "|"
    #print il.replace(",", "")
    for ss in subjid_success :
        out_str += tran_data[ss][il.replace(",", "")] + ","
    out_str = out_str[:-1]
    #print out_str
    out_file.write(out_str + "\n")

for il in inputs_list :
    stats_str = il + "|"
    #print il.replace(",", "")
    count = 0
    for ss in subjid_success :
        if (tran_data[ss][il.replace(",","")] == "A") :
            count += 1
    stats_str += str(count)
    #print out_str
    stats_file.write(stats_str + "/" +  str(len(tran_data)) + "\n")

stats_file.write("\n")

for il in inputs_list :
    stats_str = il + "|"
    #print il.replace(",", "")
    count = 0
    for ss in subjid_success :
        if (tran_data[ss][il.replace(",","")] == "A") :
            count += 1
    stats_str += str((1.0*count)/len(tran_data))
    #print out_str
    stats_file.write(stats_str + "\n")

