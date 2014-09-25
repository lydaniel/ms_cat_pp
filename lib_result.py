#!/usr/bin/python

class Result:
        def __init__(self, inp, outp):
                self.inp = inp 
                self.outp = outp
        def __repr__(self):
                return repr((self.inp, self.outp))

class Stat_result:
        def __init__(self, inp, outp_mean, outp_se, outp_var):
                self.inp = inp 
                self.outp_mean = outp_mean
                self.outp_se = outp_se
                self.outp_var = outp_var
        def __repr__(self):
                return repr((self.inp, self.outp_mean, self.outp_se, self.outp_var))

