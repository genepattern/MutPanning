#!/usr/bin/env python

#************************************************************           
#* Cancer Bayesian Selection Estimation (CBaSE):			*
#* Code accompanying Weghorn & Sunyaev, Nat. Genet. (2017).	*   
#*															*   
#* Author:		Donate Weghorn								*   
#*															*   
#* Copyright:	(C) 2018 Donate Weghorn						*   
#*															*   
#* License:		Public Domain								*   
#*															*   
#* Version:		1.1											*   
#************************************************************


from scipy.optimize import minimize
import scipy.special as sp
import sys
import itertools as it
import math
import gzip
import numpy as np
import argparse
import os
import random
import subprocess
import mpmath as mp
import scipy.stats as st
import itertools as it
import subprocess
import glob


#************************************************************************************************************
#	FUNCTION DEFINITIONS

def import_special_genes(filename):
	fin = open(filename)
	lines = fin.readlines()
	fin.close()
	c_genes = []
	for line in lines:
		c_genes.append(line.strip().split()[0])
	return c_genes
def import_syn_data(filename):
	fin = open(filename)
	lines = fin.readlines()
	fin.close()
	mut_array=[]
	for line in lines[1:]:
		field = line.strip().split("\t")
		mut_array.append({"gene": field[0], "ell_s": float(field[1]), "sobs": int(field[2])})
	return mut_array
def import_nonsyn_data(filename):
	fin = open(filename)
	lines = fin.readlines()
	fin.close()
	mut_array=[]
	for line in lines[1:]:
		field = line.strip().split("\t")
		mut_array.append({"gene": field[0], "ell_x": float(field[1]), "xobs": int(field[2])})
	return mut_array
def export_pofx_given_s(p, genes, aux):
	[modC, filename] = aux

	if [1,2].count(modC):
		a,b = p
	elif [3,4].count(modC):
		a,b,t,w = p
	elif [5,6].count(modC):
		a,b,g,d,w = p

	if modC==1:
		#*************** lambda ~ Gamma:
		def pofs(s, L):
			return (L*b)**s * (1 + L*b)**(-s-a) * math.gamma(s + a) / (math.gamma(s+1) * math.gamma(a))
		def pofx_given_s(x, s, L, r, thr):
			return np.exp( x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) / pofs(s, L)
	elif modC==2:
		#*************** lambda ~ IG:
		def pofs(s, L, thr):
			if thr:
				return 2. * mp.exp( ((s + a)/2.)*math.log(L*b) + mp.log(mp.besselk(-s + a, 2*np.sqrt(L*b))) - sp.gammaln(s+1) - sp.gammaln(a) )
			else:
				return 2. * math.exp( ((s + a)/2.)*math.log(L*b) + np.log(sp.kv(-s + a, 2*math.sqrt(L*b))) - sp.gammaln(s+1) - sp.gammaln(a) )
		def pofx_given_s(x, s, L, r, thr):
			if thr:
				return mp.exp( np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (1/2. * (-s - x + a))*np.log((L*(1 + r))/b) + a*np.log(b) + mp.log(mp.besselk(s + x - a, 2*math.sqrt(L*(1 + r)*b))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) / pofs(s, L, thr)
			else:
				return np.exp( np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (1/2. * (-s - x + a))*np.log((L*(1 + r))/b) + a*np.log(b) + np.log(sp.kv(s + x - a, 2*math.sqrt(L*(1 + r)*b))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) / pofs(s, L, thr)
	elif modC==3:
		#*************** lambda ~ w * Exp + (1-w) * Gamma:
		def pofs(s, L):
			return np.exp( np.log(w) + s*np.log(L) + np.log(t) + (-1 - s)*np.log(L + t) ) + np.exp( np.log(1.-w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) )
		def pofx_given_s(x, s, L, r, thr):
			return ( np.exp( np.log(w) + (s + x)*np.log(L) + x*np.log(r) + np.log(t) + (-1 - s - x)*np.log(L + L*r + t) + sp.gammaln(1 + s + x) - sp.gammaln(s+1) - sp.gammaln(x+1) ) + np.exp( np.log(1-w) + x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) ) / pofs(s, L)
	elif modC==4:
		#*************** lambda ~ w * Exp + (1-w) * InvGamma:
		def pofs(s, L, thr):
			if thr:
				return (w * t * mp.exp( s*np.log(L) + (-1 - s)*np.log(L + t) )) + mp.exp( np.log(1.-w) + np.log(2.) + ((s + a)/2.)*np.log(L*b) + mp.log(mp.besselk(-s + a, 2*math.sqrt(L*b))) - sp.gammaln(s+1) - sp.gammaln(a) )
			else:
				return (w * L**s * t * (L + t)**(-1 - s)) + np.exp( np.log(1.-w) + np.log(2.) + ((s + a)/2.)*np.log(L*b) + np.log(sp.kv(-s + a, 2*math.sqrt(L*b))) - sp.gammaln(s+1) - sp.gammaln(a) )
		def pofx_given_s(x, s, L, r, thr):
			if thr:
				return ( np.exp( np.log(w) + (s + x)*np.log(L) + x*np.log(r) + np.log(t) + (-1 - s - x)*np.log(L + L*r + t) + sp.gammaln(1 + s + x) - sp.gammaln(s+1) - sp.gammaln(x+1) ) + mp.exp( np.log(1.-w) + np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (0.5 * (-s - x + a))*np.log((L*(1 + r))/b) + a*np.log(b) + mp.log(mp.besselk(s + x - a, 2*np.sqrt(L*(1 + r)*b))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) ) / pofs(s, L, thr)
			else:
				return ( np.exp( np.log(w) + (s + x)*np.log(L) + x*np.log(r) + np.log(t) + (-1 - s - x)*np.log(L + L*r + t) + sp.gammaln(1 + s + x) - sp.gammaln(s+1) - sp.gammaln(x+1) ) + np.exp( np.log(1.-w) + np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (0.5 * (-s - x + a))*np.log((L*(1 + r))/b) + a*np.log(b) + np.log(sp.kv(s + x - a, 2*np.sqrt(L*(1 + r)*b))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) ) / pofs(s, L, thr)
	elif modC==5:
		#*************** lambda ~ w * Gamma + (1-w) * Gamma (Gamma mixture model):
		def pofs(s, L):
			return np.exp( np.log(w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) ) + np.exp( np.log(1.-w) + s*np.log(L*d) + (-s-g)*np.log(1 + L*d) + sp.gammaln(s + g) - sp.gammaln(s+1) - sp.gammaln(g) )
		def pofx_given_s(x, s, L, r, thr):
			return ( np.exp( np.log(w) + x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) + np.exp( np.log(1-w) + x*np.log(r) + (s + x)*np.log(L*d) + (-s - x - g)*np.log(1 + L*(1 + r)*d) + sp.gammaln(s + x + g) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(g) ) ) / pofs(s, L)
	elif modC==6:
		#*************** lambda ~ w * Gamma + (1-w) * InvGamma (mixture model):
		def pofs(s, L, thr):
			if thr:
				return np.exp( np.log(w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) ) + mp.exp( np.log(1.-w) + np.log(2.) + ((s + g)/2.)*np.log(L*d) + mp.log(mp.besselk(-s + g, 2*mp.sqrt(L*d))) - sp.gammaln(s+1) - sp.gammaln(g) )
			else:
				return np.exp( np.log(w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) ) + np.exp( np.log(1.-w) + np.log(2.) + ((s + g)/2.)*np.log(L*d) + np.log(sp.kv(-s + g, 2*np.sqrt(L*d))) - sp.gammaln(s+1) - sp.gammaln(g) )
		def pofx_given_s(x, s, L, r, thr):
			if thr:
				return ( np.exp( np.log(w) + x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) + mp.exp( np.log(1-w) + np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (0.5 * (-s - x + g))*np.log((L*(1 + r))/d) + g*np.log(d) + mp.log(mp.besselk(s + x - g, 2*np.sqrt(L*(1 + r)*d))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(g) ) ) / pofs(s, L, thr)
			else:
				return ( np.exp( np.log(w) + x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) + np.exp( np.log(1-w) + np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (0.5 * (-s - x + g))*np.log((L*(1 + r))/d) + g*np.log(d) + np.log(sp.kv(s + x - g, 2*np.sqrt(L*(1 + r)*d))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(g) ) ) / pofs(s, L, thr)

	probs=[]
	gcnt=0
	L=1.
	fout = open(filename, "w")
	for gene in genes:
		gcnt += 1
		if gcnt%1000==0:
			sys.stderr.write("%i%% done.\r" %int(float(gcnt)/len(genes)*100.))

		sobs = gene["sobs"]
		xobs = gene["xobs"]
		sexp = gene["ell_s"]
		xexp = gene["ell_x"]
		ratx = xexp/sexp
		if ratx<1e-10:
			continue

		last_p = 0.
		flag_peak = 0
		probs=[]
		for x in range(10000000):
			if last_p>1. or math.isnan(last_p):
				break
			if last_p>1e-5:
				try:
					cur_p = pofx_given_s(x, sobs, L, ratx, 0)
					if cur_p>1. or math.isnan(cur_p):
						cur_p = pofx_given_s(x, sobs, L, ratx, 1)
					probs.append([x, cur_p])
				except:
					cur_p = pofx_given_s(x, sobs, L, ratx, 1)
					probs.append([x, cur_p])
			else:
				if flag_peak and x>xobs:
					break
				else:
					try:
						cur_p = pofx_given_s(x, sobs, L, ratx, 0)
						if cur_p>1. or math.isnan(cur_p):
							cur_p = pofx_given_s(x, sobs, L, ratx, 1)
						probs.append([x, cur_p])
					except:
						cur_p = pofx_given_s(x, sobs, L, ratx, 1)
						probs.append([x, cur_p])
			if cur_p-last_p<0:
				flag_peak = 1
			last_p = cur_p

		fout.write("%s\n" %gene["gene"])
		for prob in probs:
			fout.write("%i\t%f\n" %(prob[0], prob[1]))
	fout.close()

#************************************************************************************************************
#	COMMAND LINE ARGS

infile_syn		= str(sys.argv[1])		#	synonymous somatic mutation data input file
infile_nonsyn	= str(sys.argv[2])		#	nonsynonymous somatic mutation data input file
infile_param = str(sys.argv[3])
outfile			= str(sys.argv[4])		#	output data file
#filepath		= str(sys.argv[4])		#	path to auxiliary input files folder
#mod_C			= int(sys.argv[5])		#	model choice: 1=G, 2=IG, 3=EmixG, 4=EmixIG, 5=GmixG, 6=GmixIG

#************************************************************************************************************
#	GLOBAL DEFINITIONS & AUXILIARY DATA

mod_choice	= ["", "Gamma(a,b): [a,b] =", "InverseGamma(a,b): [a,b] =", "w Exp(t) + (1-w) Gamma(a,b): [a,b,t,w] =", "w Exp(t) + (1-w) InverseGamma(a,b): [a,b,t,w] =", "w Gamma(a,b) + (1-w) Gamma(g,d): [a,b,g,d,w] =", "w Gamma(a,b) + (1-w) InverseGamma(g,d): [a,b,g,d,w] ="]		# model choice
#zero_genes = import_special_genes("%s/zero_genes.txt" %filepath)

#************************************************************************************************************
#************************************************************************************************************


p_files = glob.glob('param_estimates_*.txt')

all_models=[]
fin = open(infile_param)
lines = fin.readlines()
fin.close()
for line in lines:
	field = line.strip().split(", ")
	all_models.append([float(el) for el in field])

modC_map = [2,2,4,4,5,5]						# map model --> number of params

cur_min=1e20; cur_ind=10
for m in range(len(all_models)):
	if 2.*modC_map[int(all_models[m][-1])-1] + 2.*all_models[m][-2] < cur_min:
		cur_min = 2.*modC_map[int(all_models[m][-1])-1] + 2.*all_models[m][-2]
		cur_ind = m


mod_C	= int(all_models[cur_ind][-1])
params	= all_models[cur_ind][:-2]
AIC		= cur_min

print "Using model %s" %(mod_choice[mod_C]), params

syn_muts = import_syn_data(infile_syn)
nonsyn_muts = import_nonsyn_data(infile_nonsyn)
syn_muts = [mks for mks in syn_muts if (len(mks["gene"])>2 and mks["gene"][:2]=="OR" and ['0','1','2','3','4','5','6','7','8','9'].count(mks["gene"][2]))==0 and mks["ell_s"]>0.]
sys.stderr.write("Filtered out OR genes, leaving %i genes.\n" %len(syn_muts))
#syn_muts = [gene for gene in syn_muts if zero_genes.count(gene["gene"])==0]
#sys.stderr.write("Filtered out likely undercalled genes, leaving %i genes.\n" %len(mks_type))

nonsyn_muts = sorted([gene for gene in nonsyn_muts if gene["gene"] in [el["gene"] for el in syn_muts]], key=lambda arg: arg["gene"])
syn_muts = sorted([gene for gene in syn_muts if gene["gene"] in [el["gene"] for el in nonsyn_muts]], key=lambda arg: arg["gene"])
joint=[]
for gind in range(len(syn_muts)):
	joint.append({"gene": syn_muts[gind]["gene"], "ell_s": syn_muts[gind]["ell_s"], "sobs": syn_muts[gind]["sobs"], "ell_x": nonsyn_muts[gind]["ell_x"], "xobs": nonsyn_muts[gind]["xobs"]})

export_pofx_given_s(params, joint, [mod_C, outfile])










