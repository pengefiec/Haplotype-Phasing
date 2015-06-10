import random
import os
import sys
import time
import itertools
from itertools import izip
from data_generator import mix
from collections import defaultdict
import numpy as np, numpy.random
def trivial_phase(indivs):
	"""
	a trivial algorithm, check all combinations from a single hyplotype type.
	"""
	pool=make_pool(len(indivs[0]))

	for i in xrange(1,len(pool)+1):
		all_combi=itertools.combinations(pool,i)
		for t in all_combi:
			t+=t
			candidate_couples=list(itertools.combinations(t,2))
			geno_list=map(lambda x: mix(x[0],x[1]), candidate_couples)
	 		if check(indivs, geno_list):
	 			return list(set(t)), candidate_couples
	print "It's impossible to execute this, something must be wrong."
		
def greedy_phase(indivs):
	"""
	greedy algorithm.
	"""	
	start = time.clock()
	pool=make_pool(len(indivs[0]))
	# print time.clock()-start, "time in coombinaiton"
	s_pool=sorted(pool, key=lambda x: get_ability(x, indivs))

	s_pool=s_pool[::-1]
	s=0
	
	for i in xrange(1,len(s_pool)+1):
		hyplos=s_pool[0:i]
		hyplos+=hyplos
		
		candidate_couples=set(list(itertools.combinations(hyplos,2)))
		
		# print candidate_couples
		geno_list=map(lambda x: mix(x[0],x[1]), candidate_couples)
		# s+= time.clock()-start
		# print geno_list
 		if check(indivs, geno_list):
 			
 			return list(set(hyplos))
 	
	print "It's impossible to execute this, something must be wrong."

def clark_phase(indivs):
	hyplo_list=[]
	pairs=len(indivs)*[None]
	for i in xrange(0, len(indivs)):
		x=indivs[i]
		if check_fixed(x):
			hyplo_list.append(x)
			pairs[i]=[x,x]
	# hyplo_list=[x for x in indivs if check_fixed(x)]
	indiv_list=[x for x in indivs if x not in hyplo_list]
	candidate_couples=list(itertools.combinations(hyplo_list,2))
	geno_list=map(lambda x: mix(x[0],x[1]), candidate_couples)
	new_indiv_list=[]

	for i in xrange(0,len(indiv_list)):
		x=indiv_list[i]
		if x not in geno_list:
			new_indiv_list.append(x)
		else:
			pairs[indivs.index(x)]=candidate_couples[geno_list.index(x)]
	# print "this is unphased indivs: ",len(new_indiv_list)
	# print "this is haplotypes: ", len(hyplo_list)
	# print indiv_list
	# print len(new_indiv_list)
	for x in new_indiv_list:
		add_compl(hyplo_list, x, new_indiv_list, pairs)
	return list(set(hyplo_list)), pairs
def em_phase(indivs):
	"""
	phase using em algorithm, it return the the parsimony as well as the decompesed format of individuals.
	"""
	#get the start frequency using dirichlet distribution.
	hyplo_collection=[]
	indiv_dict=defaultdict(list)
	# hyplo_dict=defaultdict(float)
	res=[]
	res_pairs=[]
	for x in indivs:
		decom=decompose_acurate(x)
		indiv_dict[x]+=decom
		hyplo_collection+=list(itertools.chain.from_iterable(decom))
	return em(indiv_dict, hyplo_collection)

def em(indiv_dict, hyplo_collection):
	"""
	take individual dictionary and all of possible haplotypes, output a list of most possible haplotypes 
	and most possible pairs.
	"""
	hyplo_dict=defaultdict(float)
	res=[]
	res_pairs=[]
	start_fre= np.random.dirichlet(np.ones(len(hyplo_collection)),size=1)[0]
	i=0
	for x in hyplo_collection:
		hyplo_dict[x]=start_fre[i]
		i+=1
	prev_dict=hyplo_dict
	hyplo_dict=phase_cnt(indiv_dict, hyplo_dict)
	while True:
		if check_converge(prev_dict, hyplo_dict)==False:
			prev_dict=hyplo_dict
			hyplo_dict=phase_cnt(indiv_dict, hyplo_dict)
		else:
			break
	for k,v in indiv_dict.iteritems():
		pair=get_best(v, hyplo_dict)
		res+=pair
		res_pairs.append(pair)
	key_list=indiv_dict.keys()
	return list(set(res)), res_pairs

def get_best_combinations(indiv_chunks_1, indiv_chunks_2, pair_chunks_1, pair_chunks_2):
	"""
	using em to get a list of most possible combied haplotypes, and most possible combied chunk pairs.
	"""
	hyplo_collection=[]
	indiv_dict=defaultdict(list)
	for i in xrange(0, len(indiv_chunks_1)):
		decom=get_chunk_combination(indiv_chunks_1[i], indiv_chunks_2[i], pair_chunks_1[i], pair_chunks_2[i])
		indiv_combie=indiv_chunks_1[i]+ indiv_chunks_2[i]
		hyplo_collection+=list(itertools.chain.from_iterable(decom))
		indiv_dict[indiv_combie]+=decom
	return em(indiv_dict, hyplo_collection)



def get_best(pairs, hyplo_dict):
	"""
	get the best pair from all possible pairs for a single genotype, used in em.
	"""
	best_freq=0
	best_pair=[]
	# print "the pairs: ", pairs
	for p in pairs:
		freq=hyplo_dict[p[0]]*hyplo_dict[p[1]]
		if freq>best_freq:
			best_freq=freq
			best_pair=p
	# print best_pair
	# print "the best pair is: ",best_pair
	return best_pair

def check_converge(dict_1, dict_2):
	for k, v in dict_1.iteritems():
		if abs(dict_1[k]-dict_2[k])>0.0001:
			return False
	return True		
	
def phase_cnt(indiv_dict, hyplo_dict):
	new_indiv_dict=dict.fromkeys(list(hyplo_dict.keys()))
	for k,v in new_indiv_dict.iteritems():
		new_indiv_dict[k]=0.0
	for k,v in indiv_dict.iteritems():
		freq_list=cal_freq(v, hyplo_dict)
		total=sum(freq_list)
		i=0
		for pair in v:
			pair_freq=freq_list[i]/total
			new_indiv_dict[pair[0]]+=pair_freq
			new_indiv_dict[pair[1]]+=pair_freq
			i+=1
	new_total=sum(new_indiv_dict.values())
	for k,v in new_indiv_dict.iteritems():
		new_indiv_dict[k]=v/new_total
	return new_indiv_dict

def cal_freq(pairs, hyplo_dict):
	res=[]
	for p in pairs:
		res.append(hyplo_dict[p[0]]*hyplo_dict[p[1]])
	return res

def add_compl(hyplo_list, geno, indivs, pairs):
	"""
	get a complementary hyplotype if there is only one explaining hyplotype for a particular genotype and make up 2 hyplo if not.
	"""
	getted=False
	for hyplo in hyplo_list:
		if match(hyplo, geno):
			getted=True
			compl=get_compl(hyplo, geno)
			hyplo_list.append(compl)
			pairs[indivs.index(geno)]=[hyplo, compl]
			break
	if not getted:
		hyplo_list+=list(decompose(geno))
		# print indivs.index(geno)
		pairs[indivs.index(geno)]=list(decompose(geno))

def decompose_acurate(indiv):
	"""
	get all possible hyplotype for a indiv.
	"""
	res=[]
	# get the # and positions of 2.
	cnt, positions=collect_two(indiv)
	if cnt==0:
		res.append([indiv, indiv])
	else:
		pos_perm=make_pool(cnt)
		remove_compl(pos_perm)
		uncompl_list=[]
		temp=indiv
		for p in pos_perm:
			i=0
			for x in positions:
				temp=list(temp)
				temp[x]=str(p[i])
				temp="".join(temp)
				i+=1
			uncompl_list.append(temp)
		res=[[x, get_compl(x, indiv)] for x in uncompl_list]
	return res
def get_chunk_combination(indiv_1, indiv_2, pair_1, pair_2):
	"""
	output all possible chunk combinations.
	"""	
	if check_fixed(indiv_1):
		return [[indiv_1+indiv_2, indiv_1+indiv_2]] if check_fixed(indiv_2) else [[indiv_1+pair_2[0], indiv_1+pair_2[1]]]
	else:
		return [[pair_1[0]+indiv_2, pair_1[1]+indiv_2]] if check_fixed(indiv_2) else [[pair_1[0]+pair_2[0], pair_1[0]+pair_2[1]], [pair_1[1]+pair_2[0], pair_1[1]+pair_2[1]]]
	return res

def remove_compl(hyplo_list):
	"""
	remove the complementary hyplotype in the hyplotype list.
	"""
	for x in hyplo_list:
		if get_coml_s(x) in hyplo_list:
			hyplo_list.remove(x)

def get_coml_s(hyplo):
	"""
	get the complementary for a hyplotype.
	"""
	res=""
	for x in hyplo:
		if x=="1":
			res+="0"
		else:
		 	res+="1"
	return res

def collect_two(indiv):
	"""
	get the number and position of 2 for a single hyplotype.
	"""	
	cnt=0
	positions=[]
	for i in xrange(0,len(indiv)):
		if indiv[i]=="2":
			cnt+=1
			positions.append(i)
	return cnt, positions

def decompose(geno):
	"""
	decompose the geno arbitrarily.
	"""
	a=""
	b=""
	for i in xrange(len(geno)):
		if geno[i]!="2":
			a+=geno[i]
			b+=geno[i]
		else:
			a+=random.choice(["0","1"])
			b+="0" if a=="1" else "1"
			# print a, b
			# a+="0"
			# b+="1"
	return a, b

def get_compl(hyplo, geno):
	res=""
	for i in xrange(0,len(geno)):
		if geno[i]=="2" and hyplo[i]=="1":
			res+="0"
		elif geno[i]=="2" and hyplo[i]=="0":
			res+="1"
		else:
			res+=geno[i]
	return res


def get_ability(hyplo, indivs):
	"""
	get the num of individuals one hyplotype can explain
	"""
	count=0
	for indiv in indivs:
		if match(hyplo, indiv):
			count+=1
	return count

def match(hyplo, indiv):
	"""
	match one hyplo to a indiv.
	"""
	for i in xrange(0,len(hyplo)):
		if hyplo[i]=="0" and indiv[i]=="1":
			return False
		elif hyplo[i]=="1" and indiv[i]=="0":
			return False
	return True

def check_fixed(indiv):
	for i in xrange(0,len(indiv)):
		if indiv[i]=="2":
			return False
	return True
			
def check(indivs, geno_list):
	"""
	check if all the indivs in the geno_list.
	"""
	for i in xrange(0,len(indivs)):
		if indivs[i] not in geno_list:
			# print "this is not in: "+ indivs[i]
			return False
	return True
def make_pool(num_snp):
	"""
	enumerate all of the possible hyplotypes.
	"""
	c=0
	pool=[]
	for i in xrange(0,num_snp+1):
		s=make_str(i, num_snp)
		pool+=map("".join, itertools.permutations(s, num_snp))
	return list(set(pool))

def make_str(n, len):
	"""
	make 01 string base on the # of 1.
	"""
	s=""
	for i in xrange(0,n):
		s+="1"
	for i in xrange(n,len):
		s+="0"
	return s

def tune_em(indivs, n):
	"""
	vary the starting point for em algorithm.
	"""
	best_phasing_list, best_phasing_pairs=em_phase(indivs)
	best_len=len(best_phasing_list)
	for i in xrange(1, n):
		phase_list, phase_pairs=em_phase(indivs)
		if len(phase_list)<best_len:
			best_len=len(phase_list)
			best_phasing_list=phase_list
			best_phasing_pairs=phase_pairs
	print "in tune_em", len(best_phasing_list)
	return best_phasing_list, best_phasing_pairs

def get_chunks(indivs, k):
	"""
	get all of the k-mer chunks in pair.
	"""
	pair_chunk_collection=[]
	for i in xrange(0, len(indivs[0])-k+1, k):
		chunks=[]
		for x in indivs:
			chunks.append(x[i:i+k])
		partial_phase_pairs=tune_em(chunks, 5)[1]
		print partial_phase_pairs
		pair_chunk_collection.append(partial_phase_pairs)
	return pair_chunk_collection

def pairwise(iterable):
    """ 
    -> (s0,s1), (s2,s3), (s4, s5), ...
    """
    a = iter(iterable)
    return izip(a, a)

def phase_large(indivs, k):
	"""
	assemble the the chunks for each individual. 
	"""
	pair_chunks=get_chunks(indivs, 10)
	indiv_chunks=[]
	print "pair chunks ready"
	phased_indivs=[]
	res=[]
	for i in xrange(0, len(indivs[0])-k+1, k):
		indiv_chunks.append([x[i:i+k] for x in indivs])
	# for (x, y),(a,b) in pairwise()
	res_indivs=indiv_chunks[0]
	res_pairs=pair_chunks[0]
	# print pair_chunks
	# print "this is the test of # pair chunks and indivs chunks."
	# for i in xrange(1, len(indiv_chunks)):
	# 	print "the i is : ", i
	# 	print len(pair_chunks[i]), len(indiv_chunks[i])
	for i in xrange(1, len(indiv_chunks)):
		res, res_pairs=get_best_combinations(res_indivs, indiv_chunks[i], res_pairs, pair_chunks[i])
		print res_pairs
		for j in xrange(0, len(indivs)):
			res_indivs[j]+=indiv_chunks[i][j]
	return res, res_pairs




def get_end_chunck(indiv_pair, k):
	return [indiv_pair[0][-1*k:],indiv_pair[1][-1*k:]]
def get_front_chunck(indiv_pair, k):
	return [indiv_pair[0][:k],indiv_pair[1][:k]]
def get_score(h1, h2, pool):
	if check_pool(h1, pool) and check_pool(h2, pool):
		return 2
	elif check_pool(h1, pool) or check_pool(h2, pool):
		return 1
	else:
		return 0
def check_pool(h, pool):
	for x in pool:
		if h in x:
			return True
	return False
	
def switch(hyplo_1, hyplo_2, pos):
	"""
	switch the the part of hyplo_types from the pos to the end.
	"""
	hyplo_1=list(hyplo_1)
	hyplo_2=list(hyplo_2)
	temp=hyplo_1[pos:]
	hyplo_1[pos:]=hyplo_2[pos:]
	hyplo_2[pos:]=temp
	hyplo_1="".join(hyplo_1)
	hyplo_2="".join(hyplo_2)
	return hyplo_1, hyplo_2

def switch_distance(pair_1, pair_2):
	"""
	calculate the swich distance.
	"""
	h1_1=pair_1[0]
	h2_1=pair_2[0]
	h2_2=pair_2[1]
	cnt_1=0
	cnt_2=0
	# print h2_1, h1_1
	for i in xrange(0,len(pair_1[0])):
		if h2_1[i]!=h1_1[i]:
			h2_1, h2_2=switch(h2_1, h2_2, i)
			# print h2_1, h2_2, i
			cnt_1+=1
	print h2_1, h1_1
	h1_1=pair_1[0]
	h2_1=pair_2[0]
	h2_2=pair_2[1]
	for i in xrange(0,len(pair_1[0])):
		if h2_2[i]!=h1_1[i]:
			h2_1, h2_2=switch(h2_1, h2_2, i)
			# print h2_1, h2_2, i
			cnt_2+=1
	# print cnt_1, cnt_2
	return cnt_1 if cnt_1<cnt_2 else cnt_2
def avg_distance(pairs_1, pairs_2):
	# all_sw=map(lambda x: switch_distance(x[0],x[1]), pairs)
	sws=[]
	for i in xrange(0,len(pairs_1)):
		sws.append(switch_distance(pairs_1[i], pairs_2[i]))
	print sws
	print reduce(lambda x, y: x + y, sws) / float(len(sws))
def count_switch(pair_1, pair_2):
	h1_1=pair_1[0]
	h2_1=pair_2[0]
	h2_2=pair_2[1]
	cnt=0
	for i in xrange(0,len(pair_1[0])):
		if h2_1[i]!=h1_1[i]:
			h2_1, h2_2=switch(h2_1, h2_2, i)
			# print h2_1, h2_2, i
			cnt+=1
	return cnt

# def is_pow2(num):
# 	"""
# 	check if num is power of 2.
# 	"""
# 	return num!=0 and (num&(num-1)==0)

def count_flip(pair_1,pair_2):
	h1_1=pair_1[0]
	h2_1=pair_2[0]
	cnt=0
	for i in xrange(0,len(pair_1[0])):
		if h2_1[i]!=h1_1[i]:
			cnt+=1
	return cnt
if __name__ == '__main__':
	start = time.clock()
	f = open('indiv_data.txt', 'r')
	# p=open('pairs_data.txt', 'r')
	individuals=[]
	# pairs=[]
	for line in f:
		line = line.strip()
		individuals.append(line)
	# for line in p:
	# 	line = line.strip()
	# 	pairs.append([line[0:10], line[11:21]])
	
	if len(sys.argv)>1 and sys.argv[1]=="hard":
		l=phase_large(individuals, 10)[0]
		print "# of original individuals is: ", len(individuals)," which generated by 50 hyplotypes"
		print "# of hyplotypes from EM algorithm for 1000 SNP is : ", len(l)
		print time.clock() - start, 'seconds'
	else:
		print "# of original individuals is: ", len(individuals)," which generated by 10 hyplotypes"
		l1=greedy_phase(individuals)
		time1=time.clock() - start
		print "# of hyplotypes from greedy algorithm is : ", len(l1)
		print "time spend is ",time1, 'seconds'
		l2,p=clark_phase(individuals)
		time2=time.clock()-time1
		print "# of hyplotypes from clark's method is : ", len(l2)
		print "time spend for clark's is ",time2, 'seconds'
		l3, p=tune_em(individuals, 10)
		print "# of hyplotypes from EM algorithm is : ", len(l3)
		time3=time.clock()-time1-time2
		print "time spend by EM is ",time3, 'seconds'
		# le=trivial_phase(individuals)[0]
		# timee=time.clock() - time1-time2-time3
		# print "# of hyplotypes from exhaustive search is : ", len(le)
		# print "time spend is ",timee, 'seconds'
	