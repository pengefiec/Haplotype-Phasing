import random
import os
import sys
import time
from collections import defaultdict
def ge_hyplo(len):
	a=""
	for i in xrange(0,len):
		a+=str(random.randint(0,1))
	return a

def ge_hyplo_pool(num_hyplo, len):
	pool = [] 
	for x in xrange(0,num_hyplo):
		s=ge_hyplo(len)
		while s in pool:
			s=ge_hyplo(len)
		pool.append(s)
	return pool
def mix(s1, s2):
	s=""
	for i in xrange(0,len(s1)):
		if s1[i]=="1" and s2[i]=="1": 
			s+="1"
		elif s1[i]=="0" and s2[i]=="0"	:
			s+="0"
		else:
			s+="2"
	return s	
def ge_indiv(num_indiv, num_snp, num_hyplo):
	pool=ge_hyplo_pool(num_hyplo, num_snp)
	indivs=[]
	pairs=[]
	f=False
	for i in xrange(0,num_indiv):
		h1=random.choice(pool)
		h2=random.choice(pool)
		pairs.append(h1+h2)
		s=mix(h1, h2)
		while s in indivs:
			print "chunk hit"
			f=True
			h1=random.choice(pool)
			h2=random.choice(pool)
			s=mix(h1, h2)
		if f:
			pairs.pop()
			pairs.append(h1+h2)
			f=False
		indivs.append(s)
	return indivs, pairs

def ge_large_indiv(num_indiv, num_snp, num_hyplo, k):
	num_chunk=num_snp/k
	indivs=ge_indiv(num_indiv, k, num_hyplo)[0]
	for i in xrange(1, num_chunk):
		chunks=ge_indiv(num_indiv, k, num_hyplo)[0]
		while check_middle(chunks, indivs):
			print "middle hit"
			chunks=ge_indiv(num_indiv, k, num_hyplo)[0]
		for j in xrange(0, num_indiv):
			indivs[j]+=chunks[j]
	return indivs

def check_middle(chunks, indivs):
	temp=[]
	for i in xrange(0,len(indivs)):
		print chunks[i][:5]+indivs[i][-5:]
		temp.append(indivs[i][-5:]+chunks[i][:5])
	return len(temp)!=len(set(temp))



if __name__ == '__main__':
	# para: num of indiv, num of snp, num of haplotypes, (chunck size)
	output=[]
	if len(sys.argv)>1 and sys.argv[1]=="hard":
		output=ge_large_indiv(50, 1000, 10, 10)
	else:
		output=ge_indiv(20, 10, 10)[0]
	folder = './final_proj'
	output_n='indiv_data.txt'
	# output_p='pairs_data.txt'
	# output_fn=os.path.join(folder, output_n)
	with open(output_n, 'w') as output_file:
		 output_file.write('\n'.join(output))
	# print pairs
	# with open(output_p, 'w') as output_file:
	# 	 	output_file.write('\n'.join(pairs))