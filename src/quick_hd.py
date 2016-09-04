#!/usr/bin/python
# 08/15/2016
# Luis Mojica
# UTD
# Comput fast hamming distance

# from sage.all import designs
import argparse
import os.path
from utils import pa2str,hd_pairwise,parse_pa,hd_perm_list,hd,add_symbol

def plus_one(perm):
	return [x+1 for x in perm]

def inv(perm):
	min_val=min(perm)
	inverse = [0] * len(perm)
	for i, p in enumerate(perm):
		if min_val==1:
			inverse[p-1] = i+1
		else:
			inverse[p] = i
	return inverse

def mul(a,b):
	min_val=min(a)
	b.extend(list(range(len(b), len(a))))

	if min_val==1:
		perm = [b[i-1] for i in a] + b[len(a):]
	else:
		perm = [b[i] for i in a] + b[len(a):]
	return perm

''' As in the paper: pi_i^-1*pi'''
def pi_mul(p_i,p):
	return mul(inv(p_i),p)

''' Quick Check'''
def qcheck_i(pi_list,p,pa=None):
	n=len(p)

	min_hd=n

	for pi in pi_list:
		di=mul(inv(pi),p)
		hdist=hd_perm_list(di,pa)
		min_hd=min(hdist,min_hd)
	
	return min_hd

def main():

	# Parse function arguments
	parser = argparse.ArgumentParser(description='Retrieve PGL as Permuatation Arrays')
	args = parser.parse_args()
	print 'Nothing to see here'
	

if __name__ == '__main__':
	main()



