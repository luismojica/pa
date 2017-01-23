#!/usr/bin/python
# 08/15/2016
# Luis Mojica
# UTD
# Comput fast hamming distance

# from sage.all import designs
import argparse
import os.path
# from utils import write_file
from utils import pa2str,hd_pairwise,parse_pa,hd_perm_list,hd,add_symbol

''' Write to a file the given data'''
def write_file(filename, data):
	try:
		# Open, and read training file
	   	f = open(filename,'w');
	except IOError:
	   	print "Couldn't open %s" % filename
	   	exit()
	f.write(data)
	# Close file:
	f.close()

def plus_one(perm):
	return [x+1 for x in perm]

def fast_hd(n,pi):
	
	# Initialize:
	D=[n for m in xrange(n)]
	print range(n)
	print pi

	# Substract:
	for j in xrange(0,n):
		# # print j,pi[j],j-(pi[j]-j)
		# print j%n==0,(pi[j]-j)%n==0
		# print (j-(pi[j]-j))%n==0
		if (j-(pi[j]-j))%n==0:
			D[j]-=1
	print D

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

'''Compose the left coset of a permutation array and a permutation'''
def pa_p_mul(pa,p):
	out=[]
	for q in pa:
		out.append(mul(p,q))
	return out

''' As in the paper: pi_i^-1*pi'''
def pi_mul(p_i,p):
	return mul(inv(p_i),p)

# ''' Quick Check -- Orig'''
# def qcheck_i(pi_list,p,pa=None):
# 	n=len(p)
# 	if pa is None:
# 		cn=Cn(n)
# 		print 'No Pa'
# 		exit()

# 	min_hd=n

# 	for pi in pi_list:
# 		di=mul(inv(pi),p)
# 		# if di in pa:
# 		# 	min_hd=0
# 		# 	break
# 		# else:
# 		# 	min_hd=min(hd_perm_list(di,pa),min_hd)
# 		min_hd=min(hd_perm_list(di,pa),min_hd)

# 	return min_hd

''' Quick Check'''
def qcheck_i(pi_list,p,pa=None,s_hd=-1):
	n=len(p)
	# if pa is None:
	# 	cn=Cn(n)
	# 	print 'No Pa'
	# 	exit()

	min_hd=n
	# print 'min_hd',min_hd

	for pi in pi_list:
		di=mul(inv(pi),p)
		# if di in pa:
		# 	min_hd=0
		# 	break
		# else:
		# 	min_hd=min(hd_perm_list(di,pa),min_hd)
		hdist=hd_perm_list(di,pa,s_hd)
		# min_hd=min(hd_perm_list(di,pa),min_hd)
		
		min_hd=min(hdist,min_hd)

		# Eearly Termination:
		if min_hd<s_hd:
			break

	# print
	return min_hd

''' Verify that a set of cosets representatives are all pairwise at distance d from a permutation array'''
def rep_cosets_check(pi_list,pa,d):
	n=len(pa[0])
	min_hd=n

	# pi_list=pi_list[:]
	# print pi_list
	# exit()

	for i_idx in range(len(pi_list)):
		for j_idx in range(i_idx+1,len(pi_list)):
		# for j_idx in range(len(pi_list)):
			if i_idx==j_idx:
				continue
			# print i_idx,j_idx
			# print pi_list[i_idx],pi_list[j_idx]
			di=mul(inv(pi_list[i_idx]),pi_list[j_idx])
			hd=hd_perm_list(di,pa)
			min_hd=min(hd,min_hd)
			# print hd,'?',d," ;",i_idx,j_idx
			if hd<d:
				print hd,'<',d," ;",i_idx,j_idx
				# exit()

	return min_hd


def Cn(n):
	out=[]
	c=range(n)
	for x in range(n):
		c=c[1:]+[c[0]]
		out.append(c)
	return out

def main():

	# Parse function arguments
	parser = argparse.ArgumentParser(description='Retrieve PGL as Permuatation Arrays')
	parser.add_argument("--n", type=int, required=False, help='A positive Integer')
	parser.add_argument("--degree", type=int, required=False, help='The degree q, a prime power')
	parser.add_argument('--out_pa_dir', help='Directory to store the PA', required=False)
	parser.add_argument('-l', '--l_num_symbols', nargs='+', type=int)
	parser.add_argument('--in_pa_file', help='Group or PA to check against', required=False)
	parser.add_argument('--represenatives', help='Permutations Representatives of the cosets to check', required=False)
	parser.add_argument("--sought_hd", type=int, required=False, help='Desired Hamming distance')
	parser.add_argument("--add_symbol", action='store_true', required=False, help='Add a symbol at the end of each permutation of the PA')
	
	
	args = parser.parse_args()

	# # When a list of mols is asked, and the mols are save to a file:
	if args.in_pa_file and args.represenatives and args.sought_hd:
		pa=parse_pa(args.in_pa_file)

		# Add symbol:
		if args.add_symbol:
			pa=add_symbol(pa)

		# Check distances:
		reps=parse_pa(args.represenatives)

		print 'min pairwise hd:',rep_cosets_check(reps,pa,args.sought_hd)

	else:

		# 	# Process the list:
		# 	get_pgl_pa(args.n,args.degree,args.out_pa_dir)

		# Testing:
		sigma=[0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 6, 7]
		sigma=[0, 3, 2, 1]
		sigma=[1, 3, 2, 0]
		# pa=parse_pa('../data/12_10.txt')
		# print pa[100]

		# print len(pa)
		# print type(pa[10])
		# print type(sigma)

		# print pa[10]
		# print sigma

		# print 'hd:', hd_perm_list(sigma,pa)

		# Change the domain:
		# s=plus_one(sigma)
		# g1=plus_one(pa[30])
		# g2=plus_one(pa[1000])

		# # # Operations on permutations:
		# s = Permutation(s)
		# g1 = Permutation(g1)
		# g2 = Permutation(g2)

		# lhs=s.left_action_product(g1)
		# rhs=g2.left_action_product(~g1)

		# print 'hd(sg1,g2):',hd_perm_list(lhs,[g2])
		# print 'hd(s,g2g1^-1):',hd_perm_list(lhs,[g2])

		# s.left_action_product(g1)

		# ###### Fast HD:
		# fast_hd(12,sigma)
		# fast_hd(4,sigma)

		a = [3,8,5,10,9,4,6,1,7,2]
		b = [2,7,4,9,8,3,5,0,6,1]
		# Inverse of a permutation:
		# print inv([3,8,5,10,9,4,6,1,7,2])
		# print inv([2,7,4,9,8,3,5,0,6,1])

		a=[0,2,1]
		b=[2,1,0]

		print mul(a,b)
		print mul(b,inv(a))

		print pi_mul(a,b)
		sigma=pi_mul(a,b)
		print hd_perm_list(sigma,[b])

		# cn:
		# c=[1,2,3]
		# print 'Cn:',Cn(4)

		print 'qcheck_i',qcheck_i([[0,1,2,3],[2,3,1,0]],[2,3,1,0])
	

if __name__ == '__main__':
	main()



