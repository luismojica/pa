#!/usr/bin/python2.7 -S
# 06/15/2016
# Luis Mojica
# UTD
# Permutation Arrays Utilities

import numpy as np

''' Standard hamming distance on a pair of lists'''
def hd(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(bool(ch1 != ch2) for ch1, ch2 in zip(s1, s2))

''' Check the pairwise distance between the elements of the given list
	return the smallest hd'''
def hd_pairwise(a_list):
	lowest_hd=float('inf')
	for x in xrange(len(a_list)):
		for y in xrange(x+1,len(a_list)):
			hd_val=hd(a_list[x],a_list[y])
			if hd_val<lowest_hd:
				lowest_hd=hd_val

	return lowest_hd


''' Check between x permutation and the elements of the given list
	return the smallest hd'''
def hd_perm_list(perm, a_list):
	lowest_hd=float('inf')
	for p in a_list:
		hd_val=hd(perm,p)

		if hd_val<lowest_hd:
			lowest_hd=hd_val

	return lowest_hd


''' Hamming Distance between the elements of two lists'''
def hd_btw_lists(a_list,b_list,s_hd=None):
	lowest_hd=float('inf')
	for p in a_list:
		hd_val=hd_perm_list(p,b_list)
		# if s_hd is not None:
		# 	# if hd_val<s_hd:
		# 	# 	print p,b_list
		# 	# 	print hd_val
		# 	# 	exit()
		if hd_val<lowest_hd:
			lowest_hd=hd_val

	return lowest_hd

''' Parse a PA assuming that each symbol is an integer, return a list of lists'''
def parse_pa(filename,end=None):
	# Return list
	pa=[]

	with open(filename, 'r') as f:
		if end is None:
			for line in f:
				# print line.strip()
				p=[int(x) for x in line.strip().split()]
				if p:
					pa.append(p)
		else:
			for idx,line in enumerate(f):
				# print line.strip()
				pa.append([int(x) for x in line.strip().split()])
				if idx>end:
					break

    # Check the range:
	if pa:
		min_elemn=min(pa[0])
		max_elemn=max(pa[0])

	if min_elemn>0:
		print '************* Permutation Array with symbols from %i to %i insted from 0 to n-1' % min_elemn,max_elemn

	return pa

''' Parse list of the number of MOLS per of side x'''
# Parse the mols file
def parse_mols(filename):
	# Return list
	pa=[]

	with open(filename, 'r') as f:
		for line in f:
			# print line.strip()
			pa=[int(x) for x in line.split(',')]
	return pa


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

''' Append to a file the given data'''
def append_file(filename, data):
	with open(filename, "a") as myfile:
		myfile.write(data)

''' Parse the mols file'''
def parse_original_bounds(filename):
	# Return list
	out={}

	with open(filename, 'r') as f:
		for line in f:
			if len(line.strip().split()):
				# print line.strip().split()
				n,n_1=line.strip().split()
				out[int(n)]=int(n_1)
			
	return out


''' Harmonize two bounds dictionary '''
def merge_bounds(b1,b2):

	# return 'hola'

	# Return list
	out={}

	# Get the keys of b1 and b2 and find the intersection:
	b1_keys=b1.keys()
	b2_keys=b2.keys()

	# What overlaps:
	intersection=list(set(b1.keys())&set(b2.keys()))
	
	# # Symmetric difference:
	# differnce=list(set(b1.keys())^set(b2.keys()))
	# Differences:
	b1_list=list(set(b1.keys())-set(b2.keys()))
	b2_list=list(set(b2.keys())-set(b1.keys()))

	# Populate the dictionary with the difference between b1 and b2
	for key in b1_list:
		out[key]=b1[key]
	for key in b2_list:
		out[key]=b2[key]

	# Merge what overlaps:
	for key in intersection:
		val1=b1[key]
		val2=b2[key]
		out[key]=max(val1,val2)


	return out

''' Format to a string a dictionary of bounds'''
def pprint_bounds(new_bounds):

	out=[]

	for idx in sorted(new_bounds):
		# print idx,' '.join([str(x) for x in list(new_bounds[idx])[1:3]])
		out.append(str(idx)+' '+str(new_bounds[idx]))

	return '\n'.join(out)
	

'''Pretty print a pa or a block'''
def pa2str(pa):

	pretty=[]
	for perm in pa:
		pretty.append(' '.join([str(x) for x in perm]))

	return '\n'.join(pretty)


''' Check if a permuation is covered '''
def covered(perm,p,q,transform=False,new_symbol=None):
	for position in p:
		for val in q:
			# print position,val
			if perm[position]==val:

				# Transform:
				if transform:
					# Make a copy of the permutation
					perm_copy=perm[:]
					perm_copy[position]=new_symbol
					perm_copy.append(val)

					return perm_copy

				return True

	return None


''' Add new symbol to the given PA'''
def add_symbol(pa):
    num_symbol=len(pa[0])

    for perm in pa:
        perm.append(num_symbol)

    return pa

''' If the class file is not present, compile'''
def compile_java():
    from os.path import isfile,join
    import os.path,subprocess
    from subprocess import STDOUT,PIPE

    if not isfile('GroupGeneration.class'):

        cmd = ['javac', 'GroupGeneration.java']
        proc = subprocess.Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
        stdout,stderr = proc.communicate()
        if stdout is not None:
            print stdout
            return False
    else:
        return True
    
''' Call Zac's program to generate the group'''
def get_pgl_java(p, pwr):
    from os.path import isfile,join
    import os.path,subprocess
    from subprocess import STDOUT,PIPE

    # java_class,ext = os.path.splitext(java_file)
    cmd = ['java', 'GroupGeneration','pgl',str(p),str(pwr)]
    proc = subprocess.Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
    stdout,stderr = proc.communicate()
    
    if stdout is "":
        return None
    elif stdout[:3]=='Not':
        print 'aqui'
        return -1
    return parse_group_table(stdout)

''' Zach's table'''
def parse_group_table(table):
    # print table
    # exit()
    out=[]
    for line in table.split('\n'):
        if line:
            # print line.split()
            p=[int(x) for x in line.strip().split()]
            if p:
                out.append(p)
    return out

''' Parse partition file: P,Q'''
def parse_partition(filename):
	from ast import literal_eval
	# Check if the filename is a already parsed PA:
	if type(filename) is list:
		return filename

	# Return list
	P=[]; Q=[]

	with open(filename, 'r') as f:
		for line in f:
			if line[0]=='#' or line=='\n':
				continue

			if line[0]=='P':
				P=literal_eval(line[2:])
			if line[0]=='Q':
				Q=literal_eval(line[2:])

	return P,Q


''' Parse from a textfile into a list, all cosets representatives in it'''
def parse_repsentatives(filename):
	# Check if the filename is a already parsed PA:
	if type(filename) is list:
		return filename

	# Return list
	reps=[]

	with open(filename, 'r') as f:
		for line in f:
			if line[0]=='#' or line=='\n':
				continue

			# Parse:
			str_rep=line.split()
			reps.append([int(x) for x in str_rep])

	return reps

''' Rotate a list'''
def rotate(iterable):
	if iterable is None:
		return None
	if len(iterable)==1:
		return iterable
	else:
		out=[]

		for idx in xrange(len(iterable)):
			# iterable=iterable[:] ?
			shift=iterable[1:]+iterable[0:1]
			out.append(shift)
			iterable=shift

	return out

''' Compose the group/pa with the coset representatives'''
def create_cosets(g,reps,num_cosets=-1,num_new_sym=0):

    cosets=[]
    
    for idx,r in enumerate(reps):
        cosets.append(compose_list(g,r,extend=num_new_sym))
        if num_cosets>0 and idx==num_cosets-2: # -2 Because we have the group itself
            break

    # On the actual group
    if num_new_sym:
        g_ext=[]
        for p in g:
            p+=[None]*num_new_sym
    cosets.append(g)

    return cosets

''' Extend cosets with n none elements'''
def extend_cosets(cosets,num_new_sym=0):
    
    if num_new_sym:
        for coset in cosets:
            for p in coset:
                p+=[None]*num_new_sym
    return cosets

''' Get n and d from filename'''
def get_nd(fname, reobj=None):

	if reobj is None:
		import re
		reobj = re.compile(r".*/(\d+)_(\d+).*")
	found=reobj.match(fname)
	if found:
		n=int(found.group(1))
		d=int(found.group(2))
		return (n,d)
	
	return found

''' With the P extend the permutation array'''
def extend(blocks,P,Q,new_sym):
    out=[]

    # print P,list(Q),new_sym

    for idx_block,block in enumerate(blocks):
        q=Q[idx_block]
        for perm in block:
            # Only add an element to the end:
            # if P[0] is not None and new_sym==P[0][0]:
            if P[0] and new_sym==P[0][0]:
                perm.append(new_sym)
                out.append(perm)
            # Insert:
            else:
                for pos in P[idx_block]:
                    # print P[idx_block]
                    if perm[pos] in q:
                        out.append(perm[:pos]+[new_sym]+perm[pos+1:]+[perm[pos]])
                        break # To count it only once

    # Verfiy correctness:
    # assert(hd_pairwise(out)==new_sym-1)
    # print hd_pairwise(out)
    # print len(out)

    return out