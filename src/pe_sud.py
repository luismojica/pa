#!/usr/bin/python2.7 -S
# 06/24/2016
# Luis Mojica
# UTD
# Partition and Extension Algorithm proposed by Dr. Sudborough
# 11/07/2016. Added Partition and Extension if blocks are given 

import argparse
from utils import hd_pairwise,hd_perm_list,hd_btw_lists,write_file,parse_original_bounds,parse_pa,parse_repsentatives,create_cosets,extend_cosets,extend,pa2str
from utils import add_symbol as add_symbol_fun
import math
from scipy.spatial.distance import hamming as hd
from random import shuffle
import numpy as np
import os.path
import math # For the square root
from pgl_blocks import get_pgl_blocks

''' Check the max coverage'''
def check_coverage(mat_block,pos,symbols):
	
	covered_perm=[]
	delete_idx=[]

	for idx, element in enumerate(mat_block[:,pos]):
		if element in symbols:
			covered_perm.append(mat_block[idx,:].tolist())
			delete_idx.append(idx)
			# print mat_block[idx,:].tolist()

	# Delete from matrix:
	mat_block=np.delete(mat_block,delete_idx,axis=0)

	return mat_block,covered_perm

''' Partition and Extension if blocks are given 11/07/2016'''
def pe_blocks(blocks,num_blocks):

	# Number of parts:
	num_parts=num_blocks

	# New symbol: !!! Note the difference, no longe assuming square blocks! this may not be AGL
	new_symbol=len(blocks[0][0])

	# Block matrices:
	mat_blocks=[]
	
	# Make the blocks matrices:
	for x in range(num_parts):
		mat_blocks.append(np.array(blocks[x]))

	# Search for the best parition:
	best_p=None
	best_q=None
	best_pa=None
	best_pa_len=0

	# Create the partition shell
	P=[[] for x in range(num_parts)]
	Q=[[] for x in range(num_parts)]

	# Num elements in each default partition:
	def_elemnts=int(math.floor(new_symbol/float(num_parts)))
	# Check if the number of elements is an integer or not
	is_int=(new_symbol/float(num_parts)).is_integer()

	# Covered permutations:
	pa=[]

	# For when num of symbols is not integer:
	q_start=0
	q_add=0

	# Add the default part:
	for x in range(num_parts):
		P[x].append(x)
		if True and not is_int:
			if x%2==0:
				Q[x]+=range(q_start,q_start+(def_elemnts))
				q_start=Q[x][-1]+1
			else:
				Q[x]+=range(q_start,q_start+(def_elemnts))
				q_start=Q[x][-1]+1
			# print q_start
		else:
			Q[x]+=range((def_elemnts)*x,(def_elemnts)*x+(def_elemnts))
		
		# Keep track of the number of symbols in the partition Q
		q_add+=len(Q[x])
		
		mat_blocks[x],cov_perm=check_coverage(mat_blocks[x],P[x][0],Q[x])

	# Make sure tha all symbols are included in the Q partitioning
	if q_add<new_symbol:
		# get the missing symbols:
		Q[0]+=range(q_add,new_symbol)

	# Go over each of the remaining columns (symbols):
	for x in range(num_parts,new_symbol):

		# Test which block gets more coverage:
		tmp_mats=[]
		tmp_cov_perm=[]
		len_coverage=[]
		
		# Check each block
		for y in range(len(mat_blocks)):
			tmp_mat,tmp_perms=check_coverage(mat_blocks[y],x,Q[y])
			tmp_mats.append(tmp_mat)
			tmp_cov_perm.append(tmp_perms)
			len_coverage.append(len(tmp_perms))

		# Get the max block index:
		mat_idx=len_coverage.index(max(len_coverage))

		# Assign this position to P_mat_idx:
		P[mat_idx].append(x)

		# Delete the permutations in M once we found the position and values:
		mat_blocks[mat_idx]=tmp_mats[mat_idx]

		# Early termination
		if len(pa)==new_symbol*num_parts:
			break

	# Create the p partition and q partition as strings:
	p_q_str='P:'+str(P)+'\n\n'
	p_q_str+='Q:'+str(Q)+'\n'

	return P,Q

# Run the partition and extension algorithm on blocks:
def partition_extension_blocks(blocks):

    P,Q=pe_blocks(blocks[:-1],len(blocks)-1)
    freeby=blocks[-1]

    # New symbol: !!! Note the difference, no longe assuming square blocks! this may not be AGL
    new_symbol=len(blocks[0][0])

    # Extend the Permutation Array:
    pa_ext=extend(blocks[:-1],P,Q,new_symbol)

    # Extend the freeby:
    pa_ext+=extend([freeby],[[new_symbol]],[xrange(new_symbol)],new_symbol)

    # Create the p partition and q partition as strings:
    p_q_str='P:'+str(P)+'\n'
    p_q_str+='Q:'+str(Q)+'\n'

    # Convert to string:
    pa_str=pa2str(pa_ext)

    return pa_str,p_q_str,len(pa_ext)

# Run the partition and extension algorithm:
def partition_extension(pa_file,num_symbols=None,max_iter=5,num_blocks=None):

	# Num of partitions to make:
	sqr_sym=int(math.floor(math.sqrt(num_symbols)))

	# Number of parts:
	num_parts=min(sqr_sym,num_blocks-1)
	# print 'num_parts',num_parts

	# Number of lines from the PA (asumming in order):
	end=(num_parts+1)*num_symbols # The +1 is for the freeby
	pa=parse_pa(pa_file,end)


	# Find blocks
	if num_blocks is not None:
		blocks=[]
		elem_per_block=len(pa)/(num_parts+1)
		for b_idx in range((num_parts+1)):
			sough_dist_block=pa[elem_per_block*b_idx:elem_per_block*b_idx+elem_per_block]
			blocks.append(sough_dist_block)

	# New symbol:
	new_symbol=len(blocks[0])

	# Block matrices:
	mat_blocks=[]
	
	# Make the blocks matrices:
	for x in range(num_parts):
		mat_blocks.append(np.array(blocks[x]))

	# Search for the best parition:
	best_p=None
	best_q=None
	best_pa=None
	best_pa_len=0

	# Create the partition shell
	P=[[] for x in range(num_parts)]
	Q=[[] for x in range(num_parts)]

	# Num elements in each default partition:
	def_elemnts=int(math.floor(new_symbol/float(num_parts)))
	# Check if the number of elements is an integer or not
	is_int=(new_symbol/float(num_parts)).is_integer()

	# Covered permutations:
	pa=[]

	# For when num of symbols is not integer:
	q_start=0
	q_add=0

	# Add the default part:
	for x in range(num_parts):
		P[x].append(x)
		if True and not is_int:
			if x%2==0:
				Q[x]+=range(q_start,q_start+(def_elemnts)+1)
				q_start=Q[x][-1]+1
			else:
				Q[x]+=range(q_start,q_start+(def_elemnts))
				q_start=Q[x][-1]+1
			# print q_start
		else:
			Q[x]+=range((def_elemnts)*x,(def_elemnts)*x+(def_elemnts))
		
		# Keep track of the number of symbols in the partition Q
		q_add+=len(Q[x])
		
		mat_blocks[x],cov_perm=check_coverage(mat_blocks[x],P[x][0],Q[x])
		pa+=cov_perm

	# Make sure tha all symbols are included in the Q partitioning
	if q_add<new_symbol:
		# get the missing symbols:
		Q[0]+=range(q_add,new_symbol)

	# Go over each of the remaining columns (symbols):
	for x in range(num_parts,new_symbol):

		# Test which block gets more coverage:
		tmp_mats=[]
		tmp_cov_perm=[]
		len_coverage=[]
		
		# Check each block
		for y in range(len(mat_blocks)):
			tmp_mat,tmp_perms=check_coverage(mat_blocks[y],x,Q[y])
			tmp_mats.append(tmp_mat)
			tmp_cov_perm.append(tmp_perms)
			len_coverage.append(len(tmp_perms))

		# Get the max block index:
		mat_idx=len_coverage.index(max(len_coverage))

		# Assign this position to P_mat_idx:
		P[mat_idx].append(x)

		# Delete the permutations in M once we found the position and values:
		mat_blocks[mat_idx]=tmp_mats[mat_idx]

		# Add symbol:
		ext_perms=add_symbol(tmp_cov_perm[mat_idx],x,new_symbol)

		# Add the selected block permutations to the pa:
		pa+=ext_perms

		# Early termination
		if len(pa)==new_symbol*num_parts:
			break


	# Add the freeby block:
	ext_perms=add_symbol(blocks[-1],new_symbol-1,new_symbol)
	pa+=ext_perms

	# Create the partition as a string
	pa_str=pprint_pa(pa)

	# Create the p partition and q partition as strings:
	p_q_str='P:'+str(P)+'\n\n'
	p_q_str+='Q:'+str(Q)+'\n'

	# The size of the pa:
	pa_len=len(pa)

	return pa_str,p_q_str,pa_len

def main():

	# Parse function arguments
	parser = argparse.ArgumentParser(description='Find controversial comments and save them as seeds')
	parser.add_argument("--in_pa_file", required=False, help='The Permutation Array File')
	parser.add_argument("--sought_hd", type=int, required=False, help='The Permutation Array File')
	parser.add_argument("--out_pa_file", required=False, help='The filename to save the resulting permutation Array')
	parser.add_argument("--out_partition_file", required=False, help='The filename to save the patition for permutation Array ')

	# Do a list of PAs and compute new bouns
	parser.add_argument("--in_pa_folder", required=False, help='Folder to find PAs')
	parser.add_argument('-l', '--l_num_symbols', nargs='+', type=int, help='List of mols to partition and extend')
	parser.add_argument("--out_extended_pa", required=False, help='Folder to store newly created PAs')
	parser.add_argument("--out_extended_partitions", required=False, help='Folder to store newly created PAs\' partitions')
	parser.add_argument("--max_iter", type=int, required=False, help='Maximum number of random partitions made by the algorithm')
	parser.add_argument("--in_original_bounds", required=False, help='Original Bounds')
	parser.add_argument("--out_new_bounds_file", required=False, help='A filename to store the new bounds, if any')
	parser.add_argument("--num_blocks", type=int, required=False, help='The number of blocks to look for')
	parser.add_argument("--in_block_file", required=False, help='The Permutation Array File')
	parser.add_argument("--representatives", required=False, help='A file containing a list of coset representatives of the group to be extended')
	parser.add_argument("--add_symbol", action='store_true',required=False, help='Increase the number of symbols by one')
	parser.add_argument("--q", type=int, required=False, help='The prime [power] in PGL(1,q)')
	parser.add_argument("--in_pgl_file", required=False, help='The Permutation Array File from PGL')
	parser.add_argument("--use_pgl_blocks", action='store_true',required=False, help='Compute sequential parition and extension on blocks from PGL(2,q)')

	args = parser.parse_args()
	max_iter=None

	if args.l_num_symbols and args.in_pa_folder:

		# Process this list
		pa_len_dic=list_partition_extension(args.l_num_symbols,in_pa_folder=args.in_pa_folder,out_extended_pa=args.out_extended_pa,out_extended_partitions=args.out_extended_partitions,max_iter=args.max_iter)
		
		print pa_len_dic

		# Compute new bounds and save them to file if any:
		if args.in_original_bounds and args.out_new_bounds_file and pa_len_dic is not None:

			# Retrieve the original bounds:
			orig_bounds=parse_original_bounds(args.in_original_bounds)

			# Check for new bounds and create a string to save to file
			new_bounds_str=check_new_bounds(orig_bounds,pa_len_dic)

			# Save to file:
			if new_bounds_str is not None:
				# print args.out_new_bounds_file
				write_file(args.out_new_bounds_file,new_bounds_str)

		print 'Done!'


	# Get the pa from file
	elif args.in_pa_file:
		pa=parse_pa(args.in_pa_file)
		# print pa

		# Do partition and extension
		pa_str,p_q_str,pa_len=partition_extension(pa,sought_hd=args.sought_hd,num_blocks=args.num_blocks)

		# Save the extended PA
		if args.out_pa_file:

			# Save to file
			write_file(args.out_pa_file,pa_str)

		# # Save the partition:
		if args.out_partition_file:

			# Save to file
			write_file(args.out_partition_file,p_q_str)

		if args.out_pa_file is None:
			print pa_len

	# Do PE on blocks instead with coset representatives:
	elif args.in_block_file and args.representatives and args.num_blocks:
		
		# Retrieve representatives:
		reps=parse_repsentatives(args.representatives)

		# Get the block from file
		pa=parse_pa(args.in_block_file)

		if args.add_symbol:
			pa=add_symbol_fun(pa)

		# Get cosets:
		cosets=create_cosets(pa,reps,args.num_blocks)

		pa_str,p_q_str,pa_len=partition_extension_blocks(cosets)

		# Save the extended PA
		if args.out_pa_file:

			# Save to file
			write_file(args.out_pa_file,pa_str)

		# # Save the partition:
		if args.out_partition_file:

			# Save to file
			write_file(args.out_partition_file,p_q_str)

		print pa_len
		


if __name__ == '__main__':
	main()
