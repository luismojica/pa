#!/usr/bin/python
# 11/06/2016
# Luis Mojica
# UTD
# Implement parallel partition and extension using cosets and existing partitions
# 01/25/2017 Clean up

import argparse
import sys, getopt
from utils import parse_pa,covered,hd_pairwise,pa2str,write_file,parse_partition,parse_repsentatives,rotate,create_cosets,extend_cosets,get_nd
from pe_sud import pe_blocks as pe_sud_fun
from pe_ilp import pe_blocks as pe_ilp_fun
from pe_ilp import partition_extension_blocks as pe_ilp_fun_cov
from os.path import isfile,join
from pgl_blocks import get_pgl_blocks

# Get blocks coverage
def block_coverage(block,ps,qs,new_symbols=None):
    coverd_perm=[]
    for perm in block:
        for idx in xrange(len(ps)):
            if new_symbols is not None:
                perm=covered(perm,ps[idx],qs[idx],new_symbols[idx])
            else:
                perm=covered(perm,ps[idx],qs[idx])
            # If one of the symbols is not included, discard this permutation
            if perm is None:
                break
        if perm is not None:
            coverd_perm.append(perm)
    return coverd_perm

''' Extende in the right position, for the freebies'''
def extend(block,ps,new_symbols):

    for perm in block:
        for idx in xrange(len(ps)):
            perm[ps[idx]]=new_symbols[idx]

''' Single rotation of a an iterable'''
def signle_rot(iterable):
    return iterable[1:]+iterable[0:1]

''' Build blocks from a PA'''
def create_blocks(pa,num_blokcs,block_size=None):
    blocks=[]
    start=0
    num_symbols=len(pa[0])
    if len(pa)%float(num_symbols)==0:
        print 'AGL Block'
        block_size=num_symbols

    for x in xrange(num_blokcs):
        blocks.append(pa[start:start+block_size])
        start=block_size*x+1

    return blocks

''' Get a new permutation array using paralle partition and extension: '''
def partition_extension2(pa,num_new_sym,num_blocks,reps=None,P=None,Q=None,modify_pa=True,use_ILP=False):

    # New symbols:
    new_symbols=range(len(pa[0]),len(pa[0])+num_new_sym)
    print 'new_symbols',new_symbols

    # Limit the nubmer of cosets to compute based on the partitions that we have:
    num_parts=num_blocks
    needed_cosets=num_parts+num_new_sym # Total, including freebies
    print 'needed_cosets_with_freebies:',needed_cosets

    # Get cosets:
    if reps is not None:
        if needed_cosets>len(reps)+1:
            print 'Not enough cosets/representatives. %i needed, only %i available.' %(needed_cosets,len(reps))
            return 0,None,None,None
            exit()
        cosets=create_cosets(pa,reps,needed_cosets)

    else:
        # Construct 'cosets' form AGL/MOLS:'
        cosets=create_blocks(pa,needed_cosets)

    if num_blocks<num_new_sym:
        print 'Less blocks than the number of symbols'
        return 0,None,None,None

    # Parallel partition is reguluar partition and extension in this case :|
    if num_new_sym==1 and use_ILP:
        _,_,coverage=pe_ilp_fun(cosets,num_blocks=num_parts)
        return coverage,None,None,None

    elif use_ILP:
        # _,_,coverage=pe_ilp_fun_cov(cosets,num_blocks=num_parts)
        _,_,coverage=pe_ilp_fun(cosets,num_blocks=num_parts)
        coverage=coverage+num_new_sym*len(cosets[0])
        print 'coverage',coverage
        return coverage,None,None,None

    else:
        # Compute P and Q on the cosets:
        P,Q=pe_sud_fun(cosets,num_blocks=num_parts)
        print '_______'

    # If not previously extended, do it here:
    if modify_pa:
        extend_cosets(cosets,num_new_sym=num_new_sym)

    # Rotation of the new symbols to be introuduced:
    cycles=rotate(new_symbols)
    print 

    # Output container:
    out=[]

    # Extend:
    curr_part=0
    for idx in xrange(num_parts):

        # Get te corresponding parts for this block:
        PS=P[0:num_new_sym]
        QS=Q[0:num_new_sym]

        P=signle_rot(P)
        Q=signle_rot(Q)

        if modify_pa:
            extended_block=block_coverage(cosets[idx],PS,QS,new_symbols)
        else:
            extended_block=block_coverage(cosets[idx],PS,QS)

        # Print block coverage:
        cov_per=len(extended_block)/float(len(cosets[idx]))
        print idx,"%.3f" %(cov_per*100)
        out+=extended_block

    # Freebies:
    for idx in xrange(num_parts,needed_cosets):
        if modify_pa:
            if len(cycles)==1: # When extending by only one symbol is used
                cycles=[cycles]
            extend(cosets[idx],new_symbols,cycles[idx%num_new_sym])
        out+=cosets[idx]

    print "pa_size:",len(out)
    return len(out),out,P,Q

''' Do partition and extension on PGL'''
def partition_extension_pgl(pgl_file,q,num_blocks,modify_pa=True,use_ILP=False):

    # Since it is PGL one can add at most 2 new symbols:
    num_new_sym=2

    # New symbols:
    new_symbols=range(q+1,q+1+num_new_sym)

    # # Number of symbols:
    num_symbols=q+1

    # Number of parts:
    num_parts=num_blocks

    # Number of lines from the PA (asumming in order):
    end=(num_parts+2)*num_symbols # The +2 is for the freebies
    pgl_pa=parse_pa(pgl_file,end)

    # Get the pgl blocks:
    blocks=get_pgl_blocks(pgl_pa,q,num_parts+2)

    # Run parallel pe
    len_out,out,P,Q=parallel_pe_blocks(blocks,num_blocks,num_new_sym,modify_pa=True,use_ILP=use_ILP)

    print "pa_size:",len_out
    return len_out,out,P,Q


''' Do partition and extension on given blocks, probably from kronecker product'''
def parallel_pe_blocks(blocks,num_blocks,num_new_sym,modify_pa=True,use_ILP=False):

    # New symbols:
    new_symbols=range(len(blocks[0][0]),len(blocks[0][0])+num_new_sym)

    # Number of parts:
    num_parts=num_blocks

    # Compute P and Q on the blocks:
    if use_ILP:
        P,Q=pe_ilp_fun(blocks[:-num_new_sym],num_blocks=num_parts)
    else:
        P,Q=pe_sud_fun(blocks,num_blocks=num_parts)
    
    # If not previously extended, do it here:
    if modify_pa:
        extend_cosets(blocks,num_new_sym=num_new_sym)

    cycles=rotate(new_symbols)
    print 

    # Output container:
    out=[]

    # Extend:
    curr_part=0
    for idx in xrange(num_parts):

        # Get te corresponding parts for this block:
        PS=P[0:num_new_sym]
        QS=Q[0:num_new_sym]

        # # print idx,PS
        P=signle_rot(P)
        Q=signle_rot(Q)

        if modify_pa:
            extended_block=block_coverage(blocks[idx],PS,QS,new_symbols)
        else:
            extended_block=block_coverage(blocks[idx],PS,QS)

        cov_per=len(extended_block)/float(len(blocks[idx]))
        print idx,"%.3f" %(cov_per*100)
        out+=extended_block
        # print hd_pairwise(extended_block)

    # Freebies:
    for idx in xrange(num_parts,num_parts+2):
    # for idx in xrange(3):
        if modify_pa:
            if len(cycles)==1: # When extending by only one symbol is used
                cycles=[cycles]
            extend(blocks[idx],new_symbols,cycles[idx%num_new_sym])
        out+=blocks[idx]

    # print pa2str(out)
    # print "pa_size:",len(out)
    # print hd_pairwise(out)
    return len(out),out,P,Q


def main(argv):

    parser = argparse.ArgumentParser(description='Crete and ILP model for the partitioning problem')
    parser.add_argument("--in_pa_file", required=False, help='The Permutation Array File, a group')
    parser.add_argument("--sought_hd", type=int, required=False, help='The Permutation Array File')
    parser.add_argument("--num_blocks", type=int, required=False, help='The number of blocks to look for')
    parser.add_argument("--out_pa_file", required=False, help='The filename to save the resulting permutation Array')
    parser.add_argument("--out_partition_file", required=False, help='The filename to save the patition for permutation Array ')
    parser.add_argument("--in_partition", required=False, help='A file containing the partion system to be used')
    parser.add_argument("--representatives", required=False, help='A file containing a list of coset representatives of the group to be extended')
    parser.add_argument("--num_new_sym", type=int, required=False, help='The number of new symbols to extend to')
    parser.add_argument("--use_ILP", action='store_true', required=False, help='Use the ILP PE method instead of the greedy algorithm')
    parser.add_argument("--out_dir", required=False, help='A filename to keep track of current new bounds')
    parser.add_argument("--q", type=int, required=False, help='The prime [power] in PGL(1,q)')
    parser.add_argument("--in_pgl_file", required=False, help='The Permutation Array File from PGL')


    # Get the arguments:
    args = parser.parse_args()

    P=None;Q=None;
    if args.in_partition:
        P,Q=parse_partition(args.in_partition)

    # Retrieve coset representatives
    if args.representatives:
        reps=parse_repsentatives(args.representatives)
    else:
        reps=None

    # Do partition and extension on the given cosets of group:
    if args.in_pa_file and reps:

        # Try this: because we want only a part of it
        nd_tuple=get_nd(args.in_pa_file)
        if nd_tuple is not None:
            n,d=nd_tuple

        # Get the group, all of it!
        pa=parse_pa(args.in_pa_file)[:]

        # Do partition and extension on coset representatives:
        pa_pe_size,ext_pa,P,Q=partition_extension2(pa,args.num_new_sym,args.num_blocks,reps,P,Q,modify_pa=True)

    # Do PE on PGL
    elif args.in_pgl_file and args.q and args.num_blocks:
        len_out,out,P,Q=partition_extension_pgl(args.in_pgl_file,args.q,args.num_blocks,use_ILP=args.use_ILP)
    
    # If a directory is given to save the resulting files:
    if args.out_dir:
        pa_fname=join(args.out_dir,str(n+args.num_new_sym)+'_'+str(d+args.num_new_sym)+'.txt')
        part_fname=join(args.out_dir,str(n+args.num_new_sym)+'_'+str(d+args.num_new_sym)+'_partition.txt')

        pa_str=pa2str(ext_pa)
        # Save to file
        print pa_fname
        write_file(pa_fname,pa_str)

        print part_fname
        # Create the p partition and q partition as strings:
        p_q_str='P:'+str(P)+'\n'
        p_q_str+='Q:'+str(Q)+'\n'
        write_file(part_fname,p_q_str)
    
    # Save the extended PA
    if args.out_pa_file:
        pa_str=pa2str(ext_pa)
        # Save to file
        write_file(args.out_pa_file,pa_str)

    # # Save the partition:
    if args.out_partition_file:

        # Create the p partition and q partition as strings:
        p_q_str='P:'+str(P)+'\n'
        p_q_str+='Q:'+str(Q)+'\n'

        # Save to file
        write_file(args.out_partition_file,p_q_str)


if __name__ == '__main__':
    main(sys.argv[1:])
