#!/usr/bin/python
# 08/22/2016
# Luis Mojica
# UTD
# Creating the random permutation and random walk

import argparse
import sys, getopt
from utils import parse_pa,covered,hd_pairwise,pa2str,hd_perm_list,append_file
from utils import add_symbol as add_symbol_fun
from quick_hd import mul,qcheck_i,inv
from random import randint,random
import os.path

'''Create a Random permutation of n symbols'''
def rand_perm(n):
    perm=[]
    available_symbols=range(n)
    for x in range(n):
        # Draw a random number from a uniform distribution:
        syn_idx=randint(0,len(available_symbols)-1)
        
        new_symbol=available_symbols[syn_idx]
        # print new_symbol
        
        # Add to the permutation
        perm.append(new_symbol)

        # Remove the last symbol from the available symbols:
        available_symbols.pop(syn_idx)
        # print perm, available_symbols

    return perm

''' Make the transposition of two position in the permutation p'''
def transposition(p,idx1,idx2):
    p2=list(p) # Make sure not to modify p itself, it may be stored somewhere e.g. p_i
    v1=p2[idx1]
    v2=p2[idx2]
    p2[idx1]=v2
    p2[idx2]=v1
    return p2

''' A random transposition'''
def rand_transpositon(p):
    positions=range(len(p))
    idx1=randint(0,len(positions)-1)
    p1=positions[idx1]
    positions.pop(p1)
    idx2=randint(0,len(positions)-1)
    # print idx2
    p2=positions[idx2]

    # Make the transposition:
    return transposition(p,p1,p2)

'''A random permutation that is not in a PA'''
def p_not_in_pa(pa,num_symbols=None):
    # Number of symbols:
    if num_symbols is None:
        num_symbols=len(pa[0])

    # Get a random permutation:
    p=rand_perm(num_symbols)
    while p in pa:
        p=rand_perm(num_symbols)
    return p

''' Find new permutation using a random walk'''
def random_walk(pa_file, d, add_symbol, append_file_name, jump_prob, existing_representatives):
    if jump_prob == None:
        jump_prob = 0.7

    # Retrieve the Permutation Array
    pa=parse_pa(pa_file)
    # print pa

    if add_symbol:
        pa=add_symbol_fun(pa)

    # Number of symbols:
    num_symbols=len(pa[0])
    
    # Get a random permutation:
    p=p_not_in_pa(pa,num_symbols)

    keep_searching=True
    
    # Keep track of p_i:
    p_i=[range(num_symbols)]
    if existing_representatives is not None:
        p_i=existing_representatives

    while keep_searching:

        # The first pernutation to be found:
        ckc=qcheck_i(p_i,p,pa)

        if ckc>=d:
            print 'ckc',ckc
            p_i.append(p)
            print p
            if append_file_name is not None:
                append_file(append_file_name,pa2str([p])+'\n')


        # Compute a random draw from the uniform distribution
        draw=random()

        # Reject:
        if draw<=jump_prob:
            p=p_not_in_pa(pa,num_symbols)

        else:
            # print 'trans~'
            for x in range(randint(1,num_symbols/2)):
                p=rand_transpositon(p)
            while p in pa:
                # print 'trans+'
                # p=rand_transpositon(p)
                for x in range(randint(1,num_symbols/2)):
                    p=rand_transpositon(p)

    return p_i

def main(argv):

    parser = argparse.ArgumentParser(description='Crete and ILP model for the partitioning problem')
    parser.add_argument("--in_pa_file", required=False, help='The Permutation Array File')
    parser.add_argument("--sought_hd", type=int, required=False, help='Desired hamming distance for the resulting permutation Array')
    parser.add_argument("--out_pa_file", required=False, help='The filename to save the resulting permutation Array')
    parser.add_argument("--add_symbol", action='store_true',required=False, help='Increase the number of symbols by one')
    parser.add_argument("--out_append_perms", required=False, help='A filename to keep track of newly found permutations')
    parser.add_argument("--jump_prob", required=False, type=float, help='The random probability to change the current permutation and find a new one')
    parser.add_argument("--existing_representatives_list", required=False, help='Start the search including this list of cosets representatives')


    print 'Aqui mero'
    # Get the arguments:
    args = parser.parse_args()

    # The append file name:
    if args.out_append_perms:
        append_file_name=args.out_append_perms
    else:
        append_file_name=None

    # Parse the coset representatives:
    existing_representatives=None
    if args.existing_representatives_list:
        if os.path.isfile(args.existing_representatives_list):
            existing_representatives=parse_pa(args.existing_representatives_list)

    # Get the pa from file
    elif args.in_pa_file:

        # Do partition and extension
        perms=random_walk(args.in_pa_file,d=args.sought_hd,add_symbol=args.add_symbol,append_file_name=append_file_name,jump_prob=args.jump_prob,existing_representatives=existing_representatives)

        print perms

if __name__ == '__main__':
    main(sys.argv[1:])
