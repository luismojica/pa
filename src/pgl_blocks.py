#!/usr/bin/python2.7 -S
# 01/23/2017
# Luis Mojica
# UTD
# Dedicated Script to generate q-1 blocks of size q+1 from PGL(2,q)

import argparse
from utils import write_file,parse_pa,pa2str
from os import listdir
from os.path import isfile,join
from quick_hd import mul
import math

import os.path,subprocess
from subprocess import STDOUT,PIPE

''' If the class file is not present, compile'''
def compile_java():
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
    # java_class,ext = os.path.splitext(java_file)
    cmd = ['java', 'GroupGeneration','pgl',str(p),str(pwr)]
    proc = subprocess.Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
    stdout,stderr = proc.communicate()
    
    if stdout is "":
        return None
    elif stdout[:3]=='Not':
        print 'aqui'
        return -1
    return parse_table(stdout)

''' Zach's table'''
def parse_table(table):
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


def to_cycles(p):
    cycles=[]
    cycle=[0]
    prev=0
    for idx in xrange(len(p)-1):
        nxt=p[prev]
        if nxt not in cycle:
            cycle.append(nxt)
            prev=nxt
        else:
            # print cycle
            cycles.append(cycle)
            cycle=[p[idx+1]]
            prev=p[idx+1]
    cycles.append(cycle)
    return cycles

''' Check if permutation is a q+1 cycle'''
def is_q_1_cycle(p):

    cycles=to_cycles(p)
    # print cycles
    if len(cycles)==1:
        return True
    return False

'''Compute PGL blocks'''
def get_pgl_blocks(pgl_pa,q,num_blocks):
    for p in pgl_pa:
        if is_q_1_cycle(p):
            qcyc=p
            pgl_pa.remove(qcyc)

    blocks=[[]]
    counter=0

    # Precompute all qcyc powers:
    qccy_powers=[]
    power=qcyc
    for x in xrange(1,q+2):
        power=mul(qcyc,power)
        qccy_powers.append(power)
        # print power in pgl_pa
        if power in pgl_pa:
            pgl_pa.remove(power)

    while pgl_pa:
        p=pgl_pa.pop()
        if p==range(q+1):
            continue
        for x in xrange(q+1):
            comp=mul(qccy_powers[x],p)
            if comp in pgl_pa:
                pgl_pa.remove(comp)
            blocks[counter].append(comp)
        if counter<num_blocks-2:
            blocks.append([])
        else:
            blocks.append(qccy_powers)
            break
        counter+=1
    return blocks

'''Check the expected result and the size of AGL to be loaded'''
def prepare(agl_file,pgl_file,p,q):

    num_cosets=min(p-1,q*(q-1))

    # Get this number of AGL cosets:
    end=(num_cosets)*p # The +1 is for the freeby
    agl_pa=parse_pa(agl_file,end)
    agl_blocks=[]
    for i in xrange(num_cosets):
        agl_blocks.append(agl_pa[i*p:(i*p+p)])

    # Parse all pgl and find blocks
    pgl_pa=parse_pa(pgl_file,end)

    # Get the pgl blocks:
    pgl_blocks=get_pgl_blocks(pgl_pa,q,num_cosets)
    # print len(pgl_blocks),num_cosets
    # exit()

    return agl_blocks, pgl_blocks

def main():

    # Parse function arguments
    parser = argparse.ArgumentParser(description='Generate q-1 blocks of size q+1 from PGL(2,q) where q is a prime (power)')
    parser.add_argument("--in_pgl_file", required=False, help='The Permutation Array File from PGL')
    parser.add_argument("--q", type=int, required=False, help='The prime in PGL(1,q)')
    parser.add_argument("--pwr", type=int, required=False, help='Power to wich q is exponentiated')
    parser.add_argument("--out_blocks_file", required=False, help='Filename to save the resulting PGL blocks')
    parser.add_argument("--no_block_separator", dest='use_block_separator', action='store_false', required=False, help='Do NOT include a line of "###" between blocks')
    parser.set_defaults(use_block_separator=True)
    parser.set_defaults(pwr=1)
    
    args = parser.parse_args()
    pgl_blocks=None
    out_str=[]
    if args.q:
        q=args.q**args.pwr

    # # # Parralle partition and extension:
    if args.in_pgl_file is None and q is not None:
        print 'Using java group generator'
        if not compile_java():
            print 'Error compiling the group generator. Make sure that GroupGeneration.java is located in the same directory'
            exit()
        pgl_pa=get_pgl_java(args.q,args.pwr)
        if pgl_pa is None:
            print 'Error generating the group, preferably provide a file with it'
            exit()
        elif pgl_pa==-1:
            print 'Error: q is not a prime or q is a prime power. If q is a prime power, include the exponent i.e. --pwr 3'
            exit()
        
        # Compute the number of blocks to generate
        num_cosets=q*((q)-1)
        # Get the pgl blocks:
        pgl_blocks=get_pgl_blocks(pgl_pa,q,num_cosets)


    elif q:
        if not isfile(args.in_pgl_file):
            print '"'+args.in_pgl_file+'"', 'not found'
            exit()
        pgl_pa=parse_pa(args.in_pgl_file)
        num_cosets=q*((q)-1)
        
        # Get the pgl blocks:
        pgl_blocks=get_pgl_blocks(pgl_pa,q,num_cosets)

    # Generate output:
    if pgl_blocks is not None:
        for block in pgl_blocks:
            out_str.append(pa2str(block))
            # Separate blocks:
            if args.use_block_separator:
                out_str.append('###')
        if args.use_block_separator:
            out_str.pop()

    # Save the extended PA
    if args.out_blocks_file and pgl_blocks is not None:
        # Save to file
        write_file(args.out_blocks_file,'\n'.join(out_str))
    # print
    elif out_str:
        print '\n'.join(out_str)
    else:
        print 'No File processed or blocks were generated!'


if __name__ == '__main__':
    main()
