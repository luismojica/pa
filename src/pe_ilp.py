#!/usr/bin/python
# 08/04/2016
# Luis Mojica
# UTD
# Test the partitioning problem as an ILP

import gurobipy as gu
import argparse
import sys, getopt
from utils import parse_pa,pa2str,write_file,parse_repsentatives,create_cosets,extend
import math
import re

# Create an ILP problem for an instance:
def make_ilp(blocks,Q,P,timelimit=-1):

    # Model
    m = gu.Model("partitioning")
    # print m
    # The position variables:
    num_symbols=len(blocks[0][0])

    # Number of blocks:
    num_blocks=len(blocks)

    # P solution:
    P_res=P[:]

    # Make the position binary variables:
    for idx_block in xrange(len(blocks)):
        for position in xrange(num_symbols):
            P[idx_block].append(m.addVar(vtype=gu.GRB.BINARY, name='b'+str(idx_block)+'_p'+str(position)))
            P_res[idx_block]=[]

    m.update()

    # Create the objective function:
    obj = gu.LinExpr()

    # Create constrain -- Sum_P_ix >= C_ij:
    obj = gu.LinExpr()

    borrar=0

    C=[[] for x in xrange(num_blocks)]
    # Make the coverage binary variables:
    for idx_block,block in enumerate(blocks):
        q=Q[idx_block]
        for idx_perm,perm in enumerate(block):

            # Permutation Variables:
            C[idx_block].append(m.addVar(vtype=gu.GRB.BINARY, name='C_b'+str(idx_block)+'_j'+str(idx_perm)))
    m.update()


    # Test:
    # Make the coverage binary variables:
    for idx_block,block in enumerate(blocks):
        q=Q[idx_block]
        for idx_perm,perm in enumerate(block):
            tmp_P=[]
            for sym_idx,per_sym in enumerate(perm):
                if per_sym in q:
                    tmp_P.append(P[idx_block][sym_idx])
            m.addConstr(gu.quicksum(tmp_P) >= C[idx_block][idx_perm],'Sum_P_ix >= C_'+str(idx_block)+str(idx_perm))
        # Add the permutation variables to the objective function:
        obj.addTerms([1]*len(C[idx_block]),C[idx_block])

    m.update()

    # The partition constraints:
    all_P_vars=[]
    # 2. Parts are mutually exclusive:
    for position in range(num_symbols):
        # # Temporary holder for the constraint:
        tmp_P=[]
        for idx_block in range(num_blocks):
            tmp_P.append(P[idx_block][position])
            all_P_vars.append(position*P[idx_block][position])

        # Add the position/permutation constraints:
        m.addConstr(gu.quicksum(tmp_P) == 1,'Sum_Pi_'+str(position)+'==1')


    # 3. The sum of the parts is the original set:
    m.addConstr(gu.quicksum(all_P_vars) == sum(xrange(0,num_symbols)),'Sum_Pi_x==Sum_Zn')

    m.update()
    # print m

    m.setObjective(obj, gu.GRB.MAXIMIZE)
    # m.setParam( 'OutputFlag', False )

    # Set time limit:
    if timelimit>0:
        m.params.timelimit = timelimit
        m.setParam( 'OutputFlag', False)
        m.update()
        m.optimize()
        return None,m.ObjVal
    else:
        f_time=num_symbols/float(10)
        m.params.timelimit = 1000+f_time
        # m.params.timelimit = 10+f_time

    m.optimize()

    if m.status==2 or m.status==9:
        P_res=parse_ilp_sol(m,P_res)

    else:
        print 'Something went wrong with the ILP'
        print 'm.status',m.status
        return None,None

    return P_res,m.ObjVal

# Parse and infer the results from the model solution:
def parse_ilp_sol(m,P):
    # print 'P:',P
    # Create the output:
    for v in m.getVars():
        if v.x>0:
            # print v.varName
            prefix=v.varName[0]
            if prefix=='b':
                block_idx,pos=[int(s) for s in re.findall(r'\d+', v.varName)]
                P[block_idx].append(pos)
    return P


# def prepare(pa,num_blocks):
def prepare(pa_file,num_symbols=0,num_blocks=0,ceil=False):

    # Num of partitions to make:
    if ceil:
        sqr_sym=int(math.ceil(math.sqrt(num_symbols)))
    else:
        sqr_sym=int(math.floor(math.sqrt(num_symbols)))

    # Number of parts:
    num_parts=min(sqr_sym,num_blocks-1)

    # Number of lines from the PA (asumming in order):
    end=(num_parts+1)*num_symbols # The +1 is for the freeby
    pa=parse_pa(pa_file,end)

    elem_per_block=len(pa)/(num_parts+1)
    print elem_per_block
    
    # AGL blocks:
    if len(pa)%float(num_symbols)==0:
        print 'AGL blocks'
        elem_per_block=num_symbols
        print 'elem_per_block',elem_per_block
    # exit()

    # Find blocks
    if num_blocks is not None:
        blocks=[]
        for b_idx in range((num_parts+1)):
            sough_dist_block=pa[elem_per_block*b_idx:elem_per_block*b_idx+elem_per_block]
            # print len(sough_dist_block)
            # print 'b_idx',b_idx,hd_pairwise(sough_dist_block)
            blocks.append(sough_dist_block)

    # New symbol:
    new_symbol=len(blocks[0][0])

    # The freeby block:
    freeby=blocks[-1]

    # Get only the blocks that will actually be used
    blocks=blocks[0:num_parts]

    # For when num of symbols is not integer:
    q_start=0
    q_add=0

    # Num of permutations covered with the arbitrary selection:
    start_coverage=0

    # Create the partition shell
    P=[[] for x in range(num_parts)]
    Q=[[] for x in range(num_parts)]

    # Num elements in each default partition:
    def_elemnts=int(math.floor(new_symbol/float(num_parts)))
    # Check if the number of elements is an integer or not
    is_int=(new_symbol/float(num_parts)).is_integer()

    # Add the default part:
    for x in range(num_parts):
        # P[x].append(x)
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

    print P,Q

    return blocks,freeby,Q,P,new_symbol

# Run the partition and extension algorithm:
def partition_extension(pa_file,sought_hd,max_iter=5,num_blocks=None,ceil=False):

    # Retrive the blocks and freeby:
    blocks,freeby,Q,P,new_symbol=prepare(pa_file,sought_hd,num_blocks,ceil=ceil)

    # print len(pa_file)
    # print len(blocks)
    # exit()
    
    # Solve with ILP
    P,coverage=make_ilp(blocks,Q,P)

    # Something went wrong with ILP
    if P is None:
        return None,None,None

    # Extend the Permutation Array:
    pa_ext=extend(blocks,P,Q,new_symbol)

    # Extend the freeby:
    pa_ext+=extend([freeby],[[new_symbol]],[xrange(new_symbol)],new_symbol)

    # Create the p partition and q partition as strings:
    p_q_str='P:'+str(P)+'\n'
    p_q_str+='Q:'+str(Q)+'\n'

    # Convert to string:
    pa_str=pa2str(pa_ext)

    return pa_str,p_q_str,len(pa_ext)

''' Partition and Extension if blocks are given 11/07/2016'''
# def partition_extension(pa_file,sought_hd,max_iter=5,num_blocks=None,ceil=False):
def pe_blocks(blocks,num_blocks):

    # Number of parts:
    num_parts=num_blocks

    # New symbol: !!! Note the difference, no longe assuming square blocks! this may not be AGL
    new_symbol=len(blocks[0][0])

    # For when num of symbols is not integer:
    q_start=0
    q_add=0

    # Num of permutations covered with the arbitrary selection:
    start_coverage=0

    # Create the partition shell
    P=[[] for x in range(num_parts)]
    Q=[[] for x in range(num_parts)]

    # Num elements in each default partition:
    def_elemnts=int(math.floor(new_symbol/float(num_parts)))
    # Check if the number of elements is an integer or not
    is_int=(new_symbol/float(num_parts)).is_integer()

    # Add the default part:
    for x in range(num_parts):
        # P[x].append(x)
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

    # Make sure tha all symbols are included in the Q partitioning
    if q_add<new_symbol:
        # get the missing symbols:
        Q[0]+=range(q_add,new_symbol)
   
    # Solve with ILP
    P,coverage=make_ilp(blocks,Q,P)

    return P,Q

# Run the partition and extension algorithm on blocks:
def partition_extension_blocks(blocks,num_blocks=None,modify_pa=True):

    if num_blocks is None:
        num_parts=num_blocks
    else:
        num_parts=len(blocks)-1

    P,Q=pe_blocks(blocks[:-1],len(blocks)-1)
    freeby=blocks[-1]

    # New symbol: !!! Note the difference, no longe assuming square blocks! this may not be AGL
    new_symbol=len(blocks[0][0])

    # Extend the Permutation Array:
    if modify_pa:
        pa_ext=extend(blocks[:-1],P,Q,new_symbol)

        # Extend the freeby:
        pa_ext+=extend([freeby],[[new_symbol]],[xrange(new_symbol)],new_symbol)
    else:
        pa_ext=extend(blocks[:-1],P,Q)

        # Extend the freeby:
        pa_ext+=extend([freeby],[[new_symbol]],[xrange(new_symbol)])

    # Create the p partition and q partition as strings:
    p_q_str='P:'+str(P)+'\n'
    p_q_str+='Q:'+str(Q)+'\n'

    # Convert to string:
    pa_str=pa2str(pa_ext)

    return pa_str,p_q_str,len(pa_ext)


def main(argv):

    parser = argparse.ArgumentParser(description='Crete and ILP model for the partitioning problem')
    parser.add_argument("--in_pa_file", required=False, help='The Permutation Array File')
    parser.add_argument("--sought_hd", type=int, required=False, help='The Permutation Array File')
    parser.add_argument("--num_blocks", type=int, required=False, help='The number of blocks to look for')
    parser.add_argument("--out_pa_file", required=False, help='The filename to save the resulting permutation Array')
    parser.add_argument("--out_partition_file", required=False, help='The filename to save the patition for permutation Array ')
    parser.add_argument("--time_limit", required=False, help='Stop solver after these many seconds')
    parser.add_argument("--use_ceiling_function", action='store_true', required=False, help='When computing the number of blocks, use ceiling instead of floor')
    parser.add_argument("--in_block_file", required=False, help='The Permutation Array File')
    parser.add_argument("--representatives", required=False, help='A file containing a list of coset representatives of the group to be extended')


    # Get the arguments:
    args = parser.parse_args()

    # Get the pa from file
    if args.in_pa_file:
        # pa=parse_pa(args.in_pa_file)
        # print pa

        # Do partition and extension
        pa_str,p_q_str,pa_len=partition_extension(args.in_pa_file,sought_hd=args.sought_hd,num_blocks=args.num_blocks,ceil=args.use_ceiling_function)

        print p_q_str
        # exit()

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

        # Get cosets:
        cosets=create_cosets(pa,reps,args.num_blocks)

        pa_str,p_q_str,pa_len=partition_extension_blocks(cosets)

        print pa_len

if __name__ == '__main__':
    main(sys.argv[1:])
