# Required packages:

*\<numpy\>*:
Follow instructions from https://www.scipy.org/install.html

or if using Ubuntu:

sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose

if using OSX:

export PATH="$PATH:/Users/your_user/Library/Python/3.5/bin
python -m pip install --upgrade pip
pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose

*\<Gurobi License\>*:

For academic users, get a free license from gurobi.com.

Specifically, for a single machine, unlimited use, follow the instructions given in:

http://www.gurobi.com/academia/for-universities

# Permutation Array Coset Representatives Finder
Source Code for Permutation Array Research in Dr. Sudborough's Lab

From a unix type terminal call:

python src/pa_group_random.py --in_pa_file *\<a_group_file\>* --sought_hd expected_hamming_distance --out_append_perms *\<out_file\>* --jump_prob *\<a_number_between_0_1\>*

i.e.:

python src/pa_group_random.py --in_pa_file ../data/17_15.txt --sought_hd 5 --out_append_perms ../results/17_5_new_perms_RD.txt --jump_prob 0.9

The program will start searching for permutation at the given hamming distance. It will keep running until manually killed using ctrl+c

All Found permutations will be printed to standard output and saved in a file, if a file name is give in the --out_append_perms parameter.


# Find q*(q-1) blocks of size q+1 from PGL(2,q) where q is prime or prime power:

python src/pgl_blocks.py -h for options

i.e:

Give a PGL(2,n) group file to be processed:

python src/pgl_blocks.py --in_pgl_file ../data/4_2.txt --q 3

Use the internal group generator:

python src/pgl_blocks.py --q 3

Use the internal group generator when q is not a prime, but a prime power instead:

python src/pgl_blocks.py --q 2 --pwr 2

# Parallel Partition and Extension on PGL(2,q) where q is a prime or a prime power:

From a unix type terminal call:

python src/pe_parallel_ilp.py --in_pgl_file *\<a_group_file\>* --q *\<a prime or prime power\>* --out_pa_file *\<out_file_name\>* --jump_prob *\<a_number_between_0_1\>* [--use_ILP] --num_blocks *\<number of blocks to partition\>*

Note, in addition to the number of blocks to partition, two additional blocks will be added to the final permutation arrays (freebies).

Use the parameter --use_ILP to compute the partition using an Integer Linear Program. Otherwise a greedy approach is used.

i.e.:

python src/pe_parallel_ilp.py --in_pgl_file ../../74_72.txt --q 73 --num_blocks 5
