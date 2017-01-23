# Permutation Array Coset Representatives Finder
Source Code for Permutation Array Research in Dr. Sudborough's Lab

From a unix type terminal call:

python src/pa_group_random.py --in_pa_file *\<a_group_file\>* --sought_hd expected_hamming_distance --out_append_perms *\<out_file\>* --jump_prob *\<a_number_between_0_1\>*

i.e.:

python src/pa_group_random.py --in_pa_file ../data/17_15.txt --sought_hd 5 --out_append_perms ../results/17_5_new_perms_RD.txt --jump_prob 0.9

The program will start searching for permutation at the given hamming distance. It will keep running until manually killed using ctrl+c

All Found permutations will be printed to standard output and saved in a file, if a file name is give in the --out_append_perms parameter.


# Find q-1 blocks of size q+1 from PGL(2,q) where q is prime or prime power:

python src/pgl_blocks.py -h for options

i.e:

Give a PGL(2,n) group file to be processed:
python src/pgl_blocks.py --in_pgl_file ../data/4_2.txt --q 3

Use the internal group generator:
python src/pgl_blocks.py --q 3

Use the internal group generator when q is not a prime, but a prime power instead:
python src/pgl_blocks.py --q 2 --pwr 2


