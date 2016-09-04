# Permutation Array Coset Representatives Finder
Source Code for Permutation Array Research in Dr. Sudborough's Lab

From a unix type terminal call:

python src/pa_group_random.py --in_pa_file a_group_file --sought_hd expected_hamming_distance --out_append_perms out_file --jump_prob a_number_between_0_1

i.e.:

python src/pa_group_random.py --in_pa_file ../data/17_15.txt --sought_hd 5 --out_append_perms ../results/17_5_new_perms_RD.txt --jump_prob 0.9
