#!/usr/bin/env sage
# 08/12/2016
# Luis Mojica
# UTD
# Get Mathieu Groups from Sage

# from sage.all import designs
import argparse
import os.path
# from utils import write_file
from utils import pa2str

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

def pptrint(mols):
	for idx,mol in enumerate(mols):
		for row in mol:
			print ' '.join([str(x) for x in list(row)])

		if idx<len(mols)-1:
			print '###'

def mols2str(mols):
	out_str=[]
	for idx,mol in enumerate(mols):
		for row in mol:
			# print ' '.join([str(x) for x in list(row)])
			out_str.append(' '.join([str(x) for x in list(row)]))

		if idx<len(mols)-1:
			# print '###'	
			out_str.append('###')
	return out_str

def get_mathieu_pa(num_symbols,out_pa_dir=None):
		
		# # Retrieve the mols
		try:
			G = MathieuGroup(num_symbols)
		except Exception, e:
			print 'Error retrieving the group'
			return None

		X = range(num_symbols)
		act = lambda g, x: X[g(X.index(x) + 1) - 1]
		im = lambda g: [act(g,x) for x in X]

		# Get the array:
		pa=[im(g) for g in G]

		# Make a string:
		pa_str=pa2str(pa)

		# Save to file
		if pa_str and out_pa_dir is not None:
			
			# Create the filename
			out_filename=str(num_symbols)+'_m.txt'

			# Join filename with directory:
			out_filename=os.path.join(out_pa_dir,out_filename)

			print out_filename

			# Write to file
			write_file(out_filename,pa_str)

		else:
			print pa_str

def main():

	# Parse function arguments
	parser = argparse.ArgumentParser(description='Retrieve Mathieu Groups as Permuatation Arrays')
	parser.add_argument("--num_symbols", type=int, required=False, help='Number of symbols for this Mathieu grop')	
	parser.add_argument('--out_pa_dir', help='Directory to store the PA', required=False)
	parser.add_argument('-l', '--l_num_symbols', nargs='+', type=int)
	
	args = parser.parse_args()

	# When a list of mols is asked, and the mols are save to a file:
	if args.num_symbols:

		# Process the list:
		get_mathieu_pa(args.num_symbols,args.out_pa_dir)

if __name__ == '__main__':
	main()