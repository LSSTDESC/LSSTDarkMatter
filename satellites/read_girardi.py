import re
import atpy
import collections

# Code to read Parsec isochrones 

def get_columns(filename):
	prevline = None
	with open(filename) as f:
		for i, line in enumerate(f):
			if line[0] != '#':
				break
			prevline = line
			
	columns = prevline
	if columns is None:
		raise IOError('Failed to read file')
	columns = re.sub(r'#[\t ]*', '', columns)  # strip the "#	"
	columns = re.sub(r'[\t ]+[\t ]*', r'\t', columns)
	columns = columns.rstrip('\n')
	column_list = columns.split('\t')
	ret = []
	for c in column_list:
		if len(c)!=0:
			ret.append(c)

	return ret, i


def read_girardi(filename):
	columns, header_height = get_columns(filename)
	tab = atpy.Table(filename, type='ascii')

	hash = collections.OrderedDict()
	for i, cur_col in enumerate(columns):
		hash[cur_col] = tab['col%d'% (i + 1)]
	return hash
