# This script calculates the differences between sitewise log likelihoods calculated
# (using RAxML) from a set of ML trees.  The first tree is assumed to be the ML tree
# with no constraints, and the other trees inferred with constraints.  Pair-wise
# comparisons are made between the first tree and the other trees.  Sitewise differences,
# as well as differences summarized by gene are written to two ouptut files.

# specify file names
charset_fname = 'dataset.charsets' # file containing charsets block
site_like_fname = 'RAxML_perSiteLLs.set1'  # output file from RAxML
SL_fname = 'set1.site.comparisons.csv'   # sitewise likelihood differences file
GL_fname = 'set1.locus.comparisons.csv'  # gene likelihood differences file

class charset:
    def __init__(self, label, range):
        self.name = label.strip("'")
        indexes = range.strip(";").split("-")
        self.start = int(indexes[0])
        self.end = int(indexes[1])
        self.diff = []

# read charsets
charset_list = []
charset_file = open(charset_fname, 'r')
print('Reading {}...'.format(charset_fname))
for line in charset_file:
	split_line = line.split()
	if len(split_line) > 0:
		if split_line[0] == 'charset':
			ch = charset(split_line[1], split_line[3])
			charset_list.append(ch)

# read site-likelihood scores
line_num = 1
site_like_scores = []
site_like_file = open(site_like_fname, 'r')
print('Reading {}...'.format(site_like_fname))
for line in site_like_file:
	split_line = line.split()
	if line_num == 1:
		# retrieve number of trees and sites from first line
		numtrees = int(split_line[0])
		numsites = int(split_line[1])
		line_num = line_num + 1
	else:
		# check number of likelihood scores is correct
		if len(split_line) == numsites + 1:
			site_like_scores.append(split_line)
		else:
			raise Exception('Error reading line {} of the site-likelihood file'.format(line_num))

# check number of sets of likelihood scores is correct
if len(site_like_scores) <> numtrees:
	raise Exception('Missing site-likelihood scores')

# open output files
# sitewise likelihood diff file will contain values in rows per comparison
# gene likelihood diff file will contain values in columns per comparison
sl_file = open(SL_fname, 'w')
gl_file = open(GL_fname, 'w')
# write column header on gene likelihood diff file
gl_file.write('locus\t')

# calculate sitewise likelihood differences
for i in range(2, numtrees + 1):
	print('Comparing {} with {}...'.format(site_like_scores[0][0], site_like_scores[i-1][0]))
	# write row label on sitewise likelihood diff file
	sl_file.write('{}-vs-{}\t'.format(site_like_scores[0][0], site_like_scores[i-1][0]))
	# write column header on gene likelihood diff file
	gl_file.write('{}-vs-{}\t'.format(site_like_scores[0][0], site_like_scores[i-1][0]))
	SL_differences = []
	for ch in charset_list:
		locus_diff = 0
		for site in range(ch.start, ch.end + 1):
			site_diff =  float(site_like_scores[0][site]) - float(site_like_scores[i-1][site])
			SL_differences.append(str(site_diff))
			locus_diff = site_diff + locus_diff
		ch.diff.append(str(locus_diff))
	sl_file.write('\t'.join(SL_differences) + '\n')
sl_file.close()
gl_file.write('\n')
for ch in charset_list:
	gl_file.write('{}\t{}\n'.format(ch.name, '\t'.join(ch.diff)))
gl_file.close()
print('Differences in log-likelihoods written to {} and {}.'.format(SL_fname, GL_fname))
