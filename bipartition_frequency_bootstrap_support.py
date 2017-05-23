import dendropy
import os
import sys
import multiprocessing
from dendropy.utility import bitprocessing
from dendropy.datamodel import treemodel

def frequency_of_bipartition_inclusive(tree_list, taxon_labels, return_locus_list):
        is_bipartitions_updated = False
        split = tree_list.taxon_namespace.taxa_bitmask(labels=taxon_labels)
        k = len(taxon_labels)
        if bitprocessing.num_set_bits(split) != k:
            raise IndexError('Not all taxa could be mapped to bipartition (%s): %s' \
                % (tree_list.taxon_namespace.bitmask_as_bitstring(split), k))
        found = 0
        total = 0
        if return_locus_list:
            locus_list = []
        for tree in tree_list:
            tree_labels = [leaf.taxon.label for leaf in tree.leaf_nodes()]
            labels_present = [name for name in taxon_labels if name in tree_labels]
            modified_split = tree_list.taxon_namespace.taxa_bitmask(labels=labels_present)
            unnormalized_split = modified_split
            normalized_split = treemodel.Bipartition.normalize_bitmask(
                 bitmask=modified_split,
                 fill_bitmask=tree_list.taxon_namespace.all_taxa_bitmask(),
                 lowest_relevant_bit=1)
            if not is_bipartitions_updated or not tree.bipartition_encoding:
                tree.encode_bipartitions()
            bipartition_encoding = set(b.split_bitmask for b in tree.bipartition_encoding)
            total += 1
            if tree.is_unrooted and (normalized_split in bipartition_encoding):
                found += 1
                if return_locus_list:
                    locus_list.append(tree._label)
            elif (not tree.is_unrooted) and (unnormalized_split in bipartition_encoding):
                found += 1
                if return_locus_list:
                    locus_list.append(tree._label)
        if return_locus_list:
            try:
                return float(found)/total, locus_list
            except ZeroDivisionError:
                return 0, 0
        else:
            try:
                return float(found)/total
            except ZeroDivisionError:
                return 0

def frequency_of_bipartition_multiple(tree_list, taxon_sets):
        is_bipartitions_updated = False
        for taxon_labels in taxon_sets:
            split = tree_list.taxon_namespace.taxa_bitmask(labels=taxon_labels)
            k = len(taxon_labels)
            if bitprocessing.num_set_bits(split) != k:
                raise IndexError('Not all taxa could be mapped to bipartition (%s): %s' \
                    % (tree_list.taxon_namespace.bitmask_as_bitstring(split), k))
        found_all = 0
        total = 0
        num_sets = len(taxon_sets)
        for tree in tree_list:
            tree_labels = [leaf.taxon.label for leaf in tree.leaf_nodes()]
            if not is_bipartitions_updated or not tree.bipartition_encoding:
                tree.encode_bipartitions()
            bipartition_encoding = set(b.split_bitmask for b in tree.bipartition_encoding)
            found = 0
            for taxon_labels in taxon_sets:
                labels_present = [name for name in taxon_labels if name in tree_labels]
                modified_split = tree_list.taxon_namespace.taxa_bitmask(labels=labels_present)
                unnormalized_split = modified_split
                normalized_split = treemodel.Bipartition.normalize_bitmask(
                     bitmask=modified_split,
                     fill_bitmask=tree_list.taxon_namespace.all_taxa_bitmask(),
                     lowest_relevant_bit=1)
                if tree.is_unrooted and (normalized_split in bipartition_encoding):
                    found += 1
                elif (not tree.is_unrooted) and (unnormalized_split in bipartition_encoding):
                    found += 1
            if found == num_sets:
                found_all += 1
            total += 1
        try:
            return float(found_all)/total
        except ZeroDivisionError:
            return 0

class bipartition_info:
    def __init__(self, n, taxon_labels):
        self.index = n
        self.taxon_labels = taxon_labels
        self.gt_prop = 0.0
        self.ave_support = 0.0
        self.support_list = []
        self.loci_list = []

def worker(work):
    locus, bipartition_indexes = work
    support_values = []
    bs = dendropy.TreeList.get_from_path(os.path.join(bootstrap_dir,locus),'newick', taxon_namespace=tns,rooting="force-unrooted")
    for i in bipartition_indexes:
        support = frequency_of_bipartition_inclusive(bs, list[i], return_locus_list=False)
        support_values.append([i, support])
    sys.stdout.write(".")
    sys.stdout.flush()
    return support_values

cores = 11   # number of cores
bootstrap_dir = 'euth70_spades_genetree_bootreps' # directory containing 1 file for each locus with filename the same as the locus name
# with no extension.  each file should contain bootstrap replicate trees in newick format.

# make sure the same taxon name space with all taxa is used
# you can use an ml tree with all taxa in newick format
examl = dendropy.Tree.get_from_path('examl','newick')
tns = examl.taxon_namespace

print("Reading gene trees ...")
# read gene trees with appropriate taxon name space and unrooted
trees=dendropy.TreeList.get_from_path('euth70.spades.best.nexus.tre','nexus', taxon_namespace=tns, rooting="force-unrooted")

# set clade names

Xenarthra = ['bradypus variegatus 55206', 'choloepus hoffmanni', 'dasypus novemcinctus', 'tamandua tetradactyla 55944']

Afrotheria = ['chrysochloris asiatica', 'echinops telfairi', 'elephantulus edwardii', 'loxodonta africana', 'microgale drouhardi 217026', 'orycteropus afer', 'procavia capensis', 'trichechus manatus']

Scandentia = ['ptilocercus lowii 153633', 'tupaia belangeri', 'tupaia chinensis']

Primatomorpha = ['cynocephalus volans 165788', 'galeopterus variegatus', 'aotus nancymaae', 'ateles sp 165380', 'callithrix jacchus', 'cebus capucinus', 'daubentonia madagascariensis', 'eulemur macaco', 'gorilla gorilla', 'homo sapiens', 'loris tardigradus 180671', 'macaca fascicularis', 'microcebus murinus', 'nomascus leucogenys', 'otolemur garnettii', 'pan paniscus', 'propithecus coquereli', 'tarsius syrichta']

Dermoptera = ['cynocephalus volans 165788', 'galeopterus variegatus']

Glires = ['cavia porcellus', 'chinchilla lanigera', 'cricetulus griseus', 'dipodomys ordii', 'fukomys damarensis', 'heterocephalus glaber', 'jaculus jaculus', 'marmota marmota', 'mus musculus', 'nannospalax galili', 'octodon degus', 'rattus norvegicus', 'ochotona princeps', 'oryctolagus cuniculus']

Chiroptera = ['megaderma lyra', 'miniopterus natalensis', 'pteronotus parnellii', 'pteropus alecto', 'rhinolophus ferrumequinum', 'taphozous theobaldi 192565']

Ferae = ['manis javanica', 'manis pentadactyla', 'ailurus fulgens 50275', 'bdeogale nigripes 167685', 'callorhinus ursinus 111859', 'canis lupus', 'cryptoprocta ferox 7539', 'felis catus', 'leptonychotes weddellii', 'mustela putorius', 'nandinia binotata 193214', 'nasua nasua 34293', 'odobenus rosmarus', 'proteles cristata 196086', 'ursus maritimus', 'viverra tangalunga 146957']

Perissodactyla = ['ceratotherium simum', 'equus caballus', 'tapirus terrestris 41204']

Cetartiodactyla = ['antilocapra americana 156758', 'balaenoptera bonaerensis', 'bos taurus', 'camelus bactrianus', 'capreolus capreolus', 'delphinapterus leucas 168091', 'giraffa camelopardalis', 'lipotes vexillifer', 'mesoplodon stejnegeri 118820', 'orcinus orca', 'physeter catodon', 'sus scrofa', 'tayassu tajacu 142881', 'vicugna pacos']

Eulipotyphla = ['condylura cristata', 'crocidura tansaniana 149973', 'erinaceus europaeus', 'hylomys suillus 168652', 'scapanus townsendii 31414', 'sorex araneus']

Boreoeutheria = Scandentia + Primatomorpha + Glires + Chiroptera + Ferae + Perissodactyla + Cetartiodactyla + Eulipotyphla

Primates = ['aotus nancymaae', 'ateles sp 165380', 'callithrix jacchus', 'cebus capucinus', 'daubentonia madagascariensis', 'eulemur macaco', 'gorilla gorilla', 'homo sapiens', 'loris tardigradus 180671', 'macaca fascicularis', 'microcebus murinus', 'nomascus leucogenys', 'otolemur garnettii', 'pan paniscus', 'propithecus coquereli', 'tarsius syrichta']

Carnivora = ['ailurus fulgens 50275', 'bdeogale nigripes 167685', 'callorhinus ursinus 111859', 'canis lupus', 'cryptoprocta ferox 7539', 'felis catus', 'leptonychotes weddellii', 'mustela putorius', 'nandinia binotata 193214', 'nasua nasua 34293', 'odobenus rosmarus', 'proteles cristata 196086', 'ursus maritimus', 'viverra tangalunga 146957']

Atlantogenata = Xenarthra + Afrotheria

Epitheria = Afrotheria + Boreoeutheria

Exafroplacentalia = Xenarthra + Boreoeutheria

TS1 = Scandentia + Primates

TS2 = Scandentia + Glires

TS3 = Scandentia + Primatomorpha

TS4 = Scandentia + Glires + Primatomorpha

LT1 = Perissodactyla + Cetartiodactyla

LT2 = Perissodactyla + Chiroptera

LT3 = Perissodactyla + Carnivora

LT4 = Perissodactyla + Ferae

LT5 = Perissodactyla + Cetartiodactyla + Ferae

LT6 = Perissodactyla + Cetartiodactyla + Chiroptera

LT7 = Perissodactyla + Chiroptera + Ferae

LT8 = Perissodactyla + Cetartiodactyla + Chiroptera + Ferae


# list bipartitions

list = [Scandentia, Primates, Dermoptera, Glires, Primatomorpha, Perissodactyla, Cetartiodactyla, Ferae, Chiroptera, Carnivora]
print("Scandentia, Primates, Dermoptera, Glires, Primatomorpha, Perissodactyla, Cetartiodactyla, Ferae, Chiroptera, Carnivora")
# list = [Xenarthra, Afrotheria, Boreoeutheria]
# print("Xenarthra, Afrotheria, Boreoeutheria")
# list = [Atlantogenata, Epitheria, Exafroplacentalia, TS1, TS2, TS3, TS4, LT1, LT2, LT3, LT4, LT5, LT6, LT7, LT8]
# print("Atlantogenata, Epitheria, Exafroplacentalia, TS1, TS2, TS3, TS4, LT1, LT2, LT3, LT4, LT5, LT6, LT7, LT8")

# list multiple constraints
# constraints are separated by commas, 
multiple_constraints = [[Afrotheria, Xenarthra, Afrotheria + Xenarthra],
    [Afrotheria, Boreoeutheria, Afrotheria + Boreoeutheria],
    [Xenarthra, Boreoeutheria, Xenarthra + Boreoeutheria],	
    [Scandentia, Primates, Scandentia + Primates],
    [Scandentia, Primatomorpha, Scandentia + Primatomorpha],
    [Scandentia, Glires, Scandentia + Glires],
    [Scandentia, Primatomorpha + Glires, Scandentia + Primatomorpha + Glires, Primatomorpha, Glires],
    [Perissodactyla, Cetartiodactyla, Perissodactyla + Cetartiodactyla],
    [PerissodactyChiropterala + Ferae, Perissodactyla, Ferae],
    [Perissodactyla + Chiroptera, Perissodactyla, Chiroptera],
    [Perissodactyla + Cetartiodactyla + Ferae, Perissodactyla, Cetartiodactyla + Ferae, Cetartiodactyla, Ferae],
    [Perissodactyla + Ferae + Chiroptera, Perissodactyla, Ferae + Chiroptera, Ferae, Chiroptera],
    [Perissodactyla + Cetartiodactyla + Chiroptera, Perissodactyla, Cetartiodactyla + Chiroptera, Cetartiodactyla, Chiroptera],
    [Perissodactyla, Cetartiodactyla + Chiroptera + Ferae],
    [Perissodactyla, Cetartiodactyla, Chiroptera, Ferae, Chiroptera + Ferae, Cetartiodactyla + Chiroptera + Ferae]]

# Part 1.  Calculate frequency of gene trees that support a bipartition, and
# among those gene trees calculate the average bootstrap support for that bipartition.

bipartitions_list = [] # initialize

print("Generating list of loci for each bipartition ")
# generate list of loci for each bipartition
for n, bipartition in enumerate(list):
    b = bipartition_info(n, bipartition)
    b.gt_prop, b.loci_list = frequency_of_bipartition_inclusive(trees, b.taxon_labels, return_locus_list=True)
    bipartitions_list.append(b)
    sys.stdout.write(".")
    sys.stdout.flush()
print("")

# sort list of bipartitions to test for each locus into a dictionary
print("Sorting by locus ...")
work_dict = {}
for b in bipartitions_list:
    for locus in b.loci_list:
        if locus in work_dict.keys():
            work_dict[locus].append(b.index)
        else:
            work_dict[locus] = [b.index]

# multiprocess tasks per locus
work = work_dict.items()

print("Processing {} loci ...".format(len(work)))
pool = multiprocessing.Pool(cores)
bs_support = pool.map(worker, work)
print("")

# summarize average support
print("Summarizing average support ...")
for sup in bs_support:
    for s in sup:
        bipartitions_list[s[0]].support_list.append(s[1])

print("Results:")
for b in bipartitions_list:
    b.ave_support = sum(b.support_list)/float(len(b.support_list))
    print(b.index)
    print(b.taxon_labels)
    print("gene tree proportion: ", b.gt_prop)
    print("average bootstrap support: ", b.ave_support)

# Part 2: Calculate proportion of gene trees that support multiple constraints
for const in multiple_constraints:
    freq = frequency_of_bipartition_multiple(trees, const)
    print("\nFrequency of gene trees with multiple constraints")
    print(const, "frequency: {}".format(freq))
    
