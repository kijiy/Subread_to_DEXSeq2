#!/usr/bin/env python

import sys, collections, itertools, os.path, optparse, os, copy, pprint # changed here to add sed

optParser = optparse.OptionParser(

   usage = "python %prog [options] <in.gtf> <out.gff>",

   description=
      "Script to prepare annotation for DEXSeq." +
      "This script takes an annotation file in Ensembl GTF format" +
      "and outputs a 'flattened' annotation file suitable for use " +
      "with the count_in_exons.py script ",

   epilog =
      "Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology " +
      "Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General " +
      "Public License v3. Part of the 'DEXSeq' package. " +
	"Modified by Vivek Bhardwaj (just a bit!) to write featurecounts gtf as an option.")

optParser.add_option( "-r", "--aggregate", type="choice", dest="aggregate",
   choices = [ "yes" ], default = "yes",
   help = "Always aggregate genes (leaving this option argument for the legacy purpose). Indicates whether two or more genes sharing an exon should be merged into an 'aggregate gene'. If 'no', the exons that can not be assiged to a single gene are ignored." )

# add option for featurecounts output
optParser.add_option( "-f", "--featurecountsgtf", type="string", dest="fcgtf", action = "store",
   help = "gtf file to write for featurecounts." )

##

(opts, args) = optParser.parse_args()

if len( args ) != 2:
   sys.stderr.write( "Script to prepare annotation for DEXSeq.\n\n" )
   sys.stderr.write( "Usage: python %s <in.gtf> <out.gff>\n\n" % os.path.basename(sys.argv[0]) )
   sys.stderr.write( "This script takes an annotation file in Ensembl GTF format\n" )
   sys.stderr.write( "and outputs a 'flattened' annotation file suitable for use\n" )
   sys.stderr.write( "with the count_in_exons.py script.\n" )
   sys.exit(1)

try:
   import HTSeq
except ImportError:
   sys.stderr.write( "Could not import HTSeq. Please install the HTSeq Python framework\n" )
   sys.stderr.write( "available from http://www-huber.embl.de/users/anders/HTSeq\n" )
   sys.exit(1)




gtf_file = args[0]
out_file = args[1]

aggregateGenes = opts.aggregate == "yes"

# Step 1: Store all exons with their gene and transcript ID
# in a GenomicArrayOfSets
print("Step 1: Loading GTF file...", flush = True)
exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )

for f in HTSeq.GFF_Reader( gtf_file ):
   if f.type != "exon":
         continue
   # gene_now = f.attr['gene_id'].replace( ":", "_" )
   f.attr['gene_id'] = f.attr['gene_id'].replace( ":", "_" )
   exons[f.iv] += ( f.attr['gene_id'], f.attr['transcript_id'] )



# Step 1.1: Store the last exon information
last_exon_dict = dict()
for cnt, (iv, s) in enumerate(exons.steps()):
   full_set = set()
   for gene_id, transcript_id in s:
      full_set.add( gene_id )

   # Store the last exon for each gene. (Keep updating along with the gtf scanning)
   if '"ENSMUSG00000098997.1"' in full_set:
      print(s)
      print(iv, full_set)

   for gene_now in full_set:
      last_exon_dict[gene_now] = iv


   

# Step 2: Form sets of overlapping genes

# We produce the dict 'gene_sets', whose values are sets of gene IDs. Each set
# contains IDs of genes that overlap, i.e., share bases (on the same strand).
# The keys of 'gene_sets' are the IDs of all genes, and each key refers to
# the set that contains the gene.
# Each gene set forms an 'aggregate gene'.
print("\nStep 2: Storing overlapping gene information...", flush = True)
current_gene_pool = set()
intron_set = set()
is_next_region_intron = dict()
prev_gene_set = set()
prev_region = ""
flag = 0
flag2 = 0


# if aggregateGenes == True: # Regacy: aggregateGenes is always True in this script
gene_sets = dict()

for cnt, (iv, s) in enumerate(exons.steps()):
   # if cnt >= 11805:
   #    flag += 1
   #    print(cnt, iv)
   
   # Set of genes overlapping at the current location (only available for exons)
   full_set = set()
   for gene_id, transcript_id in s:
      full_set.add( gene_id )

   

   
   # Update the current gene pool set. If the current position corresponds to the last exons,
   # the gene should be removed from the set.
   current_gene_pool |= full_set
   
   # Store the current set of genes. This can be used for introns as well because ids are taken from current_pool.
   gene_sets[iv] = copy.copy(current_gene_pool)
   
   # if flag >= 1:
   #    print(full_set, current_gene_pool)
   #    print("")
   
   # if flag == 10:
   #    exit()

   # if '"ENSMUSG00000098997.1"' in full_set:
   #    print(cnt, iv, current_gene_pool, full_set, last_exon_dict['"ENSMUSG00000098997.1"'])

   # Remove gene ids if they are the last exons. By doing so, if current pool gets empty, the next non-exonic region
   # can be treated as an intergenic region.
   for gene_now in full_set:
      if iv == last_exon_dict[gene_now]:
         current_gene_pool = current_gene_pool - {gene_now}

   # Store the region to an intron pool if the current region is not exon and the size of the
   # current gene pool is more than 0 (=we are still inside at least one gene body).
   if len(s) == 0 and len(current_gene_pool) > 0:
      intron_set.add(iv)
      


   # if '"ENSMUSG00000098997.1"' in full_set:
   #    print(cnt, iv, current_gene_pool, full_set, last_exon_dict['"ENSMUSG00000098997.1"'])


   # if "ENSMUSG00000048747.14" in full_set:
   #    flag2 = 1
   
   # if flag2 == 1:
   #    print(cnt, iv, full_set, prev_gene_set, len(full_set & prev_gene_set) > 0)

   # # if "ENSMUSG00000085710.1" in full_set:
   # #    print(cnt, iv, full_set, prev_gene_set, len(full_set & prev_gene_set) > 0)

   # if cnt == 172561:
   #    exit()




   # # Judge if the previous exon was followed by intron.
   # # This can be told by checking if the current and previous
   # # gene ids have overlaps.
   # if len(s) > 0 and flag != 0: # Exon
   #    if len(full_set & prev_gene_set) > 0:
   #       is_next_region_intron[prev_region] = True
   #    else:
   #       is_next_region_intron[prev_region] = False
      
   #    # Store the prev info
   #    prev_gene_set = copy.copy(full_set)
   #    prev_region = copy.copy(iv)

   
      
   # flag = 1

   
   # if cnt == 100:
   #    break

# The final exon in the entire GTF shouldn't have introns after it
# is_next_region_intron[prev_region] = False
# print(is_next_region_intron)

# pprint.pprint(gene_sets)
# print(is_next_region_intron)



# Step 3: Go through the steps again to get the exonic sections. Each step
# becomes an 'exonic part'. The exonic part is associated with an
# aggregate gene, i.e., a gene set as determined in the previous step,
# and a transcript set, containing all transcripts that occur in the step.
# The results are stored in the dict 'aggregates', which contains, for each
# aggregate ID, a list of all its exonic_part features.

aggregates = collections.defaultdict( lambda: list() )
gene_id_prev = ""
is_intron = False
for cnt, (iv, s) in enumerate(exons.steps()):
   # Deal with non-exonic region
   if len(s) == 0: # non-exonic
      if iv in intron_set:
         is_intron = True
      else:
         is_intron = False
         continue
   else: # exonic
      is_intron = False
      gene_id = list(s)[0][0]

   ## if aggregateGenes=FALSE, ignore the exons associated to more than one gene ID
   ## This option does not exist in the current implementation

   # if aggregateGenes == False:
   #    check_set = set()
   #    for geneID, transcript_id in s:
   #       check_set.add( geneID )
   #    if( len( check_set ) > 1 ):
   #       continue
   #    else:
   #       aggregate_id = gene_id

   # # Take one of the gene IDs, find the others via gene sets, and
   # # form the aggregate ID from all of them
   # else:
   # assert set( gene_id for gene_id, transcript_id in s ) <= gene_sets[ gene_id ]
   aggregate_id = '+'.join( sorted(list(gene_sets[ iv ])) )
   
   # Make the feature and store it in 'aggregates'
   f = HTSeq.GenomicFeature( aggregate_id, "exonic_part", iv )
   f.source = os.path.basename( sys.argv[0] )
   f.attr = {}
   f.attr[ 'gene_id' ] = aggregate_id
   transcript_set = set( ( transcript_id for gene_id, transcript_id in s ) )
   f.attr[ 'transcripts' ] = '+'.join( transcript_set )
   f.attr[ 'region_type' ] = "Intron" if is_intron else "Exon"
   aggregates[ aggregate_id ].append( f )

   # if gene_id == "ENSMUSG00000102307.1":
   #    print(cnt, iv, s, aggregate_id, '+'.join( transcript_set ), f.attr[ 'region_type' ])



# Step 4: For each aggregate, number the exonic parts

aggregate_features = []
for l in aggregates.values():
   for i in range( len(l)-1 ):
      assert l[i].name == l[i+1].name, str(l[i+1]) + " has wrong name"
      assert l[i].iv.end <= l[i+1].iv.start, str(l[i+1]) + " starts too early"
      if l[i].iv.chrom != l[i+1].iv.chrom:
         raise ValueError(
            "Same name found on two chromosomes: %s, %s" % ( str(l[i]), str(l[i+1]) )
         )
      if l[i].iv.strand != l[i+1].iv.strand:
         raise ValueError(
            "Same name found on two strands: %s, %s" % ( str(l[i]), str(l[i+1]) )
         )
   aggr_feat = HTSeq.GenomicFeature( l[0].name, "aggregate_gene",
      HTSeq.GenomicInterval( l[0].iv.chrom, l[0].iv.start,
         l[-1].iv.end, l[0].iv.strand ) )
   aggr_feat.source = os.path.basename( sys.argv[0] )
   aggr_feat.attr = { 'gene_id': aggr_feat.name }
   for i in range( len(l) ):
      l[i].attr['exonic_part_number'] = "%03d" % ( i+1 )
   aggregate_features.append( aggr_feat )


# Step 5: Sort the aggregates, then write everything out

aggregate_features.sort( key = lambda f: ( f.iv.chrom, f.iv.start ) )

fout = open( out_file, "w" )
for aggr_feat in aggregate_features:
   fout.write( aggr_feat.get_gff_line() )
   for f in aggregates[ aggr_feat.name ]:
      fout.write( f.get_gff_line() )

fout.close()

## modify file to print gtf if featurecounts gtf requested
fcountgtf = opts.fcgtf

if fcountgtf :
	os.system('sed s/aggregate_gene/gene/g ' + out_file + ' > ' + fcountgtf)
	os.system('sed -i s/exonic_part/exon/g ' + fcountgtf)
	print("Done!")
else :
	print("Done!")
