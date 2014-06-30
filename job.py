#!/usr/bin/env python
'''\
create a CSV file with data from 3 separate files
'''
import sys, os, csv, traceback, time, re, urllib, argparse
sys.path.append(os.path.abspath('.'))  # for pytest
from Bio.Blast import NCBIStandalone, Record
from common import *
import tsv_blast_parser, blast2xml, parsenames
from collections import defaultdict, OrderedDict
PARSER = tsv_blast_parser
DATESTAMP = os.getenv('DATESTAMP') or time.strftime('%Y%m%d%H%M%S')
DATA = {}
CONTAMINATION_URLS = [
 'http://www.prevalentfungi.org/fungi.cfm',
 'http://www.prevalentbacteria.org/bacteria.cfm',
]
REPORT1 = 'AnnotationReport-%s.tsv' % DATESTAMP
REPORT2 = 'SNPReport-%s.tsv' % DATESTAMP
REPORT3 = 'NoSNPReport-%s.tsv' % DATESTAMP
LOGFILE = os.getenv('LOGFILE') or 'log-%s.txt' % DATESTAMP
NOHITS = os.getenv('NOHITS') or 'nohits-%s.txt' % DATESTAMP
CONTAMINATED = os.getenv('CONTAMINATED') or 'contaminated-%s.txt' % DATESTAMP
CONTAMINANTS = set()
SCRIPT_NAME = sys.argv[0]
SCRIPT_PATH = os.path.abspath(os.path.dirname(SCRIPT_NAME))
CONTAMINANTS_FILE = os.path.join(SCRIPT_PATH, 'contaminants.txt')
INFORMATIVE_DOESNT_MATTER = os.getenv('INFORMATIVE_DOESNT_MATTER') or False
VCFHEADER = []
NONZERO = .000001  # to avoid divide-by-zero errors
CUTOFF = .01  # max evalue for determination of contaminated species
INDEX = dict([[k, tsv_blast_parser.FIELDS.index(k)]
 for k in tsv_blast_parser.FIELDS])
HITS = [set(), set()]
RECORDS = defaultdict(list)
STATS = {
 'informative': [[], []],
 'uninformative': [[], []],
 'chloroplast': [0, 0],
 'fully_aligned': [[], []],
 'contaminated': [[], []],
 'species': [defaultdict(int), defaultdict(int)],
 'contaminated_species': [defaultdict(int), defaultdict(int)],
 'snps': {
  'informative': [{}, {}],
  'uninformative': [{}, {}],
 },
 'gff3': 0,
}
class Result(object):
 'storage for properties that will be needed for output and statistics'
 results = []
 resultsdict = {}
 def __init__(self, id):
  self.id = id
  self.vcfs = []
  self.plantprot = None
  self.nr = None
  self.informative = False
  self.results.append(self)
  self.resultsdict[id] = self
 def get(self, id):
  return self.resultsdict[id]
 def __str__(self):
  string = '<Result plantprot: %s, nr: %s, informative: %s>' % (
   self.plantprot, self.nr, self.informative)
  return string
def process(*ignored):
 'process data according to instructions, args from argparse'
 print >>sys.stderr, 'starting processing'
 if DEBUGGING and DEBUGGING < DEBUGLEVEL:
  print >>sys.stderr, 'Warning: debugging messages not shown by default'
 args = getargs()
 if len(args.outfiles) < 3:
  args.outfiles.extend([REPORT1, REPORT2, REPORT3][-(3 - len(args.outfiles)):])
 tsvout1, tsvout2, tsvout3, blast1, blast2, queries = init(
  args.blastfile1, args.blastfile2, args.vcf_file, args.gff3_file,
  args.fastafile, args.outfiles[0], args.outfiles[1], args.outfiles[2])
 header = file_header(args.blastfile1, args.blastfile2, args.vcf_file)
 try:
  snp_annotation_index = header.index('SNP Annotation')
 except:
  snp_annotation_index = sys.maxint
 tsvout1.writerow(header)
 RECORDS['header'] = header
 if tsvout2:
  tsvout2.writerow(header)
 debug('vcf header: %s' % VCFHEADER)
 all_queries = extract_unique_queries([blast1, blast2])
 print >>sys.stderr, 'processing %d queries\n' % len(all_queries)
 for id in all_queries:
  debug('processing %s' % id)
  result = Result(id)
  result.plantprot = blast1.get(id, Record.Blast())
  result.nr = blast2.get(id, Record.Blast())
  try:
   row = defaultlist([id])
   vcfs = vcf_data(id)
   result.plantprot.best_hit = best_hit(result.plantprot)
   update_stats(result.plantprot, result.plantprot.best_hit, vcfs, result, 0)
   result.nr.best_hit = best_hit(result.nr)
   if blast2:
    update_stats(result.nr, result.nr.best_hit, vcfs, result, 1)
   row += parts_needed(result.plantprot.best_hit) + (parts_needed(result.nr.best_hit) if blast2 else [])
   row.extend(gff3_data(id))
   if snp_annotation_index < sys.maxint:
    row[snp_annotation_index] = ''  # extends row as necessary
   rowdict = dict(zip(header, row))
   if snp_annotation_index < sys.maxint and 'POS' in rowdict and \
    not rowdict['Query_start'] <= int(rowdict['POS']) <= rowdict['Query_end']:
    row[snp_annotation_index] = 'NF'
   elif result.plantprot.best_hit.title or result.nr.best_hit.title:
    writerows(list(row), [[]], tsvout3, partial = False)
    writerows(list(row), vcfs, tsvout1, partial = True)
    if vcfs and vcfs[0] and vcfs[0][0]:
     writerows(list(row), vcfs, tsvout2, partial = True)
     result.vcfs = vcfs
   else:
    debug('unreported result %s' % result)
   RECORDS[id] = list(row)
  except Exception:
   raise
   die(result.plantprot, sys.exc_info())
 stats(blast1, blast2, queries, HITS, args)
 writeblastxml(args, queries, blast1)
def writeblastxml(args, queries, blast1):
 blast2xml.FASTADATA.update(queries)
 blastpath = os.path.basename(args.blastfile1)
 xmlout = open(blastpath + '.xml', 'w')
 sys.stdout = xmlout
 blast2xml.convert(args.blastfile1, args.fastafile, blast1.values(),
  RECORDS['contaminated'])
def getargs():
 parser = argparse.ArgumentParser()
 parser.add_argument('--blast', action = 'store', dest = 'blastfile1')
 parser.add_argument('--nrblast', action = 'store', dest = 'blastfile2')
 parser.add_argument('--vcf', action = 'store', dest = 'vcf_file')
 parser.add_argument('--gff3', action = 'store', dest = 'gff3_file')
 parser.add_argument('--fasta', action = 'store', dest = 'fastafile')
 parser.add_argument('--output', action = 'append', dest = 'outfiles',
  default = [])
 args = parser.parse_args()
 return args
def extract_unique_queries(csvdata):
 queries = []
 for spreadsheet in csvdata:
  for record in spreadsheet.values():
   if not record.query in queries:
    queries.append(record.query)
 return queries
def write_ids(filename, ids):
 outfile = open(filename, 'w')
 for id in ids:
  print >>outfile, id
 outfile.close()
def write_contaminated(filename):
 outfile = open(filename, 'w')
 csvout = csv.writer(outfile, delimiter = '\t')
 csvout.writerow(RECORDS['header'])
 for record in RECORDS['contaminated']:
  csvout.writerow(RECORDS[record])
 outfile.close()
def writefasta(filename, querydict, keys):
 fastafile = open(filename, 'w')
 for key in keys:
  print >>fastafile, '>' + key
  print >>fastafile, querydict[key]
 fastafile.close()
def get_median(lengths):
 odd, halfway = len(lengths) % 2, len(lengths) / 2
 if odd:  # odd number
  return lengths[halfway]  # e.g. for length 3, returns lengths[1]
 else:
  return 0.5 * (lengths[halfway - 1] + lengths[halfway])
def get_n50(lengths):
 '''http://en.wikipedia.org/wiki/N50_statistic

    >>> get_n50([2, 2, 2, 3, 3, 4, 8, 8])
    6.0
 '''
 n50_lengths = []
 for n in lengths:
  n50_lengths.extend([n] * n)
 return get_median(n50_lengths)
def stats(blast1, blast2, queries, hits, args):
 debug('records, blast1: %d, blast2: %d' % (len(blast1), len(blast2)))
 logfile = open(LOGFILE, 'w')
 print >>logfile, 'Total number of query sequences: %d' % len(queries)
 runs = [blast1, blast2.values() if blast2 else []]
 no_hits = [set(queries.keys()) - hits[0], set(queries.keys()) - hits[1]]
 write_ids(NOHITS, no_hits[0] & no_hits[1])
 lengths = sorted(map(len, queries.values())) or [0]
 sum_length = sum(lengths)
 average = float(sum_length) / (len(lengths) + NONZERO)
 print >>logfile, 'Average length of query sequences: %.2f' % average
 shortest, longest, median, n50 = \
  lengths[0], lengths[-1], get_median(lengths), get_n50(lengths)
 print >>logfile, \
   'Length of longest query: %d, shortest: %d, median: %.1f, n50: %.1f' % (
  longest, shortest, median, n50)
 informative, uninformative, no_hit = find_informative(queries, hits)
 print >>logfile, \
  'Total number of queries with at least 1 informative hit: %d' % \
  len(informative)
 print >>logfile, 'Total number of queries without a hit: %d' % len(no_hit)
 print >>logfile, 'Total number of queries with uninformative hits: %d' % \
  len(uninformative)
 print >>logfile, \
  'Number of query sequences without a hit, plantprot: %d, nr: %d' % (
  len(no_hits[0]), len(no_hits[1]))
 print >>logfile, \
  'Sequences with uninformative hits, plantprot: %d, nr: %d' % \
  tuple(map(len, STATS['uninformative']))
 print >>logfile, 'Chloroplast sequences, plantprot: %d, nr: %d' % \
  tuple(STATS['chloroplast'])
 if False: print >>logfile, \
  'Fully aligned sequences (full length), plantprot: %d, nr: %d' % \
  tuple(map(len, STATS['fully_aligned']))
 print >>logfile, \
  'Informative fully aligned sequences (full length), plantprot: %d, nr: %d' % \
  tuple(map(is_informative, STATS['fully_aligned']))
 if args.vcf_file:
  total = sum([len(r.vcfs) for r in Result.results])
  print >>logfile, 'Total SNPs identified: %d' % total
  average = float(total) / len(queries)
  print >>logfile, 'Average SNPs per query sequence: %.2f' % average
  snphits = [sum([STATS['snps']['informative'][0][k][1] > 0
   for k in STATS['snps']['informative'][0]]),
   sum([STATS['snps']['informative'][1][k][1] > 0
   for k in STATS['snps']['informative'][1]])]
  print >>logfile, \
   'Informative hits with at least one SNP, plantprot: %d, nr: %d' % (
   snphits[0], snphits[1])
 if args.gff3_file:
  print >>logfile, 'Number of sequences predicted to be ' + \
   'full-length (independent of similarity assessment): %d' % STATS['gff3']
 print >>logfile, 'Sequences annotated as bacterial/fungal, plantprot: ' + \
  '%d, nr: %d' % tuple(map(len, STATS['contaminated']))
 write_contaminated(CONTAMINATED)
 if not len(STATS['contaminated_species'][0]) == 0:
  print >>sys.stderr, '***** WARNING ***** plantprot hits are contaminated'
 print >>logfile, 'Hits by contaminated species,\n plantprot:\n%s nr:\n%s' % \
  (prettydict(STATS['contaminated_species'][0]),
   prettydict(STATS['contaminated_species'][1]))
 print >>logfile, 'Hits by species,\n plantprot:\n%s\n nr:\n%s' % \
  (prettydict(STATS['species'][0]), prettydict(STATS['species'][1]))
 logfile.close()
def find_informative(queries, hits):
 informative = reduce(set.__or__, map(set, STATS['informative']))
 uninformative = reduce(set.__or__, map(set, STATS['uninformative']))
 uninformative -= informative  # if informative in just one, it's "informative"
 no_hit = set(queries.keys()) - (hits[0] | hits[1])
 if hits[0] | hits[1] != informative | uninformative:
  print >>sys.stderr, 'hit discrepancy: %s' % ((hits[0] | hits[1]) ^
   (informative | uninformative))
  for query in ((hits[0] | hits[1]) ^ (informative | uninformative)):
   print >>sys.stderr, '%s in hits[0]: %s' % (query, query in hits[0])
   print >>sys.stderr, '%s in hits[1]: %s' % (query, query in hits[1])
   print >>sys.stderr, '%s in informative: %s' % (query, query in informative)
   print >>sys.stderr, '%s in uninformative: %s' % (query,
    query in uninformative)
   print >>sys.stderr, '%s in no_hit: %s' % (query, query in no_hit)
 if len(informative) + len(uninformative) + len(no_hit) != len(queries):
  raise Exception, 'problem with records: %s' % (
   set(queries.keys()) ^ (informative | uninformative | no_hit))
 return informative, uninformative, no_hit
def update_stats(record, alignment, vcfs, result, i):
 if not record.query:
  return
 if is_contaminated(getattr(alignment, 'species', alignment.title)):
  if not record.query in RECORDS['contaminated']:
   RECORDS['contaminated'].append(record.query)
  STATS['contaminated'][i].append(alignment.title)
  STATS['contaminated_species'][i][getattr(alignment, 'species', '')] += 1
 informative = not alignment.title.endswith('uninformative')
 record.informative = informative  # not same as result.informative!
 index = ['uninformative', 'informative'][informative]
 hsp = (alignment.hsps + [None])[0]
 query_length = int(record.query_letters or '0')
 query_start = int(hsp.query_start) if hsp else 0
 query_end = int(hsp.query_end) if hsp else 0
 if hsp and not record.query in RECORDS['contaminated']:
  HITS[i].add(record.query)
  STATS[index][i] += [record.query]
  if informative:
   STATS['species'][i][getattr(alignment, 'species', '')] += 1
   result.informative = True
  STATS['chloroplast'][i] += 'chloroplast' in alignment.title
  if query_length > 10 and abs(query_end - query_start) >= query_length - 10:
   STATS['fully_aligned'][i].append(record.query)
  if False: debug('query length: %d, alignment length: %d (%d)' % (
   query_length, abs(query_end - query_start), alignment.length), 1)
  STATS['snps'][index][i][record.query] = [alignment.title,
   len(filter(None, vcfs))]
 else:
  print >>sys.stderr, 'record %s has no HSP: %s' % (record.query, record)
def is_contaminated(title):
 words = title.lower().split()
 if len(words) == 2:
  return words[0] in CONTAMINANTS
 else:
  debug('is_contaminated() words: %s' % words)
  return set(words) & CONTAMINANTS
def prettydict(dictionary):
 result = ''
 keys = dictionary.keys()
 keys.sort(key = lambda k: dictionary[k], reverse = True)
 for key in keys:
  result += ' %-50s: %6d\n' % (key, dictionary[key])
 return result
def writerows(row, vcfs, tsvout, partial = False):
 for vcf in vcfs:
  if tsvout:
   tsvout.writerow(row + vcf)
  if partial:
   row = [''] * len(row)
def parts_needed(alignment = None):
 parts = ['title', 'hsp.identities', 'alignment.length', 'mismatches',
  'hsp.query_start', 'hsp.query_end', 'hsp.sbjct_start',
  'hsp.sbjct_end', 'hsp.expect', 'hsp.bits', 'description', 'species']
 if alignment and len(alignment.hsps):
  title, description = alignment.title.split(None, 1)
  hsp = alignment.hsps[0]
  mismatches = getattr(hsp, 'mismatches', '')
  species = getattr(alignment, 'species', 'unspecified')
  return map(eval, parts)
 else:
  return [''] * len(parts)
def init(blastfile1, blastfile2, vcf_file, gff3_file, fastafile,
 outfile1, outfile2, outfile3):
 CONTAMINANTS.update(contaminants())
 DATA['vcf'] = vcfread(vcf_file)
 DATA['gff3'] = gff3read(gff3_file)
 print >>sys.stderr, 'data loaded'
 output1 = open(outfile1, 'w') if outfile1 else None
 output2 = open(outfile2, 'w') if outfile2 else None
 output3 = open(outfile3, 'w') if outfile3 and vcf_file else None
 tsvout1 = csv.writer(output1, delimiter = '\t') if output1 else None
 tsvout2 = csv.writer(output2, delimiter = '\t') if output2 else None
 tsvout3 = csv.writer(output3, delimiter = '\t') if output3 else None
 print >>sys.stderr, 'loading blast file 1'
 blast = filter(lambda r: r.alignments,
  PARSER.parse_file(blastfile1, fastafile))
 blast1 = OrderedDict([[r.query, r] for r in blast])
 print >>sys.stderr, 'loading blast file 2'
 blast = filter(lambda r: r.alignments,
  PARSER.parse_file(blastfile2, fastafile))
 blast2 = OrderedDict([[r.query, r] for r in blast])
 print >>sys.stderr, 'loading fasta file'
 queries = readfasta(fastafile)
 return tsvout1, tsvout2, tsvout3, blast1, blast2, queries
def die(record, exception = sys.exc_info()):
 'print informative error message pointing to problematic record'
 if record is not None:
  print >>sys.stderr, 'error at record %s' % vars(record)
  print >>sys.stderr, 'alignments: %s' % [vars(a) for a in record.alignments]
 raise exception
def gff3_data(sequence_id, gff3_data = None):
 'fill in the blanks with GFF3 data if present'
 gff3_data = ['no', '', '', ''] if not gff3_data else gff3_data
 if DATA['gff3'].has_key(sequence_id):
  datadict = DATA['gff3'][sequence_id]
  gff3_data = ['yes', datadict['start'], datadict['stop'], datadict['orf']]
  STATS['gff3'] += 1
 return gff3_data
def vcf_data(sequence_id, vcf_data = None):
 vcf_data = [[''] * len(VCFHEADER)] if not vcf_data else vcf_data
 'fill in the blanks with VCF data if present'
 if DATA['vcf'].has_key(sequence_id):
  vcf_data = DATA['vcf'][sequence_id]
 return vcf_data
def best_hit(record, dont_care = INFORMATIVE_DOESNT_MATTER):
 alignments = record.alignments
 informative = [a for a in alignments if not a.title.endswith('uninformative')]
 pool = informative if informative and not dont_care else alignments
 if pool:
  if is_contaminated(getattr(pool[0], 'species', pool[0].title)) \
   and ((pool[0].hsps and pool[0].hsps[0].expect > CUTOFF) \
   or record.query in STATS['informative'][0] + STATS['uninformative'][0]):
   chosen = Record.Alignment()
  else:
   chosen = pool[0]
 else:
  chosen = Record.Alignment()  # empty alignment record
  debug('WARNING: no alignments found: %s' % vars(record))
  # note that WARNING means nothing if empty NR BLAST record
 return chosen
def gff3read(filename):
 data = {}
 if not filename or filename == 'None':
  return data
 debug('reading in GFF3 file "%s"' % filename)
 rawdata = tsvread(filename)
 for line in rawdata:
  if line[2] == 'gene':
   data[line[0]] = {'start': line[3], 'stop': line[4], 'orf': orf_id(line[8])}
 return data
def vcfread(filename):
 data = {}
 if not filename or filename == 'None':
  return data
 rawdata = tsvread(filename)
 for line in rawdata:
  if line[0].startswith('##'):
   continue
  elif line[0].startswith('#'):
   if not VCFHEADER:
    VCFHEADER.extend([line[0][1:]] + line[1:])
  else:
   linedict = dict(zip(VCFHEADER, line))
   if linedict['FILTER'] == 'PASS' and \
    'INDEL' not in linedict['INFO'].split(';'):
    row = [linedict[c] for c in VCFHEADER if c not in ['CHROM', 'ID']]
    if not data.has_key(linedict['CHROM']):
     data[line[0]] = [row]
    else:
     data[line[0]] += [row]
 return data
def orf_id(string):
 orf = urllib.unquote(string).split()[-1]
 if not orf.startswith('m.'):
  raise Exception, 'Bad ORF identifier: %s in %s' % (orf, string)
 return orf
def tsvread(filename):
 input = open(filename)
 tsvin = csv.reader(input, delimiter = '\t')
 data = [row for row in tsvin]
 input.close()
 return filter(None, data)
def file_header(blastfile1, blastfile2, vcf_file):
 vcfheader = [c for c in VCFHEADER if c not in ['CHROM', 'ID']]
 header = blast_header(blastfile1) + blast_header(blastfile2)[1:] + \
  gff3_header()
 if vcf_file:
  header += ['SNP Annotation']
 header += vcfheader
 return header
def blast_header(filename):
 if filename is None or filename == 'None':
  header = []
 else:
  header = tsv_blast_parser.FIELDS + ['Species']
  header.remove('Gap_openings')
 return header
def gff3_header():
 return ['has GFF3 data', 'GFF3 start', 'GFF3 stop', 'GFF3 ORF']
def clean_hit(string):
 while string[0] == '>' or string[0] == ' ':
  string = string[1:]
 return ' '.join(string.split())  # get rid of newlines and extra spaces
def readlines(filename):
 input = open(filename)
 data = map(str.strip, input.readlines())
 input.close()
 return data
def contaminants():
 names = []
 for url in CONTAMINATION_URLS:
  names.extend(firstnames(parsenames.parse(url)))
 names.extend(map(str.lower, readlines(CONTAMINANTS_FILE)))
 return set(names)
def firstnames(names):
 firstnames = []
 for name in names:
  debug('contaminant name: %s' % name)
  match = re.compile('[(](.*)[)]').match(name)
  if match:
   firstname = match.group(1).split()[0].lower()
   firstnames.append(firstname)
  else:
   if '(' in name:
    print >>sys.stderr, 'unexpected "(" in contaminant name %s' % name
   firstname = name.split()[0].lower()
  firstnames.append(firstname)
 return firstnames
def is_informative(queries):
 return len(filter(lambda q: not q.endswith('uninformative'), queries))
def _test():
 import doctest
 doctest.testmod()
if __name__ == '__main__':
 command = os.path.splitext(os.path.basename(sys.argv[0]))[0]
 job = process  # default for job.py
 print eval(command)(*sys.argv[1:])
