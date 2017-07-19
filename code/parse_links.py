# author: Felipe Zapata

import sys
from collections import defaultdict
taxa = defaultdict(list)
for line in sys.stdin:
  if line.startswith('>'):
    line = line[1:].rstrip()
    taxon, _, seq = line.partition('@')
    taxa[taxon].append(seq)
for taxon in sorted(taxa):
  print "%s:%s" % (taxon, ';'.join("%s@%s" % (taxon, seq) for seq in taxa[taxon]))
