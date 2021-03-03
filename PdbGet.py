#!/anaconda/bin/python

import urllib2
import gzip
import sys
import StringIO

if len(sys.argv) > 1:
    if type(sys.argv[1]) == type("string"):
        if len(sys.argv[1]) == 4:
            pdbcode = sys.argv[1].upper()
        else:
            sys.exit('PDB codes have 4 characters')
    else:
        sys.exit('PDB code should be a string')
else:
    sys.exit('Useage: PdbGet.py 1ABC')

url = 'http://www.rcsb.org/pdb/files/'+pdbcode+'.pdb.gz'
outfile_name = pdbcode+'.pdb'

try:
    page = urllib2.urlopen(url)
except:
    sys.exit('Cannot access PDB entry, perhaps code is invalid?')

zipped = StringIO.StringIO(page.read())
unzipped = gzip.GzipFile(fileobj=zipped)

with open(outfile_name, 'w') as outfile:
    outfile.write(unzipped.read())

print 'Normal completion. Saved '+outfile_name+'.'
    

