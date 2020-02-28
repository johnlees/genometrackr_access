#!python3

import os, sys
import pandas as pd
import re
import time

from Bio import Entrez
Entrez.email = 'j.lees@imperial.ac.uk'
Entrez.api_key = '82fc4fe1eb04a8801146737727f3b11a0b08'

def runQuery(assemblies, names, outfile):
    esearch_term = " OR ".join(assemblies)
    success = 0
    for attempt in range(5):
        try:
            handle = Entrez.esearch(db='Assembly', term=esearch_term)
            ret = Entrez.read(handle)['IdList']
            handle2 = Entrez.esummary(db='Assembly', id=",".join(ret), retmode='xml')
            records = Entrez.read(handle2)
            success = 1
            break
        except RuntimeError:
            time.sleep(3)
        except Entrez.Parser.ValidationError:
            sys.stderr.write("Entrez parse error\n")
            break

    if success == 1:
        for (record, name) in zip(records['DocumentSummarySet']['DocumentSummary'], names):
            try:
                ftp_loc = record['FtpPath_GenBank']
                rsync_match = re.search(r'genomes/all/(.+)$', ftp_loc)
                if rsync_match != None:
                    rsync_loc = rsync_match.group(1) + "/" + os.path.basename(ftp_loc) + "_genomic.fna.gz"
                    outfile.write(rsync_loc + "\n")
                else:
                    sys.stderr.write("Failed to parse FTP path " + ftp_loc + "\n")
            except (IndexError, KeyError):
                sys.stderr.write("Could not find assembly on FTP for " + sample + "\n")
                continue
    else:
        sys.stderr.write("FAILED:\n")
        sys.stderr.write("\n".join(names) + "\n")

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Download assemblies from NCBI ftp')
    parser.add_argument('clusterfile', help='File with assemblies to download')

    options = parser.parse_args()

    sys.stderr.write("Reading cluster file\n")
    cluster_def = pd.read_csv(options.clusterfile, '\t', header=0, index_col=0, comment=None, engine='c')

    row_it = cluster_def.iterrows()
    with open('rsync_files.txt', 'w') as outfile:
        while True:
            try:
                assemblies = []
                names = []
                for query_size in range(20):
                    table_row = next(row_it)
                    assembly = table_row[1]['fullasm_id']
                    name = table_row[0].split("|")[0]
                    assemblies.append(str(assembly))
                    names.append(name)
                runQuery(assemblies, names, outfile)
            except StopIteration:
                if len(assemblies) > 0:
                    runQuery(assemblies, names, outfile)
                break

    print('rsync --copy-links --times --verbose --files-from=rsync_files.txt rsync:://ftp.ncbi.nlm.nih.gov/genomes/all/')

if __name__ == '__main__':
    main()
