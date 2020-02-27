#!python3

import os, sys
import pandas as pd
import re
import time

from Bio import Entrez
Entrez.email = 'j.lees@imperial.ac.uk'
Entrez.api_key = '82fc4fe1eb04a8801146737727f3b11a0b08'

def runQuery(assemblies, names, output):
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

        if success == 1:
            for (record, name) in zip(records['DocumentSummarySet']['DocumentSummary'], names):
                try:
                    ftp_loc = record['FtpPath_GenBank']
                    rsync_loc = re.sub(r'^ftp', 'rsync', ftp_loc, 1)
                    print("rsync --times --verbose " + rsync_loc + " " + output + "/" + name + "_genomic.fna.gz")
                except (IndexError, KeyError):
                    sys.stderr.write("Could not find assembly on FTP for " + sample + "\n")
                    continue
        else:
            sys.stderr.write("FAILED:\n")
            sys.stderr.write("\n".join(names))

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Download assemblies from NCBI ftp')
    parser.add_argument('clusterfile', help='File with assemblies to download')
    parser.add_argument('--output', default="./", help='Output prefix')

    options = parser.parse_args()

    sys.stderr.write("Reading cluster file\n")
    cluster_def = pd.read_csv(options.clusterfile, '\t', header=0, index_col=0, comment=None, engine='c')

    row_it = cluster_def.iterrows()
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
            runQuery(assemblies, names, options.output)
        except StopIteration:
            if len(assemblies) > 0:
                runQuery(assemblies, names, options.output)
            break

if __name__ == '__main__':
    main()
