#!python3

import os, sys
import pandas as pd
import time

from ftplib import FTP
from Bio import Entrez
Entrez.email = 'j.lees@imperial.ac.uk'
Entrez.api_key = '82fc4fe1eb04a8801146737727f3b11a0b08'

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Download assemblies from NCBI ftp')
    parser.add_argument('clusterfile', help='File with assemblies to download')
    parser.add_argument('assemblyreport', help='Assembly report file (ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/)')
    parser.add_argument('--output', default="./", help='Output prefix')

    options = parser.parse_args()

    sys.stderr.write("Reading cluster file\n")
    cluster_def = pd.read_csv(options.clusterfile, '\t', header=0, index_col=0, comment=None, engine='c')

    sys.stderr.write("Reading assembly report file\n")
    assembly_locs = pd.read_csv(options.assemblyreport, '\t', skiprows=1, header=0, index_col=0, comment=None, engine='c')

    downloaded = 0
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('genomes/all')
    for table_row in cluster_def.iterrows():
        assembly = table_row[1]['asm_acc']
        sample = table_row[1]['biosample_acc']
        name = table_row[0].split("|")[0]
        full_name = table_row[0].split("|")[0:4]
        if pd.notna(assembly) or pd.notna(sample):
            sys.stderr.write("Downloading " + name + "\n")

            # Find file on FTP
            try:
                ftp_loc = assembly_locs.loc[assembly]['ftp_path']

            except KeyError:
                # Block for Entrez failures
                failed = 0
                for attempt in range(0,5):
                    try:
                        handle = Entrez.esearch(db='Assembly', term=sample)
                        id = Entrez.read(handle)['IdList']
                        if len(id) == 0:
                            sys.stderr.write("Could not find assembly for " + sample + "\n")
                            failed = 1
                            break
                        handle = Entrez.esummary(db='Assembly', id=id[0], retmode='xml')
                        ret = Entrez.read(handle)
                        break
                    except RuntimeError:
                        time.sleep(1)
                if failed == 1:
                    sys.stderr.write("Could not search Entrez for " + sample + "\n")
                    continue

                # Block for result format failures
                try:
                    ftp_loc = ret['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
                except (IndexError, KeyError):
                    sys.stderr.write("Could not find assembly on FTP for " + sample + "\n")
                    continue
                if ftp_loc == None or ftp_loc == '':
                    sys.stderr.write("Could not find assembly on FTP for " + sample + "\n")
                    continue

            # FTP access
            ftp_path = ftp_loc.split("/")
            for idx, path_part in enumerate(ftp_path):
                if path_part == "all":
                    ass_path = "/".join(ftp_path[idx+1:]) + "/" + ftp_path[-1] + "_genomic.fna.gz"
                    break
            if ftp.size(ass_path) == None:
                sys.stderr.write("Could not find FTP file " + ass_path + "\n")
            else:
                target_file = os.path.join(options.output, os.path.basename(name)) + "_genomic.fna.gz"
                with open(target_file ,'wb') as fhandle:
                    ftp.retrbinary('RETR %s' % ass_path, fhandle.write)
                downloaded += 1

    print("Downloaded " + int(downloaded) + " assemblies\n")


if __name__ == '__main__':
    main()
