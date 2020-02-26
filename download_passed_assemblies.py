#!python3

import os, sys
import pandas as pd
import re
import time

from ftplib import FTP, error_perm, error_temp

def ftp_reload(ftp_obj):
    ftp_obj.close()

    success = 1
    for attempt in range(5):
        time.sleep(5)
        try:
            ftp_new = FTP('ftp.ncbi.nlm.nih.gov')
            ftp_new.login()
            success = 1
            break
        except (EOFError, error_temp, error_perm):
            pass

    if success == 0:
        raise RunimeError("Could not connect to FTP after five tries")

    return ftp_new


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Download assemblies from NCBI ftp')
    parser.add_argument('clusterfile', help='File with assemblies to download')
    parser.add_argument('--output', default="./", help='Output prefix')

    options = parser.parse_args()

    sys.stderr.write("Reading cluster file\n")
    cluster_def = pd.read_csv(options.clusterfile, '\t', header=0, index_col=0, comment=None, engine='c')

    downloaded = 0
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    for table_row in cluster_def.iterrows():
        try:
            ftp.cwd('/genomes/all')
        except EOFError:
            ftp = ftp_reload(ftp)
            ftp.cwd('/genomes/all')
        assembly = table_row[1]['asm_acc']
        sample = table_row[1]['biosample_acc']
        name = table_row[0].split("|")[0]
        full_name = table_row[0].split("|")[0:4]
        try:
            if pd.notna(assembly):
                sys.stderr.write("Downloading " + name + "\n")
                m = re.match(r"^(GCA|GCF)_(\d\d\d)(\d\d\d)(\d\d\d)\..+$", assembly)
                ftp_path = "/".join([m.group(1), m.group(2), m.group(3), m.group(4)])
                # Find file on FTP
                # FTP access
                ftp.cwd(ftp_path)
                directory = ftp.nlst()
                ftp.cwd(directory[0])
                file_name = directory[0] + "_genomic.fna.gz"
                #files = ftp.nlst()
                #for file in files:
                #    fna_m = re.search(r"_genomic\.fna\.gz$", file)
                #    if fna_m:
                #        file_name = file
                #        break

                #if ftp.size(file_name) == None:
                #    sys.stderr.write("Could not find FTP file " + "/".join([ftp_path, directory, file_name]) + "\n")
                #else:
                target_file = os.path.join(options.output, os.path.basename(name)) + "_genomic.fna.gz"
                with open(target_file ,'wb') as fhandle:
                    ftp.retrbinary('RETR %s' % file_name, fhandle.write)
                downloaded += 1
        except (error_temp, error_perm):
            sys.stderr.write("Could not find " + name + " | " + assembly + "\n")

    print("Downloaded " + int(downloaded) + " assemblies\n")


if __name__ == '__main__':
    main()
