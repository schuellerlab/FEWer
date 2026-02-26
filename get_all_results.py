#! /usr/bin/env python3
#
# get_all_results.py - Retrieve all free energy calculation results from FAR or FEW
#                      MMPBSA and present a summary table with mean deltaG [kcal/mol],
#                      standard deviation, Ki [nM] (at 300 K) and pKi.
#
# AUTHOR
# Andreas Schueller <aschueller@uc.cl>
#
# HISTORY
# 2026-02-26    0.6    Andreas    Added option to print individual MD runs, not means
#                                 Also fixed a bug with RMSD_last_20 that actually
#                                 returned the mean of the last 20 RMSDs of only the
#                                 last replicate, instead of all replicates
# 2025-06-24    0.5    Andreas    Added "n" (counts) to output
# 2025-05-03    0.4.1  Andreas    Fixed bug with FEW version, output filename is
#                                 *_statistics.out 
# 2025-04-29    0.4    Andreas    Adapted to work with FEW results
# 2025-04-19    0.3    Andreas    Added reporting of ligand RMSD
# 2025-04-05    0.2    Andreas    Improved version to handle replicates in separate
#                                 subdirectories instead of duplicated ligands
# 2025-03-22    0.1    Andreas    First version
#
import subprocess
import statistics
import math
import argparse
import os

def stdev(nums):
    try:
        sd = statistics.stdev(nums)
    except:
        sd = -1
    return sd

def print_line(lig, rec, deltaGs, rmsds):
    deltaG = statistics.mean(deltaGs)
    sd = stdev(deltaGs)
    if rmsds:
        rmsds_flat = [item for sublist in rmsds for item in sublist] # flatten 2d list
        mean_rmsd_all = statistics.mean(rmsds_flat)
        sd_rmsd_all = stdev(rmsds_flat)
        print(len(rmsds_flat))
        rmsds_flat_20 = [item for sublist in rmsds for item in sublist[-20:]] # flatten 2d list, only last 20 entries
        mean_rmsd_20 = statistics.mean(rmsds_flat_20)
        sd_rmsd_20 = stdev(rmsds_flat_20)
        print(len(rmsds_flat_20))
    else:
        mean_rmsd_all = ""
        sd_rmsd_all = ""
        mean_rmsd_20 = ""
        sd_rmsd_20 = ""
    Ki = math.exp(deltaG*1000/1.985/300)
    pKi = -1 * math.log10(Ki)
    print(lig, rec, deltaG, sd, Ki * 1e9, pKi, len(deltaGs), mean_rmsd_20, sd_rmsd_20, len(rmsds_flat_20), mean_rmsd_all, sd_rmsd_all, len(rmsds_flat), sep="\t")

def print_deltaGs(deltaGs, rmsds):
    headers = ["Ligand", "Receptor", "Average(deltaG)[kcal/mol]", "StDev(deltaG)[kcal/mol]", "Ki[nM]", "pKi", "n(deltaG)", "Average(ligang_RMSD_last_20)[A]", "StDev(ligand_RMSD_last_20)[A]", "n(RMSD_last_20)", "Average(ligang_RMSD)[A]", "StDev(ligand_RMSD)[A]", "n(RMSD)"]
    print(*headers, sep="\t")
    for rec in deltaGs:
        for lig in deltaGs[rec]:
            if rmsds:
                print_line(lig, rec, deltaGs[rec][lig], rmsds[rec][lig])
            else:
                print_line(lig, rec, deltaGs[rec][lig], [])

def print_line_replicate(lig, rec, replicate, deltaG, rmsds):
    if rmsds:
        mean_rmsd_all = statistics.mean(rmsds)
        sd_rmsd_all = stdev(rmsds)
        mean_rmsd_20 = statistics.mean(rmsds[-20:])
        sd_rmsd_20 = stdev(rmsds[-20:])
    else:
        mean_rmsd_all = ""
        sd_rmsd_all = ""
        mean_rmsd_20 = ""
        sd_rmsd_20 = ""
    Ki = math.exp(deltaG*1000/1.985/300)
    pKi = -1 * math.log10(Ki)
    print(lig, rec, replicate, deltaG, Ki * 1e9, pKi, 1, mean_rmsd_20, sd_rmsd_20, len(rmsds[-20:]), mean_rmsd_all, sd_rmsd_all, len(rmsds), sep="\t")

def print_deltaGs_details(deltaGs, rmsds):
    headers = ["Ligand", "Receptor", "Replicate", "deltaG[kcal/mol]", "Ki[nM]", "pKi", "n(deltaG)", "Average(ligang_RMSD_last_20)[A]", "StDev(ligand_RMSD_last_20)[A]", "n(RMSD_last_20)", "Average(ligang_RMSD)[A]", "StDev(ligand_RMSD)[A]", "n(RMSD)"]
    print(*headers, sep="\t")
    for rec in deltaGs:
        for lig in deltaGs[rec]:
            for replicate in range(0, len(deltaGs[rec][lig])):
                if rmsds:
                    print_line_replicate(lig, rec, replicate, deltaGs[rec][lig][replicate], rmsds[rec][lig][replicate])
                else:
                    print_line_replicate(lig, rec, replicate, deltaGs[rec][lig][replicate], [])

def parse_res(res, rep):
    deltaGs = dict()
    for line in res.split("\n"):
        if line.startswith("AffinityBindingPred"):
             continue
        elif line:
            tokens = line.split("\t")
            deltaG = float(tokens[0])
            lig = tokens[1]
            if rep == "by_lig":
                lig = lig.rsplit("_", 1)[0] # Remove trailing _<num>
            rec = tokens[2]
            if not rec in deltaGs:
                 deltaGs[rec] = dict()
            if not lig in deltaGs[rec]:
                 deltaGs[rec][lig] = []
            deltaGs[rec][lig].append(deltaG)
            #print(lig,rec,deltaG)
    return deltaGs

def parse_rmsd(files, rep):
    rmsds = dict()
    for path in files.split("\n"):
        if not path:
            continue # Ignore empty lines
        # Example path: ./run_3cs7_lig_5/run_3CS7.A_4/rec_3cs7_recep/ligand_prod.rmsd
        tokens = path.split(os.sep)
        lig = tokens[2][4:] # Remove leading run_
        if rep == "by_lig":
            lig = lig.rsplit("_", 1)[0] # Remove trailing _<num>
        rec = tokens[3][4:] # Remove leading rec_
        # Read RMSD file
        with open(path, "r") as f:
            rmsd = []
            for line in f:
                if not line.startswith("#"):
                    tokens = line.split() # 2nd column is RMSD
                    rmsd.append(float(tokens[1]))
            if not rec in rmsds:
                rmsds[rec] = dict()
            rmsds[rec][lig].append(rmsd)
    return rmsds

def parse_rmsd_few(files):
    rmsds = dict()
    for path in files.split("\n"):
        if not path:
            continue # Ignore empty lines
        # Example path: ./Ccpg_225/1/rmsd/complex_all.rmsd
        tokens = path.split(os.sep)
        rec = tokens[1]
        lig = "ligand"
        # Read RMSD file
        with open(path, "r") as f:
            rmsd = []
            for line in f:
                if not line.startswith("#"):
                    tokens = line.split() # 2nd column is RMSD
                    rmsd.append(float(tokens[1]))
            if not rec in rmsds:
                rmsds[rec] = dict()
            if not lig in rmsds[rec]:
                rmsds[rec][lig] = []
            rmsds[rec][lig].append(rmsd)
    return rmsds

def parse_res_few(res):
    deltaGs = dict()
    for line in res.split("\n"):
        # Example: ./2uwp/3/calc_a_1t/ligand/s81_100_1/pb3_gb0/ligand_statistics.outPBTOT            -22.03       2.81
        if not line:
            continue # Ignore empty lines
        tokens = line.split('\x00')
        path_tokens = tokens[0].split(os.sep)
        rec = path_tokens[1]
        lig = path_tokens[4]
        energy_tokens = tokens[1].split()
        deltaG = float(energy_tokens[1])
        std = energy_tokens[2]
        if not rec in deltaGs:
             deltaGs[rec] = dict()
        if not lig in deltaGs[rec]:
             deltaGs[rec][lig] = []
        deltaGs[rec][lig].append(deltaG)
    return deltaGs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='get_all_results.py',
                    description='Retrieve all free energy calculation results from FEWer or FAR MMPBSA and present a summary table with mean deltaG [kcal/mol], standard deviation, Ki [nM] (at 300 K) and pKi. Also reports RMSD [A] of ligand atoms for the last 20 frames (as for deltaG) and for all frames.')
    parser.add_argument('--rep', default='by_lig', help='How are replicates provided? "by_lig": Ligands are replicated inside the SDF. "by_dir": Each replicate run resides in its own subdirectory. In both cases, replicate ligands need to have exactly the same name.')
    parser.add_argument('--met', default='FEW', help='FEW or FAR. Determines the MM-PBSA method.')
    parser.add_argument('-d', '--details', action='store_true', help='Report individual values (replicates) instead of summary table (means). Implemented for FEWer only.')
    args = parser.parse_args()

    if args.met == "FEW":
        res = subprocess.check_output("find . -name *_statistics.out -print0 -exec tail -1 {} \;", shell=True, universal_newlines=True)
        deltaGs = parse_res_few(res)
        #print(deltaGs)
        rmsd_files = subprocess.check_output("find . -name complex_all.rmsd -print", shell=True, universal_newlines=True)
        rmsds = parse_rmsd_few(rmsd_files)
        #import pprint
        #pprint.pprint(rmsds)
        #print(len(rmsds['2w26_083']['ligand'][0]))
        #print(len(rmsds['2w26_083']['ligand']))
        if args.details:
            print_deltaGs_details(deltaGs, rmsds)
        else:
            print_deltaGs(deltaGs, rmsds)
    else:
        print("FAR_results")
        res = subprocess.check_output("find . -name FAR_results.tsv -exec cat {} \;", shell=True, universal_newlines=True)
        deltaGs = parse_res(res, args.rep)
        print_deltaGs(deltaGs, {})
        print("MMPBSA_results")
        res = subprocess.check_output("find . -name MMPBSA_results.tsv -exec cat {} \;", shell=True, universal_newlines=True)
        deltaGs = parse_res(res, args.rep)
        rmsd_files = subprocess.check_output("find . -name ligand_prod.rmsd -print", shell=True, universal_newlines=True)
        rmsds = parse_rmsd(rmsd_files, args.rep)
        print_deltaGs(deltaGs, rmsds)
