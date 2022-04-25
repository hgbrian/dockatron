"""convert smiles to sdf

Adapted from http://bits.csb.pitt.edu/tdtCDPK1/summary.pdf for python 3, import.
"""
import sys
import tempfile
from optparse import OptionParser
from pathlib import Path

from rdkit.Chem import AllChem as Chem

# silence C++ output https://github.com/rdkit/rdkit/issues/2683
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def getRMS(mol, c1, c2):
    rms, _trans = Chem.GetAlignmentTransform(mol, mol, c1, c2)
    return rms


def rdconf(smifile_str,
        maxconfs=25, sample_multiplier=1, seed=9162006,
        rms_threshold=0, energy_window=-1, is_verbose=False,
        is_mmff=False, is_nomin=False, is_etkdg=False):

    output_fh = tempfile.NamedTemporaryFile(delete=False)

    sdwriter = Chem.SDWriter(output_fh.name)

    if sdwriter is None:
        print("Could not open SDWriter")
        sys.exit(-1)

    for line in smifile_str.splitlines():
        toks = line.split()
        smi = toks[0]
        name = ' '.join(toks[1:])

        pieces = smi.split('.')
        if len(pieces) > 1:
            smi = max(pieces, key=len) #take largest component by length
            print(f"Taking largest component: {smi}\t{name}")

        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            if is_verbose:
                print(smi)
            try:
                Chem.SanitizeMol(mol)
                mol = Chem.AddHs(mol)
                mol.SetProp("_Name", name)

                if is_etkdg:
                    cids = Chem.EmbedMultipleConfs(mol, int(sample_multiplier*maxconfs), Chem.ETKDG())
                else:
                    cids = Chem.EmbedMultipleConfs(mol, int(sample_multiplier*maxconfs),randomSeed=seed)
                if is_verbose:
                    print(len(cids),"conformers found")

                cenergy = []
                for conf in cids:
                    #not passing confID only minimizes the first conformer
                    if is_nomin:
                        cenergy.append(conf)
                    elif is_mmff:
                        converged = Chem.MMFFOptimizeMolecule(mol, confId=conf)
                        mp = Chem.MMFFGetMoleculeProperties(mol)
                        cenergy.append(Chem.MMFFGetMoleculeForceField(mol, mp, confId=conf).CalcEnergy())
                    else:
                        converged = not Chem.UFFOptimizeMolecule(mol, confId=conf)
                        cenergy.append(Chem.UFFGetMoleculeForceField(mol, confId=conf).CalcEnergy())
                    if is_verbose:
                        print("Convergence of conformer",conf,converged)

                # BN: I think this is ok / expected??
                mol = Chem.RemoveHs(mol)
                sortedcids = sorted(cids, key = lambda cid: cenergy[cid])
                if len(sortedcids) > 0:
                    mine = cenergy[sortedcids[0]]
                else:
                    mine = 0
                if rms_threshold == 0:
                    cnt = 0
                    for conf in sortedcids:
                        if cnt >= maxconfs:
                            break
                        if (energy_window < 0) or (cenergy[conf]-mine <= energy_window):
                            sdwriter.write(mol,conf)
                            cnt+=1
                else:
                    written = {}
                    for conf in sortedcids:
                        if len(written) >= maxconfs:
                            break
                        #check rmsd
                        passed = True
                        for seenconf in written.keys():
                            rms = getRMS(mol,seenconf,conf)
                            if (rms < rms_threshold) or (energy_window > 0 and cenergy[conf]-mine > energy_window):
                                passed = False
                                break
                        if passed:
                            written[conf] = True
                            sdwriter.write(mol,conf)
            except (KeyboardInterrupt, SystemExit) as err:
                raise err
            except Exception as err:
                print(f"Exception {err}")
        else:
            print(f"ERROR: {smi}")

    #sdwriter.close()
    #outf.close()

    output_fh.flush()
    output_fh.seek(0)

    # return output_fh.read()
    return Path(output_fh.name).read_text()


def rdconf_old(smifile_str,
        maxconfs=25, sample_multiplier=1, seed=9162006,
        rms_threshold=0, energy_window=-1, is_verbose=False):
    """Weird syntax because SDWriter writes to the output filehandle"""

    output_fh = tempfile.NamedTemporaryFile(delete=False)

    sdwriter = Chem.SDWriter(output_fh.name)
    if sdwriter is None:
        raise SystemExit(f"Could not open {output_fh}")

    for line in smifile_str.splitlines():
        toks = line.split()
        smi = toks[0]
        name = ' '.join(toks[1:])

        mol = Chem.MolFromSmiles(smi)

        if mol is not None:
            if is_verbose:
                print(f"smi: {smi}")

            try:
                Chem.SanitizeMol(mol)
                #mol = Chem.AddHs(mol) # added Brian 2022-04-13
                mol.SetProp("_Name", name)
                cids = Chem.EmbedMultipleConfs(mol, int(sample_multiplier*maxconfs), randomSeed=seed)

                if is_verbose:
                    print(f"{len(cids)} conformers found")

                cenergy = []
                for conf in cids:
                    # not passing confID only minimizes the first conformer
                    converged = not Chem.UFFOptimizeMolecule(mol, confId=conf)
                    if is_verbose:
                        print(f"Convergence of conformer {conf} {converged}")
                    cenergy.append(Chem.UFFGetMoleculeForceField(mol, confId=conf).CalcEnergy())

                sortedcids = sorted(cids, key=lambda cid: cenergy[cid])

                if(rms_threshold == 0):
                    cnt = 0
                    for conf in sortedcids:
                        if cnt >= maxconfs:
                            break
                        if (energy_window < 0) or (cenergy[conf]-cenergy[0] <= energy_window):
                            sdwriter.write(mol, conf)
                            cnt += 1
                else:
                    written = {}
                    for conf in sortedcids:
                        if len(written) >= maxconfs:
                            break

                        # check rmsd
                        is_passed = True
                        for seenconf in written.iterkeys():
                            rms = getRMS(mol, seenconf, conf)
                            if (rms < rms_threshold) or (energy_window > 0 and cenergy[conf]-cenergy[0] > energy_window):
                                is_passed = False
                                break
                        if is_passed:
                            written[conf] = True
                            sdwriter.write(mol, conf)

            except (KeyboardInterrupt, SystemExit) as err:
                raise err
            except Exception as err:
                print(f"Exception occurred: {sys.exc_info()[0]} {err}")
                raise err
        else:
            raise Exception(f"Unknown ERROR: {smi}")

    output_fh.flush()
    output_fh.seek(0)
    # return output_fh.read()
    return Path(output_fh.name).read_text()


if __name__ == "__main__":
    parser = OptionParser(usage="Usage: %prog [options] <input>.smi <output>.sdf")
    parser.add_option("--maxconfs", dest="maxconfs", action="store",
                      help="maximum number of conformers to generate per a molecule (default 25)", default="25", type="int", metavar="CNT")
    parser.add_option("--sample_multiplier", dest="sample", action="store",
                      help="sample N*maxconfs conformers and choose the maxconformers with lowest energy (default 1)", default="1", type="float", metavar="N")
    parser.add_option("--seed", dest="seed", action="store",
                      help="random seed (default 9162006)", default="9162006", type="int", metavar="s")
    parser.add_option("--rms_threshold", dest="rms", action="store",
                      help="filter based on rms (default 0)", default="0", type="float", metavar="R")
    parser.add_option("--energy_window", dest="energy", action="store",
                      help="filter based on energy difference with lowest energy conformer", default="-1", type="float", metavar="E")
    parser.add_option("-v", "--verbose", dest="is_verbose", action="store_true", default=False,
                      help="verbose output")
    # 3 new options
    parser.add_option("--mmff", dest="mmff",action="store_true",default=False,
                    help="use MMFF forcefield instead of UFF")
    parser.add_option("--nomin", dest="nomin",action="store_true",default=False,
                    help="don't perform energy minimization (bad idea)")
    parser.add_option("--etkdg", dest="etkdg",action="store_true",default=False,
                    help="use new ETKDG knowledge-based method instead of distance geometry")


    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("Need input and output")
        raise SystemExit

    print(f"Running with options {options}")

    input_filename = args[0]
    output_filename = args[1]

    smifile_str = open(input_filename).read()
    if options.is_verbose:
        print(f"Generating a maximum of {options.maxconfs} per a mol")

    sdfs = rdconf(smifile_str, maxconfs=options.maxconfs, sample_multiplier=options.sample,
                  seed=options.seed, rms_threshold=options.rms,
                  energy_window=options.energy, is_verbose=options.is_verbose,
                  is_mmff=options.mmff, is_nomin=options.nomin, is_etkdg=options.etkdg)

    with open(output_filename, 'w') as out:
        out.write(sdfs)
