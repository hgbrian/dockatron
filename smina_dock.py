"""This script requires the following binary and file
!curl -L https://sourceforge.net/projects/smina/files/smina.static/download --output smina && chmod +x smina
!curl -O https://raw.githubusercontent.com/openbabel/openbabel/master/data/space-groups.txt
!mkdir -p bin && mv space-groups.txt bin/ && mv smina bin/

# To install gnina, the deep learning fork of smina
!docker pull gnina/gnina

# To convert pdb to pdbqt
!apt install openbabel python-openbabel

# To install openmm and pdb fixer
mamba install -c conda-forge openmm pdbfixer

# To make a "high quality 3D model of the ligang" these guys recommend gypsum-dl. I have not tried this yet.
# https://git.durrantlab.pitt.edu/jdurrant/webina/blob/c259bd868b033a5ba32739a2db5e0532fccfdc42/README.md
!wget https://git.durrantlab.pitt.edu/jdurrant/gypsum_dl/-/archive/1.1.7/gypsum_dl-1.1.7.tar.gz
"""

import io
import os
import re
import sys
import time
import subprocess

from os.path import join as pjoin
from pathlib import Path
from random import random
from sys import platform
from tempfile import NamedTemporaryFile, gettempdir
from typing import Optional, Any, Tuple, List

import numpy as np
import pandas as pd
import requests
from tqdm.auto import tqdm

from Bio import PDB
from scipy.spatial.distance import cdist

from apps.docking.rdconf import rdconf

SMINA = "smina"
GNINA = "gnina"
DOCK_BIN_D = {SMINA: [pjoin(os.environ['HXROOT'], "apps", "docking", "bin", "smina" + ("_osx" if platform == "darwin" else ""))],
              GNINA: ["docker", "run", "-v", "/tmp/:/tmp/", "gnina/gnina", "gnina"]}
OBABEL_BIN = pjoin(os.environ['HXROOT'], "apps", "docking", "bin", "babel" + ("_osx" if platform == "darwin" else ""))
OBABEL_BIN = "obabel" # problematic binary

DEBUG = True
DEFAULT_MAX_SDF_CONFS = 16
DEFAULT_EXHAUSTIVENESS = 64
DEFAULT_SDF_CONFS = 4
DEFAULT_NUM_CPUS = 4
DEFAULT_SEED = 1
DEFAULT_SCORING = ("default", "vinardo")[0]
TMPDIR = gettempdir()
PDB_CACHE = os.path.join(TMPDIR, "apps_docking_pdb_cache")

os.makedirs(PDB_CACHE, exist_ok=True)


def flatten(alist: list) -> List[Any]:
    return [item for sublist in alist for item in sublist]


def download_smiles(pubchem_id: str, retries:int = 5) -> str:
    """Given a pubchem id, get a smiles string"""
    while True:
        req = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{pubchem_id}/property/CanonicalSMILES/CSV")
        smiles_url_csv = req.text if req.status_code == 200 else None
        if smiles_url_csv is not None:
            break
        if retries == 0:
            break
        time.sleep(5*random())
        retries -= 1

    return smiles_url_csv.splitlines()[1].split(',')[1].strip('"').strip("'") if smiles_url_csv is not None else None


def download_sdf(pubchem_id: str) -> Tuple[str, str]:
    """Download sdf_2d and sdf_3d from pubchem with pubchem id"""
    pubchem_url_sdf = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{pubchem_id}/record/SDF/"

    sdf_2d_req = requests.get(f"{pubchem_url_sdf}?record_type=2d")
    sdf_3d_req = requests.get(f"{pubchem_url_sdf}?record_type=3d")
    sdf_2d = sdf_2d_req.text if sdf_2d_req.status_code == 200 else None
    sdf_3d = sdf_3d_req.text if sdf_3d_req.status_code == 200 else None

    return sdf_2d, sdf_3d


def get_or_download_pdb(pdb_id: str) -> str:
    """Download pdb file as a string from rcsb.org"""
    # pdb file
    if pdb_id.endswith('.pdb'):
        with open(pdb_id) as f:
            return f.read()

    # url or pdb_id
    if pdb_id.startswith('http'):
        url = pdb_id
        filename = url.split('/')[-1]
    else:
        url = f"http://files.rcsb.org/view/{pdb_id}.pdb"
        filename = f'{pdb_id}.pdb'

    cache_path = os.path.join(PDB_CACHE, filename)
    if os.path.exists(cache_path):
        with open(cache_path) as cache_f:
            return cache_f.read()

    pdb_req = requests.get(url)
    pdb_req.raise_for_status()
    open(cache_path, 'w').write(pdb_req.text)
    return pdb_req.text


def link_proteome_files(from_proteome_dir:str, to_inference_path:str, sdf_file:str):
    """Copy (soft link) files to the directory structure EquiBind requires:
    {inference_path}/{protein_names}/{pdb_id}_protein.pdb,{sdf_id}_ligand.sdf
    """
    # absolute path everything because of "cp -as" peculiarities
    from_proteome_dir = Path(from_proteome_dir).resolve()
    to_inference_path = Path(to_inference_path).resolve()
    from_sdf_path = Path(sdf_file).resolve()

    if from_sdf_path.name.endswith("_ligand.sdf"):
        to_sdf_filename = from_sdf_path.name
    else:
        to_sdf_filename = from_sdf_path.stem + "_ligand.sdf"

    # copy pdb files (many to many)
    from_pdb_dirs = Path(from_proteome_dir).glob("*")
    for from_pdb_dir in tqdm(from_pdb_dirs, desc="Copying PDB files"):
        subprocess.run(["cp", "-as", from_pdb_dir.as_posix(), to_inference_path.as_posix()], check=True)

    # copy sdf (ligand) file (one to many)
    to_sdf_dirs = Path(to_inference_path).glob("*")
    for dr in tqdm(to_sdf_dirs, desc="Copying SDF files"):
        subprocess.run(["cp", "-as", from_sdf_path.as_posix(),
            Path(to_inference_path, dr, to_sdf_filename).as_posix()], check=True)


def run_smina(pdb: str,
        sdf: str,
        dock_bin:str = SMINA,
        exhaustiveness:int = DEFAULT_EXHAUSTIVENESS,
        scoring:str = DEFAULT_SCORING,
        is_pdbqt:bool = False,
        seed:int = DEFAULT_SEED,
        num_cpus:int = DEFAULT_NUM_CPUS,
        autobox_ligand:Optional[str] = None,
        smina_args:Optional[list] = None,
        progress_log:Optional[str] = None) -> str:
    """Run smina or gnina and return docked SDF files as a string"""

    with NamedTemporaryFile('w', suffix='.pdbqt' if is_pdbqt else '.pdb', delete=not DEBUG) as inf_pdb, \
            NamedTemporaryFile('w', suffix='.sdf', delete=not DEBUG) as inf_sdf, \
            NamedTemporaryFile("r+", suffix='.sdf', delete=not DEBUG) as outf_sdf:

        inf_pdb.write(pdb)
        inf_pdb.flush()

        inf_sdf.write(sdf)
        inf_sdf.flush()

        cmd = (DOCK_BIN_D[dock_bin] +
            ["--receptor", inf_pdb.name, "--ligand", inf_sdf.name, "--scoring", scoring, 
            "--out", outf_sdf.name, "--seed", seed, "--cpu", num_cpus,
            "--exhaustiveness", exhaustiveness])

        cmd += ["--no_gpu"] if dock_bin == GNINA else []
        cmd += smina_args or []

        # default autobox_ligand is the pdb file + 4 Angstroms
        if autobox_ligand is not None:
            autobox_ligand_text = get_or_download_pdb(autobox_ligand)
            with NamedTemporaryFile('w', suffix='.pdb', delete=False) as inf_al:
                inf_al.write(autobox_ligand_text)
        else:
            autobox_ligand = inf_pdb.name
        cmd += ["--autobox_ligand", inf_pdb.name, "--autobox_add", 4]

        with subprocess.Popen([str(c) for c in cmd],
                bufsize=1, universal_newlines=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p:
            if progress_log is not None:
                progress_out = open(progress_log, 'w', buffering=1)
            for line in p.stdout:
                if progress_log is not None:
                    progress_out.write(line)
                    progress_out.flush()

        outf_sdf.flush()
        outf_sdf.seek(0)
        docked_sdfs = outf_sdf.read()

    return docked_sdfs


def sdf_to_coords_scores(docked_sdf: str) -> tuple:
    """Get all the coordinates and scores from a docked sdf"""
    clines = [line for line in docked_sdf.splitlines()[1:] if len(line.split()) == 16]
    coords = np.array([(float(cl.split()[0]), float(cl.split()[1]), float(cl.split()[2]))
                       for cl in clines])

    scores_d = {k: float(v)
        for k, v in dict(re.findall(r"\<(.+?)>\n([\-\d\.]+)", docked_sdf)).items()
        if k in {"minimizedAffinity", "CNNscore", "CNNaffinity", "CNN_VS"}}

    return coords, scores_d


def get_nearest_res(docked_sdf: str, pdb_obj) -> str:
    """For a given sdf and pdb object, return the nearest chain and residue"""

    pdb_coords = {}
    for chain in pdb_obj.get_chains():
        for resi in chain.get_residues():
            for atom in resi.get_atoms():
                # e.g., "A LEU 381" where A is the chain
                full_res_name = ' '.join([chain.id, resi.get_resname(), str(resi.full_id[-1][1])])
                pdb_coords.setdefault(full_res_name, []).append(list(atom.get_coord()))

    pdb_coords = {full_res_name: np.array(pdb_coords[full_res_name]) for full_res_name in pdb_coords}

    # sdf coordinate parsing
    sdf_clines = [line for line in docked_sdf.splitlines()[1:] if len(line.split()) == 16]
    sdf_coords = np.array([(float(cl.split()[0]), float(cl.split()[1]), float(cl.split()[2]))
                           for cl in sdf_clines])

    mindist_chain = (np.inf, None)
    for chain in pdb_coords:
        dists = cdist(pdb_coords[chain], sdf_coords)
        if np.min(dists) < mindist_chain[0]:
            mindist_chain = (np.min(dists), chain)

    return mindist_chain[1]


def fix_pdb(pdb: str) -> str:
    """Fix / normalize pdb, add Hs, remove HETATMs. Necessary for EquiBind"""

    with NamedTemporaryFile('w', suffix='.pdb', delete=True) as inf_pdb, \
            NamedTemporaryFile("r+", suffix='.pdb', delete=False) as outf_pdb:
        inf_pdb.write(pdb)
        inf_pdb.flush()

        subprocess.run([OBABEL_BIN, inf_pdb.name, "-h", "-O", outf_pdb.name],
            capture_output=True, check=True)

        outf_pdb.flush()
        outf_pdb.seek(0)
        fixed_pdb = ''.join([l for l in outf_pdb.readlines() if not l.startswith("HETATM")])

    assert len(fixed_pdb) > 0, f"{OBABEL_BIN} may have had a segfault"

    return fixed_pdb


def pdb_to_pdbqt(pdb: str) -> str:
    """Convert pdb to pdbqt.
    pdbqt has additional flexibility in side-chains for better docking.
    - gnina and smina may do this by default if you just supply a pdb, with flags `-h` `-xr`
    - `-h` adds hydrogens, which seems important

    See http://openbabel.org/docs/current/Features/Radicals.html#smiles-extensions-for-radicals
    We also may need makeflex.py https://sourceforge.net/p/smina/discussion/help/thread/feb3d277/
    to recover aas from pdbqt, if necessary
    """

    with NamedTemporaryFile('w', suffix='.pdb', delete=True) as inf_pdb, \
            NamedTemporaryFile("r+", suffix='.pdbqt', delete=False) as outf_pdbqt:
        inf_pdb.write(pdb)
        inf_pdb.flush()

        subprocess.run([OBABEL_BIN, inf_pdb.name, "-h", "-xr", "-O", outf_pdbqt.name],
            capture_output=True, check=True)

        outf_pdbqt.flush()
        outf_pdbqt.seek(0)
        pdbqt = outf_pdbqt.read()

    assert len(pdbqt) > 0, f"{OBABEL_BIN} may have had a segfault"

    return pdbqt


def pdb_to_coords(pdb: str, pdb_id:Optional[str]=None, min_mol_size:float = 3.0) -> dict:
    """Get coordinates for every residue in the pdb file.
    Returns a dictionary of residue : xyz
    """
    pdb_parser = PDB.PDBParser()
    pdb_obj = pdb_parser.get_structure(file=io.StringIO(pdb), id=pdb_id)

    coords_d = {}
    for chain in pdb_obj.get_chains():
        for res in chain.get_residues():
            # res.full_id: ('5HZE_hets', 0, 'A', ('H_E62', 401, ' '))
            full_res_name = f"{res.full_id[2]}_{res.full_id[-1][0]}_{res.full_id[-1][1]}"
            coords_d[full_res_name] = np.array([atom.coord for atom in res.get_atoms()])

    return {res: coords for res, coords in coords_d.items()
            if len(coords) >= min_mol_size}


def sm_id_to_sdfs(sm_id:str, max_sdf_confs:int=1, seed:int=1) -> tuple:
    """Small molecule id (pubchem id, sdf filepath, smiles)
    to two kinds of SDF files (inferred 3d and recorded 3d)
    """
    if str(sm_id).endswith(".sdf"):
        sdf_from_smiles = None
        with open(sm_id) as sdf_3d_file:
            sdf_3d = sdf_3d_file.read()
    else:
        if str(sm_id).isdigit():  # then it is a pubchem id
            sm_smiles = download_smiles(sm_id)
            _unused_sdf_2d, sdf_3d = download_sdf(sm_id)            
        else:  # it is a smiles
            sm_smiles = sm_id
            _unused_sdf_2d, sdf_3d = None, None

        sdf_from_smiles = rdconf(sm_smiles, maxconfs=max_sdf_confs, seed=seed)

        if len(sdf_from_smiles) == 0:
            sys.stderr.write(f"rdconf of sm_smiles {sm_id} {sm_smiles} failed\n")
            sdf_from_smiles = None

        sys.stderr.write(f"### SDF {sm_id}\t2d: {'success' if _unused_sdf_2d else 'failed'}\t3d: {'success' if sdf_3d else 'failed'}\n")

    return sdf_3d, sdf_from_smiles


def dock(pdb_id:str,
        sm_id:str,
        dock_bin:str = SMINA,
        exhaustiveness:int = DEFAULT_EXHAUSTIVENESS,
        seed:int = DEFAULT_SEED,
        num_cpus:int = DEFAULT_NUM_CPUS,
        scoring:str = DEFAULT_SCORING,
        is_pdbqt:bool = True,
        is_pdb_no_hetatm:bool = True,
        autobox_ligand:Optional[str] = None,
        max_sdf_confs:int = DEFAULT_MAX_SDF_CONFS,
        smina_args:Optional[list] = None,
        out_tsv:Optional[str]=None,
        progress_log:Optional[str] = None) -> pd.DataFrame:
    """smina or gnina dock from a PDB id and a (pubchem) small molecule id
    """

    sdf_3d, sdf_from_smiles = sm_id_to_sdfs(sm_id, max_sdf_confs, seed)
    if sdf_3d is None and sdf_from_smiles is None:
        raise SystemExit(f"Fatal Error: No sdf file for {sm_id}")

    #
    # parse pdb file to obj
    # normal vs no_hets: no_hets usually makes more sense
    # because if i include HETATMs then it obstructs the binding pocket
    # use sdf_from_smiles sdf_3d
    #
    pdb = get_or_download_pdb(pdb_id)
    if pdb is None:
        raise SystemExit(f"Fatal Error: No pdb file for {pdb_id}")

    #
    # Split pdb file into amino acid residues and HETATMs
    #
    pdb_hets = "\n".join(line for line in pdb.splitlines() if line.startswith("HETATM"))
    pdb_no_hets = "\n".join(line for line in pdb.splitlines() if not line.startswith("HETATM"))
    pdb_hets_coords = pdb_to_coords(pdb_hets) if pdb_hets else []

    #print(f"Remove hetatms for {pdb_id}\n"
    #      f"lines without HETATMS:   {len([l for l in pdb_no_hets.splitlines() if not l.startswith('REMARK')])}\n"
    #      f"lines with only HETATMS: {len([l for l in pdb_hets.splitlines() if not l.startswith('REMARK')])}")

    pdb_parser = PDB.PDBParser()
    if is_pdb_no_hetatm is True:
        pdb_obj = pdb_parser.get_structure(file=io.StringIO(pdb_no_hets), id=pdb_id+"_no_hets")
        pdb_or_qt = pdb_to_pdbqt(pdb_no_hets) if is_pdbqt else pdb_no_hets
    else:
        pdb_obj = pdb_parser.get_structure(file=io.StringIO(pdb), id=pdb_id)
        pdb_or_qt = pdb_to_pdbqt(pdb) if is_pdbqt else pdb

    results = []
    for sdf_type, sdf in [("rdconf", sdf_from_smiles), ("3d", sdf_3d)]:
        if sdf is None:
            continue

        docked_sdfs = run_smina(pdb_or_qt, sdf, dock_bin=dock_bin, exhaustiveness=exhaustiveness,
            seed=seed, num_cpus=num_cpus, scoring=scoring, is_pdbqt=is_pdbqt,
            autobox_ligand=autobox_ligand, smina_args=smina_args, progress_log=progress_log)

        for n, docked_sdf in enumerate([ds for ds in docked_sdfs.split('$$$$\n') if ds.strip()]):
            #if n < 3:
            #    sys.stderr.write(f"{n+1: 2d}:nearest {pdb_id} {sm_id} {get_nearest_res(docked_sdf, pdb_obj)}\n")

            sdf_coords, scores_d = sdf_to_coords_scores(docked_sdf)
            assert "minimizedAffinity" in scores_d, f"no minimizedAffinity found in {scores_d}"

            nearest_res = get_nearest_res(docked_sdf, pdb_obj)
            distance_to_hetatms = {res: cdist(sdf_coords, pdb_hets_coords[res]).min(axis=1).mean()
                                   for res in pdb_hets_coords}

            mean_coord = [round(x, 3) for x in sdf_coords.mean(axis=0)]

            # if there is no ligand, i still want to print out results
            if len(distance_to_hetatms) == 0:
                distance_to_hetatms = {'': -1}

            for hetatm, dist in distance_to_hetatms.items():
                results.append({
                    "pdb_id": pdb_id,
                    "sm_id": sm_id,
                    **{k:round(v, 4) for k, v in sorted(scores_d.items())},
                    "ligand": hetatm,
                    "ligand_distance": round(dist, 4),
                    "mean_coord": mean_coord,
                    "dock_bin": dock_bin,
                    "scoring": scoring,
                    "smina_args": ' '.join(smina_args) if smina_args else '',
                    "is_pdbqt": is_pdbqt,
                    "is_pdb_no_hetatm": is_pdb_no_hetatm,
                    "seed": seed,
                    "exhaustiveness": exhaustiveness,
                    "autobox_ligand": autobox_ligand,
                    "nearest_res": nearest_res,
                    "sdf_type": sdf_type,
                    "max_sdf_confs": max_sdf_confs,
                    "docked_sdf": docked_sdf
                })

    if len(results) == 0:
        raise AssertionError("No results from dock")

    # fix up the results a bit and sort
    df_res = (pd.DataFrame(results)
        .assign(docked_sdf=lambda df: df['docked_sdf'].apply(lambda r: list(r.splitlines())))
        .sort_values('minimizedAffinity'))

    if out_tsv is not None:
        df_res = df_res[[c for c in df_res.columns if c != 'docked_sdf'] + ['docked_sdf']]
        df_res.to_csv(out_tsv, sep='\t', index=None)

    return df_res


def dock_all_by_all(pdb_ids_dict:dict,
        sm_ids_dict:dict,
        dock_bin:str = SMINA, 
        exhaustiveness:int = DEFAULT_EXHAUSTIVENESS,
        scoring:str = DEFAULT_SCORING,
        is_pdbqt:bool = True,
        is_pdb_no_hetatm:bool = True,
        autobox_ligand:Optional[str] = None,
        max_sdf_confs:int = DEFAULT_MAX_SDF_CONFS,
        smina_args:Optional[list]=None,
        top_n_in_sdf:int = 10,
        out_tsv:str = "out_dock.tsv") -> pd.DataFrame:
    """Dock a collection of small molecules against a collection of pdbs.
    """

    import dask

    lazy_results = []
    for pdb_name, pdb_id in sorted(pdb_ids_dict.items()):
        for sm_name, sm_id in sorted(sm_ids_dict.items()):
            print(f"\N{dog} docking {pdb_name} ({pdb_id}) with ({sm_name}) {sm_id}")
            lazy_results.append(dask.delayed(dock)(pdb_id, sm_id, dock_bin=dock_bin,
                exhaustiveness=exhaustiveness, is_pdbqt=is_pdbqt, is_pdb_no_hetatm=is_pdb_no_hetatm,
                scoring=scoring, autobox_ligand=autobox_ligand, max_sdf_confs=max_sdf_confs,
                smina_args=smina_args))

    sm_names = {v: k for k, v in sm_ids_dict.items()}
    pdb_names = {v: k for k, v in pdb_ids_dict.items()}

    df_docking_res = [item for item in dask.compute(*lazy_results, num_workers=2)]

    #
    # Combine all docking results for this run
    #
    df_res = (pd.concat(df_docking_res)
        .assign(sm_name=lambda df: df['sm_id'].apply(sm_names.get),
                pdb_name=lambda df: df['pdb_id'].apply(pdb_names.get))
        .sort_values('minimizedAffinity'))

    print("df:", df_res.head(1))

    #
    # Combine all docking results from all time
    #
    if os.path.exists(out_tsv):
        all_df_res = (pd.concat([pd.read_csv(out_tsv, sep='\t'), df_res])
            .sort_values('minimizedAffinity'))
    else:
        all_df_res = df_res

    # put docked_sdf at the end since it's a lot of text
    all_df_res = all_df_res[[c for c in all_df_res.columns if c != 'docked_sdf'] + ['docked_sdf']]
    open(out_tsv, 'w').write(all_df_res.to_csv(sep='\t', index=None))

    #
    # Output sdf file containing docked molecules
    #
    for (pdb_id, sm_id), grouped_rows in df_res.groupby(['pdb_id', 'sm_id']):
        out_sdf = f"out_{pdb_id.split('/')[-1]}_{pdb_names[pdb_id].split('/')[-1]}_{sm_names[sm_id].split('/')[-1]}.sdf".replace(' ', '_')
        with open(out_sdf, 'w') as out:
            for _n, row in grouped_rows.iterrows():
                out.write('\n'.join(row.docked_sdf))
                out.write('$$$$\n') # SDF footer
                if _n > top_n_in_sdf:
                    break

    return df_res


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_id', help="pdb id or pdb file path", nargs='?')
    parser.add_argument('sm_id', help="pubchem id or smiles or sdf file", nargs='?')
    parser.add_argument('--out_tsv', help="tsv results file", nargs='?')
    parser.add_argument('--max_sdf_confs', help="try multiple conformers", nargs='?', default=DEFAULT_SDF_CONFS, type=int)
    # smina arguments here
    parser.add_argument('--exhaustiveness', help="usually from 8 (medium) to 64 (high quality)", nargs='?', default=DEFAULT_EXHAUSTIVENESS, type=int)
    parser.add_argument('--cpu', help="number of cpus", nargs='?', default=DEFAULT_NUM_CPUS)
    args, smina_args_global = parser.parse_known_args()

    if args.pdb_id is not None and args.sm_id is not None:
        out = dock(args.pdb_id, args.sm_id,
            num_cpus=args.cpu,
            exhaustiveness=args.exhaustiveness,
            max_sdf_confs=args.max_sdf_confs,
            smina_args=smina_args_global,
            out_tsv=args.out_tsv or f"out_{args.pdb_id}_{args.sm_id}.tsv")
    else:
        sm_ids_dict_test = {
            # 'zearalenone': 5281576,
            # 'sulfozearalenone': "CC1CCCC(=O)CCC/C=C\c2cc(OS(=O)(=O)O)cc(O)c2C(=O)O1",
            # 'atpenin a5': 54676868,
            # 'atpenin b': 54676869,
            # 'mycophenolic acid': 446541,
            # 'destruxin e': 107863,
            # 'bafilomycin a1': 6436223,
            # 'alpha-trans zearalenol': 5284645,
            # 'hypothemycin': 9929643,
            # '5z-7-oxozeaenol': 9863776,
            # 'E6201': 10172827,
            "aspterric acid": 135143071,
            #"dihydroxyisovalerate": 21738659,
        }

        pdb_ids_dict_test = {
            # 'yeast vATPase state 1': '3J9T',
            # 'yeast vATPase state 2': '3J9U',
            # 'yeast vATPase state 3': '3J9V',
            # 'mammalian vATPase state 2': '6VQA',
            # 'human AURKC': '6GRA',
            # 'IMPDH MPA': '1JR1',
            # 'SARS-CoV-2 spike bound to ACE2': '6M0J'
            # 'MAP2K1 E6201': '5HZE',
            # 'porcine SDHC': '3AEE',
            # 'Ecoli SDHC': '2ACZ',
            # 'Avian SDHC': '6MYS',
            # 'Avian SDHC 2': '6MYT',
            # 'Avian SDHC 2 no AT5': '6MYT_no_AT5', # manually edited pdb!
            # 'Nematode SDHC': '3VRA',
            # 'superfolder GFP': '5B61',
            "ILVD_DHAD from TB": "6OVT",
        }

        print(f"testing code with dock_all_by_all: {pdb_ids_dict_test} {sm_ids_dict_test}")
        dock_all_by_all(pdb_ids_dict_test,
            sm_ids_dict_test,
            dock_bin=SMINA,
            is_pdbqt=True,
            is_pdb_no_hetatm=True,
            exhaustiveness=1,
            autobox_ligand=None,
            max_sdf_confs=args.max_sdf_confs,
            smina_args=smina_args_global,
            out_tsv="out_dock.test.tsv")
