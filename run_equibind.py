import os, sys, re
import dask
import subprocess
import requests
import pandas as pd
from tqdm.auto import tqdm
from tqdm.dask import TqdmCallback
from contextlib import contextmanager
from datetime import datetime
from tempfile import TemporaryDirectory, NamedTemporaryFile, gettempdir

tempdir = gettempdir()

# taken from EquiBind's default inference.yml
# could be replaced with command line args
config_yaml = """run_dirs:
  - flexible_self_docking # the resulting coordinates will be saved here as tensors in a .pt file (but also as .sdf files if you specify an "output_directory" below)
inference_path: '{input_dir}' # this should be your input file path as described in the main readme

test_names: timesplit_test
output_directory: '{output_dir}' # the predicted ligands will be saved as .sdf file here
run_corrections: True
use_rdkit_coords: False # generates the coordinates of the ligand with rdkit instead of using the provided conformer. If you already have a 3D structure that you want to use as initial conformer, then leave this as False
save_trajectories: False

num_confs: 1 # usually this should be 1
"""

@contextmanager
def using_directory(path: str):
    origin = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)

def download_sdf(sm_id: str):
    """download sdf_2d and sdf_3d from pubchem with pubchem id"""
    sm_pubchem_url_sdf = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{sm_id}/record/SDF/"

    sdf_2d = requests.get(f"{sm_pubchem_url_sdf}?record_type=2d")
    sdf_3d = requests.get(f"{sm_pubchem_url_sdf}?record_type=3d")
    sdf_2d = sdf_2d.text if sdf_2d.status_code == 200 else None
    sdf_3d = sdf_3d.text if sdf_3d.status_code == 200 else None

    return sdf_2d, sdf_3d

def smina_score(pdb_file: str, sdf_file: str, scoring="minimize"):
    """Use smina to score a docked position"""
    if scoring == "minimize":
        opt = ["--local_only", "--minimize"]
    elif scoring == "score_only":
        opt = ["--score_only"]
    else:
        raise ValueError

    cmd = ["python", f"{os.environ['HXROOT']}/apps/docking/smina_dock.py",
        pdb_file, sdf_file,
        *opt, "--out_tsv=/dev/null"]

    out = subprocess.run(cmd, capture_output=True).stdout.decode()
    rs = re.search("Affinity:\s+(\S+)", out, re.MULTILINE)
    if rs:
        return rs.group(1)
    else:
        return None


def run_equibind(yeast_or_human: str, sm_id: str, output_dir=None, friendly_id=None):
    """Run equibind against a proteome"""
    proteome_dir = os.path.abspath(f"{yeast_or_human}_proteome")

    if "." in sm_id:
        sdf = open(sm_id).read()
    else:
        sdf_2d, sdf_3d = download_sdf(sm_id)
        sdf = sdf_3d if sdf_3d else sdf_2d

    if not friendly_id:
        friendly_id = sm_id.split("/")[-1].split('.')[0]

    today = datetime.now().isoformat()[:10].replace("-","").replace("20220326", "20220325")
    if not output_dir:
        output_dir = os.path.join(tempdir, f"{today}_{yeast_or_human}_{friendly_id}")

    scores = []
    df_uniprot = pd.read_csv(f"uniprot_{yeast_or_human}.tsv", sep='\t')[["Entry", "Gene names", "Entry name"]]

    with NamedTemporaryFile('w', suffix=".yml") as config_file, \
        TemporaryDirectory() as sdf_dir, \
        TemporaryDirectory() as input_dir:

        sdf_file = open(os.path.join(sdf_dir, f"{friendly_id}_ligand.sdf"), 'w')
        sdf_file.write(sdf)
        sdf_file.flush()

        config_file.write(config_yaml.format(input_dir=input_dir, output_dir=output_dir))
        config_file.flush()
        
        print("Linking files")
        subprocess.run(f"cp -as {os.path.abspath(proteome_dir).rstrip('/')+'/*'} {input_dir}", shell=True)
        for dr in tqdm(os.listdir(input_dir)):
            subprocess.run(["cp", "-as", os.path.abspath(sdf_file.name), os.path.join(input_dir, dr)])

        print("Running EquiBind")
        with using_directory("EquiBind"):
            with subprocess.Popen(["conda", "run", "-n", "equibind",
                "python", "-u", "inference.py", f"--config={config_file.name}" ],
                bufsize=1, universal_newlines=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p:

                # TODO fix this something went wrong
                #for line in tqdm((l for l in p.stdout if l.decode().startswith("Processing")), total=len(df_uniprot)):
                #    pass
            
                print(f"wrote equibind results to {output_dir}")

            pjobs = []
            for pid in [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]:

                pdb = [f for f in os.listdir(os.path.join(proteome_dir, pid)) if f.endswith(".pdb")]
                scores.append((pid, sm_id, friendly_id))
                pjobs.append(dask.delayed(smina_score)(os.path.join(proteome_dir, pid, pdb[0]),
                    os.path.join(output_dir, pid, "lig_equibind_corrected.sdf")))

            with TqdmCallback():
                pjobs_res = dask.compute(*pjobs, num_workers=100)
            scores = [list(s) + [pj] for s, pj in zip(scores, pjobs_res)]

    out_df_file = os.path.join(output_dir, f"{today}_{yeast_or_human}_{friendly_id}_smina_affinities.tsv")
    df_scores = (pd.DataFrame(scores, columns=["pid", "sm_id", "friendly_id", "smina_kcal_mol"])
        .merge(df_uniprot.rename(columns={"Entry":"pid"}), on="pid")
        .assign(smina_kcal_mol = lambda df: df.smina_kcal_mol.astype(float))
        .sort_values("smina_kcal_mol")
        .to_csv(out_df_file, sep='\t', index=None))

    print(f"Wrote affinity results to {out_df_file}")
    return df_scores


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("yeast_or_human", choices=["yeast", "human", "test"])
    parser.add_argument("sm_id")
    parser.add_argument("--output_dir", nargs="?")
    parser.add_argument("--friendly_id", nargs="?")
    args = parser.parse_args()

    run_equibind(args.yeast_or_human, args.sm_id, args.output_dir, args.friendly_id)
