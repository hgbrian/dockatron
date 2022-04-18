"""Run EquiBind as a subprocess.
"""
import os
import subprocess

from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory, NamedTemporaryFile, gettempdir

from tqdm.auto import tqdm
from tqdm.dask import TqdmCallback

import pandas as pd
import dask
import smina_dock

tempdir = gettempdir()
NUM_WORKERS = 100

# taken from EquiBind's default inference.yml
# cannot easily be replaced with command line args for some reason (--config is required)
EB_CONFIG_YAML = """run_dirs:
  - flexible_self_docking # the resulting coordinates will be saved here as tensors in a .pt file (but also as .sdf files if you specify an "output_directory" below)
inference_path: '{input_dir}' # this should be your input file path as described in the main readme

output_directory: '{output_dir}' # the predicted ligands will be saved as .sdf file here
run_corrections: True
use_rdkit_coords: False # generates the coordinates of the ligand with rdkit instead of using the provided conformer. If you already have a 3D structure that you want to use as initial conformer, then leave this as False
save_trajectories: False
num_confs: 1 # usually this should be 1
"""


@contextmanager
def using_directory(path: str):
    """change directory context"""
    origin = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)

def smina_score(pdb_file: str, sdf_file: str, out_tsv:str = None, score_type="minimize"):
    """Use smina to score a docked position"""
    if score_type == "minimize":
        smina_args = ["--local_only", "--minimize"]
    elif score_type == "score_only":
        smina_args = ["--score_only"]
    else:
        raise ValueError(f"score_type is {score_type}")

    df_sm = smina_dock.dock(pdb_file, sdf_file, max_sdf_confs=1,
        smina_args=smina_args, out_tsv=out_tsv)

    return df_sm


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


def run_equibind(proteome_dir: str, sm_id: str, output_dir=None, out_smina_tsv=None,
        friendly_sm_id=None, friendly_proteome_id=None):
    """Run equibind against a proteome"""
    proteome_dir = Path(proteome_dir).resolve()

    sdf_3d, sdf_from_smiles = smina_dock.sm_id_to_sdfs(sm_id, max_sdf_confs=1, seed=1)
    sdf = sdf_3d if sdf_3d is not None else sdf_from_smiles

    if not friendly_sm_id:
        friendly_sm_id = sm_id.split("/")[-1].split('.')[0]

    if not friendly_proteome_id:
        friendly_proteome_id = proteome_dir.stem

    today = datetime.now().isoformat()[:10].replace("-","")
    if not output_dir:
        output_dir = Path(tempdir, f"{today}_{friendly_proteome_id}_{friendly_sm_id}")

    # ----------------------------------------------------------------------------------------------
    # Start docking
    #
    scores = []
    df_uniprot = pd.read_csv(f"uniprot_{friendly_proteome_id}.tsv", sep='\t')[["Entry", "Gene names", "Entry name"]]

    with NamedTemporaryFile('w', suffix=".yml") as config_file, \
        TemporaryDirectory() as sdf_dir, \
        TemporaryDirectory() as input_dir:

        sdf_file = open(Path(sdf_dir, f"{friendly_sm_id}_ligand.sdf"), 'w')
        sdf_file.write(sdf)
        sdf_file.flush()

        config_file.write(EB_CONFIG_YAML.format(input_dir=input_dir, output_dir=output_dir))
        config_file.flush()

        # ------------------------------------------------------------------------------------------
        # Linking files
        #
        link_proteome_files(proteome_dir, input_dir, sdf_file.name)

        with using_directory("EquiBind"):
            with subprocess.Popen(["conda", "run", "-n", "equibind", "--no-capture-output",
                "python", "-u", "inference.py", f"--config={config_file.name}" ],
                bufsize=1, universal_newlines=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p:

                for line in tqdm((l for l in iter(p.stdout) if l.startswith("Processing")),
                    desc="EquiBind", total=len(os.listdir(input_dir))):
                    pass

            # Run smina in parallel to get affinity
            pjobs = []
            for pid in [d.name for d in Path(output_dir).iterdir() if Path(output_dir, d).is_dir()]:
                pdb_file = list(Path(proteome_dir, pid).glob("*.pdb"))
                assert len(pdb_file) == 1, f"proteome directory error: {proteome_dir}/{pid}"
                pdb_file = pdb_file[0]

                scores.append((pid, sm_id, friendly_sm_id))

                pjobs.append(dask.delayed(smina_score)(
                    pdb_file = pdb_file.as_posix(), #Path(proteome_dir, pid, pdb[0]),
                    sdf_file = Path(output_dir, pid, "lig_equibind_corrected.sdf").as_posix(),
                    out_tsv = out_smina_tsv,
                    score_type = "minimize"
                    )
                )

            with TqdmCallback(desc="smina"):
                df_res = pd.concat(dask.compute(*pjobs, num_workers=NUM_WORKERS))

            print(df_res)
            #scores = [list(s) + [pj] for s, pj in zip(scores, pjobs_res)]

    # ----------------------------------------------------------------------------------------------
    # Output
    #
    out_df_file = Path(output_dir, f"{today}_{friendly_proteome_id}_{friendly_sm_id}_equibind_smina.tsv")
    df_scores = (df_res #(pd.DataFrame(scores, columns=["pid", "sm_id", "friendly_sm_id", "smina_kcal_mol"])
        .assign(pdb_id = lambda df: df.pdb_id.apply(lambda o: o.split('/')[-1]))
        .assign(sm_id = lambda df: df.sm_id.apply(lambda o: o.split('/')[-1]))
        .assign(dock_bin = lambda df: "EquiBind + smina")
        .sort_values("minimizedAffinity")
        # how to get this back?
        #.merge(df_uniprot.rename(columns={"Entry":"pid"}), on="pid")
        #.assign(smina_kcal_mol = lambda df: df.smina_kcal_mol.astype(float))
        #.sort_values("smina_kcal_mol")
        .to_csv(out_df_file, sep='\t', index=None))

    print(f"Wrote affinity results to {out_df_file}")
    return df_scores


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("proteome_dir")
    parser.add_argument("sm_id")
    parser.add_argument("--out_tsv", nargs="?")
    parser.add_argument("--output_dir", nargs="?")
    parser.add_argument("--friendly_proteome_id", nargs="?")
    parser.add_argument("--friendly_sm_id", nargs="?")
    args = parser.parse_args()

    run_equibind(args.proteome_dir, args.sm_id, args.output_dir,
        friendly_sm_id=args.friendly_sm_id,
        friendly_proteome_id=args.friendly_proteome_id)
