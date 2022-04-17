"""
# Dockatron
Dockatron relies on the following tools:

## Cite smina
- http://pubs.acs.org/doi/abs/10.1021/ci300604z

## Cite EquiBind
- https://arxiv.org/abs/2202.05146

## Cite AlphaFold:
- Jumper, J et al. Highly accurate protein structure prediction with AlphaFold. Nature (2021).
- Varadi, M et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Research (2021).
"""
import io
import os
import platform
import sys

from contextlib import contextmanager, redirect_stdout
from multiprocessing import Process
from pathlib import Path
from tempfile import NamedTemporaryFile, mkdtemp, gettempdir
from time import sleep
from typing import List

import requests
import yaml

import rich
from rich import box
from rich.align import Align
from rich.console import RenderableType
from rich.markdown import Markdown
from rich.panel import Panel
from rich.pretty import Pretty
from rich.progress import Progress, track, SpinnerColumn
from rich.progress_bar import ProgressBar
from rich.screen import Screen
from rich.status import Status
from rich.style import Style, StyleType
from rich.table import Table
from rich.text import Text, TextType
from textual import events, log
from textual.app import App
from textual.keys import Keys
from textual.message import Message
from textual.message_pump import MessagePump
from textual.view import View
from textual.widget import Reactive, Widget
from textual.widgets import Button, ButtonPressed, ScrollView, Static, TreeControl, TreeClick, TreeNode

import smina_dock
import run_equibind

# Assume EquiBind is downloaded to a folder under dockatron
EBHOME = Path("./EquiBind").resolve()
sys.path.insert(1, EBHOME.as_posix())
from EquiBind import inference

# --------------------------------------------------------------------------------------------------
# Globals
#
if platform.uname().system == "Linux":
    MAC_OR_LINUX = "linux"
elif platform.uname().system == "Darwin":
    MAC_OR_LINUX = "osx" if platform.mac_ver()[0][:2]<"12" else "osx12"
else:
    raise ValueError(f"{platform.uname()} not supported")

URLS = {
    ("equibind", "linux"): "https://github.com/HannesStark/EquiBind/archive/refs/heads/main.zip",
    ("gnina", "linux"): "https://github.com/gnina/gnina/releases/download/v1.0/gnina",
    ("smina", "linux"): "https://sourceforge.net/projects/smina/files/smina.static/download",
    ("smina", "osx"): "https://sourceforge.net/projects/smina/files/smina.osx/download",
    ("smina", "osx12"): "https://sourceforge.net/projects/smina/files/smina.osx.12/download",
}

PROTEOMES = {
    "E. coli": "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000625_83333_ECOLI_v2.tar",
    "D. melanogaster": "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000803_7227_DROME_v2.tar",
    "H. sapiens": "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v2.tar",
    "S. cerevisiae": "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v2.tar",
    "A. thaliana": "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000006548_3702_ARATH_v2.tar",
    "C. elegans": "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000001940_6239_CAEEL_v2.tar",
}

# some urls don't return a content length, so put something
DEFAULT_CONTENT_LENGTHS = {
    URLS[("equibind", "linux")]: int(4e8),
}

DEBUG = True

TMPDIR = gettempdir()
DATADIR = Path(Path.cwd(), "dockatron_files")
DATADIR.mkdir(parents=True, exist_ok=True)
OUTDIR = Path(Path.cwd(), "dockatron_files", "results")
OUTDIR.mkdir(parents=True, exist_ok=True)


# with DEBUG, all the fast options are selected
DEFAULT_MAX_SDF_CONFS = 1 if DEBUG else 4
DEFAULT_EXHAUSTIVENESS = 1 if DEBUG else 8

# FIXFIX
DEFAULT_MAX_SDF_CONFS = 10
DEFAULT_EXHAUSTIVENESS = 1

# --------------------------------------------------------------------------------------------------
# Calling other libraries
#

def gen_dl(url, chunk_size=1_048_576, out_dir=None, out_file=None):
    """Generator for downloading files that lets me show progress."""
    if out_dir is not None:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
    else:
        out_dir = Path(".").resolve()

    resp = requests.get(url, stream=True)
    total_dl = int(resp.headers.get('content-length', DEFAULT_CONTENT_LENGTHS.get(url, 1e8)))
    yield total_dl // chunk_size

    with open(Path(out_dir, out_file) if out_file
            else Path(out_dir, url.split('/')[-1]), 'wb') as out:
        for data in resp.iter_content(chunk_size=chunk_size):
            out.write(data)
            yield

@contextmanager
def using_directory(path: str):
    """use a working directory"""
    origin = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


def dock_equibind(config_dict):
    """Code taken directly from EquiBind/inference.py __main__"""
    # EquiBind expects a file, no matter what?
    placeholder_yml = Path(TMPDIR, "EquiBind_placeholder.yml")
    if not placeholder_yml.exists():
        placeholder_yml.touch()
    sys.argv.extend(["--config", placeholder_yml.as_posix()])

    # using_directory(EBHOME) maybe ??????
    with io.StringIO() as f: # redirecting stdout to f
        args = inference.parse_arguments()

        # replace yaml with my own config_dict object
        #config_dict = yaml.load(args.config, Loader=yaml.FullLoader)
        args.config.close()

        arg_dict = args.__dict__
        for key, value in config_dict.items():
            if isinstance(value, list):
                for v in value:
                    arg_dict[key].append(v)
            else:
                arg_dict[key] = value
        args.config = args.config.name

        with redirect_stdout(f):
            for run_dir in args.run_dirs:
                args.checkpoint = f'{EBHOME}/runs/{run_dir}/best_checkpoint.pt'
                config_dict['checkpoint'] = f'{EBHOME}/runs/{run_dir}/best_checkpoint.pt'
                # overwrite args with args from checkpoint except for the args that were contained in the config file
                arg_dict = args.__dict__
                with open(os.path.join(os.path.dirname(args.checkpoint), 'train_arguments.yaml'), 'r') as arg_file:
                    checkpoint_dict = yaml.load(arg_file, Loader=yaml.FullLoader)

                for key, value in checkpoint_dict.items():
                    if key not in config_dict.keys():
                        if isinstance(value, list):
                            for v in value:
                                arg_dict[key].append(v)
                        else:
                            arg_dict[key] = value
                args.model_parameters['noise_initial'] = 0
                if args.inference_path is None:
                    inference.inference(args)
                else:
                    inference.inference_from_files(args)

        return


def gen_dock_equibind(inference_path:str, output_directory:str, timeout_s:int=100_000):
    """Run EquiBind in a Process to get progress"""

    # config_dict comes straight from the default EquiBind yaml
    config_dict = dict(
        output_directory=output_directory,
        run_dirs=['flexible_self_docking'],
        inference_path=inference_path,
        run_corrections=True,
        use_rdkit_coords=False,
        save_trajectories=False,
        num_confs = 1
    )

    total_pdbs = len(list(Path(inference_path).glob("*/*.pdb")))

    p = Process(target=dock_equibind, args=(config_dict,))
    p.start()

    sleep_s = 10
    for _ in range(timeout_s // sleep_s):
        num_done = len(list(Path(output_directory).glob("*/lig_equibind_corrected.sdf")))
        yield (num_done, total_pdbs)
        sleep(sleep_s)
        if not p.is_alive():
            break
    else:
        raise Exception(f"EquiBind failed to finish before {timeout_s}")

    return

def test_equibind_gen():
    """Run EquiBind in a Process to get progress"""

    eb_iter = gen_dock_equibind(
        inference_path=f"{EBHOME}/../mycophenolic_acid_test",
        output_directory=f"{EBHOME}/data/results/output/mycophenolic_acid_yeast"
    )

    total_prots = next(eb_iter)
    for p in eb_iter:
        print(p, total_prots)

    raise SystemExit("finished testing equibind as a module")


# def gen_dock_smina_subprocess(pdb_id, sm_id, out_tsv, exhaustiveness=1, max_sdf_confs=1):
#     with subprocess.Popen(["python", Path(os.path.abspath(Path.cwd()), "smina_dock.py"), pdb_id, sm_id,
#         "--exhaustiveness", str(exhaustiveness),
#         "--max_sdf_confs", str(max_sdf_confs),
#         "--out_tsv", out_tsv],
#         bufsize=1, universal_newlines=True,
#         stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p:

#         for line in (l for l in iter(p.stdout) if "Refine time" in l):
#             yield

#         # this would output stderr to output
#         for line in (l for l in iter(p.stderr)):
#             yield line


def gen_dock_smina(pdb_id:str, sm_id:str, out_tsv:str, exhaustiveness:int=DEFAULT_EXHAUSTIVENESS,
        max_sdf_confs:int=DEFAULT_MAX_SDF_CONFS, timeout_s:int=100_000):
    """Run smina in a Process to get progress"""
    print("MAX_SDF_CONFS", max_sdf_confs)

    progress_log = NamedTemporaryFile(suffix=".log", delete=False)
    print(progress_log.name)
    sleep_s = 10

    total_sdfs = max_sdf_confs + 1

    with redirect_stdout(io.StringIO()) as f:
        p = Process(target=smina_dock.dock, kwargs=(dict(pdb_id=pdb_id, sm_id=sm_id,
            exhaustiveness=exhaustiveness, max_sdf_confs=max_sdf_confs,
            out_tsv=out_tsv, progress_log=progress_log.name)))
        p.start()

    for _ in range(timeout_s // sleep_s):
        if Path(progress_log.name).exists():
            with open(progress_log.name, 'r') as _f:
                print(_f)
                f_read = _f.read()
                print(len(f_read))
                count_iters = f_read.count("Refine")
                yield (count_iters, total_sdfs)
        sleep(sleep_s)
        if not p.is_alive():
            break
    else:
        raise Exception(f"smina failed to finish before {timeout_s}")

    return

# def gen_dock_equibind_subprocess(pdb_id, sm_id, out_tsv):
#     with subprocess.Popen(["python", "run_equibind.py", "test", "123456"],
#         bufsize=1, universal_newlines=True,
#         stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p:
#         for line in (l for l in iter(p.stdout) if "Refine time" in l):
#             yield

def prep_equibind_dir(pdb_or_proteome_id:str, sm_id:str):
    """Place PDB and SDF files in a directory for EquiBind to process"""

    inference_path = mkdtemp(prefix="equibind_")
    sdf_3d, sdf_from_smiles = smina_dock.sm_id_to_sdfs(sm_id)
    sdf = sdf_3d if sdf_3d is not None else sdf_from_smiles

    with NamedTemporaryFile('w', suffix=".sdf", delete=False    ) as sdf_file:
        sdf_file.write(sdf)
        sdf_file.flush()
        # either link from a pdb or a proteome
        proteome_dir = pdb_or_proteome_id
        print(proteome_dir, inference_path, sdf_file.name)
        run_equibind.link_proteome_files(proteome_dir, inference_path, sdf_file.name)
        print("Done linking!")

    return inference_path

def run_docking(pdb_id, sm_id, out_tsv, docking="smina"):
    """run docking with smina or equibind"""
    if docking=="smina":
        yield from gen_dock_smina(pdb_id, sm_id, out_tsv=out_tsv)
    elif docking=="equibind":
        inference_path = prep_equibind_dir("test_proteome", 123456)
        from datetime import datetime
        today = datetime.now().isoformat()[:10].replace("-","")
        output_directory = Path(gettempdir(), f"{today}_test_proteome_123456")

        yield from gen_dock_equibind(inference_path, output_directory.as_posix())

gen_dock = run_docking("6GRA", 123456, "/tmp/6GRA_123456_temp.tsv", docking="equibind")
for it in gen_dock:
    print("it", it)
print(it)
1/0

# --------------------------------------------------------------------------------------------------
# Textual stuff
#

@rich.repr.auto(angular=False)
class TextInputPanel(Widget, can_focus=True):

    has_focus: Reactive[bool] = Reactive(False)
    mouse_over: Reactive[bool] = Reactive(False)
    keypress: Reactive[str] = Reactive("")
    style: Reactive[str] = Reactive("")
    height: Reactive[int] = Reactive(None)

    def __init__(self, *, name: str = None, height: int = None, val:str = '', row:int = None, cols:List = None) -> None:
        super().__init__(name=name)
        self.height = height
        self.text = ""
        self.title = name
        self.val = val
        self.row = row
        self.cols = cols

    def __rich_repr__(self) -> rich.repr.Result:
        yield self.name

    def render(self) -> RenderableType:
        return Panel(
            Align.left(self.text),
            border_style="green" if self.mouse_over else "green",
            box=box.HEAVY if self.has_focus else box.ROUNDED,
            title=self.title,
            style=self.style,
        )

    async def on_focus(self, event: events.Focus) -> None:
        self.has_focus = True

    async def on_blur(self, event: events.Blur) -> None:
         self.has_focus = False

    async def on_key(self, event: events.Key):
        if event.key == Keys.Enter:
            await self.emit(ButtonPressed(self))
        elif event.key == Keys.ControlH:
            self.text = self.text[:-1]
            self.parent.parent.parent.vals[self.val] = self.text
            self.refresh()
        elif self.has_focus and len(event.key)==1:
            self.text = self.text + event.key
            self.parent.parent.parent.vals[self.val] = self.text
            self.refresh()

    async def on_click(self, event: events.Click):
        await self.emit(ButtonPressed(self))

@rich.repr.auto(angular=False)
class TextPanel(Widget, can_focus=True):

    has_focus: Reactive[bool] = Reactive(False)
    mouse_over: Reactive[bool] = Reactive(False)
    keypress: Reactive[str] = Reactive("")
    style: Reactive[str] = Reactive("")
    height: Reactive[int] = Reactive(None)

    def __init__(self, *, name: str = None, height: int = None, val:str = '', row:int = None, cols:List = None) -> None:
        super().__init__(name=name)
        self.height = height
        self.text = ""
        self.title = name
        self.val = val
        self.row = row
        self.cols = cols

    def __rich_repr__(self) -> rich.repr.Result:
        yield self.name

    def render(self) -> RenderableType:
        return Panel(
            Align.left(self.text),
            border_style="green" if self.mouse_over else "green",
            box=box.HEAVY if self.has_focus else box.ROUNDED,
            title=self.title,
            style=self.style,
        )

    async def on_focus(self, event: events.Focus) -> None:
        self.has_focus = True

    async def on_blur(self, event: events.Blur) -> None:
         self.has_focus = False

    async def on_key(self, event: events.Key):
        if event.key == Keys.Enter:
            await self.emit(ButtonPressed(self))

    async def on_click(self, event: events.Click):
        await self.emit(ButtonPressed(self))

@rich.repr.auto(angular=False)
class GridScrollView(ScrollView):
    def __init__(self, *, name: str = None, val:str = "", row:int = None, cols:List = None) -> None:
        super().__init__(name=name)
        self.name = name
        self.row = row
        self.cols = cols
        self.update(Markdown('\n\n'.join("self.output_md"), hyperlinks=True))


@rich.repr.auto(angular=False)
class GridButton(Button):
    def __init__(self, label:str, *, name: str = None, row:int = None, cols:List = None) -> None:
        super().__init__(label)
        self.row = row
        self.cols = cols

    # ButtonPressed does not work
    async def on_focus(self, event: events.Focus) -> None:
        self.has_focus = True
        self.label = "[ " + self.label + " ]"
        self.button_style = "white on blue"
        self.refresh()

    async def on_blur(self, event: events.Blur) -> None:
        self.has_focus = False
        self.label = self.label.lstrip("[").rstrip("]").strip()
        self.button_style = "white on dark_blue"
        self.refresh()

    # is this still necessary?
    async def on_key(self, event: events.Key):
        log(f"Button on_key {self}")
        if event.key == Keys.Enter:
            await self.emit(ButtonPressed(self))

    async def on_click(self, event: events.Click) -> None:
        log(f"XXX Button on_click {self}")

    async def on_mouse_move(self, event: events.MouseMove) -> None:
        log(f"XXX Button move {self}")

    #    #self.has_focus = True
    #    self.start_docking()
    #    self.label = "Docking..."
    #    self.button_style = "white on dark_green"

    #def start_docking(self):
    #    self.post_message_from_child_no_wait(Message(self))

    # def handle_button_pressed(self, message: ButtonPressed) -> None:
    #     print("Button pressed")
    #     self.start_docking()
    #     self.label = "Docking..."
    #     self.button_style = "white on dark_green"



class GridTest(App):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.grid = None
        self.row = 0
        self.col = 1
        self.row_proteome = 0
        self.output_md = ["# Dockatron Output"]
        self.showing_proteome_list = False
        self.vals = {}

    async def on_load(self) -> None:
        """Bind keys here."""
        await self.bind(Keys.Up, "move_up", "move up")
        await self.bind(Keys.Down, "move_down", "move down")
        await self.bind(Keys.Left, "move_left", "move left")
        await self.bind(Keys.Right, "move_right", "move right")

    async def _update_output(self):
        await self.output.update(Markdown('\n\n'.join(self.output_md), hyperlinks=True))
        self.refresh() # ???

        # another janky timer!
        def _scroll():
            self.output.target_y = self.output.max_scroll_y
            self.output.y = self.output.target_y
        self.set_timer(0.1, _scroll)

    def _get_friendly_id(self):
        def _smiles_to_id(smiles):
            return "".join(A+str(smiles.count(A)) if smiles.count(A) else "" for A in "CHONPS")

        prot = self.vals.get('pdb_id') or self.vals.get('gene_name') or self.vals.get('proteome')
        sm = self.vals.get('pubchem_id') or self.vals.get('sdf') or _smiles_to_id(self.vals.get('smiles'))

        if prot and sm:
            return f"{prot}_{sm}"
        else:
            return None

    async def _start_docking(self, message_sender):
        if not (self.vals.get("proteome") or self.vals.get("pdb_id") or self.vals.get("gene_name")):
            self.output_md.append("No proteome or pdb or gene name supplied.")
            await self._update_output()
            return
        elif self.vals.get("proteome"):
            self.output_md.append(f'Using proteome {self.vals.get("proteome")}')
            await self._update_output()
        elif self.vals.get("pdb_id"):
            self.output_md.append(f'Using PDB {self.vals.get("pdb_id")}')
            await self._update_output()
        elif self.vals.get("gene_name"):
            self.output_md.append(f'Using Gene name ID {self.vals.get("gene_name")}')
            await self._update_output()

        if not (self.vals.get("sdf") or self.vals.get("smiles") or self.vals.get("pubchem_id")):
            self.output_md.append("No SDF or SMILES or PubChem ID supplied.")
            await self._update_output()
            return
        else:
            message_sender.label = "Docking..."
            message_sender.button_style = "white on dark_green"
            self.refresh() # ???

            self.output_md.append('## Running docking with params')
            await self._update_output()

            prot = self.vals.get('pdb_id') or self.vals.get('gene_name') or self.vals.get('proteome')
            sm = self.vals.get('pubchem_id') or self.vals.get('sdf') or self.vals.get('smiles')
            out_tsv = Path(OUTDIR, f"{self._get_friendly_id()}.tsv")

            dk_gener = run_docking(prot, sm, out_tsv)

            # ---------------------------------------------------------------
            # Progress bar for docking -- max_sdf_confs is not known though?
            #
            dk_task = self.progress_bar.add_task(f"[red]smina:", total=MAX_SDF_CONFS + 1)
            self.progress_bar.update(dk_task, advance=0)
            self.progress_panel.refresh()
            self.refresh()

            err = [] # keep track of errors returned on running smina
            for l in iter(dk_gener):
                self.progress_bar.update(dk_task, advance=1)
                # Both refreshes are necessary?
                self.progress_panel.refresh()
                self.refresh()

                if l is not None:
                    err.append(l)

            message_sender.label = "Start docking"
            message_sender.button_style = "white on blue"

            if Path(out_tsv).exists():
                self.output_md.append(f"Output to {out_tsv}")
                await self._update_output()
            
            if len(err) > 0:
                self.output_md.extend(err)
                await self._update_output()


    async def _download_file(self, message_sender, url, name):
        message_sender.label = f"Downloading {name}..."
        message_sender.button_style = "white on dark_green"
        self.refresh()

        if name == "smina":
            out_file = "smina"
        else:
            out_file = None

        dl_gener = gen_dl(url, out_dir=DATADIR, out_file=out_file)
        dl_total = next(dl_gener)

        dl_task = self.progress_bar.add_task(f"[red]Downloading {name}:", total=dl_total)
        for _ in iter(dl_gener):
            self.progress_bar.update(dl_task, advance=1)
            # Both refreshes are necessary?
            self.progress_panel.refresh()
            self.refresh()

        log(f"dl_done")
        message_sender.label = f"{name} downloaded"
        self.output_md.append(f"Downloaded {name} to {DATADIR}")
        await self._update_output()

    async def handle_button_pressed(self, message: ButtonPressed) -> None:
        log("ZZZ3", message.sender.name)
        if message.sender.name == "Start docking":
            await self._start_docking(message.sender)
        elif message.sender.name == "Download EquiBind":
            await self._download_file(message.sender, URLS[("equibind", MAC_OR_LINUX)], "EquiBind")
        elif message.sender.name == "Download smina" and "smina downloaded" not in message.sender.label:
            await self._download_file(message.sender, URLS[("smina", MAC_OR_LINUX)], "smina")
        elif message.sender.name == "Download proteome":
            await self.action_show_proteome_list()
        elif message.sender.name == "Proteome":
            await self.action_show_downloaded_proteome_list()

        if message.sender.name in self.rowcol_dict:
            self.row = self.rowcol_dict[message.sender.name][0]
            self.col = self.rowcol_dict[message.sender.name][1][0]

    async def on_click(self, event: events.Click):
        log("ZZZ1", message.sender.name)

    async def handle_on_click(self, event: events.Click):
        log("ZZZ2", message.sender.name)

    # ------------------------------------------------------------------
    # Actions
    #
    async def action_hide_proteome_list(self) -> None:
        self.proteome_list.animate("layout_offset_x", -40)
        await self.bind(Keys.Up, "move_up", "move up")
        await self.bind(Keys.Down, "move_down", "move down")
        await self.bind(Keys.Left, "move_left", "move left")
        await self.bind(Keys.Right, "move_right", "move right")

        # this is nasty but otherwise it closes immediately
        def _hide_p(): self.showing_proteome_list = False
        self.set_timer(0.1, _hide_p)

    async def action_show_proteome_list(self) -> None:
        self.proteome_list.animate("layout_offset_x", 0)
        await self.bind(Keys.Up, "move_up_proteome", "move up proteome")
        await self.bind(Keys.Down, "move_down_proteome", "move down proteome")
        await self.bind(Keys.Escape, "hide_proteome_list", "hide proteome list")
        await self.bind(Keys.Enter, "enter_pressed", "enter pressed")

        # this is nasty but otherwise it closes immediately
        def _show_p(): self.showing_proteome_list = True
        self.set_timer(0.1, _show_p)

    async def action_show_downloaded_proteome_list(self) -> None:
        # TODO
        # based on action_show_proteome_list
        # with undownloaded ones removed
        for n in self.proteome_list.nodes:
            log("node1", dir(self.proteome_list.nodes[n].parent))
            log("node2", dir(self.proteome_list.nodes[n]))
            if n != 3:
                pass #self.proteome_list.nodes

        self.proteome_list.animate("layout_offset_x", 0)
        await self.bind(Keys.Up, "move_up_proteome", "move up proteome")
        await self.bind(Keys.Down, "move_down_proteome", "move down proteome")
        await self.bind(Keys.Escape, "hide_proteome_list", "hide proteome list")
        await self.bind(Keys.Enter, "enter_pressed", "enter pressed")

        # this is nasty but otherwise it closes immediately
        def _show_p(): self.showing_proteome_list = True
        self.set_timer(0.1, _show_p)

    async def show_hide_proteome_list(self) -> None:
        if self.showing_proteome_list:
            await self.action_hide_proteome_list()
        else:
            await self.action_show_proteome_list()

    async def _change_focus(self) -> None:
        for child in self.children:
            try:
                log("grandchildren found", child.children)
            except:
                pass
            if hasattr(child, "row") and hasattr(child, "cols"):
                if child.row == self.row and self.col in child.cols:
                    await self.set_focus(child)

    async def action_move_up(self) -> None:
        self.row = max(1, self.row - 1)
        await self._change_focus()

    async def action_move_down(self) -> None:
        self.row = min(4, self.row + 1)
        await self._change_focus()

    async def action_move_left(self) -> None:
        self.col = max(1, self.col - 1)
        await self._change_focus()

    async def action_move_right(self) -> None:
        self.col = min(3, self.col + 1)
        await self._change_focus()

    async def action_move_up_proteome(self) -> None:
        self.row_proteome = max(1, self.row_proteome - 1)
        self.proteome_list.hover_node = self.proteome_list.nodes[self.row_proteome].id
        log("SELECT UP", self.proteome_list.hover_node)

    async def action_move_down_proteome(self) -> None:
        self.row_proteome = min(len(PROTEOMES.keys()), self.row_proteome + 1)
        self.proteome_list.hover_node = self.proteome_list.nodes[self.row_proteome].id
        log("SELECT DOWN", self.proteome_list.hover_node)

    async def action_enter_pressed(self) -> None:
        if self.showing_proteome_list:
            log("SELECT PROTEOME", self.showing_proteome_list, self.proteome_list.hover_node)
            if self.proteome_list.hover_node:
                log("SELECT PROTEOME HOVER", self.proteome_list.hover_node)
                for child in self.children:
                    if child.name == "Proteome":
                        log("SELECT PROTEOME TEXT")
                        if self.proteome_list.hover_node != 0:
                            child.text = self.proteome_list.nodes[self.proteome_list.hover_node].label
                            child.refresh()
                await self.action_hide_proteome_list()

    async def handle_tree_click(self, message: TreeClick[dict]) -> None:
        """Called in response to a tree click."""
        log(f"Tree Click {message} {message.node.label}")
        for child in self.children:
            if child.name == "Proteome":
                if message.node.id != 0:
                    child.text = message.node.label
                    child.refresh()
        await self.action_hide_proteome_list()

    # I think this grabs messages before handle_button_pressed got to them?
    #async def on_message(self, message):
    #    if message.sender.name == "Start docking":
    #        print(message)

    async def on_mount(self) -> None:
        """Make a simple grid arrangement."""

        self.grid = await self.view.dock_grid(edge="right", name="grid")
        grid = self.grid

        grid.add_column(fraction=1, name="l1")
        grid.add_column(fraction=1, name="l2")
        grid.add_column(fraction=1, name="l3")
        grid.add_column(fraction=6, name="right")

        grid.add_row(fraction=1, name="r1")
        grid.add_row(fraction=1, name="r2")
        grid.add_row(fraction=1, name="r3")
        grid.add_row(fraction=1, name="r4")

        grid.add_areas(
            dl_equibind="l1,r1",
            dl_smina="l2,r1",
            dl_proteome="l3,r1",
            enter_pdb="l1,r2",
            enter_gene_name="l2,r2",
            enter_proteome="l3,r2",
            enter_pubchem="l1,r3",
            enter_smiles="l2,r3",
            enter_sdf="l3,r3",
            start_docking="l1-start|l3-end,r4",
            output="right,r1-start|r3-end",
            progress="right,r4"
        )

        # Right hand side panels
        self.output = ScrollView(name="Output", gutter=1)
        self.progress_bar = Progress()
        self.progress_panel = Static(name="Progress", renderable=Align.center(self.progress_bar, vertical="middle"))

        self.rowcol_dict = {
            "Download EquiBind": [1, [1]],
            "Download smina": [1, [2]],
            "Download proteome": [1, [3]],
            "PDB": [2, [1]],
            "Gene name": [2, [2]],
            "Proteome": [2, [3]],
            "PubChem ID": [3, [1]],
            "SMILES": [3, [2]],
            "SDF": [3, [3]],
            "Start docking": [4, [1,2,3]]
        }

        equibind_label = "EquiBind downloaded" if Path(DATADIR, "EquiBind").exists() else "Download EquiBind"
        smina_label = "smina downloaded" if Path(DATADIR, "smina").exists() else "Download smina"
        gnina_label = "gnina downloaded" if Path(DATADIR, "gnina").exists() else "Download gnina"

        grid.place(
            # downloads
            dl_equibind=GridButton(name="Download EquiBind", label=equibind_label, row=1, cols=[1]),
            dl_smina=GridButton(name="Download smina", label=smina_label, row=1, cols=[2]),
            dl_proteome=GridButton(name="Download proteome", label="Download proteome", row=1, cols=[3]),
            # text
            enter_pdb=TextInputPanel(name="PDB", val="pdb_id", row=2, cols=[1]),
            enter_gene_name=TextInputPanel(name="Gene name", val="gene_name", row=2, cols=[2]),
            enter_proteome=TextPanel(name="Proteome", val="proteome", row=2, cols=[3]),
            enter_pubchem=TextInputPanel(name="PubChem ID", val="pubchem_id", row=3, cols=[1]),
            enter_smiles=TextInputPanel(name="SMILES", val="smiles", row=3, cols=[2]),
            enter_sdf=TextInputPanel(name="SDF", val="sdf", row=3, cols=[3]),
            start_docking=GridButton(name="Start docking", label="Start docking", row=4, cols=[1,2,3]),
            # right panel
            output=self.output,
            progress=self.progress_panel,
        )

        # not sure what call_later does
        await self.call_later(self._update_output)

        # --------------
        # sidebar test!
        #
        self.proteome_list = TreeControl("Press Escape to dismiss\nProteomes", {})
        for pname in PROTEOMES.keys():
            await self.proteome_list.add(self.proteome_list.root.id, pname, {"pname": pname})
        await self.proteome_list.root.expand()
        self.proteome_list.layout_offset_x = -40
        await self.view.dock(self.proteome_list, edge="left", size=40, z=1)

        # TODO make this work
        #if not Path(DATADIR, "smina").exists():
        #    self._download_smina()

GridTest.run(log="textual.log")
