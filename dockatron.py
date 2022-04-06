"""
# Dockatron
Dockatron relies on the following tools.

## Cite smina
- http://pubs.acs.org/doi/abs/10.1021/ci300604z

## Cite EquiBind
- https://arxiv.org/abs/2202.05146

## Cite AlphaFold:
- Jumper, J et al. Highly accurate protein structure prediction with AlphaFold. Nature (2021).
- Varadi, M et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Research (2021).
"""

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
from textual import events
from textual.app import App
from textual.keys import Keys
from textual.message import Message
from textual.widget import Reactive, Widget
from textual.widgets import Button, ButtonPressed, ScrollView, Static

import os
import sys
import time
import requests
import platform
import subprocess

from typing import List

# -------------------------------------------------------------------------------------------------
# Globals
#
if platform.uname().system == "Linux":
    mac_or_linux = "linux"
elif platform.uname().system == "Darwin":
    mac_or_linux = "osx" if platform.mac_ver()[0][:2]<"12" else "osx12"
else:
    raise ValueError(f"{platform.uname()} not supported")

URLS = {
    ("equibind", "linux"): "https://github.com/HannesStark/EquiBind/archive/refs/heads/main.zip",
    ("gnina", "linux"): "https://github.com/gnina/gnina/releases/download/v1.0/gnina",
    ("smina", "linux"): "https://sourceforge.net/projects/smina/files/smina.static/download",
    ("smina", "osx"): "https://sourceforge.net/projects/smina/files/smina.osx/download",
    ("smina", "osx12"): "https://sourceforge.net/projects/smina/files/smina.osx.12/download",
    ("proteome", "yeast"): "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v2.tar",
    ("proteome", "human"): "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v2.tar",
}

DATADIR = os.path.join(os.getcwd(), "dockatron")

def gen_dl(url, chunk_size=1_048_576, out_dir=None, out_file=None):
    """Generator for downloading files"""
    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
    else:
        out_dir = "."

    resp = requests.get(url, stream=True)
    total = int(resp.headers.get('content-length', 0))
    yield total // chunk_size

    with open(out_file or f"{out_dir}/{url.split('/')[-1]}", 'wb') as out:
        for data in resp.iter_content(chunk_size=chunk_size):
            out.write(data)
            yield


def DEPRECATED_install(software_name):
    msg = []
    if software_name == "EquiBind":
        p = subprocess.run(["echo", "git", "clone", URLS[("equibind", mac_or_linux)]], capture_output=True)
        msg.append(p.stdout.decode())
        msg.append("`Done`")
        p = subprocess.run(["echo", "conda", "env", "create", "-f", "environment.yml"], capture_output=True)
        msg.append(p.stdout.decode())
    elif software_name == "smina":
        p = subprocess.run(["echo", "wget", URLS[("smina", mac_or_linux)]])
        msg.append(p.stdout.decode())
    elif software_name == "gnina":
        p = subprocess.run(["echo", "wget", URLS[("gnina", mac_or_linux)]])
        msg.append(p.stdout.decode())
    elif software_name == "yeast_proteome":
        p = subprocess.run(["echo", "wget", URLS[("proteome", "yeast")]], bufsize=1, universal_newlines=True)
        msg.append(p.stdout.decode())
        p = subprocess.run(["echo", "tar", "xvf", URLS[("proteome", "yeast")].split('/')[-1]])
        msg.append(p.stdout.decode())

    return '\n'.join(msg)


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
        if event.key == Keys.Enter:
            await self.emit(ButtonPressed(self))

    #async def on_click(self, event: events.Click) -> None:
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
        self.output_md = ["# Dockatron Output"]

        self.vals = {}

    async def on_load(self) -> None:
        """Bind keys here."""
        await self.bind(Keys.Up, "move_up", "move up")
        await self.bind(Keys.Down, "move_down", "move down")
        await self.bind(Keys.Left, "move_left", "move left")
        await self.bind(Keys.Right, "move_right", "move right")

    async def _update_output(self):
        await self.output.update(Markdown('\n\n'.join(self.output_md), hyperlinks=True))

    async def _start_docking(self, message_sender):
        if not (self.vals.get("proteome") or self.vals.get("uniprot_id") or self.vals.get("gene_name")):
            self.output_md.append("No proteome or uniprot_id or gene_name supplied.")
            await self._update_output()
            return
        elif self.vals.get("proteome"):
            self.output_md.append(f'Using proteome {self.vals.get("proteome")}')
            await self._update_output()
        elif self.vals.get("uniprot_id"):
            self.output_md.append(f'Using UniProt ID {self.vals.get("uniprot_id")}')
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
            # test set_timer only
            self.output_md.append('## Running docking with params')
            self.set_timer(0.1, self._update_output)


    async def _download_file(self, message_sender, url, name):
        message_sender.label = Status(f"Downloading {name}...")
        message_sender.button_style = "white on dark_green"
        self.refresh()

        dl_gener = gen_dl(url, out_dir=DATADIR)
        dl_total = next(dl_gener)

        # These two lines appears to be sometimes necessary to show the progress update???
        self.output_md.append(f"Downloading {name} to {DATADIR}")
        await self._update_output()

        dl_task = self.progress_bar.add_task(f"[red]Downloading {name}:", total=dl_total)
        for _ in iter(dl_gener):
            self.progress_bar.update(dl_task, advance=1)
            # Both refreshes are necessary!!!
            self.progress_panel.refresh()
            self.refresh()

        message_sender.label = f"{name} downloaded"
        self.output_md.append(f"Downloaded {name} to {DATADIR}")
        await self._update_output()


    async def handle_button_pressed(self, message: ButtonPressed) -> None:
        if message.sender.name == "Start docking":
            await self._start_docking(message.sender)
        elif message.sender.name == "Download EquiBind":
            await self._download_file(message.sender, URLS[("equibind", mac_or_linux)], "EquiBind")
        elif message.sender.name == "Download smina":
            await self._download_file(message.sender, URLS[("smina", mac_or_linux)], "smina")

    # ------------------------------------------------------------------
    # Selecting boxes
    #
    async def _change_focus(self) -> None:
        for child in self.children:
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
            dl_human_proteome="l3,r1",
            enter_uniprot_id="l1,r2",
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

        grid.place(
            # download
            dl_equibind=GridButton(name="Download EquiBind", label="Download EquiBind", row=1, cols=[1]),
            dl_smina=GridButton(name="Download smina", label="Download smina", row=1, cols=[2]),
            dl_human_proteome=GridButton(name="Download human proteome", label="Download human proteome", row=1, cols=[3]),
            # text entry
            enter_uniprot_id=TextInputPanel(name="UniProt ID", val="uniprot_id", row=2, cols=[1]),
            enter_gene_name=TextInputPanel(name="Gene name", val="gene_name", row=2, cols=[2]),
            enter_proteome=TextInputPanel(name="Proteome", val="proteome", row=2, cols=[3]),
            enter_pubchem=TextInputPanel(name="PubChem ID", val="pubchem_id", row=3, cols=[1]),
            enter_smiles=TextInputPanel(name="SMILES", val="smiles", row=3, cols=[2]),
            enter_sdf=TextInputPanel(name="SDF", val="sdf", row=3, cols=[3]),
            # button
            start_docking=GridButton(name="Start docking", label="Start docking", row=4, cols=[1,2,3]),
            output=self.output,
            progress=self.progress_panel,
        )

        # hmm, this has to be at the end of this class to work
        #async def init_markdown() -> None:
        #    md = Markdown(f"# Dockatron Output\n", hyperlinks=True)
        #    await self.output.update(md)
        await self.call_later(self._update_output)

GridTest.run(log="textual.log")
