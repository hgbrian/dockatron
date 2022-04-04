"""
If you find smina useful, please cite our paper:
http://pubs.acs.org/doi/abs/10.1021/ci300604z
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
from rich.style import Style, StyleType
from textual import events
from textual.message import Message
from textual.widget import Reactive, Widget
from textual.widgets import Button, ButtonPressed, ScrollView, Static
from textual.app import App
from textual.keys import Keys

import time
import subprocess

URLS = {
    "equibind": "https://github.com/HannesStark/EquiBind",
    "gnina_linux": "https://github.com/gnina/gnina/releases/download/v1.0/gnina",
    "smina_linux": "https://sourceforge.net/projects/smina/files/smina.static/download",
    "smina_mac": "https://sourceforge.net/projects/smina/files/smina.osx/download",
    "smina_mac12": "https://sourceforge.net/projects/smina/files/smina.osx.12/download",
}

def install(software_name):
    msg = []
    if software_name == "EquiBind":
        p = subprocess.run(["echo", "git", "clone", URLS["equibind"]], capture_output=True)
        msg.append(p.stdout.decode())
        subprocess.run(["echo", "conda", "env", "create", "-f", "environment.yml"], capture_output=True)
        msg.append(p.stdout.decode())
    elif software_name == "smina":
        p = subprocess.run(["echo", "wget", URLS["smina"]])
        msg.append(p.stdout.decode())
    elif software_name == "gnina":
        p = subprocess.run(["echo", "wget", URLS["gnina_linux"]])
        msg.append(p.stdout.decode())
    return ''.join(msg)


@rich.repr.auto(angular=False)
class Placeholder(Widget, can_focus=True):

    has_focus: Reactive[bool] = Reactive(False)
    mouse_over: Reactive[bool] = Reactive(False)
    style: Reactive[str] = Reactive("")
    height: Reactive = Reactive(None)

    def __init__(self, *, name: str = None, height: int = None, row:int = None, col:int = None) -> None:
        super().__init__(name=name)
        self.height = height
        self.row = row
        self.col = col

    def __rich_repr__(self) -> rich.repr.Result:
        yield "name", self.name
        yield "has_focus", self.has_focus, False
        yield "mouse_over", self.mouse_over, False

    def render(self) -> RenderableType:
        return Panel(
            Align.center(
                Pretty(self, no_wrap=True, overflow="ellipsis"), vertical="middle"
            ),
            title=self.__class__.__name__,
            border_style="green" if self.mouse_over else "blue",
            box=box.HEAVY if self.has_focus else box.ROUNDED,
            style=self.style,
            height=self.height,
        )

    async def on_focus(self, event: events.Focus) -> None:
        self.has_focus = True

    async def on_blur(self, event: events.Blur) -> None:
        self.has_focus = False

    async def on_enter(self, event: events.Enter) -> None:
        self.mouse_over = True

    async def on_leave(self, event: events.Leave) -> None:
        self.mouse_over = False



@rich.repr.auto(angular=False)
class Placeholder2(Widget, can_focus=True):

    has_focus: Reactive[bool] = Reactive(False)
    mouse_over: Reactive[bool] = Reactive(False)
    keypress: Reactive[str] = Reactive("")
    style: Reactive[str] = Reactive("")
    height: Reactive[int] = Reactive(None)

    def __init__(self, *, name: str = None, height: int = None, row:int = None, col:int = None) -> None:
        super().__init__(name=name)
        self.height = height
        self.text = ""
        self.title = name
        self.row = row
        self.col = col

    def __rich_repr__(self) -> rich.repr.Result:
        yield "???????????"

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

    # async def on_enter(self, event: events.Enter) -> None:
    #     print(event)
    #     self.mouse_over = True

    # async def on_leave(self, event: events.Leave) -> None:
    #     print(event)
    #     self.mouse_over = False

    async def on_key(self, event: events.Key):
        if event.key == Keys.Enter:
            self.emit(ButtonPressed(self))
        elif event.key == Keys.ControlH:
            self.text = self.text[:-1]
            self.refresh()
        elif self.has_focus and len(event.key)==1:
            self.text = self.text + event.key
            self.refresh()

@rich.repr.auto(angular=False)
class GridButton(Button):
    def __init__(self, label:str, *, name: str = None, row:int = None, col:int = None) -> None:
        super().__init__(label)
        self.row = row
        self.col = col

    # ButtonPressed does not work
    async def on_focus(self, event: events.Focus) -> None:
        self.has_focus = True
        self.label = "[ Start docking ]"
        self.button_style = "white on blue"
        self.refresh()

    async def on_blur(self, event: events.Blur) -> None:
        self.has_focus = False
        self.label = "Start docking"
        self.button_style = "white on dark_blue"
        self.refresh()

    async def on_key(self, event: events.Key):
        if event.key == Keys.Enter:
            await self.emit(ButtonPressed(self))

    #async def on_click(self, event: events.Click) -> None:
    #    #self.has_focus = True
    #    self.start_docking()
    #    self.label = "Docking..."
    #    self.button_style = "white on dark_green"

    def start_docking(self):
        self.post_message_from_child_no_wait(Message(self))

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
        self.output_md = ["# Output"]

    async def on_load(self) -> None:
        """Bind keys here."""
        await self.bind(Keys.Up, "move_up", "move up")
        await self.bind(Keys.Down, "move_down", "move down")
        await self.bind(Keys.Left, "move_left", "move left")
        await self.bind(Keys.Right, "move_right", "move right")

    async def _update_output(self):
        await self.output.update(Markdown('\n'.join(self.output_md), hyperlinks=True))
        # Test
        self.output_md.append(install("EquiBind"))
        self.set_timer(1, self._update_output)

    async def handle_button_pressed(self, message: ButtonPressed) -> None:
        if message.sender.name == "Start docking":
            message.sender.label = "Docking..."
            message.sender.button_style = "white on dark_green"
            self.output_md.append(install("EquiBind"))
            self.set_timer(1, self._update_output)

    async def _change_focus(self) -> None:
        for child in self.children:
            if hasattr(child, "row") and hasattr(child, "col"):
                if child.row == self.row and child.col == self.col:
                    await self.set_focus(child)

    async def action_move_up(self) -> None:
        self.row = max(1, self.row - 1)
        await self._change_focus()

    async def action_move_down(self) -> None:
        self.row = min(5, self.row + 1)
        await self._change_focus()

    async def action_move_left(self) -> None:
        self.col = max(1, self.col - 1)
        await self._change_focus()

    async def action_move_right(self) -> None:
        self.col = min(2, self.row + 1)
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
        grid.add_column(fraction=6, name="right")

        grid.add_row(fraction=1, name="r1")
        grid.add_row(fraction=1, name="r2")
        grid.add_row(fraction=1, name="r3")
        grid.add_row(fraction=1, name="r4")
        grid.add_row(fraction=1, name="r5")

        grid.add_areas(
            dl_1="l1-start|l2-end,r1",
            dl_2="l1-start|l2-end,r2",
            enter_protein="l1,r3",
            enter_proteome="l2,r3",
            enter_pubchem="l1,r4",
            enter_smiles="l2,r4",
            start_docking="l1-start|l2-end,r5",
            output="right,r1-start|r4-end",
            progress="right,r5"
        )

        self.output = ScrollView(name="Output", gutter=1)


        self.progress_bar = Progress()
        task1 = self.progress_bar.add_task("[red]Downloading:", total=1000)
        self.progress_bar.update(task1, advance=100)
        self.progress_panel = Static(name="Progess", renderable=Align.center(self.progress_bar, vertical="middle"))

        grid.place(
            dl_1=Placeholder2(name="Download stuff 1", row=1, col=1),
            dl_2=Placeholder2(name="Download stuff 2", row=2, col=1),
            enter_protein=Placeholder2(name="Enter Proteome", row=3, col=1),
            enter_proteome=Placeholder2(name="Enter UniProt ID", row=3, col=2),
            enter_pubchem=Placeholder2(name="Enter Pubchem ID", row=4, col=1),
            enter_smiles=Placeholder2(name="Enter SMILES", row=4, col=2),
            start_docking=GridButton(name="start_docking", label="Start docking", row=5, col=1),
            output=self.output,
            progress=self.progress_panel,
        )

        # hmm, this has to be at the end of the class to work
        async def get_markdown(filename: str) -> None:
            readme = Markdown(f"# Output\n{time.time()}", hyperlinks=True)
            await self.output.update(readme)

        await self.call_later(get_markdown, "richreadme.md")

GridTest.run(log="textual.log")
