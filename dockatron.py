import rich
from rich import box
from rich.align import Align
from rich.console import RenderableType
from rich.markdown import Markdown
from rich.panel import Panel
from rich.pretty import Pretty
from rich.screen import Screen
from rich.style import Style
from textual import events
from textual.message import Message
from textual.widget import Reactive, Widget
from textual.widgets import Button, ButtonPressed, ScrollView
from textual.app import App
from textual.keys import Keys

proteome = None


@rich.repr.auto(angular=False)
class Placeholder(Widget, can_focus=True):

    has_focus: Reactive[bool] = Reactive(False)
    mouse_over: Reactive[bool] = Reactive(False)
    style: Reactive[str] = Reactive("")
    height: Reactive = Reactive(None)

    def __init__(self, *, name: str = None, height: int = None) -> None:
        super().__init__(name=name)
        self.height = height

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

    def __init__(self, *, name: str = None, height: int = None) -> None:
        super().__init__(name=name)
        self.height = height
        self.text = ""
        self.title = name

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
        #print("event", event, event.key)
        if event.key == Keys.Enter:
            self.has_focus = False
        elif event.key == Keys.ControlH:
            self.text = self.text[:-1]
            self.has_focus = False
            self.has_focus = True
        if self.has_focus and len(event.key)==1:
            self.text = self.text + event.key
            self.has_focus = False
            self.has_focus = True

@rich.repr.auto(angular=False)
class GridButton(Button):
    # ButtonPressed does not work
    #async def on_button_pressed(self, message: ButtonPressed) -> None:
    async def on_focus(self, event: events.Focus) -> None:
        #print("on focus")
        self.has_focus = True
        self.start_docking()

    def start_docking(self):
        pass
        # i have to post a real Message object here
        self.post_message_from_child_no_wait(Message(self))
        #raise SystemExit(dir(self)) #.output.text = "ASDF"


class GridTest(App):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.log_verbosity = 9
        self.grid = None

    async def on_message(self, message):
        if message.sender.name =="go":
            readme = Markdown("# Equibind -> /tmp/20220402_yeast_zearalenone\nasdf", hyperlinks=True)
            await self.output.update(readme)

    async def on_mount(self) -> None:
        """Make a simple grid arrangement."""

        self.log_verbosity = 9
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
            output="right,r1-start|r5-end"
        )

        self.output = ScrollView(name="Output", gutter=1)

        grid.place(
            dl_1=Placeholder(name="Download stuff 1"),
            dl_2=Placeholder(name="Download stuff 2"),
            enter_proteome=Placeholder2(name="Enter UniProt ID"),
            enter_protein=Placeholder2(name="Enter Proteome"),
            enter_pubchem=Placeholder2(name="Enter Pubchem ID"),
            enter_smiles=Placeholder2(name="Enter SMILES"),
            start_docking=GridButton(label="Start docking", name="go"),
            output=self.output,
        )

        async def get_markdown(filename: str) -> None:
            readme = Markdown("# Output", hyperlinks=True)
            await self.output.update(readme)
        await self.call_later(get_markdown, "richreadme.md")


myapp = GridTest.run(title="Grid Test", log="textual.log")
# never gets here!
#print(myapp)
#raise SystemExit(str(myapp))