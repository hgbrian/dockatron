import rich
from rich import box
from rich.align import Align
from rich.console import RenderableType
from rich.panel import Panel
from rich.pretty import Pretty
from rich.style import Style
from textual import events
from textual.widget import Reactive, Widget
from textual.app import App
from textual.keys import Keys

#print(events.Leave.__doc__)
#1/0

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
class GridButton(Widget, can_focus=True):

    has_focus: Reactive[bool] = Reactive(False)
    mouse_over: Reactive[bool] = Reactive(False)
    keypress: Reactive[str] = Reactive("")
    style: Reactive[str] = Reactive("")
    height: Reactive[int] = Reactive(None)

    def __init__(self, *, name: str = None, height: int = None) -> None:
        super().__init__(name=name)
        self.height = height
        self.text = "Click here to dock"
        self.title = name
        self.style = Style(color="white", bgcolor="green", encircle=True)

    def render(self) -> RenderableType:
        return Panel(
            Align.center(self.text, vertical="middle"),
            border_style="red" if self.has_focus else "green",
            box=box.HEAVY if self.has_focus else box.ROUNDED,
            title=self.title,
            style=self.style,
        )

    async def on_focus(self, event: events.Focus) -> None:
        self.has_focus = True
        start_docking()

    async def on_blur(self, event: events.Blur) -> None:
         self.has_focus = False



class GridTest(App):
    async def on_mount(self) -> None:
        """Make a simple grid arrangement."""

        grid = await self.view.dock_grid(edge="right", name="grid")

        grid.add_column(fraction=1, name="left")
        grid.add_column(fraction=2, name="center")
        grid.add_column(fraction=1, name="right")

        grid.add_row(fraction=1, name="r1")
        grid.add_row(fraction=1, name="r2")
        grid.add_row(fraction=1, name="r3")
        grid.add_row(fraction=1, name="r4")
        grid.add_row(fraction=1, name="r5")
        grid.add_row(fraction=1, name="r6")

        grid.add_areas(
            area1="left,r1",
            area2="left,r2",
            area3="left,r3",
            area4="left,r4",
            area5="left,r5",
            area6="left,r6",
            progress="center-start|right-end,r1-start|r6-end"
        )

        grid.place(
            area1=Placeholder(name="Download stuff 1"),
            area2=Placeholder(name="Download stuff 2"),
            area3=Placeholder(name="Download stuff 3"),
            area4=Placeholder2(name="Enter Pubchem ID"),
            area5=Placeholder2(name="Enter SMILES"),
            area6=GridButton(name=""),
            progress=Placeholder2(name="Progress"),
        )


def start_docking():
    GridTest.post_message(GridTest, message="start")

GridTest.run(title="Grid Test", log="textual.log")
