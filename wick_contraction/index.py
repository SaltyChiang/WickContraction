class Index:
    _letter = ""
    _counter = -1
    _index = -1

    @classmethod
    def new(cls):
        name = f"{cls._letter}{cls._counter}"
        cls._counter += 1
        return cls(name)

    def __init__(self, name: str):
        self.name = name

    def __str__(self):
        return self.name

    def __repr__(self):
        return f"{type(self).__name__}({self.name!r})"

    def __eq__(self, other):
        return type(self) is type(other) and self.name == other.name

    def __hash__(self):
        return hash((type(self), self.name))


class SpinIndex(Index):
    _letter = "σ"  # "αβγδεζηθικλμνξοπρστυφχψω"
    _counter = 0
    _index = 0


class ColorIndex(Index):
    _letter = "c"  # "abcdefghijklmnopqrstuvwxyz"
    _counter = 0
    _index = 0
