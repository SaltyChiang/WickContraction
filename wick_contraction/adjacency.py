from typing import List, Optional, Tuple, Union

try:
    import sympy
except ImportError:
    sympy = None

from .index import ColorIndex, SpinIndex
from .gamma import Gamma, GAMMA_5
from .quark import QuarkPrefixType, QuarkLocal
from .tensor import (
    SpinColorTensorType,
    QuarkFieldTensor,
    SpinGammaTensor,
    SpinProjectorTensor,
    ColorDeltaTensor,
    ColorEpsilonTensor,
)


class AdjacencyEdge:
    def __init__(self) -> None:
        self.flavor: Optional[str] = None
        self.sink: Optional[str] = None
        self.source: Optional[str] = None
        self.prefix: QuarkPrefixType = QuarkLocal()
        self.suffix: QuarkPrefixType = QuarkLocal()
        self.conjugate = ""

    def set(self, quark: QuarkFieldTensor, antiquark: QuarkFieldTensor) -> None:
        if quark.flavor != antiquark.flavor:
            raise ValueError("contracted quark and antiquark should have the same flavor")
        self.flavor = quark.flavor
        self.sink, self.source = quark.location, antiquark.location
        self.prefix, self.suffix = quark.prefix, antiquark.prefix

    def adjoint(self) -> None:
        self.source, self.sink = self.sink, self.source
        self.prefix, self.suffix = self.suffix.adjoint(), self.prefix.adjoint()
        self.conjugate = "" if self.conjugate else "^*"

    def should_adjoint(self, order: List[str]) -> bool:
        return (
            self.flavor is not None
            and self.sink in order
            and self.source in order
            and order.index(self.sink) > order.index(self.source)
        )

    def degenerate(self):
        if self.flavor in ("u", "d"):
            self.flavor = "l"

    def __str__(self) -> str:
        if self.flavor is None:
            return "0"
        else:
            return f"{self.prefix}{self.flavor}{self.conjugate}({self.sink},{self.source}){self.suffix}"


class AdjacencyTerm:
    def __init__(self, num_quark: int, factor: Union[int, float, complex], tensors: List[SpinColorTensorType]) -> None:
        self.num_quark = num_quark
        self.factor = factor
        self.tensors = list(tensors)
        self.matrix: List[List[AdjacencyEdge]] = [[AdjacencyEdge() for _ in range(num_quark)] for _ in range(num_quark)]
        self.spin_indices: List[SpinIndex] = []
        self.color_indices: List[ColorIndex] = []

    def _swap_index(self, idx: str, snk: str, src: str) -> str:
        if idx == src:
            return snk
        elif idx == snk:
            return src
        else:
            return idx

    def _gamma_apply_gamma_5(self, left: str, right: str, gamma: Gamma, snk: str, src: str) -> Tuple[int, Gamma]:
        if left in (snk, src):
            gamma = GAMMA_5 @ gamma
        if right in (snk, src):
            gamma = gamma @ GAMMA_5
        factor = 1
        if gamma.factor == -1:
            factor = -1
        elif gamma.factor != 1:
            raise ValueError("unexpected factor found when swapping spin indices")
        return factor, Gamma(gamma.index)

    def _projector_apply_gamma_5(
        self, left: str, right: str, factors: List[Union[int, float, complex]], snk: str, src: str
    ) -> List[Union[int, float, complex]]:
        new_factors: List[Union[int, float, complex]] = [0 for _ in range(16)]
        if left in (snk, src):
            for i in range(16):
                if factors[i] != 0:
                    gamma = GAMMA_5 @ Gamma(i, factors[i])
                    new_factors[gamma.index] += gamma.factor
            factors = new_factors
        new_factors: List[Union[int, float, complex]] = [0 for _ in range(16)]
        if right in (snk, src):
            for i in range(16):
                if factors[i] != 0:
                    gamma = Gamma(i, factors[i]) @ GAMMA_5
                    new_factors[gamma.index] += gamma.factor
            factors = new_factors
        return factors

    def spin_swap(self, row: int, col: int) -> int:
        spin_snk = f"α{row}"
        spin_src = f"β{col}"
        factor = 1
        for i, tensor in enumerate(self.tensors):
            if isinstance(tensor, SpinGammaTensor):
                gamma, left, right = tensor
                left = self._swap_index(left, spin_snk, spin_src)
                right = self._swap_index(right, spin_snk, spin_src)
                gamma_factor, gamma = self._gamma_apply_gamma_5(left, right, gamma, spin_snk, spin_src)
                factor *= gamma_factor
                tensor = SpinGammaTensor(gamma, left, right)
                self.tensors[i] = tensor
            elif isinstance(tensor, SpinProjectorTensor):
                factors, left, right = tensor
                left = self._swap_index(left, spin_snk, spin_src)
                right = self._swap_index(right, spin_snk, spin_src)
                factors = self._projector_apply_gamma_5(left, right, factors, spin_snk, spin_src)
                tensor = SpinProjectorTensor(factors, left, right)
                self.tensors[i] = tensor
        return factor

    def color_swap(self, row: int, col: int) -> None:
        color_snk = f"a{row}"
        color_src = f"b{col}"
        for i, tensor in enumerate(self.tensors):
            if isinstance(tensor, ColorDeltaTensor):
                left, right = tensor
                left = self._swap_index(left, color_snk, color_src)
                right = self._swap_index(right, color_snk, color_src)
                self.tensors[i] = ColorDeltaTensor(left, right)
            elif isinstance(tensor, ColorEpsilonTensor):
                left, middle, right = tensor
                left = self._swap_index(left, color_snk, color_src)
                middle = self._swap_index(middle, color_snk, color_src)
                right = self._swap_index(right, color_snk, color_src)
                self.tensors[i] = ColorEpsilonTensor(left, middle, right)

    def simplify(self, degenerate: bool = False, order: Optional[List[str]] = None) -> None:
        for row, row_vector in enumerate(self.matrix):
            for col, edge in enumerate(row_vector):
                if order is not None and edge.should_adjoint(order):
                    edge.adjoint()
                    self.factor *= self.spin_swap(row, col)
                    self.color_swap(row, col)
                if degenerate:
                    edge.degenerate()

    def signature(self) -> str:
        matrix_str = "[" + ", ".join(["[" + ", ".join(str(edge) for edge in row) + "]" for row in self.matrix]) + "]"
        parts = [str(t) for t in self.tensors]
        parts.append(f"A={matrix_str}")
        content = " ".join(parts)
        return content

    def __str__(self) -> str:
        content = self.signature()
        if sympy is not None:
            factor = sympy.nsimplify(self.factor)
            if factor == -1:
                content = f"-{content}"
            elif factor != 1:
                content = f"{factor} {content}"
        else:
            if self.factor == -1:
                content = f"-{content}"
            elif self.factor != 1:
                content = f"{self.factor} {content}"
        return content
