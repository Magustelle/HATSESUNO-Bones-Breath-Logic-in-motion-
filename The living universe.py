
#!/usr/bin/env python3
"""The Living Universe — executable cumulative G0–G9 Python prototype.

Source: user's Living Universe manuscript (embedded below, species-neutral revision).
Scientific scope: educational mathematical models, NOT a validated 1:1 universe
or organism. No external dependencies. Run `python living_universe_g0_g9.py`.
Use `--source` to print the complete original manuscript, `--test` for checks.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from math import exp, log, isclose
import argparse
import json
import random
from typing import Callable

R_GAS = 8.314462618

@dataclass
class UniverseState:
    generation: int = 0
    time: float = 0.0
    symbols: set[str] = field(default_factory=lambda: {"Σ₀", "◉"})
    history: list[dict] = field(default_factory=list)

    def record(self, stage: str, **data: object) -> None:
        self.generation += 1
        self.history.append({"stage": stage, "time": self.time, **data})

# G0 + G1 — origin, information and germination

def origin() -> UniverseState:
    s = UniverseState()
    s.record("G0 origin", seed="Σ₀")
    return s

def germinate(s: UniverseState) -> None:
    s.symbols.update({"Σ₁", "state", "relation", "potential"})
    s.record("G1 germination", symbols=sorted(s.symbols))

# G2 — weighted undirected interaction graphs
@dataclass
class Graph:
    names: list[str]
    weights: list[list[float]]

    def __post_init__(self):
        n = len(self.names)
        if len(set(self.names)) != n or len(self.weights) != n:
            raise ValueError("Graph needs distinct nodes and square weights")
        if any(len(row) != n for row in self.weights):
            raise ValueError("Non-square graph")
        for i in range(n):
            for j in range(n):
                if not isclose(self.weights[i][j], self.weights[j][i], abs_tol=1e-12):
                    raise ValueError("G2 graph must be symmetric")

    def adjacency(self) -> list[list[float]]:
        return [row[:] for row in self.weights]


def branch_graph(n: int = 5) -> Graph:
    names = ["◉"] + [f"◌{i}" for i in range(1, n)]
    w = [[0.0] * n for _ in range(n)]
    for i in range(1, n):
        w[0][i] = w[i][0] = 1.0
    return Graph(names, w)

# G3 — H = αI + βA. Jacobi diagonalization for real symmetric matrices.
def quantum_hamiltonian(graph: Graph, alpha: float = 1., beta: float = .4):
    a = graph.adjacency()
    return [[(alpha if i == j else 0.) + beta * a[i][j]
             for j in range(len(a))] for i in range(len(a))]

def symmetric_eigenvalues(matrix: list[list[float]], tolerance: float = 1e-12):
    a = [row[:] for row in matrix]
    n = len(a)
    if n == 0:
        return []
    for _ in range(100 * n * n):
        p, q = max(((i, j) for i in range(n) for j in range(i+1, n)),
                   key=lambda ij: abs(a[ij[0]][ij[1]]), default=(0, 0))
        if p == q or abs(a[p][q]) < tolerance:
            break
        from math import atan2, cos, sin
        theta = .5 * atan2(2*a[p][q], a[q][q]-a[p][p])
        c, s = cos(theta), sin(theta)
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        a[p][p] = c*c*app - 2*s*c*apq + s*s*aqq
        a[q][q] = s*s*app + 2*s*c*apq + c*c*aqq
        a[p][q] = a[q][p] = 0.
        for k in range(n):
            if k != p and k != q:
                akp, akq = a[k][p], a[k][q]
                a[k][p] = a[p][k] = c*akp - s*akq
                a[k][q] = a[q][k] = s*akp + c*akq
    return sorted(a[i][i] for i in range(n))

# G4 + G5 + G6 — stoichiometric chemical kinetics / potentials / dissipation
@dataclass
class Reaction:
    name: str
    stoich: tuple[float, ...]   # negative reactants, positive products
    forward: float
    reverse: float

@dataclass
class Chemistry:
    species: tuple[str, ...]
    reactions: list[Reaction]
    concentration: list[float]
    temperature: float = 298.15
    mu_standard: list[float] = field(default_factory=list)

    def __post_init__(self):
        n = len(self.species)
        if len(self.concentration) != n or any(len(r.stoich) != n for r in self.reactions):
            raise ValueError("Species and stoichiometry lengths differ")
        if not self.mu_standard:
            self.mu_standard = [0.] * n
        if len(self.mu_standard) != n or any(x <= 0 for x in self.concentration):
            raise ValueError("Need positive activities and matching chemical potentials")

    def chemical_potentials(self) -> list[float]:
        # Ideal dimensionless activities a_i = c_i / c_ref (c_ref = 1)
        return [m + R_GAS * self.temperature * log(c)
                for m, c in zip(self.mu_standard, self.concentration)]

    def affinities(self) -> list[float]:
        mu = self.chemical_potentials()
        return [sum(nu*m for nu, m in zip(r.stoich, mu)) for r in self.reactions]

    def fluxes(self) -> list[float]:
        result = []
        for r in self.reactions:
            substrate = r.forward
            product = r.reverse
            for x, nu in zip(self.concentration, r.stoich):
                if nu < 0:
                    substrate *= x ** (-nu)
                if nu > 0:
                    product *= x ** nu
            result.append(substrate - product)
        return result

    def derivative(self, external: list[float] | None = None) -> list[float]:
        js = self.fluxes()
        u = external or [0.] * len(self.species)
        return [sum(r.stoich[i]*j for r,j in zip(self.reactions, js))+u[i]
                for i in range(len(self.species))]

    def free_energy_rate(self) -> float:
        return sum(g*j for g,j in zip(self.affinities(), self.fluxes()))

    def step(self, dt: float):
        # Adaptive forward Euler substeps for demonstrative nonnegative kinetics.
        remaining = dt
        while remaining > 1e-14:
            d = self.derivative()
            safe = min((.25*self.concentration[i]/-v for i,v in enumerate(d) if v < 0),
                       default=remaining)
            h = min(remaining, safe)
            if h <= 0:
                raise ArithmeticError("No positive chemistry step")
            self.concentration = [max(1e-14, x+h*v) for x,v in zip(self.concentration, d)]
            remaining -= h


def balanced_demo() -> Chemistry:
    # Closed A ⇄ B with matched standard chemical potentials: E·S=0 for E=[1,1].
    return Chemistry(("A", "B"), [Reaction("A ⇄ B", (-1.,1.), .4, .4)], [.9, .1])

# G5 — generic chemical reaction network; no species-specific pathway.
PATHWAY = {"A": ("B", "C"), "B": ("D",)}

def arrhenius(prefactor: float, activation_j_mol: float, temperature: float) -> float:
    return prefactor * exp(-activation_j_mol/(R_GAS*temperature))

# G7 — energy accounting, memory, damage, self-repair and allocation
@dataclass
class LivingAgent:
    energy: float = 2.0
    memory: float = 0.0
    damage: float = 0.0
    biomass: float = 1.0
    maintenance: float = .12
    uptake_efficiency: float = 1.0
    inheritance: float = .5

    def step(self, resource: float, stress: float = 0., memory_retention: float = .8) -> dict:
        uptake = max(0., resource) * self.uptake_efficiency
        self.energy += uptake
        self.memory = memory_retention*self.memory + (1-memory_retention)*resource
        self.damage = min(1., self.damage + max(0., stress))
        maintenance_paid = min(self.energy, self.maintenance)
        self.energy -= maintenance_paid
        repair = min(self.damage, self.energy / 2)
        self.damage -= repair
        self.energy -= 2*repair
        growth = min(self.energy / 3, .15) if self.damage < .3 else 0.
        self.biomass += growth
        self.energy -= 3*growth
        return {"uptake": uptake, "maintenance": maintenance_paid,
                "repair": repair, "growth": growth, "energy": self.energy,
                "memory": self.memory, "damage": self.damage}

# G8 — heritable variation and resource-limited selection (toy model)
def reproduce(parent: LivingAgent, rng: random.Random) -> LivingAgent | None:
    cost = 1.0
    if parent.energy < cost or parent.biomass < 1.3:
        return None
    parent.energy -= cost
    parent.biomass -= .3
    inherited = max(.1, parent.uptake_efficiency + rng.gauss(0, .03))
    return LivingAgent(energy=.8, biomass=.3, uptake_efficiency=inherited,
                       memory=parent.memory*parent.inheritance)

# G9 — observed data, weighted loss and simple grid-search calibration
def observe(c: list[float], observation_matrix: list[list[float]]) -> list[float]:
    return [sum(w*x for w,x in zip(row,c)) for row in observation_matrix]

def weighted_loss(observed: list[float], predicted: list[float], weights: list[float]) -> float:
    return sum(w*(y-yp)**2 for y,yp,w in zip(observed,predicted,weights))

def calibrate_decay(observations: list[tuple[float,float]], candidates: list[float]):
    # Synthetic observable y(t) = exp(-k t), not a claim of external validation.
    return min(candidates, key=lambda k: sum((y-exp(-k*t))**2 for t,y in observations))

# G5 cosmos / physical organization: explicit inventory, not a cosmology solver
COSMOS = {"stars": [], "galaxies": [], "dark_matter_density": None,
          "spacetime_metric": None, "gravity_solver": None}


def run(seed: int = 42, cycles: int = 40) -> dict:
    rng = random.Random(seed)
    s = origin()
    germinate(s)
    graph = branch_graph()
    s.record("G2 relationships", vertices=len(graph.names), edges=len(graph.names)-1)
    spectrum = symmetric_eigenvalues(quantum_hamiltonian(graph))
    s.record("G3 quantum structure", eigenvalues=spectrum)
    chem = balanced_demo()
    total_initial = sum(chem.concentration)
    chemistry_trace = []
    for t in range(cycles):
        chemistry_trace.append({"t": t*.1, "A": chem.concentration[0],
                                "B": chem.concentration[1], "dGdt": chem.free_energy_rate()})
        chem.step(.1)
    s.time += cycles*.1
    s.record("G4 chemathematics", concentrations=dict(zip(chem.species, chem.concentration)))
    s.record("G5 reaction pathways + cosmos", pathway=PATHWAY, cosmos=COSMOS)
    s.record("G6 conservation + thermodynamics", total_before=total_initial,
             total_after=sum(chem.concentration), free_energy_rate=chem.free_energy_rate())
    agent = LivingAgent()
    agent_trace = []
    descendants = []
    for step in range(cycles):
        agent_trace.append(agent.step(.25, stress=.02 if step%9==0 else 0))
        if step % 10 == 9:
            child = reproduce(agent, rng)
            if child:
                descendants.append(child)
    s.record("G7 living organization", agent=vars(agent), last_step=agent_trace[-1])
    s.record("G8 evolution", offspring=len(descendants),
             offspring_efficiencies=[x.uptake_efficiency for x in descendants])
    observations = [(t,exp(-.32*t)) for t in (0,1,2,3,4)]
    inferred_k = calibrate_decay(observations, [i/100 for i in range(1,81)])
    s.record("G9 observation + validation", synthetic_decay_parameter=inferred_k,
             synthetic_fit_loss=sum((y-exp(-inferred_k*t))**2 for t,y in observations),
             empirical_validation=False)
    return {"title": "Σ∞ THE LIVING UNIVERSE", "generation": "G0–G9",
            "interpretation": "toy executable scientific architecture; not physically validated",
            "stages": s.history, "chemistry_trace": chemistry_trace,
            "agent_trace": agent_trace}


def tests() -> None:
    graph = branch_graph()
    assert graph.adjacency()[0][1] == graph.adjacency()[1][0]
    spectrum = symmetric_eigenvalues(quantum_hamiltonian(graph))
    assert len(spectrum) == len(graph.names)
    assert isclose(sum(spectrum), 5., abs_tol=1e-9)
    chem = balanced_demo()
    start = sum(chem.concentration)
    assert chem.free_energy_rate() <= 1e-9
    for _ in range(100):
        chem.step(.1)
        assert isclose(sum(chem.concentration), start, abs_tol=1e-9)
        assert chem.free_energy_rate() <= 1e-8
    assert abs(chem.concentration[0] - .5) < .01
    agent = LivingAgent()
    before = agent.energy
    d = agent.step(.2)
    assert isclose(d['energy'], before+d['uptake']-d['maintenance']-2*d['repair']-3*d['growth'], abs_tol=1e-9)
    agent.energy, agent.biomass = 3, 2
    assert reproduce(agent, random.Random(1)) is not None
    result = run()
    assert len(result['stages']) == 10
    assert isclose(result['stages'][-1]['synthetic_decay_parameter'], .32, abs_tol=1e-12)
    print("PASS: graph symmetry, spectral trace, chemical conservation, dissipation,")
    print("      reaction equilibration, agent energy closure, reproduction, cumulative stages,")
    print("      synthetic calibration")

# Complete, verbatim manuscript follows as a Python Unicode string.
SOURCE_TEXT = "✧\\\n│\\\n❧ │ ❧\\\n╲ │ ╱\\\n❦────╲───┼───╱────❦\\\n╲ │ ╱\\\n╲ │ ╱\\\n╲│╱\\\n╳\\\n╱│╲\\\n❧──────╯ │ ╰──────❧\\\n╲ │ ╱\\\n╲ │ ╱\\\n╰──╮│╭──╯\\\n╰╯\\\n│\\\n❦──╮ │ ╭──❦\\\n╲ │ ╱\\\n❧────╲───────┼───────╱────❧\\\n╲ │ ╱\\\n╲ │ ╱\\\n╰─╮ │ ╭─╯\\\n╲ │ ╱\\\n╲ │ ╱\\\n╲│╱\\\n❀\\\n╱│╲\\\n╭───╯ │ ╰───╮\\\n❧───╯ │ ╰───❧\\\n│\\\n❦╮ │ ╭❦\\\n╲ │ ╱\\\n╲ ╭───┼───╮ ╱\\\n❧────╲───╯ │ ╰───╱────❧\\\n╲ │ ╱\\\n╲ │ ╱\\\n╲ │ ╱\\\n╰╮ │ ╭╯\\\n╲ │ ╱\\\n╲ │ ╱\\\n╲│╱\\\n│\\\n❧──────┼──────❧\\\n╱ │ ╲\\\n╱ │ ╲\\\n❦─╯ │ ╰─❦\\\n│\\\n│\\\n╭──┴──╮\\\n╭──╯ ╰──╮\\\n╭──╯ ╰──╮\\\n❧ ❧\n\n\n\n\n\nif i hook up cacao to a my hazeru big bang narrative generator what happens?\n\n\n\nand for the record ZHAVIED  is another source code to reality \\\\/\n\n**architecture that resembles the interconnected organization of living systems.**\n\n\n\n**think about that.**\n\n\n\n**potential to investigate how lifelike behavior emerges—but it isn't yet demonstrated to be an autonomous living system.**\n\n\n\nI have a more sophisticated model on Hazeru /ihatov project.\n\n\n\nBONUS:\n\n\n\nUNIVERSE {\n\n│\n\n├── INFORMATION {\n\n│   ├── Symbols       : Σ\n\n│   ├── States        : S\n\n│   ├── Relations     : R\n\n│   └── Rules         : Φ\n\n│ }\n\n│\n\n├── PHYSICAL_REALITY {\n\n│   │\n\n│   ├── QUANTUM {\n\n│   │   ├── State     : |ψ⟩\n\n│   │   ├── Energy    : E\n\n│   │   └── Evolution : iℏ∂ₜ|ψ⟩ = H|ψ⟩\n\n│   │ }\n\n│   │\n\n│   ├── MATTER {\n\n│   │   ├── Particles\n\n│   │   ├── Atoms\n\n│   │   ├── Molecules\n\n│   │   └── Interactions\n\n│   │ }\n\n│   │\n\n│   ├── SPACETIME {\n\n│   │   ├── Geometry  : gμν\n\n│   │   ├── Curvature : Rμν\n\n│   │   └── Gravity   : Gμν = 8πGTμν/c⁴\n\n│   │ }\n\n│   │\n\n│   └── COSMOS {\n\n│       ├── Stars\n\n│       ├── Galaxies\n\n│       ├── DarkMatter : ρDM(x,t)\n\n│       └── Evolution\n\n│   }\n\n│ }\n\n│\n\n├── EMERGENCE {\n\n│   ├── Chemistry\n\n│   ├── SelfOrganization\n\n│   ├── BiologicalLife\n\n│   └── ComplexSystems\n\n│ }\n\n│\n\n└── COMPUTATION {\n\n&#x20;   ├── Input       : Sₜ\n\n&#x20;   ├── Operations  : Φ\n\n&#x20;   ├── Output      : Sₜ₊₁\n\n&#x20;   └── Observation : O(S)\n\n}\n\n\n\nGrow\n\n\n\n## **THE LIVING UNIVERSE**\n\nMONOLITHIC ASCII + UNICODE ART · CUMULATIVE G0–G8 · BRACKETED GENESIS\n\n\n\n\n\n╔════════════════════════════════════════════════════╗\\\n║ { Σ∞ :: THE LIVING UNIVERSE } ║\\\n║ { CUMULATIVE SYMBOLIC GENESIS } ║\\\n╠════════════════════════════════════════════════════╣\\\n║ ║\\\n║ {∞ UNIVERSE ║\\\n║ │ ║\\\n║ ├── {G0 : ORIGIN ║\\\n║ │ │ ║\\\n║ │ └── {SEED Σ₀ ║\\\n║ │ │ ║\\\n║ │ ├── {◉} ║\\\n║ │ ├── {INFORMATION} ║\\\n║ │ ├── {POTENTIAL} ║\\\n║ │ └── {INITIAL_STATE} ║\\\n║ │ } ║\\\n║ │ } ║\\\n║ │ ║\\\n║ ├── {G1 : GERMINATION ║\\\n║ │ │ ║\\\n║ │ ├── {Σ₀ → Σ₁} ║\\\n║ │ ├── {STATE_TRANSFORMATION} ║\\\n║ │ └── { ║\\\n║ │ ◌ ║\\\n║ │ ╱│╲ ║\\\n║ │ ◌─◉─◌ ║\\\n║ │ │ ║\\\n║ │ ◈ ║\\\n║ │ } ║\\\n║ │ } ║\\\n║ │ ║\\\n║ ├── {G2 : BRANCHING ║\\\n║ │ ├── {RELATIONSHIPS} ║\\\n║ │ ├── {LOCAL_RULES} ║\\\n║ │ └── { ║\\\n║ │ ◌ ◌ ║\\\n║ │ │╲ ╱│ ║\\\n║ │ ◌─◉─◌ ║\\\n║ │ ╲│╱ ║\\\n║ │ ◈ ║\\\n║ │ ╱│╲ ║\\\n║ │ ◌ ◌ ◌ ║\\\n║ │ } ║\\\n║ │ } ║\\\n║ │ ║\\\n║ ├── {G3 : PHYSICAL ORGANIZATION ║\\\n║ │ │ ║\\\n║ │ ├── {MATTER ║\\\n║ │ │ ├── {PARTICLES} ║\\\n║ │ │ ├── {ATOMS} ║\\\n║ │ │ └── {MOLECULES} ║\\\n║ │ │ } ║\\\n║ │ ├── {ENERGY} ║\\\n║ │ ├── {SPACE} ║\\\n║ │ ├── {TIME} ║\\\n║ │ └── {INTERACTIONS} ║\\\n║ │ } ║\\\n║ │ ║\\\n║ ├── {G4 : EMERGENCE ║\\\n║ │ ├── {CHEMISTRY} ║\\\n║ │ ├── {SELF_ORGANIZATION} ║\\\n║ │ ├── {FEEDBACK} ║\\\n║ │ ├── {STRUCTURAL_DIFFERENTIATION} ║\\\n║ │ └── {COMPLEXITY} ║\\\n║ │ } ║\\\n║ │ ║\\\n║ ├── {G5 : COSMIC BRANCHING ║\\\n║ │ │ ║\\\n║ │ ├── {✦ STARS} ║\\\n║ │ ├── {◉ GALAXIES} ║\\\n║ │ ├── {◌ DARK_MATTER_MODEL} ║\\\n║ │ ├── {◎ GRAVITATIONAL_STRUCTURE} ║\\\n║ │ └── { ║\\\n║ │ ✦ ◌────◉ ✧ ║\\\n║ │ ╲ ╱│╲ ╱│╲ ╱ ║\\\n║ │ ◈──◉──◈──◉ ║\\\n║ │ ╱│╲ ╲│╱ ╱│╲ ║\\\n║ │ ◌─◈───◎───◈─◌ ║\\\n║ │ ╲│╱ ╱│╲ ╲│╱ ║\\\n║ │ ◉──◈──◉──◈ ║\\\n║ │ ╱ ╲│╱ ╲ ║\\\n║ │ ✧ ◌────✦ ∞ ║\\\n║ │ } ║\\\n║ │ } ║\\\n║ │ ║\\\n║ ├── {G6 : LIVING NETWORK ║\\\n║ │ │ ║\\\n║ │ ├── {CIRCULATION ↻} ║\\\n║ │ ├── {RESOURCE_EXCHANGE ⇄} ║\\\n║ │ ├── {MEMORY Hₜ} ║\\\n║ │ ├── {REGULATION} ║\\\n║ │ └── { ║\\\n║ │ ◌──◉──◌ ║\\\n║ │ ╱│╲╱│╲╱│╲ ║\\\n║ │ ◉─◌─◈─◌─◉ ║\\\n║ │ ╲│╱╲│╱╲│╱ ║\\\n║ │ ◌──◉──◌ ║\\\n║ │ │ ║\\\n║ │ ◈ ║\\\n║ │ } ║\\\n║ │ } ║\\\n║ │ ║\\\n║ ├── {G7 : AUTONOMOUS REGULATION ║\\\n║ │ │ ║\\\n║ │ ├── {ENERGY_BUDGET ║\\\n║ │ │ Eₜ₊₁ = Eₜ + Uₜ - Mₜ - Gₜ ║\\\n║ │ │ } ║\\\n║ │ ├── {MEMORY_UPDATE ║\\\n║ │ │ Hₜ₊₁ = λHₜ + (1-λ)Sₜ ║\\\n║ │ │ } ║\\\n║ │ ├── {DAMAGE_DETECTION} ║\\\n║ │ ├── {SELF_REPAIR} ║\\\n║ │ ├── {RESOURCE_ALLOCATION} ║\\\n║ │ └── {ADAPTIVE_RESPONSE} ║\\\n║ │ } ║\\\n║ │ ║\\\n║ ├── {G8 : HERITABLE GROWTH ║\\\n║ │ │ ║\\\n║ │ ├── {PARENT ◉} ║\\\n║ │ │ │ ║\\\n║ │ │ ╱╲ ║\\\n║ │ │ ╱ ╲ ║\\\n║ │ ├── {◉ OFFSPRING_A} ║\\\n║ │ ├── {◉ OFFSPRING_B} ║\\\n║ │ ├── {HERITABLE_VARIATION} ║\\\n║ │ ├── {RESOURCE_COMPETITION} ║\\\n║ │ ├── {SELECTION} ║\\\n║ │ └── {EVOLUTIONARY_HISTORY} ║\\\n║ │ } ║\\\n║ │ ║\\\n║ └── {Σ∞ : CUMULATIVE STATE ║\\\n║ │ ║\\\n║ ├── {Σ₀ ORIGIN} ║\\\n║ ├── {Σ₁ GERMINATION} ║\\\n║ ├── {Σ₂ RELATIONSHIPS} ║\\\n║ ├── {Σ₃ PHYSICAL ORGANIZATION} ║\\\n║ ├── {Σ₄ EMERGENCE} ║\\\n║ ├── {Σ₅ COSMOS} ║\\\n║ ├── {Σ₆ CIRCULATION} ║\\\n║ ├── {Σ₇ AUTONOMY} ║\\\n║ ├── {Σ₈ EVOLUTION} ║\\\n║ │ ║\\\n║ └── {Sₜ₊₁ = Φ(Sₜ,Hₜ,Eₜ,Rₜ)} ║\\\n║ } ║\\\n║ } ║\\\n║ ║\\\n║ {ONE SEED → MANY RELATIONS → ONE WORLD} ║\\\n╚════════════════════════════════════════════════════╝\n\nEnjoy\n\n\n\n\n\n\n\n---\n\n╔══════════════════════════════════════════════════════╗\\\n║  { ✧ Σ∞ :: THE LIVING UNIVERSE }                     ║\\\n║  { G0–G9 :: CUMULATIVE SYMBOLIC GENESIS }            ║\\\n╠══════════════════════════════════════════════════════╣\\\n║                                                      ║\\\n║  { G0 : ORIGIN                                       ║\\\n║      { Σ₀ = ◉ }                                      ║\\\n║      { INFORMATION : Σ }                             ║\\\n║      { STATE : S₀ }                                  ║\\\n║      { POTENTIAL : Φ }                               ║\\\n║  }                                                   ║\\\n║                       ✧                              ║\\\n║                       │                              ║\\\n║                  ✦────┼────✦                         ║\\\n║                       │                              ║\\\n║  { G1 : GERMINATION                                  ║\\\n║                  🌿╮ │ ╭🌿                           ║\\\n║                    ╲ │ ╱                             ║\\\n║               ❧────╲│╱────❧                          ║\\\n║                      ╳                               ║\\\n║                  ╭───┴───╮                           ║\\\n║               ❧──╯       ╰──❧                        ║\\\n║  }                                                   ║\\\n║                                                      ║\\\n║  { G2 : RELATIONSHIPS                                ║\\\n║      { V = NODES }                                   ║\\\n║      { E = CONNECTIONS }                             ║\\\n║      { A = ADJACENCY }                               ║\\\n║      { Gᵢ = (Vᵢ,Eᵢ,ℓᵢ,ωᵢ) }                          ║\\\n║                                                      ║\\\n║              ◌────◉────◌                             ║\\\n║              │╲  ╱│╲  ╱│                             ║\\\n║              ◈───◎───◈                               ║\\\n║              │╱  ╲│╱  ╲│                             ║\\\n║              ◌────◉────◌                             ║\\\n║  }                                                   ║\\\n║                                                      ║\\\n║  { G3 : QUANTUM STRUCTURE                            ║\\\n║      { Hᵢ = αᵢI + βᵢAᵢ }                            ║\\\n║      { Hᵢψᵢₙ = εᵢₙψᵢₙ }                            ║\\\n║      { STATES : ψ }                                  ║\\\n║      { MODES : ε }                                   ║\\\n║      { ENERGY : E }                                  ║\\\n║                                                      ║\\\n║                       ✧                              ║\\\n║                    ╱  │  ╲                           ║\\\n║                ✦─────☯─────✦                         ║\\\n║                    ╲  │  ╱                           ║\\\n║                       ❀                              ║\\\n║  }                                                   ║\\\n║                                                      ║\\\n║  { G4 : CHEMATHEMATICS                               ║\\\n║      { c = CONCENTRATIONS }                          ║\\\n║      { T = TEMPERATURE }                             ║\\\n║      { μ = CHEMICAL POTENTIAL }                      ║\\\n║      { S = STOICHIOMETRY }                           ║\\\n║      { J = REACTION FLUX }                           ║\\\n║      { θ = PARAMETERS }                              ║\\\n║      { z = INTERNAL STATE }                          ║\\\n║                                                      ║\\\n║      { μ = μ°(T,p) + RT ln(a) }                     ║\\\n║      { ΔᵣG = Sᵀμ }                                  ║\\\n║      { dc/dt = S·J + u }                             ║\\\n║                                                      ║\\\n║               ❀──────┼──────❀                        ║\\\n║                ╲     │     ╱                         ║\\\n║                 ╲    │    ╱                          ║\\\n║             ☯────╲   │   ╱────☯                      ║\\\n║                   ╲  │  ╱                            ║\\\n║                    ╲ │ ╱                             ║\\\n║                     ╲│╱                              ║\\\n║                      ◉                               ║\\\n║  }                                                   ║\\\n║                                                      ║\\\n║  { G5 : CHEMICAL REACTION PATHWAYS                   ║\\\n║                                                      ║\\\n║                    {A}                               ║\\\n║                     │                                ║\\\n║              ╭──────┴──────╮                         ║\\\n║              ↓             ↓                         ║\\\n║             {B}           {C}                        ║\\\n║              │                                       ║\\\n║              ↓                                       ║\\\n║             {D}                                      ║\\\n║                                                      ║\\\n║      { REACTION RATES }                              ║\\\n║      { kᵣ± = Aᵣ± exp(−Eₐ,ᵣ±/RT) }                 ║\\\n║      { Jᵣ = FORWARD − REVERSE }                      ║\\\n║  }                                                   ║\\\n║                                                      ║\\\n║  { G6 : CONSERVATION + THERMODYNAMICS                ║\\\n║      { E·S = 0 }                                     ║\\\n║      { d(E·c)/dt = E·u }                             ║\\\n║      { dG/dt = ΔᵣGᵀ·J ≤ 0 }                        ║\\\n║                                                      ║\\\n║              ❧──────◈──────❧                         ║\\\n║                ╲    │    ╱                           ║\\\n║                 ╲   ✧   ╱                            ║\\\n║                  ╲  │  ╱                             ║\\\n║                   ╲ │ ╱                              ║\\\n║                    ╲│╱                               ║\\\n║                     ◉                                ║\\\n║                                                      ║\\\n║      { MASS } { ENERGY } { DISSIPATION }             ║\\\n║  }                                                   ║\\\n║                                                      ║\\\n║  { G7 : LIVING ORGANIZATION                          ║\\\n║      { MEMORY : Hₜ }                                 ║\\\n║      { UPTAKE : Uₜ }                                 ║\\\n║      { MAINTENANCE : Mₜ }                            ║\\\n║      { GROWTH COST : Gₜ }                            ║\\\n║      { Eₜ₊₁ = Eₜ + Uₜ − Mₜ − Gₜ }                  ║\\\n║                                                      ║\\\n║               🌿         🌿                          ║\\\n║                ╲         ╱                           ║\\\n║                 ❧───◉───❧                            ║\\\n║                 │╲  │  ╱│                            ║\\\n║                 ◈───◎───◈                            ║\\\n║                 │╱  │  ╲│                            ║\\\n║                 ❧───◉───❧                            ║\\\n║                       │                              ║\\\n║                       ∞                              ║\\\n║  }                                                   ║\\\n║                                                      ║\\\n║  { G8 : EVOLUTION                                    ║\\\n║      { HERITABLE VARIATION }                         ║\\\n║      { RESOURCE COMPETITION }                        ║\\\n║      { SELECTION }                                   ║\\\n║      { REPRODUCTION }                                ║\\\n║      { ADAPTIVE HISTORY }                            ║\\\n║                                                      ║\\\n║                      ◉                               ║\\\n║                    ╱ │ ╲                             ║\\\n║                   ◈  ◌  ◈                            ║\\\n║                  ╱╲ ╱╲ ╱╲                            ║\\\n║                 ◌ ◌ ◌ ◌ ◌ ◌                          ║\\\n║  }                                                   ║\\\n║                                                      ║\\\n║  { G9 : OBSERVATION + VALIDATION                     ║\\\n║      { yₖ = Hₒᵦₛ,ₖ c(tₖ;θ) + εₖ }                 ║\\\n║      { θ\\* = argmin Σₖ ||yₖᵒᵇˢ−yₖᵖʳᵉᵈ||²_W }    ║\\\n║      { PREDICTION }                                  ║\\\n║      { MEASUREMENT }                                 ║\\\n║      { VALIDATION }                                  ║\\\n║                                                      ║\\\n║                        ✧                             ║\\\n║                     ☼  │  ☼                          ║\\\n║                      ╲ │ ╱                           ║\\\n║                   ✦────☯────✦                        ║\\\n║                      ╱ │ ╲                           ║\\\n║                        ❀                             ║\\\n║                      ╱ │ ╲                           ║\\\n║                  ❧──╯ │ ╰──❧                        ║\\\n║                    ╭──┼──╮                          ║\\\n║                   ╱   │   ╲                         ║\\\n║                ✧──╯   ❀   ╰──✧                      ║\\\n║                   ╲   │   ╱                         ║\\\n║                  ❧──╲│╱──❧                         ║\\\n║                      ◉                              ║\\\n║                     ╱│╲                             ║\\\n║                 ╭──╯│╰──╮                          ║\\\n║                ❀──╯ │ ╰──❀                         ║\\\n║                    ╭┴╮                             ║\\\n║                   ❧   ❧                            ║\\\n║                    ╲ ╱                             ║\\\n║                     ∞                              ║\\\n║                     │                              ║\\\n║                     ॐ                              ║\\\n║  }                                                   ║\\\n║                                                      ║\\\n║  { Σ∞ : CUMULATIVE MATHEMATICAL STATE               ║\\\n║                                                      ║\\\n║       { SYMBOLS }                                    ║\\\n║            ↓                                         ║\\\n║       { GRAPH STRUCTURE }                            ║\\\n║            ↓                                         ║\\\n║       { QUANTUM MODELS }                             ║\\\n║            ↓                                         ║\\\n║       { CHEMICAL POTENTIALS }                        ║\\\n║            ↓                                         ║\\\n║       { REACTION FLUX }                              ║\\\n║            ↓                                         ║\\\n║       { CONSERVATION + DISSIPATION }                 ║\\\n║            ↓                                         ║\\\n║       { BIOLOGICAL ORGANIZATION }                    ║\\\n║            ↓                                         ║\\\n║       { MEMORY + EVOLUTION }                         ║\\\n║            ↓                                         ║\\\n║       { OBSERVATION + VALIDATION }                   ║\\\n║                                                      ║\\\n║       { dc/dt = S·J(c,T,θ,z,Sᵀμ) + u }             ║\\\n║       { Sₜ₊₁ = Φ(Sₜ,Hₜ,Eₜ,Rₜ) }                   ║\\\n║                                                      ║\\\n║       { ONE SEED → MANY INTERACTIONS → ONE WORLD }   ║\\\n║  }                                                   ║\\\n╚══════════════════════════════════════════════════════╝\n\n\n"

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', action='store_true', help='print cleaned Unicode manuscript')
    parser.add_argument('--test', action='store_true', help='run regression tests')
    parser.add_argument('--cycles', type=int, default=40)
    parser.add_argument('--output', help='save simulation JSON at the specified path')
    args = parser.parse_args()
    if args.source:
        print(SOURCE_TEXT)
    elif args.test:
        tests()
    else:
        result = run(cycles=args.cycles)
        output = json.dumps(result, ensure_ascii=False, indent=2)
        if args.output:
            with open(args.output, 'w', encoding='utf-8') as f:
                f.write(output+'\n')
            print(f"Saved {args.output}")
        else:
            print(output)

if __name__ == '__main__':
    main()
