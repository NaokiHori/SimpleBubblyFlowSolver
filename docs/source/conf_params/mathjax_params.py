from mathjax_macro.general import add as add_general
from mathjax_macro.momentum import add as add_momentum
from mathjax_macro.energy import add as add_energy
from mathjax_macro.discrete import add as add_discrete

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

macros = dict()

add_general(macros)
add_momentum(macros)
add_energy(macros)
add_discrete(macros)

mathjax3_config = {"TeX": {"Macros": macros}}
