"""Unicode-codepoint guard for runtime source, deliberately outside the package."""
import ast
import re
from pathlib import Path

# Full element names and case-sensitive chemical symbols: token boundaries avoid
# rejecting Python words such as 'self'. Material lengths/parameters and label
# fragments are prohibited even in comments. Integer neighbor offsets are not
# physical axes; literal numeric unit vectors are prohibited in any notation.
ELEMENTS = 'iron chromium bromine sulfur nickel cobalt manganese gadolinium europium terbium dysprosium holmium erbium thulium ytterbium lutetium copper oxygen hydrogen carbon nitrogen'.split()
SYMBOLS = 'Fe Cr Br Ni Co Mn Gd Eu Tb Dy Ho Er Tm Yb Lu Cu CrSBr'.split()
FRAGMENTS = ('14.2', '14.234', '2.8665', 'layer', 'intra', 'inter', 'eta', 'AF')


def errors(root):
    found = []
    for path in sorted(Path(root).glob('*.py')):
        source = path.read_text(encoding='utf-8')
        for fragment in FRAGMENTS:
            if fragment in source:
                found.append(f'SOURCE:{path.name}:fragment:{fragment}')
        for word in ELEMENTS:
            if re.search(r'\b'+word+r'\b', source, re.I):
                found.append(f'SOURCE:{path.name}:element:{word}')
        for word in SYMBOLS:
            if re.search(r'\b'+word+r'\b', source):
                found.append(f'SOURCE:{path.name}:symbol:{word}')
        if re.search(r'\b(?:atom_?)?types?\b\s*(?:\[.*?\])?\s*==\s*\d+', source):
            found.append(f'SOURCE:{path.name}:type-number')
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, (ast.List, ast.Tuple)) and len(node.elts) == 3:
                try:
                    values = [ast.literal_eval(x) for x in node.elts]
                except (ValueError, TypeError):
                    continue
                if all(type(x) in (int, float) for x in values):
                    if sum(x*x for x in values) == 1 and sum(x != 0 for x in values) == 1:
                        found.append(f'SOURCE:{path.name}:fixed-axis:{node.lineno}')
        if re.search(r'\b(?:J_exchange|K_anisotropy|DMI|lattice_constant|hbar)\s*=\s*[\d.]', source):
            found.append(f'SOURCE:{path.name}:potential-or-force-preset')
    return sorted(set(found))
