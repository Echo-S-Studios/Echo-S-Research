# Put this folder on sys.path so the sibling helper module `_opalg_ops`
# is importable under pytest's importlib mode. The helper carries a unique
# name (`_opalg_ops`) so it cannot collide with other agents' helpers.
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
