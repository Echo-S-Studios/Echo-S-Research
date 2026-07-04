"""Make the shared helper module importable under pytest importlib mode.

pytest runs with --import-mode=importlib (root pytest.ini), which does not put
each test file's directory on sys.path.  This conftest, always collected for
this folder, adds the folder so `import _relcharge` works from every test file.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
