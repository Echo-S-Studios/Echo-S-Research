import os
import sys

# Ensure this test folder is importable (pytest importlib mode does not add it to
# sys.path automatically), so `import emgap_util` resolves.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
