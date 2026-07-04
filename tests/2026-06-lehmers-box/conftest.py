"""Make the local _helpers module importable under pytest importlib mode."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
