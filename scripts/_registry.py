"""Load the member registry (.github/members.yml) — the single source of truth.

Shared by scripts/new_contribution.py and scripts/check_structure.py so there is
exactly one place that parses the registry.
"""

import os

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MEMBERS_YML = os.path.join(REPO_ROOT, ".github", "members.yml")


def load_members():
    """Return {github_handle: display_name} from .github/members.yml.

    Commented reserved slots (member_slot_4/5) are YAML comments and are ignored.
    """
    with open(MEMBERS_YML, encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    members = data.get("members")
    if not isinstance(members, dict) or not members:
        raise SystemExit(
            "ERROR: .github/members.yml must contain a non-empty `members:` mapping "
            "of github_handle -> display_name."
        )
    return {str(k): str(v) for k, v in members.items()}
