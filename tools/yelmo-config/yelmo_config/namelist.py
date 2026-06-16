"""A small, comment-preserving Fortran namelist parser/writer.

`f90nml` is deliberately *not* used: the inline comments in
``yelmo_defaults.nml`` are the parameter documentation, and we need to keep
section sub-headers (``! === Grounding line ===``) and ordering intact so the
``write`` command can emit a clean, documented, complete par file.

The model is intentionally line-oriented so a file round-trips faithfully:

    Namelist
      preamble : list[str]          # comment/blank lines before the first group
      groups   : "ordered" dict[str, Group]

    Group
      name  : str
      items : list[Param | Comment] # interleaved, in source order

Values are kept verbatim as ``raw`` strings (so ``1e3`` is never rewritten as
``1000.0``); :func:`normalize` provides a typed, comparable view for diff/check.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


# --------------------------------------------------------------------------- #
# Data model
# --------------------------------------------------------------------------- #
@dataclass
class Param:
    name: str
    raw: str                    # value exactly as written, e.g. '"impl-lis"' or '1e3'
    comment: str | None = None  # inline doc text (without the leading '!'), or None

    @property
    def value(self):
        return normalize(self.raw)


@dataclass
class Comment:
    """A standalone comment or blank line inside a group (preserved verbatim)."""

    text: str   # full line text without trailing newline ('' for a blank line)


@dataclass
class Group:
    name: str
    items: list = field(default_factory=list)   # list[Param | Comment]

    def params(self) -> "dict[str, Param]":
        return {it.name: it for it in self.items if isinstance(it, Param)}

    def get(self, name: str) -> Param | None:
        for it in self.items:
            if isinstance(it, Param) and it.name == name:
                return it
        return None

    def set_value(self, name: str, raw: str) -> None:
        p = self.get(name)
        if p is not None:
            p.raw = raw
        else:
            self.items.append(Param(name=name, raw=raw))


@dataclass
class Namelist:
    groups: dict = field(default_factory=dict)   # ordered dict[str, Group]
    preamble: list = field(default_factory=list)  # list[str]

    def __contains__(self, group: str) -> bool:
        return group in self.groups

    def __getitem__(self, group: str) -> Group:
        return self.groups[group]

    def get(self, group: str) -> Group | None:
        return self.groups.get(group)


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
_GROUP_RE = re.compile(r"^\s*&(\w+)")
_END_RE = re.compile(r"^\s*/\s*$")
_PARAM_RE = re.compile(r"^\s*([A-Za-z]\w*)\s*=\s*(.*)$")


def _split_comment(text: str) -> tuple[str, str | None]:
    """Split a line at the first '!' that is not inside a quoted string."""
    in_s = in_d = False
    for i, ch in enumerate(text):
        if ch == "'" and not in_d:
            in_s = not in_s
        elif ch == '"' and not in_s:
            in_d = not in_d
        elif ch == "!" and not in_s and not in_d:
            return text[:i], text[i + 1:]
    return text, None


def parse_string(text: str) -> Namelist:
    nml = Namelist()
    current: Group | None = None
    seen_group = False
    pending: list[str] = []   # comment/blank lines sitting between groups
    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")

        if current is None:
            m = _GROUP_RE.match(line)
            if m:
                current = Group(name=m.group(1))
                nml.groups[current.name] = current
                seen_group = True
                # carry over any non-blank comment lines as group lead-in
                for pl in pending:
                    if pl.strip():
                        current.items.append(Comment(text=pl))
                pending = []
            elif not seen_group:
                nml.preamble.append(line)
            else:
                pending.append(line)   # between groups; re-emitted by dump
            continue

        if _END_RE.match(line):
            current = None
            continue

        code, comment = _split_comment(line)
        m = _PARAM_RE.match(code)
        if m and code.strip():
            name = m.group(1)
            value = m.group(2).strip().rstrip(",").strip()
            current.items.append(
                Param(name=name, raw=value,
                      comment=comment.strip() if comment is not None else None)
            )
        else:
            # blank line or standalone comment (e.g. a section sub-header)
            current.items.append(Comment(text=line))
    return nml


def parse_file(path: str | Path) -> Namelist:
    return parse_string(Path(path).read_text())


# --------------------------------------------------------------------------- #
# Value normalization (typed, comparable view)
# --------------------------------------------------------------------------- #
_BOOLS_TRUE = {"true", ".true.", "t"}
_BOOLS_FALSE = {"false", ".false.", "f"}


def _token_value(tok: str):
    tok = tok.strip()
    if not tok:
        return ""
    if (tok[0] == '"' and tok[-1] == '"') or (tok[0] == "'" and tok[-1] == "'"):
        return tok[1:-1]
    low = tok.lower()
    if low in _BOOLS_TRUE:
        return True
    if low in _BOOLS_FALSE:
        return False
    # Fortran double-precision exponent 'd' -> 'e'
    num = tok.replace("D", "e").replace("d", "e")
    try:
        if re.fullmatch(r"[+-]?\d+", num):
            return int(num)
        return float(num)
    except ValueError:
        return tok  # unquoted bareword


def _split_tokens(raw: str) -> list[str]:
    """Split a namelist value into tokens, respecting quotes.

    Handles both comma-separated (``1.0, 2.0``) and space-separated
    (``"a" "b"``) Fortran array syntaxes.
    """
    toks, buf = [], []
    in_s = in_d = False
    for ch in raw:
        if ch == "'" and not in_d:
            in_s = not in_s
            buf.append(ch)
        elif ch == '"' and not in_s:
            in_d = not in_d
            buf.append(ch)
        elif (ch == "," or ch.isspace()) and not in_s and not in_d:
            if buf:
                toks.append("".join(buf))
                buf = []
        else:
            buf.append(ch)
    if buf:
        toks.append("".join(buf))
    return toks


def normalize(raw: str):
    """Return a typed, comparable value for ``raw``.

    Scalars return a scalar; multi-element values return a tuple. This makes
    ``1e3`` == ``1000.0`` and ``"x"`` == ``x`` for comparison purposes.
    """
    toks = _split_tokens(raw)
    if len(toks) <= 1:
        return _token_value(raw)
    return tuple(_token_value(t) for t in toks)


def infer_type(raw: str) -> str:
    """Human-readable type label for a raw value."""
    toks = _split_tokens(raw)
    if len(toks) > 1:
        inner = {infer_type(t) for t in toks}
        return f"{inner.pop() if len(inner) == 1 else 'mixed'}[]"
    v = _token_value(raw)
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, int):
        return "int"
    if isinstance(v, float):
        return "real"
    if (raw[:1], raw[-1:]) in {('"', '"'), ("'", "'")}:
        return "string"
    return "string"


# --------------------------------------------------------------------------- #
# Writing
# --------------------------------------------------------------------------- #
def dump(nml: Namelist, *, align: bool = True) -> str:
    """Render a namelist back to text, preserving comments and ordering.

    Within each group, parameter names / ``=`` / inline comments are aligned to
    columns for readability.
    """
    out: list[str] = []
    out.extend(nml.preamble)
    if nml.preamble and nml.preamble[-1] != "":
        out.append("")

    for gi, group in enumerate(nml.groups.values()):
        if gi > 0:
            out.append("")
        out.append(f"&{group.name}")

        params = [it for it in group.items if isinstance(it, Param)]
        name_w = max((len(p.name) for p in params), default=0) if align else 0
        val_w = max((len(p.raw) for p in params if p.comment), default=0) if align else 0

        for it in group.items:
            if isinstance(it, Comment):
                out.append(it.text)
                continue
            line = f"    {it.name.ljust(name_w)} = "
            if it.comment is not None:
                line += f"{it.raw.ljust(val_w)}   ! {it.comment}"
            else:
                line += it.raw
            out.append(line.rstrip())
        out.append("/")

    return "\n".join(out) + "\n"
