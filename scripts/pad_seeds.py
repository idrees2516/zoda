#!/usr/bin/env python3
"""Pad `*b"..."` seed literals to exactly 32 bytes (developer helper)."""
import re, sys
for p in sys.argv[1:]:
    src = open(p).read()
    for m in re.finditer(r'\*b"([^"]+)"', src):
        s = m.group(1)
        if len(s) != 32:
            src = src.replace('*b"' + s + '"', '*b"' + (s + '0'*32)[:32] + '"')
    open(p, 'w').write(src)
    print("padded", p)
