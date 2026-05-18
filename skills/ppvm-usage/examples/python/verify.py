"""
Skill quick-verification snippet.

Mirrors the "Verifying you got the API right" block in
``skills/ppvm-usage/SKILL.md``. A failing run means the install or the
method names are off; the rest of the skill is wasted effort until this
passes.
"""

from ppvm import PauliSum

ps = PauliSum.new(2, "ZZ")
ps.cnot(0, 1)
ps.h(0)
assert ps.overlap_with_zero() == 1.0, ps.overlap_with_zero()
print("ok")
