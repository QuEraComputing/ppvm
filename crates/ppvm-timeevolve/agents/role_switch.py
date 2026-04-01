#!/usr/bin/env python3
"""
Auto role-switcher for the ppvm-timeevolve agent workflow.

Called by the Claude Code Stop hook. Reads the session transcript to find
the last assistant message, detects role-switch triggers, and updates
agents/state.md frontmatter accordingly.

Triggers (from CLAUDE.md):
  developer/perf-developer + handoff language  → reviewer (same task)
  reviewer + "Task N approved."                → perf-developer (task + 1)
  reviewer + feedback without approval         → perf-developer (same task)
"""
import json
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
STATE_FILE = PROJECT_ROOT / "agents" / "state.md"
PROJECTS_DIR = Path.home() / ".claude" / "projects"


def sanitize_path(p: Path) -> str:
    """Reproduce Claude Code's project-directory naming: replace / with -."""
    return str(p).replace("/", "-")


def read_state():
    text = STATE_FILE.read_text()
    role = re.search(r"active_role:\s*(\S+)", text)
    task = re.search(r"current_task:\s*(\d+)", text)
    return text, role.group(1) if role else None, int(task.group(1)) if task else None


def write_state(text: str, new_role: str, new_task: int) -> None:
    text = re.sub(r"(active_role:\s*)\S+", rf"\g<1>{new_role}", text)
    text = re.sub(r"(current_task:\s*)\d+", rf"\g<1>{new_task}", text)
    STATE_FILE.write_text(text)


def get_last_assistant_text(session_id: str) -> str:
    """
    Scan the session transcript JSONL and return the text of the last
    assistant message. Returns '' if the file is missing or unreadable.
    """
    sanitized = sanitize_path(PROJECT_ROOT)
    transcript = PROJECTS_DIR / sanitized / f"{session_id}.jsonl"

    if not transcript.exists():
        return ""

    last = ""
    with open(transcript) as f:
        for line in f:
            try:
                entry = json.loads(line)
                # Transcript lines have type="assistant"; the actual message
                # is nested under entry["message"]["content"].
                if entry.get("type") != "assistant":
                    continue
                content = entry.get("message", {}).get("content", "")
                if isinstance(content, list):
                    parts = [
                        c.get("text", "")
                        for c in content
                        if isinstance(c, dict) and c.get("type") == "text"
                    ]
                    last = " ".join(parts)
                elif isinstance(content, str):
                    last = content
            except Exception:
                pass
    return last


def detect_transition(role: str, task: int, text: str):
    """
    Return (new_role, new_task) if a trigger fired, else (role, task).
    """
    if role in ("developer", "perf-developer"):
        # Developer handoff: summary posted and review requested.
        if re.search(
            r"(hand.?off|please review|for the reviewer|ready for review|"
            r"reviewer to begin|review checklist|explicitly prompt the reviewer)",
            text,
            re.IGNORECASE,
        ):
            return "reviewer", task

    elif role == "reviewer":
        # Approval trigger (exact phrase from CLAUDE.md).
        if re.search(r"\bTask\s+\d+\s+approved\b", text, re.IGNORECASE):
            return "perf-developer", task + 1

        # Feedback-without-approval: reviewer returned changes to developer.
        if re.search(
            r"(please fix|fix forward|return.{0,20}task|resubmit|"
            r"checklist item fail|specific.{0,10}actionable|"
            r"return it to the developer)",
            text,
            re.IGNORECASE,
        ):
            return "perf-developer", task

    return role, task


def main() -> None:
    try:
        data = json.load(sys.stdin)
    except Exception:
        return

    session_id = data.get("session_id", "")
    if not session_id:
        return

    state_text, current_role, current_task = read_state()
    if current_role is None or current_task is None:
        return

    last_msg = get_last_assistant_text(session_id)
    if not last_msg:
        return

    new_role, new_task = detect_transition(current_role, current_task, last_msg)

    if new_role != current_role or new_task != current_task:
        write_state(state_text, new_role, new_task)


if __name__ == "__main__":
    main()
