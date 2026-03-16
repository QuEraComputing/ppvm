#!/usr/bin/env bash
# Toggle Claude autopilot permissions for ppvm-timeevolve development.
# Usage: .claude/autopilot.sh [on|off|status]

DIR="$(cd "$(dirname "$0")" && pwd)"
TARGET="$DIR/settings.local.json"

case "${1:-status}" in
  on)
    cp "$DIR/settings.autopilot.json" "$TARGET"
    echo "Autopilot ON — expanded permissions active for ppvm-timeevolve."
    ;;
  off)
    cp "$DIR/settings.default.json" "$TARGET"
    echo "Autopilot OFF — permissions restored to default."
    ;;
  status)
    if diff -q "$TARGET" "$DIR/settings.autopilot.json" &>/dev/null; then
      echo "Autopilot: ON"
    elif diff -q "$TARGET" "$DIR/settings.default.json" &>/dev/null; then
      echo "Autopilot: OFF"
    else
      echo "Autopilot: UNKNOWN (settings.local.json has been modified manually)"
    fi
    ;;
  *)
    echo "Usage: $0 [on|off|status]"
    exit 1
    ;;
esac
