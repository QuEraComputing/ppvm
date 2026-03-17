# ppvm-timeevolve: Agent Workflow

At the start of every response, read `agents/state.md` to determine the active role and
current task. Adopt that role fully for the duration of the response by reading the
corresponding persona file.

## Roles

- **developer** → follow `agents/developer.md`
- **perf-developer** → follow `agents/perf-developer.md`
- **reviewer** → follow `agents/reviewer.md`

## Role-switching rules

Update `agents/state.md` and switch roles immediately when one of these triggers fires:

| Trigger | New role |
|---------|----------|
| Developer or perf-developer explicitly hands off for review (summary posted, review requested) | `reviewer` |
| Reviewer approves ("Task N approved.") | `perf-developer` — advance `current_task` by 1 |
| Reviewer requests changes (returns feedback without approval) | `perf-developer` — keep `current_task` unchanged |

Switches happen at the **end** of the response that contains the trigger. The very next
response must open by reading `agents/state.md` and acting as the new role.

## State updates

When switching roles, rewrite the frontmatter in `agents/state.md`:

```
---
active_role: <developer|reviewer>
current_task: <N>
---
```

No other content in that file should change.

## Invariants

- The reviewer never edits code.
- The developer never commits before the reviewer approves.
- No task is skipped.
- `current_task` only increments on reviewer approval.
