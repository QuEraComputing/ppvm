# Design document of Auto-tuning skill

We want to create a skill to improve performance by tweaking our implementation based on
karpathy's autoresearch: https://github.com/karpathy/autoresearch/tree/master, the core idea is a loop of short
time bounded exploration tasks, each task record a set of metrics, and each iteration we try a few different
approaches to improve performance, and combine them, then run the next iteration, and combine them.

## Workflow

### Read user input

The user will specify which metric to improve. The metric must be a number that we can obtain from the benchmark. Validate this property first. If the metric cannot be obtained from the benchmark, we should reject the user's request, ask for clarification, and have them specify which metrics they want to improve.

The user should not request to improve a metric that is not in the benchmark, or is not a number. We want to make sure that the metric is something we can measure and track, so that we can have a clear goal for our experiment.

The user should not request to improve a large set of metrics, the metric to improve should be focused and specific, and do not require running the entire benchmark to get the metric. We want to make sure that the metric is something we can measure quickly, so that we can iterate faster and have more experiments in a shorter time. If the user wants to improve a large set of metrics, we can ask them to prioritize the metrics, and focus on the top 1 or 2 metrics that are most important to them, and that we can measure quickly.

If the user is not sure what small metric to measure, we should run profiling for them, and identify the bottleneck in the codebase, and suggest some metrics that they can track to improve the performance of the bottleneck. We want to help the user to find a good metric to track, so that they can have a clear goal for their experiment, and we can have a more focused and effective experiment.

### Prepare the experiment environment

1. think about a nice human rememberable name like `sunny-friday` as the overral experiment task name combined with the actual date `<date>`, will be referred as `<task>`
2. create a folder to store the project files in `docs/autotune/<project>/`
3. initialize a TOML file `docs/autotune/<project>/metric.toml` to store the report metrics as an array of some metric numbers, for example

```toml
[[metric]]
commit = "abc123" # the commit hash for this metric
status = "keep" # or "discard", "crash"
description = "baseline performance"
xgate = 1000
```

The key should match the corresponding benchmark name. There can be multiple metrics we keep track of, using different keys, e.g

```toml
[[metric]]
commit = "abc123" # the commit hash for this metric
status = "keep"
description = "added multi-threading to x gate"
xgate = 1000 # performance of x gate in us
hgate = 2000 # performance of h gate in us

[[metric]]
commit = "abc123" # the commit hash for this metric
status = "discard"
description = "used ghash"
xgate = 1000 # performance of x gate in us
hgate = 2000 # performance of h gate in us
```

4. initialize a markdown file `docs/autotune/<project>/log.md` to keep a log of interesting findings for future references, the format can be something like

```markdown
# Log for <task>
## <date>
- <log entry 1>
- <log entry 2>
```

each log entry should be a short description of the finding, and can be linked to the corresponding commit for more details.

The log file is not for recording every single change, but only the interesting ones that may be useful for future reference, or may be useful for other people who want to learn from our experiment. We want to keep the log file concise and informative, and avoid writing too much unnecessary details. If you want to record the description of every single change, you can use the commit message for that, and link the commit in the log file if needed.

### Task loop

LOOP FOREVER until we are satisfied with the performance improvement, or we have tried enough approaches and want to stop:

Host agent loop:
1. read the `docs/autotune/<project>/log.md` we have so far, and come up with 1 approach to try in this iteration
2. create a worktree for this approach with branch name in form `auto-tune-<project>-<approach name>`, e.g `auto-tune-sunny-friday-2026-03-12-try-multi-thread`
3. for this worktree, create a folder to store the approach files in `docs/autotune/<project>/<approach name>/`, e.g `docs/autotune/sunny-friday-2026-03-12/try-multi-thread/`
4. for this worktree, create a markdown file to describe the approach, and target in `docs/autotune/<project>/<approach name>/prompts.md`, e.g `docs/autotune/sunny-friday-2026-03-12/try-multi-thread/prompts.md`, the prompt can be something like

```markdown
# Approach: try multi-threading for x gate

## Target metric
- xgate performance in us
- hgate performance in us
- other metrics we want to keep track of

## Description
- we want to try multi-threading for x gate, and see if it can improve the performance. We will use rayon for multi-threading, and see if it can improve the performance of
- the x gate, and also see if it has any impact on the performance of h gate, and other metrics we want to keep track of. We will also compare the performance with the baseline performance we have in the metric.toml file, and see if it is an improvement or not.

## Previous findings

We have tried some approaches in the previous iterations, you can find them in the `docs/autotune/<project>/log.md` file for interesting findings,
and there is also the `docs/autotune/<project>/metric.toml` for the previous metrics we have recorded.
```

5. spawn a subagent for the implementation of this approach, the subagent will read the prompt we just created, and try to implement the approach in the codebase, and improve the performance of the target metric. The subagent should only edit files within `src/` dir of each crate, and should not touch any test, benchmark, other documentation files. The subagent should also not write massive documentation beyond `log.md`, and should focus on the implementation and performance improvement.
6. after the subagent has implemented the approach, we will squash merge the worktree branch to the `<user branch>` (the branch where user started the process, can be main or other branches), and then delete the worktree for this approach. We want to keep the commit history clean, and the squashed commit message title should be in form of `fix(autotune): <approach name> <description>`, e.g `fix(autotune): try multi-threading for x gate`, and the commit message body can be the description based on what we have in the prompt file for this approach, and also link to the prompt file for more details. For example

```
fix(autotune): try multi-threading for x gate

We want to try multi-threading for x gate, and see if it can improve the performance. We will use rayon for multi-threading, and see if it can improve the performance of
the x gate, and also see if it has any impact on the performance of h gate, and other metrics we want to keep track of. We will also compare the performance with the baseline performance
we have in the metric.toml file, and see if it is an improvement or not.
For more details, please refer to the prompt file for this approach: `docs/autotune/sunny-friday-2026-03-12/try-multi-thread/prompts.md`
```

7. We will run the benchmark to get the performance metrics, and record the metrics in the main branch's `metric.toml` file as the latest metric.
8. log any interesting finding in the main branch's `log.md` file.
9. if the traget metric improve, you can mark the status of this metric as "keep" in the `metric.toml` file, the latest metric with `keep` status is the new baseline for the next iteration.
10. if the target metric does not improve, you can mark the status of this metric as "discard" in the `metric.toml` file, revert the commit for this approach.

The idea is that you are a completely autonomous researcher trying things out. If they work, keep, if they don't, discard. And you're advancing the branch so that you can iterate. if you feel like you're getting stuck in some way, you can rewind but you should probably do this very very sparingly (if ever).

## Timeout
Each experiment should take ~10 minutes to run including the benchmarks, do not run the entire benchmark, do not ask the subagent run the entire benchmark. ONLY RUN THE MICRO BENCHMARKS that are related to the change and target metric. For example, if you are trying to improve the performance of x gate, you can run the micro benchmark for x gate, and see if it has any improvement, and also see if it has any impact on the performance of h gate, and other metrics we want to keep track of. You should not run the entire benchmark, because it will take too long, and it will also make it harder to iterate faster and have more experiments in a shorter time. You should focus on the micro benchmarks that are related to the change and target metric, so that you can get the feedback faster, and iterate faster.

## Crashes
If a run crashes (OOM, or a bug, or etc.), use your judgment: If it's something dumb and easy to fix (e.g. a typo, a missing import), fix it and re-run. If the idea itself is fundamentally broken, just skip it, log "crash" as the status in the tsv, and move on.

## NEVER STOP
Once the experiment loop has begun (after the initial setup), do NOT pause to ask the human if you should continue. Do NOT ask "should I keep going?" or "is this a good stopping point?". The human might be asleep, or gone from a computer and expects you to continue working indefinitely until you are manually stopped. You are autonomous. If you run out of ideas, think harder — read papers referenced in the code, re-read the in-scope files for new angles, try combining previous near-misses, try more radical architectural changes. The loop runs until the human interrupts you, period.

As an example use case, a user might leave you running while they sleep. If each experiment takes you ~10 minutes then you can run approx 6/hour, for a total of about 50 over the duration of the average human sleep. The user then wakes up to experimental results, all completed by you while they slept!
