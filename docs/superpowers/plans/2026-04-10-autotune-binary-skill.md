# Autotune Binary Skill Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone Rust binary skill at `/Users/roger/Code/rust/autotune` that gives agents a deterministic command surface for experiment-driven performance tuning via `ion run autotune ...`.

**Architecture:** Start from Ion's binary-skill scaffold (`ion init --bin`) and replace the generic template with a focused command-oriented CLI. Keep repository-specific behavior in `.autotune.toml`, implement deterministic git/worktree and benchmark operations in the binary, and persist experiment state in structured repo-local files. Terminal-native reporting and JSON output share the same internal model.

**Tech Stack:** Rust 2024, clap, ionem, serde, toml, anyhow/thiserror, tempfile, std::process::Command, ASCII terminal rendering

---

### Task 1: Scaffold the standalone binary skill crate

**Files:**
- Create: `/Users/roger/Code/rust/autotune/Cargo.toml`
- Create: `/Users/roger/Code/rust/autotune/build.rs`
- Create: `/Users/roger/Code/rust/autotune/SKILL.md`
- Create: `/Users/roger/Code/rust/autotune/src/main.rs`

- [ ] **Step 1: Scaffold the binary skill with Ion**

Run:

```bash
ion init --bin /Users/roger/Code/rust/autotune
```

Expected: Ion prints `Created binary skill project in /Users/roger/Code/rust/autotune`.

- [ ] **Step 2: Verify the scaffolded files exist**

Run:

```bash
find /Users/roger/Code/rust/autotune -maxdepth 2 -type f | sort
```

Expected: output includes:

```text
/Users/roger/Code/rust/autotune/Cargo.toml
/Users/roger/Code/rust/autotune/SKILL.md
/Users/roger/Code/rust/autotune/build.rs
/Users/roger/Code/rust/autotune/src/main.rs
```

- [ ] **Step 3: Build the generated scaffold before modifying it**

Run:

```bash
cargo build --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml
cargo run --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml -- self skill
```

Expected: the crate compiles and `self skill` prints the generated binary SKILL.md to stdout.

- [ ] **Step 4: Commit the raw scaffold**

```bash
cd /Users/roger/Code/rust/autotune
git init
git add Cargo.toml build.rs SKILL.md src/main.rs .gitignore
git commit -m "feat: scaffold autotune binary skill"
```

---

### Task 2: Replace the template with the real CLI surface and core models

**Files:**
- Modify: `/Users/roger/Code/rust/autotune/Cargo.toml`
- Modify: `/Users/roger/Code/rust/autotune/SKILL.md`
- Modify: `/Users/roger/Code/rust/autotune/src/main.rs`
- Create: `/Users/roger/Code/rust/autotune/src/cli.rs`
- Create: `/Users/roger/Code/rust/autotune/src/error.rs`
- Create: `/Users/roger/Code/rust/autotune/src/config.rs`
- Create: `/Users/roger/Code/rust/autotune/src/state.rs`
- Create: `/Users/roger/Code/rust/autotune/src/model.rs`

- [ ] **Step 1: Add the crate dependencies**

Update `/Users/roger/Code/rust/autotune/Cargo.toml` so the dependency section contains:

```toml
[dependencies]
anyhow = "1"
clap = { version = "4", features = ["derive"] }
ionem = "0.2"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
thiserror = "2"
toml = "0.8"

[dev-dependencies]
tempfile = "3"
```

- [ ] **Step 2: Replace the generated SKILL with the binary-specific trigger text**

Replace `/Users/roger/Code/rust/autotune/SKILL.md` with:

```md
---
name: autotune
description: Deterministic CLI for experiment-driven performance tuning. Use when an agent needs validated commands for tuning workflows such as experiment initialization, worktree preparation, benchmark execution, metric recording, keep/discard application, and terminal-native reporting. Invoke with `ion run autotune ...`.
metadata:
  binary: autotune
---

# autotune

Use `ion run autotune <command>` for deterministic tuning operations.

## Commands

- `init` — validate `.autotune.toml` and initialize experiment state
- `prepare-approach` — create approach metadata, branch, and worktree
- `benchmark` — run configured benchmark commands and parse metrics
- `record` — validate metrics, score them, and append them to state
- `apply` — integrate winning code or preserve losing attempts without integrating code
- `report` — show terminal-native summaries and ASCII charts

Use `ion run autotune self info`, `self check`, and `self update` for binary self-management.
```

- [ ] **Step 3: Create the core CLI model files**

Create `/Users/roger/Code/rust/autotune/src/model.rs` with the core serializable types:

```rust
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Direction {
    Maximize,
    Minimize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimaryMetric {
    pub name: String,
    pub direction: Direction,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuardrailMetric {
    pub name: String,
    pub direction: Direction,
    pub max_regression: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetricValue {
    pub name: String,
    pub value: f64,
}
```

- [ ] **Step 4: Create the error type**

Create `/Users/roger/Code/rust/autotune/src/error.rs`:

```rust
use thiserror::Error;

#[derive(Debug, Error)]
pub enum AutotuneError {
    #[error("config error: {0}")]
    Config(String),
    #[error("state error: {0}")]
    State(String),
    #[error("git error: {0}")]
    Git(String),
    #[error("benchmark error: {0}")]
    Benchmark(String),
}
```

- [ ] **Step 5: Define the config and state file formats**

Create `/Users/roger/Code/rust/autotune/src/config.rs`:

```rust
use serde::{Deserialize, Serialize};

use crate::model::{GuardrailMetric, PrimaryMetric};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkCommand {
    pub name: String,
    pub command: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetricExtractor {
    pub metric: String,
    pub pattern: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AutotuneConfig {
    pub experiment_dir: String,
    pub prompt_dir: String,
    pub canonical_branch: Option<String>,
    pub primary_metrics: Vec<PrimaryMetric>,
    pub guardrail_metrics: Vec<GuardrailMetric>,
    pub benchmark_commands: Vec<BenchmarkCommand>,
    pub metric_extractors: Vec<MetricExtractor>,
}
```

Create `/Users/roger/Code/rust/autotune/src/state.rs`:

```rust
use serde::{Deserialize, Serialize};

use crate::model::MetricValue;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IterationRecord {
    pub iteration: usize,
    pub approach: String,
    pub status: String,
    pub score: f64,
    pub metrics: Vec<MetricValue>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExperimentState {
    pub name: String,
    pub canonical_branch: String,
    pub baseline_iteration: Option<usize>,
    pub iterations: Vec<IterationRecord>,
}
```

- [ ] **Step 6: Replace `src/main.rs` with a modular entry point**

Replace `/Users/roger/Code/rust/autotune/src/main.rs` with:

```rust
mod cli;
mod config;
mod error;
mod model;
mod state;

use clap::Parser;
use ionem::self_update::SelfManager;

use crate::cli::{Cli, Commands, SelfCommands};

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::SelfCmd { command } => {
            let manager = SelfManager::new(env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
            match command {
                SelfCommands::Skill => print!("{}", include_str!("../SKILL.md")),
                SelfCommands::Info => manager.info()?,
                SelfCommands::Check => manager.check()?,
                SelfCommands::Update => manager.update()?,
            }
        }
        _ => anyhow::bail!("subcommands not implemented yet"),
    }
    Ok(())
}
```

- [ ] **Step 7: Add the clap command definitions**

Create `/Users/roger/Code/rust/autotune/src/cli.rs`:

```rust
use clap::{Parser, Subcommand};

#[derive(Debug, Parser)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    Init,
    PrepareApproach,
    Benchmark,
    Record,
    Apply,
    Report,
    SelfCmd {
        #[command(subcommand)]
        command: SelfCommands,
    },
}

#[derive(Debug, Subcommand)]
pub enum SelfCommands {
    Skill,
    Info,
    Check,
    Update,
}
```

- [ ] **Step 8: Verify the real command skeleton compiles**

Run:

```bash
cargo build --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml
cargo run --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml -- self skill
cargo run --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml -- init
```

Expected:
- first two commands succeed
- `init` exits with the temporary `subcommands not implemented yet` error

- [ ] **Step 9: Commit the command skeleton**

```bash
cd /Users/roger/Code/rust/autotune
git add Cargo.toml SKILL.md src/main.rs src/cli.rs src/error.rs src/config.rs src/state.rs src/model.rs
git commit -m "feat: define autotune binary command surface"
```

---

### Task 3: Implement config parsing, state persistence, and scoring

**Files:**
- Modify: `/Users/roger/Code/rust/autotune/src/config.rs`
- Modify: `/Users/roger/Code/rust/autotune/src/state.rs`
- Create: `/Users/roger/Code/rust/autotune/src/score.rs`
- Create: `/Users/roger/Code/rust/autotune/tests/config_state.rs`

- [ ] **Step 1: Write the failing tests for config and score behavior**

Create `/Users/roger/Code/rust/autotune/tests/config_state.rs`:

```rust
use autotune::config::AutotuneConfig;
use autotune::model::{Direction, MetricValue, PrimaryMetric};

#[test]
fn parses_primary_and_guardrail_metrics() {
    let raw = r#"
experiment_dir = "docs/autotune"
prompt_dir = "docs/autotune/prompts"

[[primary_metrics]]
name = "throughput"
direction = "Maximize"

[[guardrail_metrics]]
name = "latency"
direction = "Minimize"
max_regression = 0.05

[[benchmark_commands]]
name = "bench"
command = ["cargo", "bench", "--bench", "micro"]

[[metric_extractors]]
metric = "throughput"
pattern = "throughput=(?P<value>[0-9.]+)"
"#;

    let parsed: AutotuneConfig = toml::from_str(raw).unwrap();
    assert_eq!(parsed.primary_metrics.len(), 1);
    assert_eq!(parsed.guardrail_metrics.len(), 1);
}
```

- [ ] **Step 2: Run the tests to verify the crate layout fails**

Run:

```bash
cargo test --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml --test config_state
```

Expected: FAIL because the crate does not yet export the tested modules as a library.

- [ ] **Step 3: Add a library entry point and scoring logic**

Create `/Users/roger/Code/rust/autotune/src/lib.rs`:

```rust
pub mod cli;
pub mod config;
pub mod error;
pub mod model;
pub mod score;
pub mod state;
```

Create `/Users/roger/Code/rust/autotune/src/score.rs`:

```rust
use crate::model::{Direction, MetricValue, PrimaryMetric};

pub fn rank_score(primary: &[PrimaryMetric], baseline: &[MetricValue], candidate: &[MetricValue]) -> f64 {
    primary.iter().map(|metric| {
        let base = baseline.iter().find(|m| m.name == metric.name).unwrap().value;
        let next = candidate.iter().find(|m| m.name == metric.name).unwrap().value;
        match metric.direction {
            Direction::Maximize => (next - base) / base,
            Direction::Minimize => (base - next) / base,
        }
    }).sum()
}
```

- [ ] **Step 4: Add state load/save helpers**

Append to `/Users/roger/Code/rust/autotune/src/state.rs`:

```rust
impl ExperimentState {
    pub fn load(path: &std::path::Path) -> anyhow::Result<Self> {
        Ok(serde_json::from_str(&std::fs::read_to_string(path)?)?)
    }

    pub fn save(&self, path: &std::path::Path) -> anyhow::Result<()> {
        std::fs::write(path, serde_json::to_string_pretty(self)?)?;
        Ok(())
    }
}
```

- [ ] **Step 5: Expand the tests to cover ranking and state round-trips**

Append to `/Users/roger/Code/rust/autotune/tests/config_state.rs`:

```rust
use autotune::score::rank_score;

#[test]
fn sums_multi_primary_metric_improvements() {
    let primary = vec![
        PrimaryMetric { name: "throughput".into(), direction: Direction::Maximize },
        PrimaryMetric { name: "latency".into(), direction: Direction::Minimize },
    ];
    let baseline = vec![
        MetricValue { name: "throughput".into(), value: 100.0 },
        MetricValue { name: "latency".into(), value: 10.0 },
    ];
    let candidate = vec![
        MetricValue { name: "throughput".into(), value: 110.0 },
        MetricValue { name: "latency".into(), value: 9.0 },
    ];
    assert!(rank_score(&primary, &baseline, &candidate) > 0.0);
}
```

- [ ] **Step 6: Verify config and state tests pass**

Run:

```bash
cargo test --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml --test config_state
```

Expected: PASS.

- [ ] **Step 7: Commit the config/state layer**

```bash
cd /Users/roger/Code/rust/autotune
git add src/lib.rs src/config.rs src/state.rs src/score.rs tests/config_state.rs
git commit -m "feat: add autotune config parsing and score model"
```

---

### Task 4: Implement deterministic git, benchmark, and record/apply commands

**Files:**
- Create: `/Users/roger/Code/rust/autotune/src/git_ops.rs`
- Create: `/Users/roger/Code/rust/autotune/src/benchmark.rs`
- Create: `/Users/roger/Code/rust/autotune/src/commands/init.rs`
- Create: `/Users/roger/Code/rust/autotune/src/commands/prepare_approach.rs`
- Create: `/Users/roger/Code/rust/autotune/src/commands/benchmark.rs`
- Create: `/Users/roger/Code/rust/autotune/src/commands/record.rs`
- Create: `/Users/roger/Code/rust/autotune/src/commands/apply.rs`
- Create: `/Users/roger/Code/rust/autotune/tests/workflow.rs`
- Modify: `/Users/roger/Code/rust/autotune/src/main.rs`

- [ ] **Step 1: Write one end-to-end failing workflow test**

Create `/Users/roger/Code/rust/autotune/tests/workflow.rs`:

```rust
#[test]
fn prepare_benchmark_record_apply_round_trip() {
    let repo = tempfile::tempdir().unwrap();
    let config = repo.path().join(".autotune.toml");
    std::fs::write(&config, r#"
experiment_dir = "docs/autotune"
prompt_dir = "docs/autotune/prompts"
canonical_branch = "main"

[[primary_metrics]]
name = "score"
direction = "Maximize"

[[benchmark_commands]]
name = "bench"
command = ["sh", "-c", "printf 'score=11\\n'"]

[[metric_extractors]]
metric = "score"
pattern = "score=(?P<value>[0-9.]+)"
"#).unwrap();

    assert!(repo.path().join(".autotune.toml").exists());
}
```

- [ ] **Step 2: Add the git and benchmark helpers**

Create `/Users/roger/Code/rust/autotune/src/git_ops.rs`:

```rust
use std::path::Path;
use std::process::Command;

pub fn git(repo: &Path, args: &[&str]) -> anyhow::Result<String> {
    let output = Command::new("git").args(args).current_dir(repo).output()?;
    if !output.status.success() {
        anyhow::bail!("{}", String::from_utf8_lossy(&output.stderr));
    }
    Ok(String::from_utf8_lossy(&output.stdout).trim().to_string())
}
```

Create `/Users/roger/Code/rust/autotune/src/benchmark.rs`:

```rust
use std::process::Command;

use crate::config::{BenchmarkCommand, MetricExtractor};
use crate::model::MetricValue;

pub fn run_and_parse(cmd: &BenchmarkCommand, extractors: &[MetricExtractor]) -> anyhow::Result<Vec<MetricValue>> {
    let (program, args) = cmd.command.split_first().ok_or_else(|| anyhow::anyhow!("empty command"))?;
    let output = Command::new(program).args(args).output()?;
    if !output.status.success() {
        anyhow::bail!("{}", String::from_utf8_lossy(&output.stderr));
    }
    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut values = Vec::new();
    for extractor in extractors {
        let (_, rhs) = stdout
            .lines()
            .find_map(|line| line.split_once('='))
            .ok_or_else(|| anyhow::anyhow!("missing metric output"))?;
        values.push(MetricValue { name: extractor.metric.clone(), value: rhs.parse()? });
    }
    Ok(values)
}
```

- [ ] **Step 3: Implement concrete subcommand modules**

Create one file per command in `/Users/roger/Code/rust/autotune/src/commands/` with `pub fn run(...) -> anyhow::Result<()>` entry points. The command modules should:

- `init.rs`: read `.autotune.toml`, infer or validate canonical branch, create state dir, write initial state file
- `prepare_approach.rs`: create approach branch and worktree using `git worktree add -b`
- `benchmark.rs`: run one named benchmark command and print metrics as JSON
- `record.rs`: load baseline, score primary metrics, validate guardrails, append iteration to state
- `apply.rs`: integrate winning commit onto canonical branch or preserve a losing attempt without integration

- [ ] **Step 4: Wire the command modules into `src/main.rs`**

Replace the non-self branch in `/Users/roger/Code/rust/autotune/src/main.rs` with dispatch like:

```rust
match cli.command {
    Commands::Init => commands::init::run()?,
    Commands::PrepareApproach => commands::prepare_approach::run()?,
    Commands::Benchmark => commands::benchmark::run()?,
    Commands::Record => commands::record::run()?,
    Commands::Apply => commands::apply::run()?,
    Commands::Report => commands::report::run()?,
    Commands::SelfCmd { command } => { /* existing self manager */ }
}
```

- [ ] **Step 5: Verify the deterministic workflow compiles and the tests run**

Run:

```bash
cargo test --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml --test workflow
cargo test --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml
```

Expected: PASS.

- [ ] **Step 6: Commit the deterministic command layer**

```bash
cd /Users/roger/Code/rust/autotune
git add src/main.rs src/git_ops.rs src/benchmark.rs src/commands tests/workflow.rs
git commit -m "feat: implement deterministic autotune command flow"
```

---

### Task 5: Implement terminal-native reporting and ppvm development fixture

**Files:**
- Create: `/Users/roger/Code/rust/autotune/src/report.rs`
- Create: `/Users/roger/Code/rust/autotune/src/commands/report.rs`
- Create: `/Users/roger/Code/rust/autotune/tests/report.rs`
- Create: `/Users/roger/Code/rust/ppvm/.autotune.toml`

- [ ] **Step 1: Write the failing report test**

Create `/Users/roger/Code/rust/autotune/tests/report.rs`:

```rust
#[test]
fn report_renders_ascii_summary() {
    let rendered = autotune::report::render_chart(
        &[("try-a".to_string(), 0.1), ("try-b".to_string(), 0.3)],
        "score",
    );
    assert!(rendered.contains("score"));
    assert!(rendered.contains("try-a"));
    assert!(rendered.contains("try-b"));
}
```

- [ ] **Step 2: Implement the report renderer**

Create `/Users/roger/Code/rust/autotune/src/report.rs`:

```rust
pub fn render_chart(points: &[(String, f64)], metric: &str) -> String {
    let mut out = format!("metric: {metric}\n");
    for (name, value) in points {
        let bars = (value.max(0.0) * 20.0).round() as usize;
        out.push_str(&format!("{name:16} | {}\n", "#".repeat(bars.max(1))));
    }
    out
}
```

- [ ] **Step 3: Implement the `report` command**

Create `/Users/roger/Code/rust/autotune/src/commands/report.rs`:

```rust
use crate::report::render_chart;

pub fn run() -> anyhow::Result<()> {
    let output = render_chart(&[("baseline".into(), 1.0)], "score");
    print!("{output}");
    Ok(())
}
```

- [ ] **Step 4: Add a real ppvm config fixture for development**

Create `/Users/roger/Code/rust/ppvm/.autotune.toml` with an initial config like:

```toml
experiment_dir = "docs/autotune"
prompt_dir = "docs/autotune/prompts"
canonical_branch = "main"

[[primary_metrics]]
name = "xgate"
direction = "Minimize"

[[guardrail_metrics]]
name = "hgate"
direction = "Minimize"
max_regression = 0.05

[[benchmark_commands]]
name = "micro-xgate"
command = ["cargo", "bench", "--bench", "micro", "--", "gates/single-qubit/x"]

[[metric_extractors]]
metric = "xgate"
pattern = "xgate=(?P<value>[0-9.]+)"
```

- [ ] **Step 5: Verify the report and full test suite**

Run:

```bash
cargo test --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml --test report
cargo test --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml
cargo run --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml -- report
```

Expected:
- tests pass
- `report` prints terminal-native text with an ASCII bar chart

- [ ] **Step 6: Commit the reporting layer and ppvm fixture**

```bash
cd /Users/roger/Code/rust/autotune
git add src/report.rs src/commands/report.rs tests/report.rs
git commit -m "feat: add terminal-native autotune reporting"

cd /Users/roger/Code/rust/ppvm
git add .autotune.toml
git commit -m "chore(autotune): add ppvm autotune config fixture"
```

---

### Task 6: Validate binary-skill behavior and repo integration

**Files:**
- Modify: `/Users/roger/Code/rust/autotune/SKILL.md` (only if validation requires fixes)
- Modify: `/Users/roger/Code/rust/autotune/src/main.rs` (only if validation requires fixes)

- [ ] **Step 1: Verify the binary-skill self interface**

Run:

```bash
cargo run --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml -- self skill
cargo run --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml -- self info
```

Expected:
- `self skill` prints valid SKILL.md with `metadata.binary: autotune`
- `self info` prints version/build information without errors

- [ ] **Step 2: Validate the generated SKILL**

Run:

```bash
tmpdir="$(mktemp -d)"
cargo run --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml -- self skill > "$tmpdir/SKILL.md"
ion --json skill validate "$tmpdir/SKILL.md"
```

Expected: validation succeeds, or if Ion only validates directories in this version, create a temp skill directory containing that SKILL and validate the directory.

- [ ] **Step 3: Run formatting, lint, and full tests**

Run:

```bash
cargo fmt --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml --all
cargo clippy --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml --all-targets --all-features -- -D warnings
cargo test --manifest-path /Users/roger/Code/rust/autotune/Cargo.toml
```

Expected: all commands succeed without warnings or failures.

- [ ] **Step 4: Commit validation-driven fixes if any**

```bash
cd /Users/roger/Code/rust/autotune
git add .
git commit -m "fix: resolve autotune validation issues"
```
