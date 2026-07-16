// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashMap;
use std::hash::Hash;
use std::time::{SystemTime, UNIX_EPOCH};

/// Runtime configuration for trajectory continuation caching.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CacheConfig {
    pub enabled: bool,
    pub max_nodes: Option<usize>,
}

impl CacheConfig {
    pub const fn disabled() -> Self {
        Self {
            enabled: false,
            max_nodes: None,
        }
    }

    pub const fn unbounded() -> Self {
        Self {
            enabled: true,
            max_nodes: None,
        }
    }

    pub const fn bounded(max_nodes: usize) -> Self {
        Self {
            enabled: true,
            max_nodes: Some(max_nodes),
        }
    }
}

impl Default for CacheConfig {
    fn default() -> Self {
        Self::disabled()
    }
}

/// Counters collected by the shared cache runner.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CacheStats {
    pub nodes: usize,
    pub hits: u64,
    pub misses: u64,
    pub evictions: u64,
    pub terminal_hits: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CachedRun<O> {
    pub output: O,
    pub cache_stats: CacheStats,
}

/// Result of running deterministic work up to the next stochastic boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TrajectoryEvent<O> {
    Boundary,
    Terminal(O),
}

/// Adapter implemented by concrete program runners.
///
/// The contract is deliberately narrow: `run_until_boundary` must execute only
/// deterministic work and stop immediately before the next stochastic operation.
/// `execute_boundary` then executes exactly that operation and returns a choice
/// that fully determines the continuation state for cache-key purposes.
pub trait TrajectoryProgram {
    type Snapshot: Clone;
    type Choice: Eq + Hash + Clone;
    type Output: Clone;
    type Error;

    fn reset_for_shot(&mut self, shot: usize) -> Result<(), Self::Error>;
    fn snapshot(&self) -> Self::Snapshot;
    fn restore(&mut self, snapshot: &Self::Snapshot) -> Result<(), Self::Error>;
    fn reseed(&mut self, seed: u64) -> Result<(), Self::Error>;
    fn run_until_boundary(&mut self) -> Result<TrajectoryEvent<Self::Output>, Self::Error>;
    fn execute_boundary(&mut self) -> Result<Self::Choice, Self::Error>;
}

#[derive(Clone)]
enum CacheEntry<S, O> {
    Boundary(S),
    Terminal(O),
}

struct Node<S, C, O> {
    entry: Option<CacheEntry<S, O>>,
    children: HashMap<C, usize>,
    parent: Option<usize>,
    parent_choice: Option<C>,
    hits: u64,
    last_used: u64,
}

impl<S, C, O> Node<S, C, O> {
    fn root() -> Self {
        Self {
            entry: None,
            children: HashMap::new(),
            parent: None,
            parent_choice: None,
            hits: 0,
            last_used: 0,
        }
    }

    fn child(entry: CacheEntry<S, O>, parent: usize, parent_choice: C, tick: u64) -> Self {
        Self {
            entry: Some(entry),
            children: HashMap::new(),
            parent: Some(parent),
            parent_choice: Some(parent_choice),
            hits: 0,
            last_used: tick,
        }
    }
}

pub fn random_base_seed() -> u64 {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos() as u64;
    splitmix64(now ^ (&now as *const u64 as usize as u64))
}

pub fn continuation_seed(base_seed: u64, shot: usize, depth: usize) -> u64 {
    let x = base_seed
        ^ (shot as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ (depth as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    splitmix64(x)
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

pub fn run_cached_shots<P>(
    program: &mut P,
    shots: usize,
    config: CacheConfig,
    base_seed: u64,
) -> Result<CachedRun<Vec<P::Output>>, P::Error>
where
    P: TrajectoryProgram,
{
    if !config.enabled {
        return run_uncached_shots(program, shots);
    }

    let mut cache = Cache::<P::Snapshot, P::Choice, P::Output>::new(config);
    let mut output = Vec::with_capacity(shots);
    for shot in 0..shots {
        output.push(cache.run_one(program, shot, base_seed)?);
    }

    Ok(CachedRun {
        output,
        cache_stats: cache.stats,
    })
}

pub fn run_uncached_shots<P>(
    program: &mut P,
    shots: usize,
) -> Result<CachedRun<Vec<P::Output>>, P::Error>
where
    P: TrajectoryProgram,
{
    let mut output = Vec::with_capacity(shots);
    for shot in 0..shots {
        program.reset_for_shot(shot)?;
        loop {
            match program.run_until_boundary()? {
                TrajectoryEvent::Boundary => {
                    let _ = program.execute_boundary()?;
                }
                TrajectoryEvent::Terminal(result) => {
                    output.push(result);
                    break;
                }
            }
        }
    }
    Ok(CachedRun {
        output,
        cache_stats: CacheStats::default(),
    })
}

struct Cache<S, C, O> {
    nodes: Vec<Node<S, C, O>>,
    config: CacheConfig,
    stats: CacheStats,
    tick: u64,
}

impl<S, C, O> Cache<S, C, O>
where
    S: Clone,
    C: Eq + Hash + Clone,
    O: Clone,
{
    fn new(config: CacheConfig) -> Self {
        Self {
            nodes: vec![Node::root()],
            config,
            stats: CacheStats::default(),
            tick: 0,
        }
    }

    fn run_one<P>(&mut self, program: &mut P, shot: usize, base_seed: u64) -> Result<O, P::Error>
    where
        P: TrajectoryProgram<Snapshot = S, Choice = C, Output = O>,
    {
        program.reset_for_shot(shot)?;
        let mut node = self.ensure_root(program)?;
        let mut depth = 0usize;

        loop {
            self.touch(node);
            match self.nodes[node].entry.as_ref().expect("initialized node") {
                CacheEntry::Terminal(output) => {
                    self.stats.terminal_hits += 1;
                    return Ok(output.clone());
                }
                CacheEntry::Boundary(snapshot) => {
                    program.restore(snapshot)?;
                    program.reseed(continuation_seed(base_seed, shot, depth))?;
                }
            }

            let choice = program.execute_boundary()?;
            depth += 1;

            if let Some(&child) = self.nodes[node].children.get(&choice)
                && self.nodes[child].entry.is_some()
            {
                self.stats.hits += 1;
                self.nodes[child].hits += 1;
                node = child;
                continue;
            }

            self.stats.misses += 1;
            program.reseed(continuation_seed(base_seed, shot, depth))?;
            let (entry, terminal_output) = match program.run_until_boundary()? {
                TrajectoryEvent::Boundary => (CacheEntry::Boundary(program.snapshot()), None),
                TrajectoryEvent::Terminal(output) => {
                    (CacheEntry::Terminal(output.clone()), Some(output))
                }
            };
            node = self.insert_child(node, choice, entry);

            if let Some(output) = terminal_output {
                return Ok(output);
            }
        }
    }

    fn ensure_root<P>(&mut self, program: &mut P) -> Result<usize, P::Error>
    where
        P: TrajectoryProgram<Snapshot = S, Choice = C, Output = O>,
    {
        if self.nodes[0].entry.is_none() {
            self.stats.misses += 1;
            let entry = match program.run_until_boundary()? {
                TrajectoryEvent::Boundary => CacheEntry::Boundary(program.snapshot()),
                TrajectoryEvent::Terminal(output) => CacheEntry::Terminal(output),
            };
            self.nodes[0].entry = Some(entry);
            self.stats.nodes += 1;
            self.enforce_limit();
        } else {
            self.stats.hits += 1;
        }
        Ok(0)
    }

    fn insert_child(&mut self, parent: usize, choice: C, entry: CacheEntry<S, O>) -> usize {
        if let Some(&existing) = self.nodes[parent].children.get(&choice) {
            self.nodes[existing].entry = Some(entry);
            self.touch(existing);
            return existing;
        }

        self.tick += 1;
        let idx = self.nodes.len();
        self.nodes
            .push(Node::child(entry, parent, choice.clone(), self.tick));
        self.nodes[parent].children.insert(choice, idx);
        self.stats.nodes += 1;
        self.enforce_limit();
        idx
    }

    fn touch(&mut self, idx: usize) {
        self.tick += 1;
        self.nodes[idx].last_used = self.tick;
    }

    fn enforce_limit(&mut self) {
        let Some(limit) = self.config.max_nodes else {
            return;
        };

        while self.stats.nodes > limit {
            let Some(victim) = self
                .nodes
                .iter()
                .enumerate()
                .filter(|(idx, node)| *idx != 0 && node.entry.is_some() && node.children.is_empty())
                .min_by_key(|(_, node)| (node.hits, node.last_used))
                .map(|(idx, _)| idx)
            else {
                break;
            };

            if let (Some(parent), Some(choice)) = (
                self.nodes[victim].parent,
                self.nodes[victim].parent_choice.clone(),
            ) {
                self.nodes[parent].children.remove(&choice);
            }
            self.nodes[victim].entry = None;
            self.stats.nodes -= 1;
            self.stats.evictions += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::VecDeque;

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct Snapshot {
        pc: usize,
        output: Vec<u8>,
    }

    struct FakeProgram {
        pc: usize,
        output: Vec<u8>,
        choices: VecDeque<u8>,
        runs_to_boundary: usize,
    }

    impl FakeProgram {
        fn new(choices: impl IntoIterator<Item = u8>) -> Self {
            Self {
                pc: 0,
                output: Vec::new(),
                choices: choices.into_iter().collect(),
                runs_to_boundary: 0,
            }
        }
    }

    impl TrajectoryProgram for FakeProgram {
        type Snapshot = Snapshot;
        type Choice = u8;
        type Output = Vec<u8>;
        type Error = &'static str;

        fn reset_for_shot(&mut self, _shot: usize) -> Result<(), Self::Error> {
            self.pc = 0;
            self.output.clear();
            Ok(())
        }

        fn snapshot(&self) -> Self::Snapshot {
            Snapshot {
                pc: self.pc,
                output: self.output.clone(),
            }
        }

        fn restore(&mut self, snapshot: &Self::Snapshot) -> Result<(), Self::Error> {
            self.pc = snapshot.pc;
            self.output = snapshot.output.clone();
            Ok(())
        }

        fn reseed(&mut self, _seed: u64) -> Result<(), Self::Error> {
            Ok(())
        }

        fn run_until_boundary(&mut self) -> Result<TrajectoryEvent<Self::Output>, Self::Error> {
            self.runs_to_boundary += 1;
            match self.pc {
                0 => {
                    self.output.push(9);
                    self.pc = 1;
                    Ok(TrajectoryEvent::Boundary)
                }
                2 => {
                    self.output.push(7);
                    self.pc = 3;
                    Ok(TrajectoryEvent::Terminal(self.output.clone()))
                }
                _ => Err("unexpected pc"),
            }
        }

        fn execute_boundary(&mut self) -> Result<Self::Choice, Self::Error> {
            let choice = self.choices.pop_front().ok_or("missing choice")?;
            self.output.push(choice);
            self.pc = 2;
            Ok(choice)
        }
    }

    #[test]
    fn cached_runner_reuses_deterministic_segments_for_repeated_choices() {
        let mut program = FakeProgram::new([1, 1, 1]);
        let report = run_cached_shots(&mut program, 3, CacheConfig::bounded(16), 123).unwrap();

        assert_eq!(report.output, vec![vec![9, 1, 7]; 3]);
        assert!(report.cache_stats.hits >= 2);
        assert!(
            program.runs_to_boundary < 6,
            "cache should skip at least one deterministic segment"
        );
    }

    #[test]
    fn bounded_cache_evicts_cold_leaf_states() {
        let mut program = FakeProgram::new([0, 1, 0, 1]);
        let report = run_cached_shots(&mut program, 4, CacheConfig::bounded(1), 456).unwrap();

        assert_eq!(report.output.len(), 4);
        assert!(report.cache_stats.evictions > 0);
        assert!(report.cache_stats.nodes <= 1);
    }
}
