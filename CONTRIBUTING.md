# Contributing to PPVM

Thanks for your interest in contributing. This document explains how
contributions are licensed and how to get a change merged.

## Licensing of contributions

By opening a pull request against this repository, **you agree that your
contribution is licensed under the [Apache License, Version 2.0](LICENSE)
and that you accept the terms of the [PPVM Contributor License
Agreement](CLA.md)**, which grants QuEra Computing Inc. and downstream
recipients the rights described in that document.

You retain copyright on your contribution; the CLA grants licenses, it
does not assign ownership. See [`CLA.md`](CLA.md) for the full text. If
you are contributing on behalf of an employer, please make sure your
employer is aware of and authorizes the contribution — the CLA's
representations cover this in §4.

If you do not agree to those terms, please do not open a pull request.

## Where to start

- **Code overview**: read [`docs/src/pages/develop.astro`](docs/src/pages/develop.astro)
  first (rendered at `/develop/`). It is the canonical developer guide.
- **Build & test**: `cargo test --workspace` for Rust;
  `uv run --project ppvm-python --group dev pytest …` for Python;
  `cd docs && npm run build` for the docs site.
- **Style**: respect existing patterns in each crate. `prek run --all-files`
  must pass before pushing (formatting, lints, license headers).

## Workflow

1. Fork the repository and create a topic branch off `main`.
2. Make your change. Keep commits focused and use
   [Conventional Commits](https://www.conventionalcommits.org/)
   (e.g. `feat(runtime): add new gate`).
3. Add or update tests for the behavior you changed.
4. Run `prek run --all-files` and the relevant test command(s) locally.
5. Open a pull request. CI runs Rust + Python + license-header checks.
   First-time contributors will see an automated welcome comment with
   a pointer to this document.
6. A maintainer will review. Be prepared to iterate.

## Reporting issues

Use GitHub Issues for bugs and feature requests. For security-sensitive
reports, please email the maintainers privately rather than filing a
public issue.
