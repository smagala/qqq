# Branch & Merge Discipline — cdcent/oamd-bio-qqq

## Branch Structure

| Branch     | Purpose                                      |
|------------|----------------------------------------------|
| `main`     | Stable releases only — never commit directly |
| `dev`      | Integration branch — all PRs target here     |
| `TEMPLATE` | nf-core template tracking — automated sync   |

## Creating Feature Branches

Always branch off `dev`:

```bash
git checkout dev
git pull origin dev
git checkout -b feature/<ticket-id>-short-description
```

Use conventional prefixes:
- `feature/` — new functionality
- `fix/` — bug fixes
- `task/` — maintenance, docs, CI changes

## Submitting Pull Requests

- PRs **must** target `dev`, never `main`
- Title should follow conventional commits: `feat:`, `fix:`, `chore:`, etc.
- Link the relevant Multica issue in the PR description
- Ensure `nf-test` passes locally before opening PR

## Merging to main

Only `dev` -> `main` merges are allowed, triggered by a release:

```bash
git checkout main
git merge --no-ff dev
git tag v<version>
git push origin main --tags
```

## nf-core TEMPLATE Sync

When nf-core/tools releases a new version, the `TEMPLATE` branch is updated via `nf-core sync`. Conflicts are resolved in a dedicated PR from `TEMPLATE` -> `dev`. Do not merge TEMPLATE into main directly.

## Agent Rules (AgentAMD)

- Pick up **one ticket at a time**, mark `in_progress` before starting
- Create branch off `dev`, open PR against `dev` when done
- Use `gh` CLI for all GitHub operations
- Commit scope matches the ticket — no drive-by refactors
- Report blockers to the board; do not guess past them
- Conventional commit format required: `type(scope): message`
- Use Fibonacci points to estimate ticket complexity
