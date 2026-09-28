# Onboarding: machine setup and Git workflow

This guide takes you from a fresh laptop to your first merged pull request. Once set up, follow
[GETTING_STARTED.md](GETTING_STARTED.md) to work on a package with the Claude Code agent. No prior Git or
Linux experience is assumed. Follow it in order and copy the commands exactly. Lines starting
with `#` are comments; you don't type them.

Budget about an hour for part B, most of it downloads.

- **Part A** is for the repository owner (a few clicks).
- **Part B** is a one-time setup for each teammate.
- **Part C** is the daily workflow.
- **Part D** covers common problems.

---

## Part A: repository owner (once per teammate)

1. The teammate creates a free account at <https://github.com/signup> and sends you their
   **GitHub username**.
2. Open <https://github.com/harsha31bs07-blip/samaya/settings/access>. Choose **Add people**,
   enter the username, and give the **Write** role.
3. The teammate accepts the invite from the email or from
   <https://github.com/harsha31bs07-blip/samaya/invitations>.
4. Recommended, so nobody pushes to `main` directly:
   - Go to **Settings → Branches → Add branch ruleset** (or "Add rule") for `main`.
   - Tick **Require a pull request before merging** and **Require status checks to pass**,
     then choose the CI checks.
5. Tell the teammate which **work package** they own; the list is in
   [ARCHITECTURE.md §6](ARCHITECTURE.md#6-work-packages).

---

## Part B: teammate setup (once)

### B1. Get a Linux shell

The code builds with Linux tools. On Windows, use WSL (Linux inside Windows).

**Windows 10/11:**
1. Open **PowerShell as Administrator** (Start → type "PowerShell" → right-click → *Run as
   administrator*) and run:
   ```powershell
   wsl --install -d Ubuntu-24.04
   ```
2. Restart when asked. An **Ubuntu** window opens, or you can start it from the Start menu.
   Choose a Linux username and password; the password is needed for `sudo`.
3. Everything below is typed in this Ubuntu window.

**macOS:** open **Terminal** and install Homebrew from <https://brew.sh>. Where the steps below
use `apt`, use the brew commands given instead.

**Linux (Ubuntu/Debian):** just open a terminal.

### B2. Install the tools

Ubuntu / WSL:

```sh
sudo apt update
sudo apt install -y git cmake ninja-build g++ clang python3 python3-venv curl unzip bc
```

macOS:

```sh
xcode-select --install          # Apple compiler
brew install git cmake ninja python
```

Check the versions. You need CMake 3.22 or newer and g++ 11 or newer (or clang 14 or newer):

```sh
cmake --version && g++ --version | head -1 && git --version
```

### B3. Install the GitHub command-line tool and log in

Ubuntu / WSL (these are GitHub's official install commands; paste the whole block):

```sh
sudo mkdir -p -m 755 /etc/apt/keyrings
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg \
  | sudo tee /etc/apt/keyrings/githubcli-archive-keyring.gpg > /dev/null
sudo chmod go+r /etc/apt/keyrings/githubcli-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" \
  | sudo tee /etc/apt/sources.list.d/github-cli.list > /dev/null
sudo apt update && sudo apt install -y gh
```

macOS: `brew install gh`.

Then log in:

```sh
gh auth login
#   GitHub.com → HTTPS → "Yes" to authenticate Git → "Login with a web browser"
#   Copy the one-time code, press Enter, paste it in the browser page, approve.
gh auth setup-git
```

Tell Git who you are. Use the same email as your GitHub account:

```sh
git config --global user.name  "Your Name"
git config --global user.email "you@example.com"
git config --global pull.rebase true
```

### B4. Get the code and build it

```sh
mkdir -p ~/code && cd ~/code
gh repo clone harsha31bs07-blip/samaya
cd samaya

cmake --preset release            # configure (once)
cmake --build --preset release    # compile, about 1 minute
ctest --preset release            # run the tests: expect "100% tests passed"
```

Try the solver on the small instance in the repo:

```sh
build/release/apps/cli/samaya tests/instances/tiny_lp.mps
```

You should see `Optimal objective 11 ... (verified)`.

### B5. Benchmark tools (Python)

```sh
python3 -m venv ~/venv
source ~/venv/bin/activate            # do this in every new terminal before benchmarking
pip install highspy                   # comparison solver, only used by the harness

bench/fetch_instances.sh netlib       # downloads about 95 public LPs (a few minutes)
python3 bench/harness.py bench/instances/netlib --baseline highspy --time-limit 60
```

The last line should report `samaya solved 93/93` and `disagreements with baselines: 0`.

### B6. Editor (recommended: VS Code)

1. Install VS Code on Windows or macOS from <https://code.visualstudio.com>.
2. Install the extensions **WSL** (Windows only), **C/C++** or **clangd**, and **CMake Tools**.
3. On Windows, run `code .` from the Ubuntu window inside `~/code/mrpl-pro`. VS Code opens
   *inside* WSL, which is what you want. Do not edit the files from Windows Explorer.

### B7. Read before coding

1. [README.md](../README.md): what the project is.
2. [docs/ARCHITECTURE.md](ARCHITECTURE.md): sections 1–5, then **your work package** in
   section 6.
3. The code your package plugs into, listed in the package's "Edits" column.

---

## Part C: daily workflow

Each work package lives on its own **branch**. You never commit to `main` directly; changes
reach `main` through a **pull request (PR)** that CI checks and a teammate reviews.

### C1. Start your work package (once)

```sh
cd ~/code/mrpl-pro
git checkout main && git pull            # latest code
git checkout -b wp1-pdlp                 # your branch; use your package's name
```

### C2. The edit → build → test loop

```sh
cmake --build --preset debug             # debug build: warnings are errors, like CI
ctest --preset debug                     # all tests
build/debug/tests/samaya_tests 2>&1 | grep -E "FAIL|passed"   # quicker summary
```

The first time, run `cmake --preset debug` once before building.

Before pushing, also run the sanitizer build, which catches memory bugs:

```sh
cmake --preset asan && cmake --build --preset asan && ctest --preset asan
```

Rules that CI enforces or reviewers check are listed in
[ARCHITECTURE.md §5](ARCHITECTURE.md#5-engineering-rules). The important ones:
- No compiler warnings.
- Every new algorithm is tested against a reference solver.
- Plant a bug once to prove your tests catch it.
- Never skip or disable a failing test.

### C3. Save your work (commit) and upload it (push)

Commit small steps, each with a message saying what changed and why:

```sh
git status                               # what changed
git add src/lp/pdlp.cpp src/lp/pdlp.hpp tests/test_pdlp.cpp CMakeLists.txt tests/CMakeLists.txt
git commit -m "PDLP: plain PDHG iteration with projection onto bounds"
git push -u origin wp1-pdlp              # the first push; later pushes are just: git push
```

Push at least once a day, even unfinished work. It's your backup, and others can see progress.

### C4. Open a pull request when a milestone works

```sh
gh pr create --base main --fill --draft   # draft PR while still working
gh pr checks                              # watch CI (6 jobs: gcc/clang × debug/release/asan)
gh pr ready                               # when the milestone is done and CI is green
```

In the PR description, write:
- what works;
- test results (counts, mismatches);
- how you proved the tests can fail;
- known gaps.

The repository owner reviews and merges.

If a CI job fails, open it (`gh pr checks` shows links), fix the cause locally, commit, and
`git push`. CI reruns automatically.

### C5. Stay up to date with `main`

Other packages merge into `main` while you work. Pull their changes into your branch every few
days:

```sh
git fetch origin
git rebase origin/main          # replay your commits on top of the latest main
#   On a conflict: edit the files it names, then  git add <file> && git rebase --continue
cmake --build --preset debug && ctest --preset debug
git push --force-with-lease     # needed after a rebase; only ever on YOUR branch
```

`src/core/solver.cpp` is shared by every package. Keep your edits there small.

### C6. After your PR is merged

```sh
git checkout main && git pull
git checkout -b wp1-pdlp-m3     # a fresh branch for the next milestone
```

---

## Part D: troubleshooting

| Problem | Fix |
|---|---|
| `cmake: command not found` or a version below 3.22 | Ubuntu 22.04 or newer ships a new-enough CMake; otherwise `pip install cmake` inside the venv |
| `Could not find a configuration file for package ... Ninja` | `sudo apt install ninja-build` |
| Build error `-Werror` in the debug preset | That's a warning. Fix it; don't turn it off. The release preset shows the same warning without failing |
| `ctest` shows a FAIL | Run `build/debug/tests/samaya_tests` directly: failing checks print file:line |
| `git push` asks for a password | Run `gh auth setup-git` again (B3) |
| `permission denied` on push | Ask the owner for Write access (Part A) and accept the invite |
| `rejected (non-fast-forward)` on push | `git pull --rebase`, then push again |
| `pip install` says "externally managed environment" | Activate the venv first: `source ~/venv/bin/activate` |
| WSL is slow | Keep the code inside the Linux home (`~/code`), not under `/mnt/c/...` |
| Out of memory while building | `cmake --build --preset release -j2` |
| GPU work (WP2) | Needs an NVIDIA GPU. On Windows install the NVIDIA driver for WSL, then the *CUDA Toolkit for WSL-Ubuntu* from <https://developer.nvidia.com/cuda-downloads>, then `cmake --preset cuda`. Write and test the CPU version (WP1) first; it needs no GPU |

## Quick reference

```sh
git checkout main && git pull                          # update
git checkout -b <branch>                               # new work
cmake --build --preset debug && ctest --preset debug   # build + test
git add <files> && git commit -m "message"             # save
git push                                               # upload
gh pr create --base main --fill                        # propose for merge
gh pr checks                                           # CI status
git fetch origin && git rebase origin/main             # catch up with main
```
