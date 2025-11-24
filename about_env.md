
# Understanding Python Environments, `.venv`, `pyproject.toml`, and uv

This document gives a clear explanation of how Python environments work, what `.venv` contains, how `pyproject.toml` fits in, what roles lock files play, and how uv manages everything.

---

## What a Python Environment Actually Is

A Python environment is:

1. A specific Python interpreter
2. A private set of installed packages for that interpreter

This lets different projects use different versions of libraries without interfering with each other. It also makes analyses reproducible.

---

## What `.venv` Is

`.venv` is simply a folder.
It is the directory that **contains the virtual environment** for your project.

Inside, you will find:

* A Python interpreter
* A `site-packages` folder containing all installed packages
* Scripts such as `python`, `pip`, `ipython`, and others

When you activate it:

```bash
source .venv/bin/activate
```

your shell starts using:

* `.venv/bin/python`
* `.venv/bin/pip`

instead of the system Python.

You install packages into that isolated environment instead of globally.

---

## What `pyproject.toml` Is

`pyproject.toml` is **a configuration file**, not an environment.

It describes:

* The project (name, version, metadata)
* Required Python version
* Dependencies your project needs
* Which build system is used (for example hatch, setuptools, poetry)

Example:

```toml
[project]
name = "vcc"
version = "0.1.0"
requires-python = ">=3.10"

dependencies = [
  "pandas",
  "numpy",
  "matplotlib",
]

[project.optional-dependencies]
dev = [
  "pytest",
  "ruff",
]
```

Think of it like a recipe:

* `pyproject.toml` defines *what* should be installed
* `.venv` is *where* it gets installed

---

## What Lock Files Are

Lock files freeze exact versions of everything installed, including sub-dependencies.

Examples:

* `uv.lock`
* `poetry.lock`
* `Pipfile.lock`
* `requirements.txt` (when used as a freeze list)

Purpose:

* Ensures someone else can recreate **exactly** the same environment
* Ensures repeatability across machines and over time

Your project uses `uv.lock`.

---

## Other Common Environment Files

* `.python-version`
  Used by tools like pyenv to select a specific Python version for this folder.

* `.env`
  Stores environment variables for configuration.

These are optional.

---

## What uv Is (Versus a Python Environment)

uv is **not** an environment.

uv is a *tool* that:

* Creates and manages environments
* Installs packages
* Resolves dependency versions
* Updates `pyproject.toml` and `uv.lock`
* Runs commands inside the correct environment

It wraps and improves on the traditional `python -m venv` + `pip install` workflow.

### Commands and what they actually do

#### `uv venv`

* Creates `.venv`
* Ensures correct Python version is available

Equivalent idea to:

```bash
python -m venv .venv
```

but faster and more reliable.

---

#### `uv add PACKAGE_NAME`

uv does several things:

1. Adds the package to `pyproject.toml`
2. Resolves and chooses exact versions of dependencies
3. Updates `uv.lock`
4. Installs packages into `.venv`

---

#### `uv pip list`

Shows which packages are installed in the uv-managed environment.
Equivalent to running `pip list` inside an activated `.venv`.

---

#### `uv run python script.py`

* Ensures dependencies are installed
* Executes Python inside the environment
* No need to activate `.venv` manually

---

## How Everything Fits Together in a Project

**Workflow:**

1. Define dependencies in `pyproject.toml`
2. Create environment with `uv venv`
3. Install dependencies with `uv add ...`
4. Run code using the environment (`uv run python script.py`)
5. Commit `pyproject.toml` and `uv.lock`
6. `.venv` can always be deleted and recreated

Someone else cloning the project does:

```bash
uv sync
```

uv reads `pyproject.toml` + `uv.lock`, creates `.venv`, installs the exact same package versions.

---

## Short Summary

* A Python environment is just a Python interpreter plus its installed packages.
* `.venv` is the folder holding that environment.
* `pyproject.toml` is the project configuration and dependency declaration.
* `uv.lock` stores exact versions to make the environment reproducible.
* uv manages `.venv`, dependencies, locking, and running commands inside the environment.

This setup makes your project clean, reproducible, and easy to share.

---

If you want, I can also generate a companion file explaining **how to use uv daily**, including a minimal checklist for adding dependencies, syncing environments, and running scripts.
