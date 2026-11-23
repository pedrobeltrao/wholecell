"""Legacy Jurkat structure-library builders, ported and path-parameterised. Since B3 (2026-11-17) their worker
functions run inside the core library stage (wholecell/library/jurkat.py sets their module paths to the build, makes
them offline and replaces their `main` drivers); the `__main__` blocks are kept for reference only. Scripts import
each other by bare name (run with this folder on sys.path)."""
