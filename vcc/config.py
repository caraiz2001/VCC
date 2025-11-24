
from pathlib import Path
import os


# Resolve project root directory from this file's location
PROJECT_ROOT_DIRECTORY_PATH = Path(__file__).resolve().parent.parent


DATA_DIRECTORY_PATH = PROJECT_ROOT_DIRECTORY_PATH / "data"
RESULTS_DIRECTORY_PATH = PROJECT_ROOT_DIRECTORY_PATH / "results"
PLOTS_DIRECTORY_PATH = PROJECT_ROOT_DIRECTORY_PATH / "plots"

# Optional environment variable overrides
CUSTOM_RESULTS_DIRECTORY = os.getenv("VCC_RESULTS_DIRECTORY")
CUSTOM_PLOTS_DIRECTORY = os.getenv("VCC_PLOTS_DIRECTORY")

if CUSTOM_RESULTS_DIRECTORY is not None:
    RESULTS_DIRECTORY_PATH = Path(CUSTOM_RESULTS_DIRECTORY)

if CUSTOM_PLOTS_DIRECTORY is not None:
    PLOTS_DIRECTORY_PATH = Path(CUSTOM_PLOTS_DIRECTORY)

# Ensure directories exist
for path in [
    DATA_DIRECTORY_PATH,
    RESULTS_DIRECTORY_PATH,
    PLOTS_DIRECTORY_PATH,
]:
    path.mkdir(parents=True, exist_ok=True)

DEFAULT_RANDOM_SEED_VALUE = 42
DEFAULT_NUMBER_OF_PARALLEL_WORKERS = 4
