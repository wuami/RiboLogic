import os

file_path = os.path.realpath(__file__)
spl = file_path.split("/")
base_dir = "/".join(spl[:-1])

RESOURCE_DIR = os.path.join(base_dir, "resources")
PUZZLE_DIR = os.path.join(RESOURCE_DIR, "puzzles")
VIENNA_DIR = os.path.join(RESOURCE_DIR, "vienna", "linux")

