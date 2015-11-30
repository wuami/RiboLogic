import os

file_path = os.path.realpath(__file__)
spl = file_path.split('/')
base_dir = '/'.join(spl[:-1])

os = 'linux'

RESOURCE_DIR = os.path.join(base_dir, 'resources')
STRATEGY_DIR = os.path.join(base_dir, 'strategies')
VIENNA_DIR = os.path.join(RESOURCE_DIR, 'vienna', os)
NUPACK_DIR = os.path.join(RESOURCE_DIR, 'nupack', os)

