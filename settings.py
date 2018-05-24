import os, sys, platform

file_path = os.path.realpath(__file__)
spl = file_path.split('/')
base_dir = '/'.join(spl[:-1])

RESOURCE_DIR = os.path.join(base_dir, 'resources')
STRATEGY_DIR = os.path.join(base_dir, 'strategies')
TEMP_DIR = os.path.join(base_dir, 'tmp')
VIENNA_DIR = '%s/bin' % os.environ['VIENNAHOME']
NUPACK_DIR = '%s/bin' % os.environ['NUPACKHOME']

