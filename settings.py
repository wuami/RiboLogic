import os, sys, platform

file_path = os.path.realpath(__file__)
spl = file_path.split('/')
base_dir = '/'.join(spl[:-1])

pf = platform.system()
if pf == 'Linux':
    os_dir = 'linux'
elif pf == 'Darwin':
    os_dir = 'osx'
else:
    print "%s platform not supported"
    sys.exit()

RESOURCE_DIR = os.path.join(base_dir, 'resources')
STRATEGY_DIR = os.path.join(base_dir, 'strategies')
VIENNA_DIR = os.path.join(RESOURCE_DIR, 'vienna', os_dir)
NUPACK_DIR = os.path.join(RESOURCE_DIR, 'nupack', os_dir)

