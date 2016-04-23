import os, sys, platform

file_path = os.path.realpath(__file__)
spl = file_path.split('/')
base_dir = '/'.join(spl[:-1])

pf = platform.system()
if pf == 'Linux':
    os_ = 'linux'
elif pf == 'Darwin':
    os_ = 'osx'
else:
    print "%s platform not supported"
    sys.exit()

vienna_version = "2.1.9"

RESOURCE_DIR = os.path.join(base_dir, 'resources')
STRATEGY_DIR = os.path.join(base_dir, 'strategies')
TEMP_DIR = os.path.join(base_dir, 'tmp')
VIENNA_DIR = os.path.join(RESOURCE_DIR, 'vienna', os_)
NUPACK_DIR = os.path.join(RESOURCE_DIR, 'nupack', os_)

