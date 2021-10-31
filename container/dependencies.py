import pymaple
import os

basedir = os.getenv('PWD')

# checkout dependencies
#os.system('git submodule update --init {0}/tools/bubblebox && \
#           cd {0}/tools/bubblebox && \
#           python3 dependencies.py'.format(basedir))


# Build container
flowx = pymaple.Maple(image='akashdhruv/flowx:latest',container='flowx_latest',
                      source=basedir,target='/home/mount/flowx')

flowx.build()
flowx.pour()
flowx.execute('"cd tools/bubblebox && \
                ./setup develop && ./setup clean"')
flowx.commit()
flowx.rinse()
