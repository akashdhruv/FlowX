import pymaple
import os

basedir = os.getenv('PWD')+'/..'

# Build container
flowx = pymaple.Maple(image='akashdhruv/flowx:latest',container='flowx_latest',
                      source=basedir,target='/home/mount/flowx')

flowx.build()
flowx.pour()
flowx.execute('"cd tools/bubblebox && \
                ./setup develop && ./setup clean"')
flowx.commit()
flowx.rinse()
