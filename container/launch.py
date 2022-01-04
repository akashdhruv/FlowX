import pymaple
import os

if __name__ == "__main__":
    """
    """
    # Set basedir
    basedir = os.getenv('PWD') + '/..'

    # Create a container for the application
    flowx = pymaple.Maple(image='akashdhruv/flowx:latest',container='flowx',
                          source=basedir,target='/home/mount/flowx')
    flowx.build()
    flowx.pour()
    flowx.execute('cd tools/bubblebox && ./setup develop && ./setup clean')

    flowx.execute('"./setup develop && ./setup clean && python3 tests/container && \
                     jupyter notebook --port=8888 --no-browser --ip=0.0.0.0"')
    flowx.rinse()
    flowx.clean()
    flowx.remove()
