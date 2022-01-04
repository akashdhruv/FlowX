import pymaple
import os

if __name__ == "__main__":
    """ 
    """
    # Set base directory
    basedir = os.getenv('PWD')+'/..'

    # Create a container for the application
    flowx = pymaple.Maple(image='akashdhruv/flowx:debug',container='flowx',
                          source=basedir,target='/home/mount/flowx')

    flowx.build()
    flowx.pour()
    flowx.execute('cd tools/bubblebox && ./setup develop && ./setup clean')
    flowx.commit()
    flowx.rinse()
