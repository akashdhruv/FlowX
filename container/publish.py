import pymaple
import os

if __name__ == "__main__":
    basedir = os.getenv('PWD')+'/..'

    # Publish container for use with other applications
    flowx = pymaple.Maple(image='akashdhruv/flowx:latest',container='flowx',
                          source=basedir,target='/home/mount/flowx')

    flowx.build()
    flowx.pour()
    flowx.execute('cd tools/bubblebox && ./setup develop && ./setup build && \
                                      ./setup install && ./setup clean')

    flowx.execute('./setup develop && ./setup build && ./setup install && ./setup clean')
    flowx.commit()
    flowx.rinse()
    flowx.push('akashdhruv/flowx:publish')
    flowx.clean()
    flowx.remove()
