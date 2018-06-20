import numpy as np
import subprocess as sp
import signal

FFMPEG="ffmpeg"
BUFFER_SIZE=10**8


class Video_Reader(object):
    """docstring for Video_Reader"""
    command=None
    pipe=None
    resolution=None

    def __init__(self, filename,resolution):
        super(Video_Reader, self).__init__()
        self.command = [ FFMPEG,
                    '-i', filename,
                    '-f', 'image2pipe',
                    '-pix_fmt', 'rgb24',
                    '-vcodec', 'rawvideo', '-']
        # note: stdout pipes the data we want
        # stderr captures the status updates ffmpeg would print to the screen
        # stdin prevents ffmpeg from capturing keyboard input.
        self.pipe = sp.Popen(self.command, bufsize=BUFFER_SIZE,
                             stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)


        self.resolution=resolution
        self.npixels=resolution[0]*resolution[1]*resolution[2]

    def __iter__(self):
        """make the reader object iterable"""
        return self

    def __next__(self):
        try:
            raw_image = self.pipe.stdout.read(self.npixels)
        except:
            # if something breaks here, the pipe may already be closed
            raise StopIteration

        if (len(raw_image)<self.resolution[0]*self.resolution[1]*self.resolution[2]):
            self.close()
            raise StopIteration

        image = np.fromstring(raw_image, dtype='uint8')

        return image.reshape(self.resolution)

    next=__next__

    def close(self):
        """docstring for close"""
        self.pipe.stdout.close()
        self.pipe.stderr.close()
        self.pipe.send_signal(signal.SIGINT)
        self.pipe.wait()


def main():
    """Print purpose of library"""
    print("This is a library for reading video sequences into python via ffmpeg. ")
    print("Provides the 'Video_Reader' iterator class. ")
    print("Requires ffmpeg be installed. ")

if __name__ == '__main__':
    main()
