# to make all print statements in a python file automatically flush output
# from flushprint import Flushfile
# import sys
# sys.stdout=Flushfile(sys.stdout)
class Flushfile(object):
    def __init__(self, fd):
        self.fd = fd

    def write(self, x):
        ret=self.fd.write(x)
        self.fd.flush()
        return ret

    def writelines(self, lines):
        ret=self.writelines(line)
        self.fd.flush()
        return ret

    def flush(self):
        return self.fd.flush

    def close(self):
        return self.fd.close()

    def fileno(self):
        return self.fd.fileno()