import sys
import time


class ResidualVis(object):
    def __init__(self):
        print "new residual"


if __name__ == '__main__':
    modif = None
    logfile = open(sys.argv[1], 'r')
    lastpos = logfile.tell()
    for line in logfile:
        if line.startswith("=========="):
            lastpos = logfile.tell() - len(line)
            print lastpos
    logfile.seek(lastpos)

    while True:
        for line in logfile.readlines():
            if line.strip() == "==========":
                vis = ResidualVis()
        time.sleep(2)

