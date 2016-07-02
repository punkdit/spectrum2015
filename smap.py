#!/usr/bin/env python

import sys
import fcntl, termios, struct, os


class SMap(object):
    def __init__(self):
        self.data = {}
        self.rows = 0
        self.cols = 0

    def __setitem__(self, key, s):
        i, j = key
        for idx, c in enumerate(s):
            self.data[i, j+idx] = c
        self.rows = max(self.rows, i)
        self.cols = max(self.cols, j+idx)

    def __getitem__(self, key):
        return self.data.get(key, ' ')

    def __str__(self):
        lines = []
        for row in range(self.rows+1):
            line = ''
            for col in range(self.cols+1):
                line += self.data.get((row, col), ' ')
            lines.append(line)
        return '\n'.join(lines)


def get_cols():
    fd = 0 # stdin
    cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
    return cr[1]


class FlowPrint(object):
    def __init__(self, stdout=sys.stdout, space=True):
        self.lines = []
        self.col = 0
        self.stdout = stdout
        self.space = space
        #self.blank_line = False

    def output(self, s):
        #if self.blank_line and not s.strip() and s.endswith('\n'):
        #    return
        #self.blank_line = False
        self.stdout.write(s)

    def write(self, *args):
        s = ' '.join(str(s) for s in args)
        #nl = s.endswith('\n')
        nl = s == '\n' or not s
        if nl:
            s = s[:-1]
        lines = s.split('\n')
        cols = max(len(line) for line in lines)+1 # add a space...

        #debug = "%s:%s"%(self.lines, lines)
        #print >>self.stdout, debug
        #print >>self.stdout, "write(%r)"%s

        col = self.col + cols

        if col >= get_cols():
            for line in self.lines:
                self.output(line+'\n')
            if self.space:
                self.output('\n')
            self.lines = []
            self.col = 0

        idx = 0
        while idx < len(self.lines) or idx < len(lines):
            if idx < len(self.lines):
                line = self.lines[idx]
                assert len(line) == self.col
            else:
                line = ' '*self.col
                self.lines.append(line)
            extra = lines[idx] if idx < len(lines) else ''
            extra = extra + ' '*(cols-len(extra))
            self.lines[idx] = line + extra
            idx += 1

        self.col += cols

        if nl:
            for line in self.lines:
                self.output(line+'\n')
            if self.space:
                self.output('\n')
            self.lines = []
            self.col = 0

        for idx, line in enumerate(self.lines):
            assert len(line) == self.col
            assert '\n' not in line

        #print self.lines

flprint = FlowPrint().write

def install_flprint():
    sys.stdout = FlowPrint(sys.stdout)


def test():
    install_flprint()

    for i in range(10):
        #flprint('hi\nthere...\n')
        print "hi\nthere!",
        print "hi\nthere!\nworld!",
        if i == 2:
            print


if __name__=="__main__":
    test()







