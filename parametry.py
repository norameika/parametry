
# coding:utf-8

# Script Name   : parametry.py
# Author        : YY
# Created       : 2017 1228
# Last Modified : 2017 1228
# Modifications : 

# Description   : pragmatical 3d model creation

import numpy as np


class mesh:
    def __init__(self, name, ext="stl", cls=1, dnt=0):
        self.nodes = np.array([])
        self.csecs = np.array([])
        self.cells = []
        self.name = name
        self.ext = ext
        self.csec_close = cls
        self.donut = dnt

    def emboy(self, func, k=1):
        """
        fucs: list or generator to genenrate cross sections
        """
        self.csecs = np.array([np.array(i) for i in func])
        self.csec = self.csecs[0]
        self.nodes = self.csecs.reshape((len(self.csecs)*len(self.csec), 3))

        self.n = len(self.csec)
        self.m = len(self.csecs)
        if self.csec_close:
            ns = list(zip(range(self.n), list(range(1, self.n)) + [0]))
        else:
            ns = list(zip(range(self.n), range(1, self.n)))

        if self.donut:
            mlist = list(range(len(self.csecs)))
        else:
            mlist= list(range(len(self.csecs) -1))

        if self.ext == "stl":
            for l in mlist:
                for e0, e1 in ns:
                    e0 += self.n * l
                    e1 += self.n * l
                    e2 = (e1 + self.n) % (self.n * self.m)
                    e3 = (e0 + self.n) % (self.n * self.m)
                    self.cells += [[e0, e1, e2], [e2, e3, e0]]

        elif self.ext == "obj":
            for l in mlist:
                for e0, e1 in ns:
                    e0 += self.n * l
                    e1 += self.n * l
                    e2 = (e1 + self.n) % (self.n * self.m)
                    e3 = (e0 + self.n) % (self.n * self.m)
                    self.cells += [[e0, e1, e2, e3]]

        if not self.donut:
            if self.ext == "stl":
                nimax = -1
                nimax = len(self.csecs[0]) // 2
                a = list(range(self.n))[:nimax]
                b = list(range(self.n))[nimax:][::-1]
                nilst = list(sum(zip(a, b), ()))
                for ni0, ni1, ni2 in zip(nilst, nilst[1:], nilst[2:]):
                    if ni0 < ni1:
                        self.cells += [[ni0, ni1, ni2]]
                    else:
                        self.cells += [[ni1, ni0, ni2]]

                a = list(range(self.n * (self.m - 1), self.n * (self.m - 0)))[:nimax]
                b = list(range(self.n * (self.m - 1), self.n * (self.m - 0)))[nimax:][::-1]
                nilst = list(sum(zip(a, b), ()))
                for ni0, ni1, ni2 in zip(nilst, nilst[1:], nilst[2:]):
                    if ni0 < ni1:
                        self.cells += [[ni1, ni0, ni2]]
                    else:
                        self.cells += [[ni0, ni1, ni2]]

                if not len(self.csecs[0]) % 2 == 0:
                    self.cells += [[self.n // 2 + 1, self.n // 2, self.n // 2 - 1]]
                    self.cells += [[self.n * (self.m - 1) + self.n // 2, self.n * (self.m-1) + self.n // 2 + 1, self.n * (self.m - 1) + self.n // 2 - 1]]

            elif self.ext == "obj":
                nimax = len(self.csecs[0]) // 2
                a = list(range(self.n))[:nimax]
                b = list(range(self.n))[nimax:][::-1]
                nilst = list(sum(zip(a, b), ()))
                for ni0, ni1, ni2, ni3 in list(zip(nilst, nilst[1:], nilst[2:], nilst[3:]))[::2]:
                    self.cells += [[ni0, ni1, ni3, ni2]]

                a = list(range(self.n * (self.m - 1), self.n * (self.m - 0)))[:nimax]
                b = list(range(self.n * (self.m - 1), self.n * (self.m - 0)))[nimax:][::-1]
                nilst = list(sum(zip(a, b), ()))
                for ni0, ni1, ni2, ni3 in list(zip(nilst, nilst[1:], nilst[2:], nilst[3:]))[::2]:
                    self.cells += [[ni0, ni1, ni3, ni2][::-1]]
                if not len(self.csecs[0]) % 2 == 0:
                    print("cross section must have even N of nodes")
                    raise Exception

    def merge(self, mesh):
        self.cells = self.cells + (np.array(mesh.cells) + len(self.nodes)).tolist()
        self.nodes = np.append(self.nodes, mesh.nodes, axis=0)

    def rotate(self, theta, axis):
        self.nodes = funcs.rotate(theta, self.nodes, axis=axis)

    def shift(self, shift):
        self.nodes =self.nodes + np.array(shift)

    def describe(self):
        print('file name %s' % self.name)
        print('# of nodes %s' % len(self.nodes))
        print('# of cells %s' % len(self.cells))

    def export_stl(self):
        f = open('./stl/%s.stl' % (self.name), 'w')
        f.write('solid %s\n' % self.name)
        for cell in self.cells:
            p0, p1, p2 = [self.nodes[n] for n in cell]
            norm = np.cross(p1 - p0, p2 - p0)
            norm = norm / np.linalg.norm(norm)
            f.write('facet normal %f %f %f\n' % (norm[0], norm[1], norm[2]))
            f.write('outer loop\n')
            for n in cell:
                p = np.array(self.nodes[n])
                f.write('vertex %f %f %f\n' % (p[0], p[1], p[2]))
            f.write('endloop\n')
            f.write('endfacet\n')
        f.close()

    def export_obj(self):
        vn = []
        f = []
        for cell in self.cells:
            p0, p1, p2, p3 = [self.nodes[n] for n in cell]
            norm = np.cross(p1 - p0, p2 - p0)
            norm = (norm / np.linalg.norm(norm)).tolist()
            f += [[p0, p1, p2, p3]]
            vn += [norm]
        f = open('./obj/%s.obj' % (self.name), 'w')
        f.write('o %s\n' % self.name)
        f.write('usemtl None\n')
        f.write('s off\n')
        for (p0, p1, p2) in self.nodes.tolist():
            f.write('v %s %s %s\n' %(p0, p1, p2))
        for (p0, p1, p2) in vn:
            f.write('vn %s %s %s\n' %(p0, p1, p2))
        for cnt, (n0, n1, n2, n3) in enumerate(self.cells):
            f.write('f %s//%s %s//%s %s//%s %s//%s\n' % (n0+1, cnt+1, n1+1, cnt+1, n2+1, cnt+1, n3+1, cnt+1))
        f.close()

    def export(self):
        if self.ext == "stl":
            self.export_stl()
        elif self.ext == "obj":
            self.export_obj()

    def load_stl(self, fp, m="r"):
        f = open(fp, m)
        nodes = {(-999, -999, -999): -1}
        for l in f:
            if l.strip().startswith('solid'):
                self.name = l.strip('\n').split(' ')[-1]
            elif l.strip().startswith('facet normal'):
                next(f)
                n0 = next(f)
                n1 = next(f)
                n2 = next(f)
                n0 = map(float, n0.strip().strip('\n').split(' ')[-3:])
                n1 = map(float, n1.strip().strip('\n').split(' ')[-3:])
                n2 = map(float, n2.strip().strip('\n').split(' ')[-3:])

                keys = []
                for n in list(map(tuple, [n0, n1, n2])):
                    if n not in nodes:
                        key = max(nodes.values()) + 1
                        nodes.update({n: key})
                        keys.append(key)
                    else:
                        keys.append(nodes[n])

                self.cells.append(keys)
        del nodes[(-999, -999, -999)]
        self.nodes = np.array([list(n[0]) for n in sorted(nodes.items(), key=lambda x: x[1])])

if __name__ == '__main__':
    pass
