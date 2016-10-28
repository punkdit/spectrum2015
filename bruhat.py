#!/usr/bin/env python

import sys, os

from argv import argv


def pairs(gen):
    for s in gen:
        for t in gen:
            yield (s, t)

def distinct_pairs(gen):
    for s in gen:
        for t in gen:
            if s != t:
                yield (s, t)

def rewrite(word, src, tgt):
    if src not in word:
        return
    assert src != tgt
    stem = ''
    while src in word:
        idx = word.index(src)
        n = len(src)
        yield stem + word[:idx] + tgt + word[idx+n:]
        stem, word = word[:idx+1], word[idx+1:]

assert list(rewrite("zxxxw", 'xx', 'y')) == ['zyxw', 'zxyw']


def cancel_pairs(word):
    for i in range(len(word)-1):
        if word[i]==word[i+1]:
            yield word[:i] + word[i+2:]

assert list(cancel_pairs('xyzzz')) == ['xyz', 'xyz']


class Coxeter(object):
    """ Coxeter group
    """

    def __init__(self, gen, rel, identity=''):
        """
        gen : list of generators
        rel : map pairs (i,j) of generators to m_{ij} (this is the Coxeter matrix)
        """
        for i in gen:
          for j in gen:
            if rel.get((i, j)) and not rel.get((j, i)):
                rel[j, i] = rel[i, j]
        for i in gen:
          for j in gen:
            m = rel.get((i, j))
            if m is None:
                rel[i, j] = 2 # default 
        for i in gen:
          for j in gen:
            assert rel[i, j] == rel[j, i]
            assert rel[i, j] in (2, 3, 4, 6)
#        rules = {} # map word -> word
#        for i in gen:
#            rules[i+i] = identity
#        for i in gen:
#          for j in gen:
#            m = rel[i, j]
#            word = (i+j)*m
#            rules[word] = identity
        self.derived = {}
        self.gen = gen
        self.rel = rel

    def get_derived(self, word):
        #items = self.derived.setdefault(word, set())
        #if items:
        #    return items
        gen = self.gen
        rel = self.rel
        items = set()
        items.add(word)
        while 1:
            n = len(items)
            for s, t in distinct_pairs(gen): # only need s<t
                m = rel[s, t]
                stm = ''.join(([s, t]*m)[:m])
                tsm = ''.join(([t, s]*m)[:m])
                for word in list(items):
                    if stm in word:
                        for word1 in rewrite(word, stm, tsm):
                            items.add(word1)
            for word in list(items):
                for word1 in cancel_pairs(word):
                    items.add(word1)
            if len(items)==n:
                break
        return items

    def is_equal(self, lhs, rhs):
        if lhs==rhs:
            return True
        left = self.get_derived(lhs)
        right = self.get_derived(rhs)
        if left.intersection(right):
            return True
        return False

    def mul(self, a, b):
        w = a + b
        w = self.rewrite(w)
        return w

    def build(self, max_size=None):
        group = set([])
        canonical = {} # map group element -> derived set
        group.add('')
        canonical[''] = self.get_derived('')
        for g in self.gen:
            group.add(g)
            canonical[g] = self.get_derived(g)
        done = False
        while not done:
            done = True
            for a, b in pairs(list(group)):
                c = a+b
                if c in group:
                    continue
                derived = self.get_derived(c)
                for other in canonical.values():
                    if derived.intersection(other):
                        break
                else:
                    group.add(c)
                    if max_size and len(group)>=max_size:
                        return group
                    canonical[c] = derived
                    done = False
            #print len(group)
        return group


def main():

    A_2 = Coxeter("LP", {("L", "P") : 3})

    print A_2.build()

    A_3 = Coxeter("LPS", {("L", "P") : 3, ("L", "S"):3})

    #print A_3.build()

    A_4 = Coxeter("LPSH", {("L", "P") : 3, ("L", "S"):3, ("S", "H"):3})

    print A_4.build(max_size=120)




if __name__ == "__main__":

     if argv.profile:
        import cProfile as profile
        profile.run("main()")

     else:
        main()


