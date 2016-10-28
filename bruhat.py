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
        self.gen = gen
        self.rel = rel

        reduced = {'':('',)} # map word -> sorted tuple of equivalent reduced words
        lookup = {('',):''} # map sorted tuple -> word
        for g in gen:
            reduced[g] = (g,)
            lookup[(g,)] = g
            reduced[g+g] = reduced['']
            for h in gen:
                if g==h:
                    continue
                m = rel[g, h]
                ghm = ''.join(([g, h]*m)[:m])
                hgm = ''.join(([h, g]*m)[:m])
                r = [ghm, hgm]
                r.sort()
                r = tuple(r)
                reduced[ghm] = reduced[hgm] = r
                lookup[r] = ghm
        self.reduced = reduced
        self.lookup = lookup

    def get_canonical(self, orig):
        reduced = self.reduced
        if orig in reduced:
            return reduced[orig]
        gen = self.gen
        rel = self.rel
        items = set([orig])
        done = False
        #print "E"
        while not done:
            done = True
            for word in list(items):
                for word1 in cancel_pairs(word):
                    items = set([word1])
                    assert len(word1)<len(word)
                    #print "X"
                    done = False
                    break
                else:
                    continue
                break

            for s, t in distinct_pairs(gen):
                m = rel[s, t]
                stm = ''.join(([s, t]*m)[:m])
                tsm = ''.join(([t, s]*m)[:m])
                for word in list(items):
                    if stm in word:
                        for word1 in rewrite(word, stm, tsm):
                            if word1 not in items:
                                items.add(word1)
                                #print "Y"
                                done = False
        #print "F"
        items = list(items)
        items.sort()
        items = tuple(items)
        reduced[orig] = items
        lookup = self.lookup
        if items not in lookup:
            lookup[items] = orig
        elif len(lookup[items]) > orig:
            lookup[items] = orig
        return items

    def is_equal(self, lhs, rhs):
        if lhs==rhs:
            return True
        left = self.get_canonical(lhs)
        right = self.get_canonical(rhs)
        return left == right

    def build(self, max_size=None):
        lookup = self.lookup # map canonical -> word
        group = set(lookup.keys()) # canonical forms
        #print "group:", group
        pairs = [(i, j) for i in group for j in group]
        mul = {} # map (i,j) -> i*j
        while 1:
            newpairs = []
            #print "pairs:", pairs
            for i, j in pairs:
                g = lookup[i]
                h = lookup[j]
                gh = g+h # multiply words
                k = self.get_canonical(gh) # updates lookup aswell
                gh = lookup[k]
                mul[g, h] = gh
                if k not in group:
                    newpairs += [(k, g) for g in group]
                    newpairs += [(g, k) for g in group]
                    group.add(k)
            if not newpairs:
                break
            pairs = newpairs
        self.mul = mul
        self.group = lookup.values()
        return self.group



def main():

    A_2 = Coxeter("LP", {("L", "P") : 3})

    A_2.build()
    print A_2.group

    A_3 = Coxeter("LPS", {("L", "P") : 3, ("L", "S"):3})

    print len(A_3.build())

    A_4 = Coxeter("LPSH", {("L", "P") : 3, ("L", "S"):3, ("S", "H"):3})

    print len(A_4.build(max_size=120))




if __name__ == "__main__":

     if argv.profile:
        import cProfile as profile
        profile.run("main()")

     else:
        main()


