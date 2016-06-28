
def orbit(char *src):
    assert len(src)==9
    #print "orbit", src
    items = []

    cdef char target[9+1]
    cdef char target1[9+1]
    cdef char tmp
    target[9] = 0
    target1[9] = 0
    
    # di, dj = (0, 0)
    target[0] = src[0]
    target[1] = src[1]
    target[2] = src[2]
    target[3] = src[3]
    target[4] = src[4]
    target[5] = src[5]
    target[6] = src[6]
    target[7] = src[7]
    target[8] = src[8]
    # idx =  (0, 0, 0)
    target1[0] = target[0]
    target1[1] = target[1]
    target1[2] = target[2]
    target1[3] = target[3]
    target1[4] = target[4]
    target1[5] = target[5]
    target1[6] = target[6]
    target1[7] = target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (0, 1, 1)
    target1[0] = target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 0, 1)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 1, 0)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # di, dj = (0, 1)
    target[1] = src[0]
    target[2] = src[1]
    target[0] = src[2]
    target[4] = src[3]
    target[5] = src[4]
    target[3] = src[5]
    target[7] = src[6]
    target[8] = src[7]
    target[6] = src[8]
    # idx =  (0, 0, 0)
    target1[0] = target[0]
    target1[1] = target[1]
    target1[2] = target[2]
    target1[3] = target[3]
    target1[4] = target[4]
    target1[5] = target[5]
    target1[6] = target[6]
    target1[7] = target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (0, 1, 1)
    target1[0] = target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 0, 1)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 1, 0)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # di, dj = (0, 2)
    target[2] = src[0]
    target[0] = src[1]
    target[1] = src[2]
    target[5] = src[3]
    target[3] = src[4]
    target[4] = src[5]
    target[8] = src[6]
    target[6] = src[7]
    target[7] = src[8]
    # idx =  (0, 0, 0)
    target1[0] = target[0]
    target1[1] = target[1]
    target1[2] = target[2]
    target1[3] = target[3]
    target1[4] = target[4]
    target1[5] = target[5]
    target1[6] = target[6]
    target1[7] = target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (0, 1, 1)
    target1[0] = target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 0, 1)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 1, 0)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # di, dj = (1, 0)
    target[3] = src[0]
    target[4] = src[1]
    target[5] = src[2]
    target[6] = src[3]
    target[7] = src[4]
    target[8] = src[5]
    target[0] = src[6]
    target[1] = src[7]
    target[2] = src[8]
    # idx =  (0, 0, 0)
    target1[0] = target[0]
    target1[1] = target[1]
    target1[2] = target[2]
    target1[3] = target[3]
    target1[4] = target[4]
    target1[5] = target[5]
    target1[6] = target[6]
    target1[7] = target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (0, 1, 1)
    target1[0] = target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 0, 1)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 1, 0)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # di, dj = (1, 1)
    target[4] = src[0]
    target[5] = src[1]
    target[3] = src[2]
    target[7] = src[3]
    target[8] = src[4]
    target[6] = src[5]
    target[1] = src[6]
    target[2] = src[7]
    target[0] = src[8]
    # idx =  (0, 0, 0)
    target1[0] = target[0]
    target1[1] = target[1]
    target1[2] = target[2]
    target1[3] = target[3]
    target1[4] = target[4]
    target1[5] = target[5]
    target1[6] = target[6]
    target1[7] = target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (0, 1, 1)
    target1[0] = target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 0, 1)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 1, 0)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # di, dj = (1, 2)
    target[5] = src[0]
    target[3] = src[1]
    target[4] = src[2]
    target[8] = src[3]
    target[6] = src[4]
    target[7] = src[5]
    target[2] = src[6]
    target[0] = src[7]
    target[1] = src[8]
    # idx =  (0, 0, 0)
    target1[0] = target[0]
    target1[1] = target[1]
    target1[2] = target[2]
    target1[3] = target[3]
    target1[4] = target[4]
    target1[5] = target[5]
    target1[6] = target[6]
    target1[7] = target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (0, 1, 1)
    target1[0] = target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 0, 1)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 1, 0)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # di, dj = (2, 0)
    target[6] = src[0]
    target[7] = src[1]
    target[8] = src[2]
    target[0] = src[3]
    target[1] = src[4]
    target[2] = src[5]
    target[3] = src[6]
    target[4] = src[7]
    target[5] = src[8]
    # idx =  (0, 0, 0)
    target1[0] = target[0]
    target1[1] = target[1]
    target1[2] = target[2]
    target1[3] = target[3]
    target1[4] = target[4]
    target1[5] = target[5]
    target1[6] = target[6]
    target1[7] = target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (0, 1, 1)
    target1[0] = target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 0, 1)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 1, 0)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # di, dj = (2, 1)
    target[7] = src[0]
    target[8] = src[1]
    target[6] = src[2]
    target[1] = src[3]
    target[2] = src[4]
    target[0] = src[5]
    target[4] = src[6]
    target[5] = src[7]
    target[3] = src[8]
    # idx =  (0, 0, 0)
    target1[0] = target[0]
    target1[1] = target[1]
    target1[2] = target[2]
    target1[3] = target[3]
    target1[4] = target[4]
    target1[5] = target[5]
    target1[6] = target[6]
    target1[7] = target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (0, 1, 1)
    target1[0] = target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 0, 1)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 1, 0)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # di, dj = (2, 2)
    target[8] = src[0]
    target[6] = src[1]
    target[7] = src[2]
    target[2] = src[3]
    target[0] = src[4]
    target[1] = src[5]
    target[5] = src[6]
    target[3] = src[7]
    target[4] = src[8]
    # idx =  (0, 0, 0)
    target1[0] = target[0]
    target1[1] = target[1]
    target1[2] = target[2]
    target1[3] = target[3]
    target1[4] = target[4]
    target1[5] = target[5]
    target1[6] = target[6]
    target1[7] = target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (0, 1, 1)
    target1[0] = target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 0, 1)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = target[1]
    target1[2] = ord('0')+ord('1')-target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = target[4]
    target1[5] = ord('0')+ord('1')-target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = target[7]
    target1[8] = ord('0')+ord('1')-target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    # idx =  (1, 1, 0)
    target1[0] = ord('0')+ord('1')-target[0]
    target1[1] = ord('0')+ord('1')-target[1]
    target1[2] = target[2]
    target1[3] = ord('0')+ord('1')-target[3]
    target1[4] = ord('0')+ord('1')-target[4]
    target1[5] = target[5]
    target1[6] = ord('0')+ord('1')-target[6]
    target1[7] = ord('0')+ord('1')-target[7]
    target1[8] = target[8]
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[2]
    target1[2] = tmp
    tmp = target1[3]
    target1[3] = target1[5]
    target1[5] = tmp
    tmp = target1[6]
    target1[6] = target1[8]
    target1[8] = tmp
    items.append(target1)
    tmp = target1[0]
    target1[0] = target1[6]
    target1[6] = tmp
    tmp = target1[1]
    target1[1] = target1[7]
    target1[7] = tmp
    tmp = target1[2]
    target1[2] = target1[8]
    target1[8] = tmp
    items.append(target1)
    return items

def gauge(char *src):
    assert len(src)==9
    #print "orbit", src
    items = []

    cdef char target[9+1]
    target[9] = 0
    
    # (0, 0)
    target[0] = ord('0')+ord('1')-src[0]
    target[1] = ord('0')+ord('1')-src[1]
    target[2] = src[2]
    target[3] = src[3]
    target[4] = src[4]
    target[5] = src[5]
    target[6] = src[6]
    target[7] = src[7]
    target[8] = src[8]
    items.append(target)
    # (0, 1)
    target[0] = src[0]
    target[1] = ord('0')+ord('1')-src[1]
    target[2] = ord('0')+ord('1')-src[2]
    target[3] = src[3]
    target[4] = src[4]
    target[5] = src[5]
    target[6] = src[6]
    target[7] = src[7]
    target[8] = src[8]
    items.append(target)
    # (0, 2)
    target[0] = ord('0')+ord('1')-src[0]
    target[1] = src[1]
    target[2] = ord('0')+ord('1')-src[2]
    target[3] = src[3]
    target[4] = src[4]
    target[5] = src[5]
    target[6] = src[6]
    target[7] = src[7]
    target[8] = src[8]
    items.append(target)
    # (1, 0)
    target[0] = src[0]
    target[1] = src[1]
    target[2] = src[2]
    target[3] = ord('0')+ord('1')-src[3]
    target[4] = ord('0')+ord('1')-src[4]
    target[5] = src[5]
    target[6] = src[6]
    target[7] = src[7]
    target[8] = src[8]
    items.append(target)
    # (1, 1)
    target[0] = src[0]
    target[1] = src[1]
    target[2] = src[2]
    target[3] = src[3]
    target[4] = ord('0')+ord('1')-src[4]
    target[5] = ord('0')+ord('1')-src[5]
    target[6] = src[6]
    target[7] = src[7]
    target[8] = src[8]
    items.append(target)
    # (1, 2)
    target[0] = src[0]
    target[1] = src[1]
    target[2] = src[2]
    target[3] = ord('0')+ord('1')-src[3]
    target[4] = src[4]
    target[5] = ord('0')+ord('1')-src[5]
    target[6] = src[6]
    target[7] = src[7]
    target[8] = src[8]
    items.append(target)
    # (2, 0)
    target[0] = src[0]
    target[1] = src[1]
    target[2] = src[2]
    target[3] = src[3]
    target[4] = src[4]
    target[5] = src[5]
    target[6] = ord('0')+ord('1')-src[6]
    target[7] = ord('0')+ord('1')-src[7]
    target[8] = src[8]
    items.append(target)
    # (2, 1)
    target[0] = src[0]
    target[1] = src[1]
    target[2] = src[2]
    target[3] = src[3]
    target[4] = src[4]
    target[5] = src[5]
    target[6] = src[6]
    target[7] = ord('0')+ord('1')-src[7]
    target[8] = ord('0')+ord('1')-src[8]
    items.append(target)
    # (2, 2)
    target[0] = src[0]
    target[1] = src[1]
    target[2] = src[2]
    target[3] = src[3]
    target[4] = src[4]
    target[5] = src[5]
    target[6] = ord('0')+ord('1')-src[6]
    target[7] = src[7]
    target[8] = ord('0')+ord('1')-src[8]
    items.append(target)
    return items
