def bitWiseMult(a,b):
    s = 0b00000000
    while b:
        if b & 0b00000001:
            s = bitWiseAdd(s, a)
        a <<= 1
        b >>= 1
    return s

def bitWiseAdd(a, b):
    if a&b:
        return bitWiseAdd(a^b, (a&b)<<1)
    return b^a


    

print(bitWiseMult(0b00001110, 0b00001011))
