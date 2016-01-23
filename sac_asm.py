from pycca.asm import *
import numpy as np


def sac(arrays, max_delta):
    histogram = np.zeros(max_delta*2 + 1, dtype=int)
    # Python equivalent
    #for i in range(len(arrays)):
        #for j in range(i+1, len(arrays)):
            #arr1 = arrays[i]
            #arr2 = arrays[j]
            #for m in range(len(arr1)):
                #for n in range(len(arr2)):
                    #dt = arr1[m] - arr2[n]
                    #if dt < -max_delta:
                        #break
                    #if dt <= max_delta:
                        #histogram[dt + max_delta] += 1
    
    addrs = np.array([arr.ctypes.data for arr in arrays])
    sizes = np.array([len(arr) for arr in arrays])
    
    fn = mkfunction([
        
        #------------------------------------ for i in range(len(arrays)):
        mov(rsi, -1),  # rsi === i
        label('for_i'),
        add(rsi, 1),
        cmp(rsi, len(addrs)),
        jge('end_for_i'),
        
        #------------------------------------     for j in range(i+1, len(arrays)):
        mov(rdi, -1), # rsi),  # rdi === j
        label('for_j'),
        add(rdi, 1),
        cmp(rdi, rsi),
        je('for_j'),
        cmp(rdi, len(addrs)),
        jge('end_for_j'),

        # set up array pointers for i, j
        mov(r8, addrs.ctypes.data),
        mov(r8, qword([rsi*8 + r8])),
        mov(r9, addrs.ctypes.data),
        mov(r9, qword([rdi*8 + r9])),
        
        #------------------------------------ for m in range(len(arr1)):
        mov(rcx, -1),  # rcx holds arr1 index
        label('for_m'),
        add(rcx, 1),
        mov(r10, sizes.ctypes.data),
        cmp(rcx, qword([rsi*8 + r10])),
        jge('end_for_m'),
        
        #------------------------------------ for n in range(len(arr2)):
        mov(rdx, -1),  # rdx holds arr2 index
        label('for_n'),
        add(rdx, 1),
        mov(r10, sizes.ctypes.data),
        cmp(rdx, qword([rdi*8 + r10])),
        jge('end_for_n'),
        
        #------------------------------------ dt = arr1[m] - arr2[n]
        mov(r11, qword([r8+rcx*8])),
        sub(r11, qword([r9+rdx*8])),
        
        #------------------------------------ if dt < -max_delta:
        #------------------------------------    break
        cmp(r11, -max_delta),
        jl('end_for_n'),
        
        #------------------------------------ if dt <= max_delta:
        #------------------------------------     histogram[dt + max_delta] += 1
        cmp(r11, max_delta),
        jg('after_increment'),
        mov(r10, histogram.ctypes.data),
        add(qword([r10+r11*8+max_delta*8]), 1),
        label('after_increment'),
        
        jmp('for_n'),
        label('end_for_n'),
        
        jmp('for_m'),
        label('end_for_m'),

        jmp('for_j'),
        label('end_for_j'),
        
        jmp('for_i'),
        label('end_for_i'),
        
        ret(),
    ])
    
    fn()

    return histogram

if __name__ == '__main__':
    arr1 = np.array([0.1, 0.5, 1.2, 1.3, 4.5, 4.7, 7.8])
    arr2 = np.array([0.13, 0.49, 1.3, 1.33, 4.8, 4.6, 7.9, 8.5, 9.0])
    
    binw = 0.05
    i1 = (arr1/binw).astype(int)
    i2 = (arr2/binw).astype(int)
    max_dif = int(1.0/binw)
    hist = sac([i1, i2], max_dif)
    import pyqtgraph as pg
    bins = np.arange(-max_dif-0.5, max_dif+1.5)
    pg.plot(bins, hist, stepMode=True, fillBrush='y', fillLevel=0)
