import multiprocessing
from multiprocessing import Pool
import PlotFunc as pf

leap = 250
to = 50001
if __name__ == '__main__':  
    try:
        nProc = int(multiprocessing.cpu_count() / 2)
        with Pool(nProc) as p:
            results = [p.apply_async(pf.qc_qr_th_u_w, (t, )) for t in range(0, to, leap)]
            final = [result.get() for result in results]
    except:
        print("finish1")
        
    try:
        nProc = int(multiprocessing.cpu_count() / 2)
        with Pool(nProc) as p:
            results = [p.apply_async(pf.zeta, (t, )) for t in range(0, to, leap)]
            final = [result.get() for result in results]
    except:
        print("finish2")