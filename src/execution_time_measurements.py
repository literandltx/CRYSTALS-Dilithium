import random
import time
import string
from tqdm import tqdm

import main

def measure_execution_time(iterations: int=40):
    measurements_Gen = 0
    measurements_Sign = 0
    measurements_Verify = 0

    for i in tqdm(range(iterations), bar_format="{l_bar}{bar} [ time left: {remaining}, time spent: {elapsed}]"):
        M = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(256))
        
        start = time.time()
        pk, sk = main.Gen()
        end = time.time()
        measurements_Gen += end - start

        start = time.time()
        sigma = main.Sign(sk=sk, M=M)
        end = time.time()
        measurements_Sign += end - start

        start = time.time()
        main.Verify(pk=pk, M=M, sigma=sigma)
        end = time.time()
        measurements_Verify += end - start

    print("Gen:")
    print(f"Average: {measurements_Gen / iterations:0.3f} s.")

    print("Sign:")
    print(f"Average: {measurements_Sign / iterations:0.3f} s.")

    print("Verify:")
    print(f"Average: {measurements_Verify / iterations:0.3f} s.")

if __name__ == "__main__":
    measure_execution_time()