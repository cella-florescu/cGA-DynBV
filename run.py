import subprocess
import matplotlib.pyplot as plt
from io import StringIO
import numpy as np
import multiprocessing


def run_algorithm(n,k, dbv, plotting):
    command = ["./a.out"]
    process = subprocess.Popen(
        command,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True  # This assumes you are working with text, not binary data
    )
    input_string = str(n) + " " + str(k) + " " + str(dbv) + " " + str(plotting)

    stdout, stderr = process.communicate(input=input_string)
    return_code = process.wait()
    if (plotting != 1):
        return int(stdout)

    numbers_1 =  [int(x) for x in stdout.strip().split("\n")[0].strip().split()]
    numbers_2 =  [int(x) for x in stdout.strip().split("\n")[1].strip().split()]
    numbers_3 =  [int(x) for x in stdout.strip().split("\n")[2].strip().split()]

    return [numbers_1, numbers_2, numbers_3]

def get_runs(n, k, num_runs, mode):
    results = []

    for i in range(num_runs):
        results.append(run_algorithm(n, k, 1, mode))

    print(results)
    np_array = np.array(results)
    return np.median(np_array)



K_range = list(range(6, 21, 1))
num_runs = 50
n = 300

medians = []

mode = 0

for k in K_range:
    medians.append(get_runs(n, k, num_runs, mode))

file = open("safe.txt", "w+")

content = str(medians)
file.write(content)
file.close()



plt.plot(K_range, medians)
plt.xlabel('K')
if mode == 0:
    plt.ylabel('Number of Iterations')
    plt.title(f"Median Time Until Convergence for N = {n}")
elif mode == 2:
    plt.ylabel('Number of Bits')
    plt.title(f"Median Number of Bits at Lower Border for N = {n}")
else:
    plt.ylabel('Number of Bits')
    plt.title(f"Median Number of Bits from Upper to Lower Border for N = {n}")
plt.show()

