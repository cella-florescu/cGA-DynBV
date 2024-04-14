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
    if (plotting == 0):
        return int(stdout)

    numbers_1 =  [int(x) for x in stdout.strip().split("\n")[0].strip().split()]
    numbers_2 =  [int(x) for x in stdout.strip().split("\n")[1].strip().split()]
    numbers_3 =  [int(x) for x in stdout.strip().split("\n")[2].strip().split()]
    numbers_4 =  [float(x) for x in stdout.strip().split("\n")[3].strip().split()]

    return [numbers_1, numbers_2, numbers_3, numbers_4]

n = 500
k = 1000
numbers_1, numbers_2, numbers_3, numbers_4 = run_algorithm(n, k, 1, 1)



plt.plot(numbers_1, label="low")
plt.plot(numbers_2, label="mid")
plt.plot(numbers_3, label="high")
plt.plot(numbers_4, label="potential")
plt.legend()
plt.xlabel('Iteration')
plt.ylabel('Number of Marginal Frequencies')
plt.title(f'Evolution of Marginal Frequencies over Time for N = {n} and K = {k}')
plt.show()

