import logging
import sys
from datetime import datetime
load("Arion_permutation_density.sage")

# Set up the log file
now = datetime.now()
timestamp = datetime.timestamp(now)
log_file_name = "Arion_permutation_density_estimation_results_" + str(timestamp) + ".log"
file_handler = logging.FileHandler(filename=log_file_name)
stdout_handler = logging.StreamHandler(sys.stdout)
handlers = [file_handler, stdout_handler]
logging.basicConfig(format="%(asctime)-15s %(levelname)-8s %(message)s", level=logging.DEBUG, handlers=handlers)
print = logging.info

def density_estimate(prime, branches, rounds):
    K = GF(prime)
    arion = ArionDensity(K, branches=branches, rounds=rounds)
    arion.density_of_Arion_permutation()

print("Density estimation for Arion-pi.")
primes = [11, 13, 17, 19, 23, 29, 31, 37]
branches = [3, 4, 5]
rounds = [6] #[3, 4, 5, 6]

for p in primes:
    for b in branches:
        for r in rounds:
            density_estimate(p, b, r)
            print("\n")
