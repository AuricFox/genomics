# GUI terminal
# Runs main program that calls other scripts

import sys
import os
import sequence as sq
import alignment as al
import de_bruijn as db
import utils

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "temp")


def getCodons(file):
    data = utils.get_data(os.path.join(path, file))
    seq = sq.sequence(data[0][0], data[1])
    return seq.get_data()

# ==============================================================================================================
def main():
    getCodons("./input/sars_spike_protein_assembled.fna")
    return


if __name__ == "__main__":
    main()
