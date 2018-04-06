import numpy as np
import analyse_cube
import glob
from tqdm import tqdm
from click.testing import CliRunner


def main():
    runner = CliRunner()

    for w in np.linspace(4, 20, 5):
        print(w)
        result = runner.invoke(analyse_cube.main, ['build/n10_s60_t2_cube.hdf5', '-w', str(int(w))])
        if result.exit_code != 0:
            print(result.output)


if __name__ == '__main__':
    main()
