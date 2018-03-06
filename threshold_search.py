import transient_alert
import glob
from tqdm import tqdm
from click.testing import CliRunner


def main():
    files = glob.glob('build/n200_s60_t*_trigger.hdf5')
    runner = CliRunner()

    for template in tqdm(files):
        for th in range(1, 26):
            result = runner.invoke(transient_alert.main, [template, '--output_path', 'build/threshold_studies/grid_search', '-t', str(th)])
            if result.exit_code != 0:
                print(result.output)


if __name__ == '__main__':
    main()
