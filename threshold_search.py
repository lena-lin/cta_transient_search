import transient_alert
import glob
from tqdm import tqdm
from click.testing import CliRunner


def run_threshold_transients():
    files = glob.glob('build/wavelet_studies/n200_s60_t*_trigger.hdf5')
    print(files)
    runner = CliRunner()

    for template_trigger in tqdm(files):
        for th in range(1, 26):
            result = runner.invoke(transient_alert.main, [template_trigger, '--output_path', 'build/wavelet_studies/grid_search', '-t', str(th)])
            if result.exit_code != 0:
                print(result.exit_code)


def run_threshold_background():
    background_trigger = 'build/wavelet_studies/nNone_s60_tNone_trigger.hdf5'
    runner = CliRunner()
    for th in tqdm(range(1, 26)):
        result = runner.invoke(transient_alert.main, [background_trigger, '--output_path', 'build/wavelet_studies/grid_search', '-t', str(th)])
        if result.exit_code != 0:
            print(result.exit_code)

def main():
    run_threshold_transients()
    run_threshold_background()

if __name__ == '__main__':
    main()
