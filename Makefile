irf_path=''
n_transient=10
num_slices_per_part=20
num_slices=60
percentage_transient_noise = 2# as a percentage 


all: build/evaluation_score.txt

build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_cube.hdf5 : simulate_random_cube.py | build
	python simulate_random_cube.py -n $(n_transient) -s $(num_slices_per_part) -p $(percentage_transient_noise)

build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_trans.hdf5 : build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_cube.hdf5 | build

build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_denoised.hdf5 : analyse_random_cube.py build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_cube.hdf5 | build
	python analyse_random_cube.py build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_cube.hdf5

build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_alert.hdf5 : random_transient_alert.py build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_denoised.hdf5 | build
	python random_transient_alert.py build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_denoised.hdf5

build/evaluation_score.txt : evaluation.py build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_trans.hdf5 build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_alert.hdf5 | build
	python evaluation.py build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_trans.hdf5 build/n$(n_transient)_s$(num_slices)_p$(percentage_transient_noise)_alert.hdf5 evaluation_score.txt

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean
