#irf_path='/home/lena/Documents/CTA' ## Lena
#irf_path= '/home/jana/Schreibtisch/Projekt_Master' ## Jana
irf_path= '..' ## Jana auf Vollmond
n_transient = 4
num_slices_per_part=20
num_slices=60
transient_template_filename=random
random_flag=-r
cu_bool = 1



all: build/evaluation_score.txt

build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5 : simulate_cube.py | build
	python simulate_cube.py 3 -f $(irf_path) -n $(n_transient) -s $(num_slices_per_part) $(random_flag)


build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trans.hdf5 : build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5 | build

build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5 : analyse_cube.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5 | build
	python analyse_cube.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5


build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trigger.hdf5 : build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5 | build

build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_alert.hdf5 : transient_alert.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5 | build
	python transient_alert.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5

build/evaluation_score.txt: build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trans.hdf5

build/evaluation_score.txt: build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_alert.hdf5

build/evaluation_score.txt: evaluation.py  | build
	python evaluation.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trans.hdf5 \
		build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_alert.hdf5 \ 3

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean
