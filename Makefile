irf_path='/home/lena/Documents/CTA'
n_transient=2
num_slices_per_part=20
num_slices=60
random_transient_template=True
transient_template_filename=random

all: build/evaluation_score.txt

build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5 : simulate_cube.py | build
	python simulate_cube.py -n $(n_transient) -s $(num_slices_per_part) -rand_temp $(random_transient_template) --irf_path $(irf_path)

build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trans.hdf5 : build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5 | build

build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5 : analyse_cube.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5 | build
	python analyse_cube.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5

build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_alert.hdf5 : transient_alert.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5 | build
	python transient_alert.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5

build/evaluation_score.txt: build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trans.hdf5

build/evaluation_score.txt: build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_alert.hdf5

build/evaluation_score.txt: evaluation.py  | build
	python evaluation.py \
	  build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trans.hdf5 \
		build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_alert.hdf5 \
		evaluation_score.txt

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean
