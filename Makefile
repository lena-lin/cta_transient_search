irf_path='/home/lena/Dokumente/CTA'
n_transient=2
num_slices_per_part=20
num_slices=60
transient_template_index=2

all: build/evaluation_score.txt

build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_cube.hdf5 : simulate_cube.py | build
	python simulate_cube.py -n $(n_transient) -s $(num_slices_per_part) -temp $(transient_template_index)

build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_trans.hdf5 : build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_cube.hdf5 | build

build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_denoised.hdf5 : analyse_cube.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_cube.hdf5 | build
	python analyse_cube.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_cube.hdf5

build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_alert.hdf5 : transient_alert.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_denoised.hdf5 | build
	python transient_alert.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_denoised.hdf5

build/evaluation_score.txt : evaluation.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_trans.hdf5 build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_alert.hdf5 | build
	python evaluation.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_trans.hdf5 build/n$(n_transient)_s$(num_slices)_t$(transient_template_index)_alert.hdf5 evaluation_score.txt

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean
