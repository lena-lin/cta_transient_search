#irf_path='/home/lena/Documents/CTA' ## Lena
#irf_path= '/home/jana/Schreibtisch/Projekt_Master' ## Jana
irf_path= '..' ## Jana auf Vollmond
n_transient = 20
num_slices_per_part=20
num_slices=60
transient_template_filename=random
random_flag=-r
threshold = 5
Ra = 4
Dec = 5

all: build/evaluation_20_$(threshold)_Ra$(Ra)_Dec$(Dec).txt

build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5 : simulate_cube.py | build
	python simulate_cube.py -f $(irf_path) -n $(n_transient) -s $(num_slices_per_part) $(random_flag) -x $(Ra) -y $(Dec) -p False


build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trans.hdf5 : build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5 | build

build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5 : analyse_cube.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5 | build
	python analyse_cube.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_cube.hdf5


build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trigger.hdf5 : build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5 | build

build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_thr$(threshold)_alert.hdf5 : transient_alert.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5 | build
	python transient_alert.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_denoised.hdf5 --output_path build -t $(threshold)

build/evaluation_20_$(threshold)_Ra$(Ra)_Dec$(Dec).txt: build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trans.hdf5

build/evaluation_20_$(threshold)_Ra$(Ra)_Dec$(Dec).txt: build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_thr$(threshold)_alert.hdf5

build/evaluation_20_$(threshold)_Ra$(Ra)_Dec$(Dec).txt: evaluation.py  | build
	python evaluation.py build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_trans.hdf5 \
		build/n$(n_transient)_s$(num_slices)_t$(transient_template_filename)_thr$(threshold)_alert.hdf5

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean
