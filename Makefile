irf_path='/home/lena/Documents/CTA'
n_transient=2
num_slices_per_part=20
num_slices=60
transient_template_filename=random
random_flag=-r



build/threshold_studies/n200_s60_t2_th5_alert.hdf5 : transient_alert.py build/n200_s60_t2_denoised.hdf5 | build
	python transient_alert.py build/n200_s60_t2_denoised.hdf5 --output_path build/threshold_studies -t 5

build/threshold_studies/n200_s60_t2_th10_alert.hdf5 : transient_alert.py build/n200_s60_t2_denoised.hdf5 | build
	python transient_alert.py build/n200_s60_t2_denoised.hdf5 --output_path build/threshold_studies -t 10

build/threshold_studies/n200_s60_t2_th15_alert.hdf5 : transient_alert.py build/n200_s60_t2_denoised.hdf5 | build
	python transient_alert.py build/n200_s60_t2_denoised.hdf5 --output_path build/threshold_studies -t 15

build/threshold_studies/n200_s60_t2_th20_alert.hdf5 : transient_alert.py build/n200_s60_t2_denoised.hdf5 | build
	python transient_alert.py build/n200_s60_t2_denoised.hdf5 --output_path build/threshold_studies -t 20


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean
