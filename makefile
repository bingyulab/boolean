install:
	pip install --upgrade pip
# 	conda install -c potassco clingo
	pip install --force-reinstall git+https://github.com/hklarner/pyboolnet
	pip install --no-cache-dir numpy==1.25.2 pandas==1.4.0
	pip install -r requirements.txt
	Rscript Install/install.R

test:
	rm network_analysis.log
	python Step_03_Performance.py 

toy: 
	rm network_analysis.log
	python Step_03_Performance.py -p -d toy
	python Step_04_Plot.py -d toy

dream:
	rm network_analysis.log
	python Step_03_Performance.py -d dream -p
	python Step_04_Plot.py -d dream

.PHONY: install