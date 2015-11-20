all: pylmm_results.csv lmmlite_results.csv

pylmm_results.csv: try_pylmm.py lmm.py kinship.csv
	$< > $@

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

kinship.csv: write_data_to_files.R
	R $(R_OPTS) -e "source('$<')"

lmmlite_results.csv: try_lmmlite.R
	R $(R_OPTS) -e "source('$<')"

lmm.py:
	curl -O https://raw.githubusercontent.com/nickFurlotte/pylmm/master/pylmm/lmm.py
	autopep8 -i lmm.py    # install with "pip install lmm.py" (was getting error about tabs)

clean:
	rm *.csv lmm.py
