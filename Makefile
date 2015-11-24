R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

all: assets/compare2pylmm.html         \
	 assets/performance.html           \
	 assets/lmmlite.html               \
	 compare2pylmm/pylmm_results.csv   \
	 compare2pylmm/pylmm_llvals.csv    \
	 compare2pylmm/lmmlite_results.csv \
	 compare2pylmm/lmmlite_llvals.csv

assets/compare2pylmm.html: assets/compare2pylmm.Rmd assets/lmm.py
	cd $(<D);R $(R_OPTS) -e "rmarkdown::render('$(<F)')"

assets/performance.html: assets/performance.Rmd assets/lmm.py
	cd $(<D);R $(R_OPTS) -e "rmarkdown::render('$(<F)')"

compare2pylmm/pylmm_results.csv: compare2pylmm/try_pylmm.py compare2pylmm/lmm.py compare2pylmm/kinship.csv
	cd $(<D);$(<F) > $(@F)

compare2pylmm/kinship.csv: compare2pylmm/write_data_to_files.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

compare2pylmm/lmmlite_results.csv: compare2pylmm/try_lmmlite.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

assets/lmm.py: compare2pylmm/lmm.py
	cp $< $@

compare2pylmm/lmm.py:
	cd $(@D);curl -O https://raw.githubusercontent.com/nickFurlotte/pylmm/master/pylmm/lmm.py
	cd $(@D);autopep8 -i lmm.py    # all with "pip install lmm.py" (was getting error about tabs)

compare2pylmm/pylmm_llvals.csv: compare2pylmm/pylmm_calcll.py compare2pylmm/lmm.py compare2pylmm/kinship.csv
	cd $(<D);$(<F) > $(@F)

compare2pylmm/lmmlite_llvals.csv: compare2pylmm/lmmlite_calcll.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

assets/lmmlite.html: assets/lmmlite.Rmd
	cd $(<D);R $(R_OPTS) -e "rmarkdown::render('$(<F)')"

assets/lmmlite.Rmd: ../lmmlite/vignettes/lmmlite.Rmd
	cp $< $@
