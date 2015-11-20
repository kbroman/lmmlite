R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

all: doc example_data vignette build/vignette.rds
.PHONY: doc example_data

# build package documentation
doc:
	R -e 'devtools::document()'

example_data: data/recla.RData

data/recla.RData: inst/extdata/create_example_data.R
	cd $(@D);R $(R_OPTS) -e "source('../$<')"

VIGNETTE=inst/doc/compare2pylmm.html
vignette: $(VIGNETTE)

inst/doc/%.html: vignettes/%.Rmd inst/compare2pylmm/pylmm_results.csv
	cp $< $(@D)
	cd $(<D);R -e "rmarkdown::render('$(<F)')"
	mv $(<D)/$(@F) $(@D)

build/vignette.rds: vignettes/make_vignette_index.R $(VIGNETTE)
	cd $(<D);R -e "source('$(<F)')"

inst/compare2pylmm/pylmm_results.csv: inst/compare2pylmm/try_pylmm.py inst/compare2pylmm/lmm.py inst/compare2pylmm/kinship.csv
	cd $(<D);$(<F) > ($@F)

inst/compare2pylmm/kinship.csv: inst/compare2pylmm/write_data_to_files.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

inst/compare2pylmm/lmmlite_results.csv: inst/compare2pylmm/try_lmmlite.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

inst/compare2pylmm/lmm.py:
	cd $(@D);curl -O https://raw.githubusercontent.com/nickFurlotte/pylmm/master/pylmm/lmm.py
	cd $(@D);autopep8 -i lmm.py    # install with "pip install lmm.py" (was getting error about tabs)
