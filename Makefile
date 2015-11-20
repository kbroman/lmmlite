R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

all: doc example_data
.PHONY: doc example_data

# build package documentation
doc:
	R -e 'devtools::document()'

example_data: data/recla.RData

data/recla.RData: inst/extdata/create_example_data.R
	cd $(@D);R $(R_OPTS) -e "source('../$<')"
