.PHONY: environment remove-env clean build_sys # .PHONY is something we can add when our target dependencies are not files.

ENVIRONMENT=discussion-9

environment: remove-env
	conda env create -f environment.yaml

remove-env:
	conda remove --name $(ENVIRONMENT) --all --yes