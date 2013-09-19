# recursively process `Makefile.am` in these subdirectories
SUBDIRS = src

.PHONY: tests

tests:
	make
	nosetests -v

clean-tests:
	rm -rf tests/test_sample_run/*/output_* tests/test_sample_run/*/_SUCCESSFUL_RUN* tests/test_sample_run/*/geotop.log
