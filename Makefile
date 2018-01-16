.PHONY: all clean ready

all: data/globom

data/globom: sources/tinycompo.hpp
	@cd sources && make --no-print-directory -j8

sources/tinycompo.hpp:
	@curl https://raw.githubusercontent.com/vlanore/tinycompo/master/tinycompo.hpp > $@

clean:
	@cd sources && make --no-print-directory clean

ready:
	@cd sources && make format
	@git status
