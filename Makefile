.PHONY: all clean ready globom

all: globom

globom: sources/tinycompo.hpp
	@cd sources && make --no-print-directory -j8

sources/tinycompo.hpp:
	@curl https://raw.githubusercontent.com/vlanore/tinycompo/master/tinycompo.hpp > $@

clean:
	@cd sources && make --no-print-directory clean

format:
	@cd sources && make --no-print-directory format

ready: all
	@cd sources && make format
	@git status
