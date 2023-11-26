all: section2_1 section2_6

section2_1: section2_1.m plots
	matlab -batch "run('./section2_1.m')"

section2_6: section2_6.m plots
	matlab -batch "run('./section2_6.m')"

plots:
	mkdir plots

clean:
	rm plots/*.pdf
