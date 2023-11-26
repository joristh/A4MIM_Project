all: section2_2 section2_5

section2_2: section2_2.m plots
	matlab -batch "run('./section2_2.m')"

section2_5: section2_5.m plots
	matlab -batch "run('./section2_5.m')"

plots:
	mkdir plots

clean:
	rm plots/*.pdf
