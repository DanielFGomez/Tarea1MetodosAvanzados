all : Plots

Plots : UpwindGodunov_finalstate.dat plot_Profiles.py Analytic.dat
		python plot_Profiles.py

Analytic.dat : Analytic.py
		python Analytic.py

UpwindGodunov_finalstate.dat : sodtube.x plot_Profiles.py
		./sodtube.x

sodtube.x : sodtube.c 
		gcc -o sodtube.x sodtube.c -lm

clean : 
	rm sodtube.x
	rm UpwindGodunov_finalstate.dat
	rm Analytic.dat
	rm -r Plots
