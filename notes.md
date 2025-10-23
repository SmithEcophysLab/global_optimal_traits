# global optimal traits notes

## june 16, 2025
- recently updated optimal vcmax code to include Rd
- now want to run model globally for deciduous and evergreen C3 plants and produce
interpretable PCA plots
	- problems
		1. LMA not set up to do evergreen option
			- **FIXED!**
		2. no input for "f" for deciduous LMA
			- **FIXED!**
- next steps
	- make plots with deciduous and evergreen and correct traits
		- need to think about best traits to look at

## june 25, 2025
- fixed package dependencies
- received PCA plot script from Alissar (see email)
- next steps
	- filter out non-vegetated areas and redo plots
	
## june 27, 2025
- created figures with only vegetated land areas from modis
- next steps
	- figure out what we can glean from all this!
		- chi and photosynthetic capacity are orthogonal
			- nuances LCT
			- photosynthesis is in the middle (of course?)
		- n does different things for deciduous and evergreen
	- added c4

## aug 29, 2025
- interest from LT group about this
- for within site analyses
	- randomly sample from parameter distributions using mean and sd of data?
		- cheaib beta values
		- others?

## sep 5, 2025
- main points and corresponding figures/tables
	1. EEO can be used to predict global distributions of optimal photosynthetic strategies
		- figure: PCA plot
		- analysis: explain the predicted strategy/spatial variation
		- stage: mostly done
	2. Traits separate between PFTs in some cases, but not others
		- figure: trait histograms by PFT
		- analysis: confirmation of PFT strategies across environments
		- stage: simulations done, need to make figures
	3. within-site variability can be predicted through simulating variability in costs and environment
		- figure: PCA for a select number of NEON sites representing:
			- c4 grassland
			- temp deciduous
			- temp evergreem
			- temp mixed
		- analysis: TBD, where is there the most variability?
		- stage: simulations needed
	4. future condition impact on strategies
		- figures: similar plots with future temperature and CO2 conditions (overlay PCAs??)
		- analysis: similar strategies, but maybe in different part of axis?
		- stage: simulations still needed
- to do: same simulations with added within site variability
	- do this using mean and SD random draw for beta distributions from Cheaib
		- for a select few sites: neon sites?

## oct 23, 2025
- main points revisited with status update
	1. EEO can be used to predict global distributions of optimal photosynthetic strategies
		- figure: global PCA plot with 3 PFTs for Anet, gsw, nphoto, chi, lma, nue, wue, rd25
		- analysis: explain the predicted strategy/spatial variation
		- stage: plots complete
			- [PC1PC2 plot](results/plots/global_optimal_traits_all_pca_plot_PC1PC2.jpeg)
			- [PC1PC3 plot](global_optimal_traits_all_pca_plot_PC1PC3.jpeg)
	2. Traits separate between PFTs in some cases, but not others
		- figure: trait histograms by PFT for Anet, gsw, nphoto, chi, lma, nue, wue, rd25
		- analysis: confirmation of PFT strategies across environments
		- stage: plots complete
			- [global trait historgram plot](results/plots/global_optimal_traits_hist_all.jpeg)
- to do: work on points 3 and 4 (start line 688)












	