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