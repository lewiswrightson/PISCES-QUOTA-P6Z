# PISCES-QUOTA-P6Z

### Project Title:

PISCES QUOTA P6Z includes explicit diazotrophy for the standard P5Z PISCES QUOTA framework

Able to switch between two prevalent marine diazotrophs, _Trichodesmium_ and _Crocosphaera_

Able to account for the temperature dependence of the nitrogen fixation elemental use efficiencies of both Iron (Fe) and Phosphorus (P) for both diazotrophs

Nitrogen fixation is facultative within the model

Diazotroph thermal performance curves for growth and elemental use efficiencies based on Jiang et al. (2018) and Yang et al. (2021) for _Trichodesmium_ and _Crocosphaera_, respectively

Jiang, H.-B., Fu, F.-X., Rivero-Calle, S., Levine, N. M., Sañudo-Wilhelmy, S. A., Qu, P.-P., Wang, X.-W., Pinedo-Gonzalez, P., Zhu, Z., & Hutchins, D. A. (2018). Ocean warming alleviates iron limitation of marine nitrogen fixation. Nature climate change, 8(8), 709-712. doi:10.1038/s41558-018-0216-8


Yang, N., Merkel, C. A., Lin, Y.-A., Levine, N. M., Hawco, N. J., Jiang, H.-B., Qu, P.-P., DeMers, M. A., Webb, E. A., Fu, F.-X., & Hutchins, D. A. (2021). Warming iron-limited oceans enhance nitrogen fixation and drive biogeographic specialization of the globally important cyanobacterium Crocosphaera. Frontiers in Marine Science, 8(118). doi:10.3389/fmars.2021.628363

### Papers:


### Prerequisites / Getting Started:

svn: http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0 r11143

Apply the below to the /cfg/ directory:

MY_SRC directory

EXPREF directory (namelists, file_def, field_def and other xml files)


_In namelist_pisces_ref_p6z:_

Select flag ln_p6z to select explicit diazotrophy

Flag ln_tricho chooses diazotroph setup: _Trichodesmium_ (True) and _Crocosphaera_ (False)

Flag ln_tiue controls the temperature dependence of the nitrogen fixation Fe use efficiency 

Flag ln_tiue controls the temperature dependence of the nitrogen fixation P use efficiency 


### Installing: 

The configuration is compiled using ./makenemo -r P6Z

### Authors:

Lewis Wrightson and Alessandro Tagliabue 

### Acknowledgements:

Conducted during a Ph.D. study supported by the Natural Environment Research Council (NERC) EAO Doctoral Training Partnership funded by NERC (NE/L002469/1)

Funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation program (grant agreement no. 724289)


