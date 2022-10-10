## README ## 

Code and data for the manuscript:

Kristen Mehalko, Minhoo Kim, Sanjana Paye, Kelly Koh, Ryan J. Lu & Bérénice A. Benayoun
"Lack of accelerated ovarian aging in a follicle-stimulating hormone receptor haploinsufficiency model"
         (biorXiv preprint at: XXXXXXXXXXXXXXXXXXX)

The code is broadly arranged by the RNA species being analyzed:

	# 1_Genotyping      : 
			- All gels and annotation
	
	# 2_Fertility_assay : 
			- 5_Table_S1.Fertility_assay_rawdata.xlsx: Litter size and latency calculation
			- Fertility_assay_age_of_dams.txt: age distribution of dams
			- Fertility_assay_counts.txt, Fertility_assay_latency_V2.txt: parsed breeding information
			- 2022-10-09_Fertility_assay_analysis.R: R script for statistics and plotting
			
	# 3_Ovary_HE        : 
			- Original_images: microscopy images for quantification
			- Fshr_AC_follicle_count_final.txt, Fshr_AC_follicle_count_final_average.txt: counts and aaverags from microscopy
			- 2022-10-09_Fshr_ovary_HE_follicle_count_V2.R: R script for statistics and plotting
			
	# 4_Hormone_data    : 
			- Rawdata: Raw measurements from UVA on serum
			- Plots/2022-10-09_MeMo_Fshr_hormone_quantification_V2.R: R script for statistics and plotting
	
	# 5_RT-qPCR         : 
			- Fshr_RT-qPCR_calculations.xlsx: Raw Ct and Delta-Delta-Ct calculations
			- Fshr_RT-qPCR.txt: parsed expression by sample for analysis
			- Fshr_RT-qPCR.R: R script for statistics and plotting
