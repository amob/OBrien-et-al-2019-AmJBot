# OBrien-et-al-2019-AmJBot
Repository includes data and analyses for O'Brien et al 2019, American Journal of Botany.

R Script: "BZS1analysis_final_annotated.R"
Goal: Analyze experimental data
Requires input files:
  "StartDatBZS1Nov16.csv" (raw data file from image analysis)
  "StartMapBZS1Nov16.csv" (raw data file from image analysis)
  "Nov20datBZS1.csv" (raw data file from image analysis)
  "Nov20mapBZS1.csv" (raw data file from image analysis)
  "Nov23datBZS1.csv" (raw data file from image analysis)
  "Nov23mapBZS1.csv" (raw data file from image analysis)
  "EndDatBZS1Nov26.csv" (raw data file from image analysis)
  "EndMapBZS1Nov26.csv" (raw data file from image analysis)
  "SaltBTduckweeddesTRTS.csv" (experimental design)
  "AO BZS.1 12mar2019_sorted.csv" (optical density data)
  "BZT in BZS1 for trts pooled within geno.csv" (benzotriazole concentration data)
  "ORB Samples all genos.csv" (treatment mapping file for wells sampled for benzotriazole concentration)
Creates output files:
  "BZS1_Trait_results_plot.pdf" (Figure 2)
  "qqPix_liveplantsAS.pdf" (Supplmenentary figure, Appendix S4)
  "saltsurvdatonly.pdf" (Figure 3)
  "summarygrowthres.pdf" (Figure 4)
  "timeseriesPixFitAS.pdf" (Supplementary figure, Appendix S6)
  "BZTinBZS1_trts_within_geno_plusgrowth.csv" (dataframe formatted in useful way)
  "treatmenteffsBZTloss.pdf" (Figure 5)
  *also creates r objects (linear models, primarily), or output to terminal used in manuscript tables or text*
Depends on libraries: MCMCglmm


R Script:
Goal: Analyze microbiome data
Requires input files:
Creates output files:
