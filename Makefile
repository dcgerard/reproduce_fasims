# ADJUST THESE VARIABLES AS NEEDED TO SUIT YOUR COMPUTING ENVIRONMENT
# -------------------------------------------------------------------
# This variable specifies the number of threads to use for the
# parallelization. This could also be specified automatically using
# environment variables. For example, in SLURM, SLURM_CPUS_PER_TASK
# specifies the number of CPUs allocated for each task.
nc = 12

# R scripting front-end. Note that makeCluster sometimes fails to
# connect to a socker when using Rscript, so we are using the "R CMD
# BATCH" interface instead.
rexec = R CMD BATCH --no-save --no-restore

# AVOID EDITING ANYTHING BELOW THIS LINE
# --------------------------------------
# These variables specify the RDS files created by ./code/format_gtex.R.
tissue_dir = ./output/tissue_data
tissue_dat = $(tissue_dir)/adipose_subcutaneous.RDS \
	     $(tissue_dir)/adipose_visceral_omentum.RDS \
	     $(tissue_dir)/adrenal_gland.RDS \
	     $(tissue_dir)/artery_aorta.RDS \
	     $(tissue_dir)/artery_coronary.RDS \
	     $(tissue_dir)/artery_tibial.RDS \
	     $(tissue_dir)/bladder.RDS \
	     $(tissue_dir)/brain_amygdala.RDS \
	     $(tissue_dir)/brain_anterior_cingulate_cortex_ba24.RDS \
	     $(tissue_dir)/brain_caudate_basal_ganglia.RDS \
	     $(tissue_dir)/brain_cerebellar_hemisphere.RDS \
	     $(tissue_dir)/brain_cerebellum.RDS \
	     $(tissue_dir)/brain_cortex.RDS \
	     $(tissue_dir)/brain_frontal_cortex_ba9.RDS \
	     $(tissue_dir)/brain_hippocampus.RDS \
	     $(tissue_dir)/brain_hypothalamus.RDS \
	     $(tissue_dir)/brain_nucleus_accumbens_basal_ganglia.RDS \
	     $(tissue_dir)/brain_putamen_basal_ganglia.RDS \
	     $(tissue_dir)/brain_spinal_cord_cervical_c1.RDS \
	     $(tissue_dir)/brain_substantia_nigra.RDS \
	     $(tissue_dir)/breast_mammary_tissue.RDS \
	     $(tissue_dir)/cells_ebvtransformed_lymphocytes.RDS \
	     $(tissue_dir)/cells_transformed_fibroblasts.RDS \
	     $(tissue_dir)/cervix_ectocervix.RDS \
	     $(tissue_dir)/cervix_endocervix.RDS \
	     $(tissue_dir)/colon_sigmoid.RDS \
	     $(tissue_dir)/colon_transverse.RDS \
	     $(tissue_dir)/esophagus_gastroesophageal_junction.RDS \
	     $(tissue_dir)/esophagus_mucosa.RDS \
	     $(tissue_dir)/esophagus_muscularis.RDS \
	     $(tissue_dir)/fallopian_tube.RDS \
	     $(tissue_dir)/heart_atrial_appendage.RDS \
	     $(tissue_dir)/heart_left_ventricle.RDS \
	     $(tissue_dir)/kidney_cortex.RDS \
	     $(tissue_dir)/liver.RDS \
	     $(tissue_dir)/lung.RDS \
	     $(tissue_dir)/minor_salivary_gland.RDS \
	     $(tissue_dir)/muscle_skeletal.RDS \
	     $(tissue_dir)/nerve_tibial.RDS \
	     $(tissue_dir)/ovary.RDS \
	     $(tissue_dir)/pancreas.RDS \
	     $(tissue_dir)/pituitary.RDS \
	     $(tissue_dir)/prostate.RDS \
	     $(tissue_dir)/skin_not_sun_exposed_suprapubic.RDS \
	     $(tissue_dir)/skin_sun_exposed_lower_leg.RDS \
	     $(tissue_dir)/small_intestine_terminal_ileum.RDS \
	     $(tissue_dir)/spleen.RDS \
	     $(tissue_dir)/stomach.RDS \
	     $(tissue_dir)/testis.RDS \
	     $(tissue_dir)/thyroid.RDS \
	     $(tissue_dir)/uterus.RDS \
	     $(tissue_dir)/vagina.RDS \
	     $(tissue_dir)/whole_blood.RDS

# Output of comparing theoretical and seqgendiff datasets from ./code/data_features.R
feature_dir = ./output/figures/powsimr_vs_seqgendiff
feature_plots = $(feature_dir)/count_dist.pdf \
                $(feature_dir)/mean_gene_depth.pdf \
		$(feature_dir)/pc1_pc4.pdf \
                $(feature_dir)/pc_plot.pdf \
                $(feature_dir)/scree.pdf \
                $(feature_dir)/seqgendiff_truevsfits.pdf \
                $(feature_dir)/voom_plots.pdf

# Output of factor analysis simulation plots from ./code/fasim_plots.R
fasimplots_dir = ./output/figures/fasim_plots
fasimplots = $(fasimplots_dir)/minloadmse_0_0.pdf \
             $(fasimplots_dir)/minloadmse_5_0.pdf \
             $(fasimplots_dir)/angle_0_0.pdf \
             $(fasimplots_dir)/angle_5_0.pdf \
             $(fasimplots_dir)/minmse_0_0.pdf \
             $(fasimplots_dir)/minmse_5_0.pdf

# Output of the differential expression analysis simulations from ./code/diff_exp_plots.R
diff_exp_dir = ./output/figures/diff_exp
diff_exp_plots = $(diff_exp_dir)/fdp.pdf \
		 $(diff_exp_dir)/mpve.pdf \
		 $(diff_exp_dir)/mse.pdf \
		 $(diff_exp_dir)/power.pdf

# Simulation output of ./code/run_sims.R
fasims_out = ./output/fa_sims/fa_results.csv

# These variables specify the raw GTEx data, which we will clean and format into the tissue_dat
gtex_dir = ./data/gtex
gtex_dat = $(gtex_dir)/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct \
           $(gtex_dir)/GTEx_v7_Annotations_SubjectPhenotypesDS.txt \
           $(gtex_dir)/GTEx_v7_Annotations_SampleAttributesDS.txt

# The raw PBMC data ----------------------------------------------------
pbmc_dat = ./data/pbmc/barcodes.tsv \
           ./data/pbmc/genes.tsv \
           ./data/pbmc/matrix.mtx

# single cell plots ----------------------------------------------------
sc_plots = ./output/figures/sc_plots/sc_angle.pdf \
           ./output/figures/sc_plots/sc_mse.pdf \
           ./output/figures/sc_plots/sc_loadmse.pdf

all: NBplots FAsims powsimr corest diff_exp sc_fa

# Extract tissue data --------------------------------------------------
$(tissue_dat) : ./code/format_gtex.R $(gtex_dat)
	mkdir -p ./output/tissue_data
	mkdir -p ./output/rout
	$(rexec) $< output/rout/$(basename $(notdir $<)).Rout

# Mixtures of binomial and negative binomial plot ----------------------
./output/figures/mix_dists.pdf : ./code/flexible_nb.R
	mkdir -p ./output/figures
	mkdir -p ./output/rout
	$(rexec) $< output/rout/$(basename $(notdir $<)).Rout

.PHONY : NBplots
NBplots : ./output/figures/mix_dists.pdf

# FA Simulations -------------------------------------------------------

# factor analysis simulations -------------------------------------------
$(fasims_out) : ./code/run_sims.R ./code/fa_methods.R ./code/signal_funs.R $(tissue_dat)
	mkdir -p ./output/fa_sims
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< output/rout/$(basename $(notdir $<)).Rout

$(fasimplots) : ./code/fasim_plots.R $(fasims_out)
	mkdir -p ./output/figures
	mkdir -p ./output/figures/fasim_plots
	mkdir -p ./output/rout
	$(rexec) $< output/rout/$(basename $(notdir $<)).Rout

$(fasimplots_dir)/mpve_fasims.pdf : ./code/mpve_plots.R $(fasims_out)
	mkdir -p ./output/figures
	mkdir -p ./output/figures/fasim_plots
	mkdir -p ./output/rout
	$(rexec) $< output/rout/$(basename $(notdir $<)).Rout

.PHONY : FAsims
FAsims : $(fasimplots) $(fasimplots_dir)/mpve_fasims.pdf


# powsimR simulations ---------------------------------------------------
./output/compare_powsimR/powsim_params.RDS : ./code/get_powsimr_params.R $(tissue_dat)
	mkdir -p ./output/compare_powsimR
	mkdir -p ./output/rout
	$(rexec) $< output/rout/$(basename $(notdir $<)).Rout

$(feature_plots) : ./code/data_features.R ./output/compare_powsimR/powsim_params.RDS
	mkdir -p ./output/figures
	mkdir -p ./output/figures/powsimr_vs_seqgendiff
	mkdir -p ./output/rout
	$(rexec) $< output/rout/$(basename $(notdir $<)).Rout

.Phony : powsimr
powsimr : $(feature_plots)


# Correlation estimator simulations -------------------------------------
./output/est_cor/corsimout.RDS : ./code/est_cor.R
	mkdir -p ./output/est_cor
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< output/rout/$(basename $(notdir $<)).Rout

./output/figures/cor_est.pdf : ./code/corsim_plots.R ./output/est_cor/corsimout.RDS
	mkdir -p ./output/figures
	mkdir -p ./output/rout
	$(rexec) $< output/rout/$(basename $(notdir $<)).Rout

.Phony : corest
corest : ./output/figures/cor_est.pdf

# Differential expression simulations ----------------------------------
./output/diff_exp_out/powsimr_sims.RDS : ./code/run_powsimr_sims.R ./output/compare_powsimR/powsim_params.RDS ./code/de_methods.R
	mkdir -p ./output/diff_exp_out
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< output/rout/$(basename $(notdir $<)).Rout

./output/diff_exp_out/seqgendiff_sims.RDS : ./code/run_seqgendiff_sims.R $(tissue_dat) ./code/de_methods.R
	mkdir -p ./output/diff_exp_out
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< output/rout/$(basename $(notdir $<)).Rout

./output/diff_exp_out/simseq_sims.RDS : ./code/run_simseq_sims.R $(tissue_dat) ./code/de_methods.R
	mkdir -p ./output/diff_exp_out
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< output/rout/$(basename $(notdir $<)).Rout

./output/simseq_vs_seqgendiff/ash_ests.csv : ./code/prop_nonnull.R ./output/tissue_data/muscle_skeletal.RDS
	mkdir -p ./output/simseq_vs_seqgendiff
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< output/rout/$(basename $(notdir $<)).Rout

./output/diff_exp_out/seqgendiff_small_mpve_sims.RDS : ./code/run_seqgendiff_small_mpve_sims.R $(tissue_dat) ./code/de_methods.R ./output/simseq_vs_seqgendiff/ash_ests.csv
	mkdir -p ./output/diff_exp_out
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< output/rout/$(basename $(notdir $<)).Rout

$(diff_exp_plots) : ./code/diff_exp_plots.R ./output/diff_exp_out/powsimr_sims.RDS ./output/diff_exp_out/seqgendiff_sims.RDS
	mkdir -p ./output/figures
	mkdir -p ./output/figures/diff_exp
	mkdir -p ./output/rout
	$(rexec) $< output/rout/$(basename $(notdir $<)).Rout

./output/figures/simseq_plots/simseq_gtex.pdf ./output/figures/simseq_plots/simseq_time.pdf: ./code/plot_simseq_seqgendiff_gtex.R ./output/diff_exp_out/seqgendiff_small_mpve_sims.RDS ./output/diff_exp_out/simseq_sims.RDS
	mkdir -p ./output/diff_exp_out
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< output/rout/$(basename $(notdir $<)).Rout

.Phony : diff_exp
diff_exp : $(diff_exp_plots) ./output/figures/simseq_plots/simseq_gtex.pdf ./output/figures/simseq_plots/simseq_time.pdf

# Single cell factor analysis simulations ---------------------------------------

./output/sc/pbmc_cleaned.RDS : ./code/extract_pbmc.R $(pbmc_dat)
	mkdir -p ./output/sc
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< output/rout/$(basename $(notdir $<)).Rout

./output/sc/sc_fa_sims.csv : ./code/run_sims_sc.R ./output/sc/pbmc_cleaned.RDS ./code/fa_methods.R ./code/signal_funs.R
	mkdir -p ./output/sc
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< output/rout/$(basename $(notdir $<)).Rout

$(sc_plots) : ./code/sc_plots.R ./output/sc/sc_fa_sims.csv
	mkdir -p ./output/figures
	mkdir -p ./output/figures/sc_plots
	mkdir -p ./output/rout
	$(rexec) $< output/rout/$(basename $(notdir $<)).Rout

.Phony : sc_fa
sc_fa : $(sc_plots)

# remove tissue data to start over ---------------------------------------------
clean:
	rm -f $(tissue_dat)
