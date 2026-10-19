## ResolveOME single-cell DNA workflow -- Sprocket port.
##
## This is a Sprocket-native adaptation of ../main.wdl (the Cromwell-oriented
## general-purpose WDL workflow), built for St. Jude's own next-generation
## WDL engine (https://github.com/stjude-rust-labs/sprocket). Same pipeline,
## same bug fixes, same stages -- only the container-handling mechanics
## differ. See ../README.md ("this sprocket/ folder") for the full
## explanation of every difference and a step-by-step guide to actually
## running this (not done as part of writing this port).
##
## =============================================================================
## THE ONE ARCHITECTURAL DIFFERENCE THAT MATTERS: ONE CONTAINER PER TASK
## =============================================================================
## Cromwell's "Local"/"LSF" backends (../main.wdl) can run a task's command
## block directly on the host (conda-activating a named environment), and
## only specific tasks manually `apptainer exec` a second, separate
## container mid-command for tools that need one (DeepVariant/VEP/Manta).
## Sprocket's execution model is different and simpler: EVERY task runs
## inside exactly one container, chosen by that task's own `requirements {
## container: ... }` value -- there is no "run bare on the host" option, and
## no way to invoke a second container from inside a task's command block.
##
## This has two concrete consequences, both implemented in ./tasks/*.wdl:
##   1. Every tool-stack this workflow used to reach via `conda activate`
##      (resolveome-align, freec, mutsig) now needs an actual container
##      image instead -- see ./containers/*.Dockerfile. These are NEW
##      images this port needs that did not exist before; building and
##      pushing them (e.g. to ghcr.io) is a required one-time setup step
##      this port does NOT perform for you (see ../README.md).
##   2. The handful of tasks that used to run TWO toolsets in one command
##      block (manually `apptainer exec` one container, then `conda
##      activate` a second environment just for a trailing `bcftools
##      index`) had to be restructured: `MantaCall` now produces its own
##      output un-conditionally, and a separate, shared `IndexVcf` task
##      (./tasks/common.wdl) is called below ONLY if Manta didn't already
##      produce its own index -- see the `if (!defined(...))` block in
##      Stage 3. DeepVariant/VEP's chunk-level `bcftools index` steps were
##      dropped outright, because tracing ../main.wdl's own data flow shows
##      that index was never actually consumed downstream (MergeVcfChunks
##      re-sorts and re-indexes every chunk itself regardless).
##
## =============================================================================
## HOW INPUTS ARE SUPPLIED -- ONE inputs.json FOR THE WHOLE WORKFLOW
## =============================================================================
## Same convention as ../main.wdl: supply exactly one inputs JSON (or YAML)
## file, with every key prefixed "ResolveOmeDnaWorkflow.<field>", e.g.
## `"ResolveOmeDnaWorkflow.sample_ids": [...]`. See ./example_inputs.json
## for a complete, annotated starting point -- note that it differs from
## ../wdl/example_inputs.json in exactly the ways described above (no more
## `*_sif` File inputs or `conda_env`/`freec_conda_env`/`r_conda_env` String
## inputs; new `*_container` String inputs instead).
## =============================================================================
version 1.2

import "tasks/align_qc.wdl" as align_qc_t
import "tasks/cell_qc_gate.wdl" as qc_gate_t
import "tasks/variant_calling.wdl" as variants_t
import "tasks/cnv_sv.wdl" as cnvsv_t
import "tasks/qc_aggregate.wdl" as qc_t
import "tasks/phasing.wdl" as phasing_t
import "tasks/loh_ploidy.wdl" as loh_t
import "tasks/lineage.wdl" as lineage_t
import "tasks/extra_qc_dashboard.wdl" as extra_t
import "tasks/common.wdl" as common_t

workflow ResolveOmeDnaWorkflow {
  # ===========================================================================
  # WORKFLOW INPUTS -- every field below is one key in your single inputs.json
  # (prefixed "ResolveOmeDnaWorkflow.<field>"); fields with "=" have a default
  # and may be omitted from inputs.json. Grouped by the pipeline stage that
  # consumes them; see that stage's `tasks/*.wdl` file for how each is used.
  # ===========================================================================
  input {
    # ---- cohort: identifies the donor/cells and points at their FASTQs ----
    String donor_id              # free-text cohort/donor label, used in lineage tree output naming
    Array[String] sample_ids     # one cell/sample ID per cell, e.g. ["cell1","cell2",...]
    Array[File] r1_fastqs       # same order as sample_ids -- r1_fastqs[i]/r2_fastqs[i] belong to sample_ids[i]
    Array[File] r2_fastqs
    String gender = "male"      # "male" or "female" -- controls chrX phasing (see phasing.wdl)

    # ---- reference: GRCh38 FASTA + its indices, used by nearly every stage ----
    # (produced by ../setup/setup.sh's "Reference genome" section)
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File ref_bwamem2_0123
    File ref_bwamem2_bwt2bit64
    File ref_bwamem2_ann
    File ref_bwamem2_amb
    File ref_bwamem2_pac
    File metrics_intervals_bed   # BED of regions Picard/samtools depth QC should summarize coverage over

    # ---- DeepVariant: small-variant calling, consumed by tasks/variant_calling.wdl ----
    String deepvariant_container = "docker://google/deepvariant:1.6.1"   # free, stock, public
    File? deepvariant_model_dir_tar              # optional: custom checkpoint (e.g. BioSkryb's private
                                                   # single-cell sc_af model); leave unset for the free,
                                                   # general-purpose default (DeepVariant's own bundled
                                                   # --model_type=WGS model, no extra download needed)
    String deepvariant_checkpoint_filename = "sc_af.model.ckpt"   # only used if deepvariant_model_dir_tar is set
    Boolean dv_use_allele_frequency = false     # true only if deepvariant_population_vcfs_tar is BioSkryb's private AF bundle
    File? deepvariant_population_vcfs_tar        # optional: tar of per-chromosome population VCFs (BioSkryb-private)
    String chroms = "all"             # which chromosomes PlanChunks should split into DeepVariant regions
    Int chunk_size_bp = 20000000      # region chunk size in bp, for parallelizing DeepVariant across each chromosome

    # ---- VEP / SnpEff: variant annotation, consumed by tasks/variant_calling.wdl ----
    String vep_container = "docker://ensemblorg/ensembl-vep:release_111.0"   # free, stock, public
    File vep_cache_dir_tar         # tarball of the Ensembl VEP species cache, downloaded by ../setup/setup.sh
    Int vep_cache_version = 111
    String snpeff_db = "GRCh38.111"   # snpEff database name -- MUST match ../setup/setup.sh's
                                        # ENSEMBL_RELEASE (default 111), since setup.sh downloads
                                        # "GRCh38.${ENSEMBL_RELEASE}"

    # ---- Per-cell QC gate (BaseJumper-derived composite score): tasks/cell_qc_gate.wdl ----
    # Runs right after Stage 1, before any downstream analysis. A cell that
    # scores <=2/5 (FAIL) has every downstream stage below (DeepVariant/VEP/
    # SnpEff, Manta, phasing, LOH/ploidy, lineage) skipped entirely -- see
    # extra_methods/cell_qc_gate.py's own docstring for the exact formula,
    # ported from BaseJumper's basej-dnaqc. Cutoffs default to BaseJumper's
    # own published values; override only if you have a specific reason to.
    File cell_qc_gate_script
    Float qc_gate_cutoff_preseq = 3500000000
    Float qc_gate_cutoff_chimeras = 0.20
    Float qc_gate_cutoff_chrm = 0.20
    Float qc_gate_cutoff_mapd = 0.25
    Float qc_gate_cutoff_skew = 0.25

    # ---- CNVkit / Manta: copy-number + structural-variant calling, tasks/cnv_sv.wdl ----
    File? cnv_reference_cnn       # optional pre-built CNVkit flat/pooled reference .cnn; omit to build one from this run's own BAMs
    String manta_container = "docker://quay.io/biocontainers/manta:1.6.0--py27h9948957_6"   # free, stock, public;
                                                                                               # same exact build this project's own manta_env conda env pins

    # ---- Phasing (Beagle): autosome-only haplotype phasing, tasks/phasing.wdl ----
    File beagle_jar                               # Beagle 5.5 jar, downloaded by ../setup/setup.sh
    Map[String, File] phasing_ref_panel_by_chrom     # chromosome "1".."22" -> matching 1000G GRCh38 panel VCF
    Map[String, File] phasing_ref_panel_tbi_by_chrom # same keys -> that VCF's .tbi index
    Map[String, File] genetic_map_by_chrom           # same keys -> that chromosome's Beagle/PLINK genetic map

    # ---- LOH / ploidy (Control-FREEC): tasks/loh_ploidy.wdl -- run TWICE per cell ----
    # (once at a fixed ploidy=2 for a simple baseline, once letting FREEC pick
    # the best-fitting ploidy from freec_candidate_ploidies -- see that task's
    # docstring for the model-selection caveat this produces)
    File chrom_fasta_dir_tar       # tarball of per-chromosome FASTAs (chr1-22+X ONLY), built by ../setup/setup.sh
    File chrom_len_file            # samtools .fai restricted to chr1-22+X ONLY -- see loh_ploidy.wdl
    File snp_positions_vcf         # site-only 1000G VCF (can reuse the phasing panel, stripped of genotypes)
    File snp_positions_vcf_tbi
    String freec_fixed_ploidy = "2"
    String freec_candidate_ploidies = "2,3,4,5,6,7,8"

    # ---- Lineage (Sequoia, standalone GPLv3): tasks/lineage.wdl ----
    File sequoia_build_phylogeny_script   # Sequoia's build_phylogeny.R, from ../setup/setup.sh's git clone
    File verify_lineage_r_script          # this project's own bootstrap/resample verification R script
    File mpboot_binary                    # MPBoot tree-building binary, downloaded by ../setup/setup.sh
    Int mpboot_bootstrap = 1000            # number of MPBoot bootstrap replicates for tree support values

    # ---- extra_methods scripts: Python/R calculation tasks, tasks/extra_qc_dashboard.wdl ----
    # (these are the exact scripts under ../extra_methods/, passed in as WDL
    # `File` inputs rather than baked into a container -- see each script's
    # own docstring for what it computes)
    File make_region_chunks_script       # used by PlanChunks (tasks/variant_calling.wdl), not extra_qc_dashboard.wdl
    File mapd_qc_script
    File fdr_estimate_script
    File artifact_panel_script
    File phasing_summary_script
    File loh_summary_script
    File ploidy_summary_script
    File loh_by_region_script
    File qc_summary_script
    File mutational_signatures_r_script
    File lineage_track_script
    File verify_lineage_script
    File build_dashboard_script

    # ---- container image references, built by ../setup/ (see ../README.md) ----
    # These replace the Cromwell port's `conda_env`/`freec_conda_env`/
    # `r_conda_env` String inputs -- every task that used to `conda activate`
    # one of those three environments now instead runs inside the matching
    # container image named here.
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
    String freec_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-freec:latest"
    String lineage_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-lineage:latest"
    String mutsig_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-mutsig:latest"
  }

  # ===========================================================================
  # STAGE 1 -- Per-cell: align + QC  (tasks/align_qc.wdl, task AlignTrimQC)
  # ===========================================================================
  scatter (i in range(length(sample_ids))) {
    call align_qc_t.AlignTrimQC as align {
      input:
        sample_id = sample_ids[i],
        r1 = r1_fastqs[i],
        r2 = r2_fastqs[i],
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        ref_fasta_dict = ref_fasta_dict,
        ref_bwamem2_0123 = ref_bwamem2_0123,
        ref_bwamem2_bwt2bit64 = ref_bwamem2_bwt2bit64,
        ref_bwamem2_ann = ref_bwamem2_ann,
        ref_bwamem2_amb = ref_bwamem2_amb,
        ref_bwamem2_pac = ref_bwamem2_pac,
        metrics_intervals_bed = metrics_intervals_bed,
        resolveome_align_container = resolveome_align_container,
    }
  }

  # ===========================================================================
  # STAGE 2 -- Small variants: DeepVariant (chunked) + VEP + SnpEff annotation
  # ===========================================================================
  call variants_t.PlanChunks as plan_chunks {
    input:
      ref_fasta_fai = ref_fasta_fai,
      chroms = chroms,
      chunk_size_bp = chunk_size_bp,
      make_region_chunks_script = make_region_chunks_script,
      resolveome_align_container = resolveome_align_container,
  }

  scatter (i in range(length(sample_ids))) {
    # ===========================================================================
    # STAGE 3 (CNVkit half only) -- run BEFORE the QC gate below, deliberately.
    # CnvkitCall only needs Stage 1's dedup_bam, so it's moved up from its
    # original position (right alongside MantaCall, further down) to double as
    # the gate's own MAPD/SKEW source -- this mirrors BaseJumper's own
    # architecture, where `basej-dnaqc` runs its own dedicated CNV call purely
    # for QC scoring, separate from whatever deeper CNV/SV pipeline runs
    # afterwards. This also means CNVkit runs exactly ONCE per cell (not
    # duplicated for gating vs. the real analysis) -- the gate's own
    # `CellQcGate` task (next, below) reads THIS call's .cnr/.call_cns.
    # ===========================================================================
    call cnvsv_t.CnvkitCall as cnvkit {
      input:
        sample_id = sample_ids[i],
        bam = align.dedup_bam[i],
        bam_bai = align.dedup_bam_bai[i],
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        cnv_reference_cnn = cnv_reference_cnn,
        resolveome_align_container = resolveome_align_container,
    }

    # ===========================================================================
    # QC GATE -- BaseJumper-derived composite score (tasks/cell_qc_gate.wdl).
    # Runs unconditionally (it IS the gate), using Stage 1's Picard/idxstats/
    # preseq outputs plus the CnvkitCall just above. Every stage below this
    # point in the SAME scatter iteration -- and Stage 4/5 further down --
    # is wrapped in `if (gate.qc_pass)` and is SKIPPED for a FAILing cell. See
    # ../extra_methods/cell_qc_gate.py's docstring for the exact formula and
    # gate.qc_gate_explanation for a human-readable per-cell breakdown of
    # which metric(s) failed and why.
    # ===========================================================================
    call qc_gate_t.CellQcGate as gate {
      input:
        sample_id = sample_ids[i],
        bam = align.dedup_bam[i],
        bam_bai = align.dedup_bam_bai[i],
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        cnv_reference_cnn = cnv_reference_cnn,
        idxstats = align.idxstats[i],
        alignment_summary_metrics = align.alignment_summary_metrics[i],
        preseq_result = align.preseq_result[i],
        cell_qc_gate_script = cell_qc_gate_script,
        cutoff_preseq = qc_gate_cutoff_preseq,
        cutoff_chimeras = qc_gate_cutoff_chimeras,
        cutoff_chrm = qc_gate_cutoff_chrm,
        cutoff_mapd = qc_gate_cutoff_mapd,
        cutoff_skew = qc_gate_cutoff_skew,
        resolveome_align_container = resolveome_align_container,
    }

    if (gate.qc_pass) {
      # STAGE 2 -- Small variants: DeepVariant (chunked) + VEP + SnpEff
      # annotation. plan_chunks.chunks (computed ONCE above, outside this
      # scatter) is reused by every PASSING cell.
      scatter (chunk in plan_chunks.chunks) {
        # chunk = [fai_name, canonical_chrom, start, end, region, chunk_id]
        call variants_t.DeepVariantCall as dv_chunk {
          input:
            sample_id = sample_ids[i],
            bam = align.dedup_bam[i],
            bam_bai = align.dedup_bam_bai[i],
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            region = chunk[4],
            chunk_id = chunk[5],
            canonical_chrom = chunk[1],
            deepvariant_model_dir_tar = deepvariant_model_dir_tar,
            deepvariant_checkpoint_filename = deepvariant_checkpoint_filename,
            population_vcfs_tar = deepvariant_population_vcfs_tar,
            use_allele_frequency = dv_use_allele_frequency,
            deepvariant_container = deepvariant_container,
        }
        call variants_t.VepAnnotateCall as vep_chunk {
          input:
            sample_id = sample_ids[i],
            chunk_vcf = dv_chunk.chunk_vcf,
            chunk_id = chunk[5],
            vep_cache_dir_tar = vep_cache_dir_tar,
            ref_fasta = ref_fasta,
            vep_cache_version = vep_cache_version,
            vep_container = vep_container,
        }
      }

      call variants_t.MergeVcfChunks as merge_raw {
        input:
          sample_id = sample_ids[i],
          chunk_vcfs = dv_chunk.chunk_vcf,
          out_suffix = ".deepvariant.raw.vcf.gz",
          resolveome_align_container = resolveome_align_container,
      }
      call variants_t.MergeVcfChunks as merge_vep {
        input:
          sample_id = sample_ids[i],
          chunk_vcfs = vep_chunk.vep_chunk_vcf,
          out_suffix = ".deepvariant.raw.vep.vcf.gz",
          require_vep_ann = true,
          resolveome_align_container = resolveome_align_container,
      }
      call variants_t.SnpEffAnnotate as snpeff {
        input:
          sample_id = sample_ids[i],
          merged_vep_vcf = merge_vep.merged_vcf,
          snpeff_db = snpeff_db,
          resolveome_align_container = resolveome_align_container,
      }

      # ===========================================================================
      # STAGE 3 (Manta half) -- SV calling, per cell, reading the SAME dedup_bam
      # as CnvkitCall above (CnvkitCall itself already ran, unconditionally,
      # before the gate -- see the comment at the top of this scatter body).
      # ---------------------------------------------------------------------------
      # MantaCall's own `.tbi` output is `File?` -- if it's missing (Manta
      # didn't produce one, a real possibility the Cromwell version also
      # guarded against), fall back to the shared IndexVcf helper task. This
      # `if` block is how Sprocket's one-container-per-task model reproduces
      # the Cromwell version's "apptainer exec Manta, then conda-activate
      # bcftools as a fallback" pattern without ever running two toolsets in
      # one task.
      # ===========================================================================
      call cnvsv_t.MantaCall as manta {
        input:
          sample_id = sample_ids[i],
          bam = align.dedup_bam[i],
          bam_bai = align.dedup_bam_bai[i],
          ref_fasta = ref_fasta,
          ref_fasta_fai = ref_fasta_fai,
          manta_container = manta_container,
      }
      if (!defined(manta.manta_vcf_tbi)) {
        call common_t.IndexVcf as manta_reindex {
          input:
            vcf = manta.manta_vcf,
            out_name = "~{sample_ids[i]}.manta.vcf.gz",
            resolveome_align_container = resolveome_align_container,
        }
      }

      # Pass-throughs of already-unconditional Stage 1 values, gated by the
      # SAME `if (gate.qc_pass)` condition as everything else in this block --
      # cohort-wide tasks further down (Stage 7) `select_all()` these
      # alongside snpeff/manta's own (already-optional-by-construction)
      # outputs to build a consistent "passing cells only" view, in the same
      # relative order, without a separate manual filter step.
      String sample_id_passed = sample_ids[i]
      File wgs_metrics_passed = align.wgs_metrics[i]
    }
  }

  # ===========================================================================
  # STAGE 4 -- Phasing: Beagle 5.5, autosomes ONLY
  # ===========================================================================
  Array[String] autosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
                              "13","14","15","16","17","18","19","20","21","22"]

  scatter (i in range(length(sample_ids))) {
    # Skip phasing entirely for a FAILing cell (same gate as Stage 2/3 above) --
    # merge_raw.merged_vcf[i]/merged_vcf_csi[i] are themselves optional (only
    # produced when gate.qc_pass was true for this same cell).
    if (gate.qc_pass[i]) {
      scatter (chrom in autosomes) {
        call phasing_t.PhaseChromosome as phase_chrom {
          input:
            sample_id = sample_ids[i],
            raw_vcf = select_first([merge_raw.merged_vcf[i]]),
            raw_vcf_csi = select_first([merge_raw.merged_vcf_csi[i]]),
            chrom = chrom,
            ref_panel_vcf = phasing_ref_panel_by_chrom[chrom],
            ref_panel_vcf_tbi = phasing_ref_panel_tbi_by_chrom[chrom],
            genetic_map = genetic_map_by_chrom[chrom],
            beagle_jar = beagle_jar,
            resolveome_align_container = resolveome_align_container,
        }
      }
      call phasing_t.ConcatPhasedChromosomes as concat_phased {
        input:
          sample_id = sample_ids[i],
          phased_chrom_vcfs = phase_chrom.phased_chrom_vcf,
          resolveome_align_container = resolveome_align_container,
      }
    }
  }

  # ===========================================================================
  # STAGE 5 -- LOH + ploidy: Control-FREEC, per cell, run TWICE
  # ===========================================================================
  scatter (i in range(length(sample_ids))) {
    # Stage 5 reads align.dedup_bam directly (not a gated Stage 2/3 output),
    # so it needs its own explicit qc_pass gate.
    if (gate.qc_pass[i]) {
      call loh_t.ControlFreec as freec_fixed {
        input:
          sample_id = sample_ids[i],
          bam = align.dedup_bam[i],
          bam_bai = align.dedup_bam_bai[i],
          ref_fasta = ref_fasta,
          chrom_fasta_dir_tar = chrom_fasta_dir_tar,
          chr_len_file = chrom_len_file,
          snp_positions_vcf = snp_positions_vcf,
          snp_positions_vcf_tbi = snp_positions_vcf_tbi,
          ploidy = freec_fixed_ploidy,
          freec_container = freec_container,
      }
      call loh_t.ControlFreec as freec_ploidy {
        input:
          sample_id = sample_ids[i],
          bam = align.dedup_bam[i],
          bam_bai = align.dedup_bam_bai[i],
          ref_fasta = ref_fasta,
          chrom_fasta_dir_tar = chrom_fasta_dir_tar,
          chr_len_file = chrom_len_file,
          snp_positions_vcf = snp_positions_vcf,
          snp_positions_vcf_tbi = snp_positions_vcf_tbi,
          ploidy = freec_candidate_ploidies,
          freec_container = freec_container,
      }
    }
  }

  # ===========================================================================
  # STAGE 6 -- Cohort-wide QC rollup: MultiQC
  # ===========================================================================
  call qc_t.MultiQC as multiqc {
    input:
      input_files = flatten([flatten(align.fastp_out), flatten(align.fastqc_raw_out),
                              flatten(align.fastqc_trimmed_out), flatten(align.picard_out),
                              select_all(snpeff.snpeff_html), select_all(snpeff.snpeff_genes_csv)]),
      resolveome_align_container = resolveome_align_container,
  }

  # ===========================================================================
  # STAGE 7 -- Cohort-wide extra_methods calculations + the final dashboard
  # ---------------------------------------------------------------------------
  # `passed_sample_ids` is the gate-filtered sample_ids list (see Stage 2/3's
  # `if (gate.qc_pass)` block above for the `sample_id_passed`/
  # `wgs_metrics_passed` pass-throughs it's built from) -- every select_all()'d
  # array below stays index-aligned with it, since all of them were produced
  # by scatters gated on that SAME per-cell condition.
  # ===========================================================================
  Array[String] passed_sample_ids = select_all(sample_id_passed)

  call extra_t.MapdQc as mapd {
    input: cnr_files = cnvkit.cnr, sample_ids = sample_ids,
           mapd_qc_script = mapd_qc_script, resolveome_align_container = resolveome_align_container,
  }
  call extra_t.FdrEstimate as fdr {
    input: final_annotated_vcfs = select_all(snpeff.final_annotated_vcf),
           wgs_metrics_files = select_all(wgs_metrics_passed),
           sample_ids = passed_sample_ids, fdr_estimate_script = fdr_estimate_script,
           resolveome_align_container = resolveome_align_container,
  }
  call extra_t.ArtifactPanel as artifact {
    input: final_annotated_vcfs = select_all(snpeff.final_annotated_vcf), sample_ids = passed_sample_ids,
           artifact_panel_script = artifact_panel_script, resolveome_align_container = resolveome_align_container,
  }
  call extra_t.PhasingSummary as phasing_summary {
    input: phased_vcfs = select_all(concat_phased.phased_vcf), sample_ids = passed_sample_ids,
           phasing_summary_script = phasing_summary_script, resolveome_align_container = resolveome_align_container,
  }
  call extra_t.LohAndPloidySummary as loh_summary {
    input: cnvs_files = select_all(freec_ploidy.cnvs), info_files = select_all(freec_ploidy.info),
           sample_ids = passed_sample_ids,
           loh_summary_script = loh_summary_script, ploidy_summary_script = ploidy_summary_script,
           loh_by_region_script = loh_by_region_script, chrom_len_file = chrom_len_file,
           resolveome_align_container = resolveome_align_container,
  }
  call extra_t.QcSummary as qc_summary {
    input: multiqc_general_stats = multiqc.data[0], sample_ids = sample_ids,
           qc_summary_script = qc_summary_script, resolveome_align_container = resolveome_align_container,
  }
  call extra_t.MutationalSignatures as mutsig {
    input: artifact_panel_sites_tsv = artifact.sites_tsv,
           mutational_signatures_r_script = mutational_signatures_r_script, mutsig_container = mutsig_container,
  }
  call lineage_t.BuildLineageTree as lineage {
    input: final_annotated_vcfs = select_all(snpeff.final_annotated_vcf), sample_ids = passed_sample_ids,
           lineage_track_script = lineage_track_script,
           sequoia_build_phylogeny_script = sequoia_build_phylogeny_script,
           mpboot_binary = mpboot_binary, donor_id = donor_id, gender = gender,
           mpboot_bootstrap = mpboot_bootstrap, lineage_container = lineage_container,
  }
  call lineage_t.VerifyLineage as lineage_verify {
    input: final_annotated_vcfs = select_all(snpeff.final_annotated_vcf), sample_ids = passed_sample_ids,
           lineage_track_script = lineage_track_script, verify_lineage_script = verify_lineage_script,
           verify_lineage_r_script = verify_lineage_r_script,
           sequoia_build_phylogeny_script = sequoia_build_phylogeny_script,
           mpboot_binary = mpboot_binary, primary_lineage_outputs = lineage.lineage_outputs,
           donor_id = donor_id, gender = gender, mpboot_bootstrap = mpboot_bootstrap,
           lineage_container = lineage_container,
  }
  call extra_t.BuildDashboard as dashboard {
    input:
      mapd_qc_tsv = mapd.mapd_qc_tsv,
      fdr_estimate_tsv = fdr.fdr_estimate_tsv,
      artifact_panel_summary_tsv = artifact.summary_tsv,
      phasing_summary_tsv = phasing_summary.phasing_summary_tsv,
      loh_summary_tsv = loh_summary.loh_summary_tsv,
      ploidy_summary_tsv = loh_summary.ploidy_summary_tsv,
      loh_by_region_tsv = loh_summary.loh_by_region_tsv,
      loh_by_region_png = loh_summary.loh_by_region_png,
      qc_summary_tsv = qc_summary.qc_summary_tsv,
      qc_overview_png = qc_summary.qc_overview_png,
      signature_contributions_tsv = mutsig.signature_contributions_tsv,
      cosine_similarity_tsv = mutsig.cosine_similarity_tsv,
      mutation_spectrum_png = mutsig.mutation_spectrum_png,
      signature_barplot_png = mutsig.signature_barplot_png,
      lineage_outputs = lineage.lineage_outputs,
      bootstrap_support_report = lineage_verify.bootstrap_support_report,
      verification_summary_json = lineage_verify.verification_summary_json,
      build_dashboard_script = build_dashboard_script,
      donor_id = donor_id,
      resolveome_align_container = resolveome_align_container,
  }

  output {
    # =========================================================================
    # WORKFLOW OUTPUTS -- this is what ends up in your outputs.json (prefixed
    # "ResolveOmeDnaWorkflow.<field>"). Feed this straight into
    # ../wdl/scripts/organize_outputs.py (works unchanged against Sprocket's
    # outputs.json too, which uses the same flat key-prefix convention) to
    # lay all of it out into a clean, browsable results folder by stage.
    # =========================================================================
    # ---- per-cell: alignment ----
    Array[File] dedup_bams = align.dedup_bam
    Array[File] dedup_bam_bais = align.dedup_bam_bai
    Array[File] flagstats = align.flagstat
    Array[File] idxstats = align.idxstats

    # ---- per-cell: QC (raw per-tool outputs; see qc_summary/multiqc for the aggregated view) ----
    Array[File] wgs_metrics = align.wgs_metrics
    Array[Array[File]] fastp_out = align.fastp_out
    Array[Array[File]] fastqc_raw_out = align.fastqc_raw_out
    Array[Array[File]] fastqc_trimmed_out = align.fastqc_trimmed_out
    Array[Array[File]] picard_out = align.picard_out

    # ---- per-cell: small variants -- PASSED CELLS ONLY (snpeff is skipped
    # for a FAILing cell; see qc_passed_sample_ids below for which these are) ----
    Array[File] deepvariant_vep_snpeff_vcfs = select_all(snpeff.final_annotated_vcf)
    Array[File] deepvariant_vep_snpeff_vcf_csis = select_all(snpeff.final_annotated_vcf_csi)
    Array[File] snpeff_htmls = select_all(snpeff.snpeff_html)
    Array[File] snpeff_genes_csvs = select_all(snpeff.snpeff_genes_csv)

    # ---- per-cell: CNV / SV ----
    Array[File] cnvkit_cns = cnvkit.cns
    Array[File] cnvkit_cnr = cnvkit.cnr
    Array[File] cnvkit_call_cns = cnvkit.call_cns
    Array[File] cnvkit_genemetrics = cnvkit.genemetrics
    Array[File] manta_vcfs = select_all(manta.manta_vcf)
    # Manta's own .tbi where it produced one, else the IndexVcf fallback's .csi
    Array[File] manta_vcf_indexes = select_all(flatten([manta.manta_vcf_tbi, manta_reindex.vcf_csi]))

    # ---- per-cell: QC gate (BaseJumper-derived composite score; ALL cells,
    # pass or fail -- see gate.qc_gate_explanation for why a specific cell
    # failed, and ../extra_methods/cell_qc_gate.py's docstring for the formula) ----
    Array[String] qc_statuses = gate.qc_status
    Array[File] qc_gate_tsvs = gate.qc_gate_tsv
    Array[File] qc_gate_explanations = gate.qc_gate_explanation
    Array[String] qc_passed_sample_ids = passed_sample_ids

    # ---- per-cell: phasing (Beagle; autosomes only, see phasing.wdl) --
    # PASSED CELLS ONLY (see qc_passed_sample_ids above for which these are) ----
    Array[File] phased_vcfs = select_all(concat_phased.phased_vcf)
    Array[File] phased_vcf_csis = select_all(concat_phased.phased_vcf_csi)

    # ---- per-cell: LOH / ploidy (Control-FREEC) -- PASSED CELLS ONLY ----
    Array[File] freec_cnvs_fixed_ploidy = select_all(freec_fixed.cnvs)
    Array[File] freec_ratio_fixed_ploidy = select_all(freec_fixed.ratio)
    Array[File] freec_info_fixed_ploidy = select_all(freec_fixed.info)
    Array[File] freec_cnvs_candidate_ploidy = select_all(freec_ploidy.cnvs)
    Array[File] freec_ratio_candidate_ploidy = select_all(freec_ploidy.ratio)
    Array[File] freec_info_candidate_ploidy = select_all(freec_ploidy.info)

    # ---- cohort: aggregated QC ----
    File multiqc_report = multiqc.report_html
    Array[File] multiqc_data = multiqc.data
    File qc_summary_tsv = qc_summary.qc_summary_tsv
    File qc_overview_png = qc_summary.qc_overview_png
    File mapd_qc_tsv = mapd.mapd_qc_tsv
    File fdr_estimate_tsv = fdr.fdr_estimate_tsv
    File artifact_panel_summary_tsv = artifact.summary_tsv
    File artifact_panel_sites_tsv = artifact.sites_tsv

    # ---- cohort: phasing / LOH summaries ----
    File phasing_summary_tsv = phasing_summary.phasing_summary_tsv
    File loh_summary_tsv = loh_summary.loh_summary_tsv
    File ploidy_summary_tsv = loh_summary.ploidy_summary_tsv
    File loh_by_region_tsv = loh_summary.loh_by_region_tsv
    File loh_by_region_png = loh_summary.loh_by_region_png

    # ---- cohort: mutational signatures ----
    File signature_contributions_tsv = mutsig.signature_contributions_tsv
    File signature_cosine_similarity_tsv = mutsig.cosine_similarity_tsv
    File mutation_spectrum_png = mutsig.mutation_spectrum_png
    File signature_barplot_png = mutsig.signature_barplot_png

    # ---- cohort: lineage + verification ----
    Array[File] lineage_outputs = lineage.lineage_outputs
    File lineage_nv_matrix = lineage.nv_matrix
    File lineage_nr_matrix = lineage.nr_matrix
    File bootstrap_support_report = lineage_verify.bootstrap_support_report
    File verification_summary_json = lineage_verify.verification_summary_json

    # ---- final dashboard ----
    File dashboard_html = dashboard.dashboard_html
  }
}
