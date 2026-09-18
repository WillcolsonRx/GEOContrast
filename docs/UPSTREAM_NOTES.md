# GEO2R — Multi-Species v0.4.8


## v0.4.8 — optional multi-platform sample view and RNA-seq analysis

The platform selector now keeps the existing **single-platform default**, but when a selected organism has more than one GPL the user can enable **Select multiple platforms** and choose two or more platforms at once.

### Sample-table behavior

- **Default / unchecked:** one GPL is selected and only that platform's samples are shown, exactly as before.
- **Multi-platform enabled:** the Samples table shows the union of samples from every selected GPL for the current organism.
- Group assignments remain stored by GSM accession, so assignments made on one platform are preserved when platforms are hidden and reappear when that GPL is selected again.
- Table filtering, bulk assignment, auto-grouping, sample counts, and the pre-analysis Excel export all work on the currently selected platform set.
- Excel `Series Info` now records all selected GPLs.

### Analysis safeguards

Multi-platform selection does **not** mean that arbitrary expression matrices are silently merged.

- **RNA-seq:** when one validated raw-count matrix contains samples from all selected GPLs, the app can analyze the selected samples together. `Platform` is automatically added to the DESeq2 design (`~ Platform + condition`) to account for platform-associated differences. Because this is an extended model, exact official-GEO2R parity is not claimed in multi-platform mode.
- The app checks that every selected GPL is represented in the active count matrix and stops with a clear message if any selected platform is missing.
- The app also checks for a rank-deficient `Platform + condition` design. If group and platform are confounded (for example all cases are on GPL-A and all controls are on GPL-B), joint analysis is blocked because the biological group effect cannot be separated from platform.
- **Microarray:** samples from multiple GPLs can be viewed, filtered, grouped, and exported together, but platform-specific probe/expression matrices are not directly merged. Switch to one GPL before running limma.
- **Mixed assay types:** combined sample viewing/grouping is allowed, but joint differential-expression analysis is disabled.

## v0.4.6 — GEO2R RNA-seq parity fix

This release fixes a statistical mismatch that was most visible in RNA-seq `padj` values when comparing the app with official GEO2R.

### What was wrong

The previous app used a fixed low-count rule of **count >= 10 in at least 2 samples**. GEO2R's RNA-seq workflow instead pre-filters genes using **count >= 10 in at least N samples, where N is the size of the smallest assigned group**. Because Benjamini-Hochberg/FDR adjustment depends on the set of hypotheses entering `DESeq2::results()` (and DESeq2 also performs independent filtering), testing a different number of genes changes `padj` even when raw P-values are close.

The two-group pipeline also used DESeq2's default size-factor estimator (`ratio`). GEO2R's two-group Wald path uses the positive-count estimator (`sfType = "poscounts"`), which can cause small differences in `baseMean`, `log2FoldChange`, `stat`, and raw `pvalue`.

Finally, the app previously treated the second defined group as the numerator. GEO2R group creation order is meaningful: for two groups, the **first defined group is compared with the second defined group**.

### New default: GEO2R parity mode

For RNA-seq, **GEO2R parity (recommended)** is now the default. It:

- fixes the pre-filter at raw count >=10 in at least the smallest-group sample count;
- disables optional covariates, so the model is the same simple group comparison used by GEO2R;
- uses a two-group DESeq2 Wald test with `sfType = "poscounts"`;
- uses first-defined group vs second-defined group;
- uses `fdr` for the Benjamini-Hochberg option (the R alias used by GEO2R; numerically equivalent to `BH`);
- leaves DESeq2 independent filtering enabled and uses the selected significance alpha;
- reports the exact number of features before/after pre-filtering, number of non-NA adjusted P-values, independent-filter threshold, size-factor estimator, and local DESeq2 version.

An **Extended/custom** mode remains available if you intentionally want covariates or a custom low-count filter. Results from that mode are not expected to match official GEO2R exactly.

The Excel Results export now records these parity diagnostics in the **Options** sheet.

For a bit-for-bit match, the raw-count file and DESeq2/Bioconductor version must also match the version running on NCBI. The app now records the local DESeq2 version and count source so any residual difference can be diagnosed rather than hidden.


## v0.4.4 DataTables JSON hotfix

This release fixes the `DataTables warning: Invalid JSON response` dialog that could appear after defining groups or changing sample assignments. The Samples table now runs in DT server-side mode, which is compatible with `replaceData()` / `reloadData()`. When group definitions themselves change, the table is rebuilt so categorical filter levels are regenerated safely.

## v0.4.3 button/state fix

This release fixes the sample-table action buttons when filters/sorting are active. Group changes are now pushed into the existing DT widget with `DT::replaceData()` instead of rebuilding the whole table. As a result, column filters, sort order, pagination and row selections stay in place after assigning/unassigning samples. The buttons are now explicit: **Assign selected**, **Unassign selected**, **Select filtered rows**, **Clear selected rows**, plus **Assign filtered samples** and **Unassign filtered samples**. Each action reports how many samples were affected.


A GEO2R-style R/Shiny application that supports GEO microarray Series and RNA-seq Series across organisms. Microarrays are analyzed with **limma** and RNA-seq raw counts with **DESeq2**.

## v0.4.2: simplified sample filtering

The separate Cohort Builder has been removed. Large GEO Series are now filtered directly inside the main sample table, closer to the way researchers expect to work with a spreadsheet or GEO2R sample list.

### New sample-selection workflow

1. Load a GSE and choose the organism/platform. Optionally enable multi-platform mode when more than one GPL is available for the organism.
2. Click **Define groups** and create groups such as `Control`, `MDD`, `Treated`, etc.
3. Use the filter row directly under the sample-table column headings.
4. Combine filters simply by filtering more than one column. For example:
   - `Phenotype = MDD`
   - `Gender = Female`
   - `Medication = No`
   - `Smoking = No`
5. The table immediately shows only rows satisfying the active column filters.
6. Choose a group in **Assign selected to**.
7. Click **Assign filtered samples** to assign every currently filtered row at once.
8. Change the filters and repeat for the next group.
9. Click **Clear all filters** whenever you want to return to the full sample table.

### Filter types

The app prepares metadata types automatically so DT can provide appropriate column filters:

- **Categorical metadata** with a manageable number of values (phenotype, sex/gender, tissue, brain region, medication, smoking, alcohol, treatment, genotype, etc.) are converted to factors and use a multi-select dropdown filter.
- **Numeric metadata** such as age, PMI, pH, BMI, dose, duration, and similar fields use numeric/range filtering when the values are genuinely numeric.
- **High-cardinality text** such as sample titles remains a normal text-search filter.
- The existing global search box is still available for quick searches across all columns.

Multiple column filters are combined automatically, so researchers do not need to build AND/OR rules manually.

### Useful buttons

- **Assign** — assign manually selected rows.
- **Unassign** — unassign manually selected rows.
- **Select filtered** — select every row that currently passes the table filters.
- **Assign filtered samples** — directly assign every filtered row to the group selected in `Assign selected to`.
- **Unassign filtered samples** — clear group assignments for the current filtered rows.
- **Clear all filters** — restore the complete sample table.
- **Auto-group from metadata…** — still available when one metadata field already contains the desired group labels.

The header reports how many samples are currently shown, e.g. `Showing 23 of 263 samples`, so it is easy to verify a filtered sample set before assigning it.

## Excel sample export

`Download sample table Excel` is available before expression analysis and contains:

- **Samples** — complete organism/platform sample table with current group assignments.
- **Filtered Samples** — the current filtered subset when table filters are active.
- **Selected Samples** — manually selected rows when present.
- **Series Info** — GSE, title, organism, platform, defined groups, sample count, and export timestamp.

## GEO2R-style analysis

### Microarray

- Detects microarray Series Matrix data.
- Uses `limma`.
- Supports GEO2R-like log-transform / normalization / vooma options.
- Produces GEO2R-style columns such as `adj.P.Val`, `P.Value`, `t`, `B`, `logFC`, and `F` plus platform annotation.

### RNA-seq

- Uses validated non-negative integer raw counts only.
- Tries NCBI-computed counts where available and supports supplementary/uploaded count matrices otherwise.
- Uses `DESeq2`.
- Produces GEO2R-style columns such as `padj`, `pvalue`, `lfcSE`, `stat`, `log2FoldChange`, and `baseMean` plus available gene annotation.

## Install and run

Open the project in RStudio and run:

```r
source("install_packages.R")
shiny::runApp()
```

## Notes

- Group membership is stored by GSM accession, so sorting/filtering the table does not alter which sample belongs to which group.
- Table filters only determine which samples are displayed/assigned. Differential-expression analysis still uses the explicit group assignments.
- Do not use TPM/FPKM/CPM/log-normalized RNA-seq values as DESeq2 raw counts.
- Cross-species samples are not combined in one expression model; select one organism. Multiple GPLs from that organism may be displayed together, with the assay-specific analysis safeguards described above.

## Version history

- **v0.4.8** — added optional multi-platform selection, combined sample viewing/grouping/export, platform-adjusted joint RNA-seq analysis when valid, and safeguards against incompatible microarray/mixed-assay merging.
- **v0.4.6** — corrected GEO2R RNA-seq pre-filtering, two-group size-factor estimation, group-order contrast direction, and added parity diagnostics.

- **v0.4.2** — removed the Cohort Builder and replaced it with direct table-header dropdown/range/text filters plus one-click assignment of filtered samples.
- **v0.4.1** — fixed the startup parse error in the previous Cohort Builder implementation.
- **v0.3** — added automatic microarray/limma vs RNA-seq/DESeq2 routing.


## v0.4.5 — gene annotation recovery

RNA-seq result annotation is now checked automatically before DESeq2 results are built. The app:

- preserves annotation already present in NCBI/GEO or submitter count files;
- repairs feature-to-annotation matching using `FeatureID`, `GeneID`, and version-stripped Ensembl IDs;
- for human NCBI-count datasets, validates the cached `Human.GRCh38.p13.annot.tsv.gz` file, deletes stale/HTML/truncated cache files, and downloads a clean copy when needed;
- can recover annotation from `GEOquery::getRNASeqData()` row metadata when available;
- optionally uses `org.Hs.eg.db`, `org.Mm.eg.db`, or `org.Rn.eg.db` if those packages are installed and the primary annotation remains sparse;
- displays live Symbol/Description/GeneType coverage so missing annotation is visible rather than silently producing blank columns.

Blank values can still occur for genuinely unannotated/retired features; the app does not invent gene symbols.

To enable the optional local OrgDb fallback for a common organism, install the relevant package once, for example:

```r
BiocManager::install("org.Hs.eg.db")   # human
BiocManager::install("org.Mm.eg.db")   # mouse
BiocManager::install("org.Rn.eg.db")   # rat
```


## v0.4.8 — microarray annotation normalization

Microarray result annotation is now normalized across platform schemas. The app creates canonical `GeneID`, `Symbol`, and `Description` columns before retaining the original GPL annotation columns.

This specifically fixes alternative-CDF platforms such as **GPL17027** in **GSE92538**, where GEO provides `SPOT_ID` as the Gene ID and `DESCRIPTION` but no `SYMBOL` column. When `SPOT_ID` is overwhelmingly numeric, it is treated as an Entrez/NCBI Gene ID and blank symbols are recovered from NCBI Gene annotation for human datasets. Existing GEO symbols/descriptions are never overwritten.

The results page now reports microarray symbol coverage and whether annotation was recovered from an external authoritative annotation source.
