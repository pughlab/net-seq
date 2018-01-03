# LOH Suite of Scripts for Visualizations in the NET-Seq LOH project #
**formatSqzaToAbs.R** : Takes WES Sequenza output and converts it to an ABSOLUTE style output that will be used in this suite of tools

**plotRawLoh.R** : Plots a chromosomal depiction of all LOH in our samples against the pan cancer samples with each row being 1 sample track.

**aggregateLoh.R** : Takes a group of diseases and plots a Charles Swanton style chromosomal LOH visualization.  Added a Totals Row track to find areas of overlap between the samples in the given disease

**compareSeqAbs_NET.R** : Tailored to NET-seq project.  Will create aggregateLoh.R plots (Charles Swanton style) to look at overlaps in LOH between Sequenza and ABSOLUTE segmentation
