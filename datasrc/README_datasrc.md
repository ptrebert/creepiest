# Selected source files

This README merely lists the selection criteria. Not all of the files matching these criteria are actually processed.

## Source: ENCODE portal (encodeproject.org)

Date: 2016-02-11

* Assay
    * ChIP-seq
    * DNase-seq
    * RNA-seq
* Project
    * ENCODE
* Genome assembly
    * hg19
    * mm9
* Biosample type
    * immortalized cell line
    * stem cell
* Available data
    * bigWig
    * bigBed narrowPeak
    * bigBed broadPeak
    * gtf
* Run type
    * single-ended
* Read length (nt)
    * 36
* Replication type
    * isogenic
    

Note 2016-04-29

* Some histone peak files contain spurious peaks (peak size as low as 1bp)
    * Example: ENCFF001MXB [entry: 1 -- chr8 - 11007712 - 11007713]
    * Filter out these peaks, i.e. require peaks to be at least "nucleosome-sized"