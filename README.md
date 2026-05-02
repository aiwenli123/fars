Raw FARS data are not included in this repository due to file size. 
To reproduce the analysis, download the FARS National CSV files from https://www.nhtsa.gov/file-downloads?p=nhtsa/downloads/FARS/ and place them in:

data/raw/

Alcohol variable decision:
For 2015–2024 FARS, person$DRINKING and vehicle$DR_DRINK do not match perfectly among driver rows. Mismatches are concentrated among cases where person$DRINKING is No/Unknown/Not Reported. vehicle$DR_DRINK always classifies drivers as Yes when observed BAC is positive, suggesting DR_DRINK incorporates BAC evidence. We therefore treat DRINKING, DR_DRINK, and BAC-threshold definitions as separate alcohol classifications in replication analyses.