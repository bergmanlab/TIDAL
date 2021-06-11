# TIDAL
Transposon Insertion and Depletion AnaLyzer

Detailed documentation for installing and running TIDAL1.2 is available here:
https://tidal12.readthedocs.io/en/latest/

Please use Github issues to bring up any errors that occur with software.

## Installation
```
git clone https://github.com/bergmanlab/TIDAL
cd TIDAL/
conda env create -f TIDAL.yml  --name TIDAL
conda activate TIDAL
cd CODE
make
```
## Run
```
python /path/to/TIDAL/CODE/TIDAL.py \
  -f reads.fastq \
  -r reference.fasta \
  -m repeatmasked_reference.fasta \
  -o /path/to/output
```

