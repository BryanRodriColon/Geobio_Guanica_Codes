#!/bin/bash
#SBATCH --job-name=dada2_run_BRC
#SBATCH --partition=sixhour
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=b324r112@ku.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4               # Increase threads
#SBATCH --mem=64G                       # Increase memory
#SBATCH --time=6:00:00                  # Increase runtime


# Define timestamp
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)

# Redirect output and error logs to a file with a timestamp
LOGFILE="dada2_run_${TIMESTAMP}.log"
exec > $LOGFILE 2>&1



mv $LOGFILE logs/


# Load Conda environment
module load conda
conda activate qiime2-amplicon-2024.10

# Define variables from command-line arguments
qiimedatapath="/kuhpc/work/sturm/b324r112/PROJECTS/Puerto_Rico/Amplicon/Sci_Rep_data_2019"
trimmedfile="SciRep_Data_demux_trimmed.qza"
trimmedfolder="02_trimmed_files"
dada2folder="03_denoise_files"
truncf=240
truncr=190

echo "QIIME data path: $qiimedatapath"
echo "Trimmed file: $trimmedfile"
echo "DADA2 output folder: $dada2folder"
echo "Truncation lengths: forward=$truncf, reverse=$truncr"

# Change directory to QIIME data path
cd "$qiimedatapath" || { echo "Failed to change directory to $qiimedatapath"; exit 1; }

# Create DADA2 output folder if it doesn't exist
mkdir -p "$dada2folder"
echo "Created/verified folder: $dada2folder"


# Run DADA2 denoise-paired
echo "Running DADA2 denoise-paired..."
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$trimmedfolder/$trimmedfile" \
  --p-trunc-len-f "$truncf" \
  --p-trunc-len-r "$truncr" \
  --p-chimera-method "consensus" \
  --p-trim-left-f 10 \
  --p-trim-left-r 10 \
  --p-n-threads 4 \
  --p-n-reads-learn 400000 \
  --o-table "$dada2folder/SciRep-Data-dada2-table.qza" \
  --o-representative-sequences "$dada2folder/SciRep-Data-dada2-rep-seqs.qza" \
  --o-denoising-stats "$dada2folder/SciRep-Data-dada2-stats.qza" \
  --verbose
echo "DADA2 denoise-paired completed."


# Generate visualization for denoising stats
echo "Generating visualization for DADA2 stats..."
qiime metadata tabulate \
  --m-input-file "$dada2folder/SciRep-Data-dada2-stats.qza" \
  --o-visualization "$dada2folder/SciRep-Data-dada2-stats.qzv"
echo "DADA2 stats visualization generated."

#Summary
qiime feature-table summarize \
  --i-table  "$dada2folder/SciRep-Data-dada2-table.qza" \
  --o-visualization   "$dada2folder/SciRep-Data-dada2-table.qzv" \
  --m-sample-metadata-file "/kuhpc/work/sturm/b324r112/PROJECTS/Puerto_Rico/Amplicon/Sci_Rep_data_2019/00_ref_files/SciRep_metadata.tsv"

qiime feature-table tabulate-seqs \
  --i-data "$dada2folder/SciRep-Data-dada2-rep-seqs.qza" \
  --o-visualization SciRep-Data-dada2-rep-seqs.qzv


echo "DONE!"
