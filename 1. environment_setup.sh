## Step 0: Create environment
mamba create -n eTAM -c bioconda umi_tools samtools cutadapt
mamba activate eTAM

# Download and compile HISAT-3N
git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
cd hisat-3n
git checkout -b hisat-3n origin/hisat-3n
make

# Download and extract pileup2var for variant calling
wget https://github.com/shunliubio/pileup2var/archive/refs/tags/v1.1.0.tar.gz
tar -xvzf v1.1.0.tar.gz
