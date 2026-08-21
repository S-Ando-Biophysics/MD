# Codes for MD simulation

    git clone https://github.com/S-Ando-Biophysics/MD.git
    cd MD
    
    bash 01_AMBER-1.sh
    bash 02_AMBER-2.sh
    bash 03_ACPYPE.sh
    bash 04_GROMACS.sh
    
    bash 06_PROC-INDEX.sh
    bash 07_PROC-CONV.sh
    bash 08_PROC-OUTPUT.sh

    bash ANALYSIS-CODES/A01_RMSD.sh
    bash ANALYSIS-CODES/A02_RMSF.sh
    bash ANALYSIS-CODES/A03_CLUSTER-GROMOS.sh

## Installation of required software

### Miniconda

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash ~/Miniconda3-latest-Linux-x86_64.sh
    source ~/.bashrc

### AmberTools

    tar xvfj ambertools26.tar.bz2
    cd ambertools26_src
    ./update_amber --update
    cd build
    ./run_cmake
    make install
    echo "source ~/ambertools26/amber.sh" >> ~/.bashrc
    source ~/.bashrc

### ACPYPE

    conda create -n acpype
    conda activate acpype
    conda install -c conda-forge acpype
    conda deactivate

### GROMACS

    wget https://ftp.gromacs.org/gromacs/gromacs-2026.2.tar.gz
    tar xfz gromacs-2026.2.tar.gz
    cd gromacs-2026.2
    mkdir build
    cd build
    cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA
    make -j"$(nproc)"
    sudo make install
    echo 'source /usr/local/gromacs/bin/GMXRC' >> ~/.bashrc
    source ~/.bashrc

### 3DNA

    sudo apt update
    sudo apt install ruby
    sudo su
    cd /usr/local
    tar -pzxvf x3dna-v2.4-linux-64bit.tar.gz
    cd x3dna-v2.4/bin
    ./x3dna_setup
    exit
    echo 'export X3DNA=/usr/local/x3dna-v2.4' >> ~/.bashrc
    echo 'export PATH="$X3DNA/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc

### Curves+

    tar -xvf curves+_v3.0nc.tar
    make


