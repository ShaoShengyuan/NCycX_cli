

conda create -n ncycx -y
conda activate ncycx
conda install -c conda-forge -c bioconda diamond blast -y
conda install -c conda-forge -c bioconda hmmer -y



conda install -c conda-forge -c bioconda seqkit -y
seqkit stat -T test_pro/*.faa | tail -n +2 | sed 's/\.faa//g' | cut -f 1,4 | xargs -I{} basename {} > SI.txt
