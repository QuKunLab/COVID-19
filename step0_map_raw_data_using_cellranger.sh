cellranger count --id=P1-1r1_outs --transcriptome=~/refdata-cellranger-GRCh38-3.0.0/ --fastqs=P1-1r1/
cellranger count --id=P1-1r2_outs --transcriptome=~/refdata-cellranger-GRCh38-3.0.0/ --fastqs=P1-1r2/
cellranger count --id=P1-2r1_outs --transcriptome=~/refdata-cellranger-GRCh38-3.0.0/ --fastqs=P1-2r1/
cellranger count --id=P1-2r2_outs --transcriptome=~/refdata-cellranger-GRCh38-3.0.0/ --fastqs=P1-2r2/
cellranger count --id=P2-1_outs --transcriptome=~/refdata-cellranger-GRCh38-3.0.0/ --fastqs=P2-1/
cellranger count --id=P2-2_outs --transcriptome=~/refdata-cellranger-GRCh38-3.0.0/ --fastqs=P2-2/
cellranger count --id=P2-3_outs --transcriptome=~/refdata-cellranger-GRCh38-3.0.0/ --fastqs=P2-3/
cellranger aggr --id=COVID19 --csv=COVID19.csv --normalize=mapped