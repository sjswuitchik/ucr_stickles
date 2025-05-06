conda activate plink2

plink --bfile gasAcu.plink19.maxCranElev --extract mce_sig_snps.txt --r2 square --allow-extra-chr --out mce_ld_matrix
plink --bfile gasAcu.plink19.maxDecel --extract md_sig_snps.txt --r2 square --allow-extra-chr --out md_ld_matrix
