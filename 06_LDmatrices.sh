conda activate plink2

plink --bfile gasAcu.plink19.maxCranElev --extract mce_sig_snps.txt --r2 square --allow-extra-chr --out mce_ld_matrix
plink --bfile gasAcu.plink19.maxDecel --extract md_sig_snps.txt --r2 square --allow-extra-chr --out md_ld_matrix
plink --bfile gasAcu.plink19.time_HDvMG --extract time_hdvmg_sugg_sig_snps.txt --r2 square --allow-extra-chr --out time_hdvmg_ld_matrix
plink --bfile gasAcu.plink.time_maxDecelvMG --extract time_mdvmg_sugg_sig_snps.txt --r2 square --allow-extra-chr --out time_mdvmg_ld_matrix
plink --bfile gasAcu.plink19.ttpg --extract ttpg_sig_snps.txt --r2 square --allow-extra-chr --out ttpg_ld_matrix
