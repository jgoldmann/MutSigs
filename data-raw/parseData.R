
url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
cancer_signatures = read.table(url, sep = "\t", header = T)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
cancer_signatures = as.matrix(cancer_signatures[,4:33])

devtools::use_data(cancer_signatures)
