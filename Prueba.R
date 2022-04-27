library(devtools)
library(limma)
library(DT)

SourceURL <- "https://raw.github.com/christophergandrud/christophergandrud.github.com/master/SourceCode/CarsScatterExample.R"
source_url(SourceURL)

res_de <- expresion_diferencial(aux2, "Group",NULL,mi_names ,"record_id")
DE <- res_de[[2]][which((res_de[[2]]$FC>1.25 | res_de[[2]]$FC<0.8) & res_de[[2]]$p.value <0.05),"Names"]

