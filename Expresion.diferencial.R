# Expresion diferencial 



expresion_diferencial <- function(datos, var_grupo, vars_ajuste, noms_proteines,names_id) {
  
  ## Table: Models have been adjusted for age, sex and BMI.
  datos <- datos[complete.cases(datos[,c( var_grupo, vars_ajuste,names_id)]),]
  
  # "unadjusted"
  if(is.null(vars_ajuste)){
    design <- model.matrix(as.formula(paste0("~  0 + ", var_grupo)), data = datos)
    colnames(design) <- c("gr1","gr2")
    rownames(design) <- datos[,names_id]
    cont.matrix <- makeContrasts(DIF = gr2-gr1, levels=design) 
  }
  # if "adjusted"
  else{
    design <- model.matrix(as.formula(paste0("~  ", var_grupo, "+ ", paste0(vars_ajuste, collapse="+"))), data = datos)
    colnames(design)[1:2] <- c("INTERCEPT", "gr2")
    rownames(design) <- datos[,names_id]
    cont.matrix <- makeContrasts(gr2, levels=design) 
  }
  
  aux <- t(datos[, noms_proteines])
  rownames(aux) <- noms_proteines
  lmF <- lmFit(aux, design)
  
  ## definir los contrastes ##
  cont.fit <- contrasts.fit(lmF, cont.matrix)
  eBa <- eBayes(cont.fit)
  
  # Obtencion de la tabla de resultados con los 2^-ddCt y los p-values ##
  limma <- topTable(eBa, adjust.method = "fdr", number=length(noms_proteines), sort.by="p")
  # logFC: estimate of the log2-fold-change corresponding to the effect or contrast
  #head(limma) # logFC coincide con el slope que reporta Olink en el informe
  #head(limma[order(limma$P.Value, decreasing=FALSE),])
  #limma[row.names(limma)=="Infl_IL_20",]
  
  ddCt <- data.frame(rownames(limma),10^(limma$logFC),limma$P.Value,limma$adj.P.Val)
  colnames(ddCt) <- c("Names","FC","p.value","FDR")
  
  pval_FDR0.20 <- ifelse(length(ddCt$p.value[ddCt$FDR<0.20])==0, NA, max(ddCt$p.value[ddCt$FDR<0.20]))
  
  # Reordenacion de la tabla de resultados ddCt
  ddCt <- ddCt[order(ddCt$p.value, decreasing=FALSE),]
  #DT::datatable(ddCt, extensions = "Buttons", options = list(pageLength = 20, dom = "Bftip", buttons = c("copy", "csv", "excel", "pdf")), caption = "Differentially expressed") %>% DT::formatRound(names(ddCt)[-1], 3)
  #DE <- ddCt[which((ddCt$FC>1.25 | ddCt$FC<0.8) & ddCt$p.value <0.05),"Names"]
  
  
  ### Volcano plot (adjusted Filtrado dicotomizado)
  
  # Puntos de corte:
  #  -ln(pvalue): (-1)*log(0.05)
  #             : (-1)*log(max(ddCt$p.value[ddCt$FDR<0.20]))
  #  log2(Fold Change): (-1)*log2(1.15)
  #                   : log2(1.15)
  
  cols <- c("P-value<0.05" = "blue", "FDR<0.20" = "darkgreen","Non-significant" = "Gray")
  # x = log2(FC) = log2(ddCt$FC)
  # y = (-1)log(pval) = (-1)log(ddCt$p.value)
  
  ss <- data.frame(x=eBa$coefficients[,1], y=c(-1)*log(eBa$p.value[,1]),names = rownames(eBa$lods))
  ss$color <- ifelse(ss$y >= (-1)*log(pval_FDR0.20) & !is.na(pval_FDR0.20) & (ss$x < -(1)*log10(1.15) | ss$x >= log10(1.15)),"FDR<0.20",
                     ifelse(ss$y >= (-1)*log(0.05) & (ss$x < -(1)*log10(1.15) | ss$x >= log10(1.15)),"P-value<0.05","Non-significant"))
  ss$names[ss$y < (-1)*log(0.05) & ss$x >= -(1)*log10(1.15)] <- NA
  ss$names[ss$y < (-1)*log(0.05) & ss$x < log10(1.15)] <- NA
  
  #cairo_pdf("Olink post-covid/09_results/Volcano.pdf",width = 12,height = 10)
  bp <- ggplot(ss, aes(x = x, y = y, labels=names,color = color)) + geom_point(size = 2.5, alpha = 1, na.rm = T)  +
    scale_colour_manual(values = cols) +
    ggtitle(label = "", subtitle = "") +
    geom_point(size = 2.5, alpha = 1, na.rm = T) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          axis.text=element_text(size=23),
          axis.title=element_text(size=23,face="bold"),
          legend.text = element_text(size=23)) +
    xlab(expression(log10("Fold Change"))) +
    ylab(expression(-ln("p value"))) +
    geom_hline(yintercept = (-1)*log(0.05), colour="blue", linetype="dashed") +
    geom_vline(xintercept = (-1)*log10(1.15), colour="#990000", linetype="dashed") +
    geom_vline(xintercept = log10(1.15), colour="#990000", linetype="dashed") +
    scale_y_continuous(trans = "log1p") #+ xlim(-0.1,0.1)
  # + ggrepel::geom_text_repel(aes(x = x, y = y, label=names)) 
  if(!is.na(pval_FDR0.20)) bp  <- bp + geom_hline(yintercept = (-1)*log(pval_FDR0.20), colour="darkgreen", linetype="dashed")
  
  #bp
  #dev.off()
  list(limma %>% mutate_if(is.numeric, round, digits=5),
       ddCt  %>% mutate_if(is.numeric, round, digits=3), 
       bp,
       eBa)
}

