# Expresion diferencial 
# https://raw.githubusercontent.com/idbenitez/Code_lab/main/Expresion.diferencial.R


expresion_diferencial_logbase <- function(datos, var_grupo, vars_ajuste, noms_proteines, names_id, logbase = 2,
                                          pval = 0.05, FC_limit = 1.15) {
  
  ## Parametros de entrada en la función:
  ##   datos: dataframe con los datos
  ##   var_grupo: Variable para hacer la comparativa
  ##   vars_ajuste: listado de variables de ajuste o NULL si el modelo es crudo
  ##   noms_proteines: listado con los nombres de las proteinas
  ##   names_id: nombre de la variable que contiene el id del paciente
  ##   logbase: base del logaritmo con la que trabajamos
  ##   pval: pvalor utilizado para la comparativa
  ##   FC_limit: Valor de Fold Change utilizado para la comparativa
  
  
  ## Check parametre logbase  ##
  if (!logbase%in%c(2,10)) stop("'logbase' must be 2 or 10")
  
  ## Data without missing  ##
  datos <- datos[complete.cases(datos[, c(var_grupo, vars_ajuste,names_id)]),]
  
  ## "unadjusted"  ##
  if(is.null(vars_ajuste)){
    design <- model.matrix(as.formula(paste0("~  0 + ", var_grupo)), data = datos)
    colnames(design) <- c("gr1","gr2")
    rownames(design) <- datos[,names_id]
    cont.matrix <- makeContrasts(DIF = gr2-gr1, levels = design) 
  }
  ## "adjusted"  ##
  else{
    design <- model.matrix(as.formula(paste0("~  ", var_grupo, "+ ", paste0(vars_ajuste, collapse = "+"))), data = datos)
    colnames(design)[1:2] <- c("INTERCEPT", "gr2")
    rownames(design) <- datos[, names_id]
    cont.matrix <- makeContrasts(gr2, levels = design) 
  }
  
  ## Protein data (with rows corresponding to genes/proteins and columns to samples)
  aux <- t(datos[, noms_proteines])
  rownames(aux) <- noms_proteines
  
  ## lmFit: Linear Model ##
  ## (data and design matrix of the comparison groups, with rows corresponding to samples and columns to coefficients to be estimated) 
  lmF <- lmFit(aux, design)
  
  ## Definir los contrastes junto con el modelo ##
  cont.fit <- contrasts.fit(lmF, cont.matrix)
  
  ## ebayes: Empirical Bayes Statistics for Differential Expression  ##
  eBa <- eBayes(cont.fit)
    # eBa$coefficients: contrast difference by groups
    # eBa$p.value: pvalue of the contrast diffference by groups

  ## Obtencion de la tabla de resultados con los ddCt y los p-values ##
  limma <- topTable(eBa, adjust.method = "fdr", number = length(noms_proteines), sort.by = "p")
    # logFC: estimate of the log2-fold-change corresponding to the effect or contrast (for topTableF there may be several columns of log-fold-changes)
    # AveExpr: average log2-expression for the probe over all arrays and channels, same as Amean in the MarrayLM object
    # t: moderated t-statistic
    # P.Value: raw p-value
    # adj.P.Value: adjusted p-value or q-value
    # B: log-odds that the gene/protein is differentially expressed

  ## Seleccion de informacion a reportar [Delta-Delta-Ct algorithm (ddCt)]
  ddCt <- data.frame(rownames(limma), logbase^(limma$logFC), limma$P.Value, limma$adj.P.Val)
  colnames(ddCt) <- c("Names", "FC", "p.value", "FDR")
  
  ## Determinar el pvalor máximo para un FDR < 0.20
  pval_FDR0.20 <- ifelse(length(ddCt$p.value[ddCt$FDR<0.20]) == 0, NA, max(ddCt$p.value[ddCt$FDR<0.20]))
  
  ## Reordenacion de la tabla de resultados ddCt según pvalor
  ddCt <- ddCt[order(ddCt$p.value, decreasing=FALSE),]
  
  ## Volcano plot
  
  cols <- c("blue", "darkgreen","Gray")
  names(cols) <- c(paste0("P-value<",pval), "FDR<0.20", "Non-significant")
  
  
  if (logbase==10){
    # Puntos de corte:
    #  -ln(pvalue): (-1)*log(0.05)
    #             : (-1)*log(max(ddCt$p.value[ddCt$FDR<0.20]))
    #  log10(Fold Change): (-1)*log10(1.15)
    #                   : log10(1.15)
    
    ss <- data.frame(x = eBa$coefficients[,1], y = c(-1)*log(eBa$p.value[,1]),names = rownames(eBa$lods))
    ss$color <- ifelse(ss$y >= (-1)*log(pval_FDR0.20) & !is.na(pval_FDR0.20) & (ss$x < -(1)*log10(FC_limit) | ss$x >= log10(FC_limit)), "FDR<0.20",
                       ifelse(ss$y >= (-1)*log(pval) & (ss$x < -(1)*log10(FC_limit) | ss$x >= log10(FC_limit)), paste0("P-value<",pval), "Non-significant"))
    ss$names[ss$y < (-1)*log(pval) & ss$x >= -(1)*log10(FC_limit)] <- NA
    ss$names[ss$y < (-1)*log(pval) & ss$x < log10(FC_limit)] <- NA
    
    bp <- ggplot(ss, aes(x = x, y = y, labels = names, color = color)) + geom_point(size = 2.5, alpha = 1, na.rm = T)  +
      scale_colour_manual(values = cols) +
      ggtitle(label = "", subtitle = "") +
      geom_point(size = 2.5, alpha = 1, na.rm = T) +
      theme_bw(base_size = 14) +
      theme(legend.position = "top",
            legend.title = element_blank(),
            axis.text=element_text(size = 23),
            axis.title=element_text(size = 23, face = "bold"),
            legend.text = element_text(size = 23)) +
      xlab(expression(log10("Fold Change"))) +
      ylab(expression(-ln("p value"))) +
      geom_hline(yintercept = (-1)*log(pval), colour = "blue", linetype = "dashed") +
      geom_vline(xintercept = (-1)*log10(FC_limit), colour = "#990000", linetype = "dashed") +
      geom_vline(xintercept = log10(FC_limit), colour = "#990000", linetype = "dashed") +
      scale_y_continuous(trans = "log1p") + 
      annotate("text", x = log10(FC_limit), y = 0, label = paste0("FC = ", FC_limit), hjust = 0.15, vjust = 1, angle=90, colour = "#990000")+ 
      annotate("text", x = -log10(FC_limit), y = 0, label = paste0("FC = ", round(1/FC_limit,2)), hjust = 0.15, vjust = -0.25, angle=90, colour = "#990000")
    
    if(!is.na(pval_FDR0.20)) bp  <- bp + geom_hline(yintercept = (-1)*log(pval_FDR0.20), colour = "darkgreen", linetype = "dashed")
  }
  else if (logbase==2){
    # Puntos de corte:
    #  -ln(pvalue): (-1)*log(0.05)
    #             : (-1)*log(max(ddCt$p.value[ddCt$FDR<0.20]))
    #  log2(Fold Change): (-1)*log2(1.15)
    #                   : log2(1.15)

    ss <- data.frame(x = eBa$coefficients[,1], y = c(-1)*log(eBa$p.value[,1]), names = rownames(eBa$lods))
    ss$color <- ifelse(ss$y >= (-1)*log(pval_FDR0.20) & !is.na(pval_FDR0.20) & (ss$x < -(1)*log2(FC_limit) | ss$x >= log2(FC_limit)),"FDR<0.20",
                       ifelse(ss$y >= (-1)*log(pval) & (ss$x < -(1)*log2(FC_limit) | ss$x >= log2(FC_limit)),paste0("P-value<",pval),"Non-significant"))
    ss$names[ss$y < (-1)*log(pval) & ss$x >= -(1)*log2(FC_limit)] <- NA
    ss$names[ss$y < (-1)*log(pval) & ss$x < log2(FC_limit)] <- NA
    
    bp <- ggplot(ss, aes(x = x, y = y, labels = names,color = color)) + geom_point(size = 2.5, alpha = 1, na.rm = T)  +
      scale_colour_manual(values = cols) +
      ggtitle(label = "", subtitle = "") +
      geom_point(size = 2.5, alpha = 1, na.rm = T) +
      theme_bw(base_size = 14) +
      theme(legend.position = "top",
            legend.title = element_blank(),
            axis.text=element_text(size = 23),
            axis.title=element_text(size = 23, face = "bold"),
            legend.text = element_text(size = 23)) +
      xlab(expression(log2("Fold Change"))) +
      ylab(expression(-ln("p value"))) +
      geom_hline(yintercept = (-1)*log(pval), colour = "blue", linetype = "dashed") +
      geom_vline(xintercept = (-1)*log2(FC_limit), colour = "#990000", linetype = "dashed") +
      geom_vline(xintercept = log2(FC_limit), colour = "#990000", linetype = "dashed") +
      scale_y_continuous(trans = "log1p")  + 
      annotate("text", x = log2(FC_limit), y = 0, label = paste0("FC = ", FC_limit), hjust = 0.15, vjust = 1, angle=90, colour = "#990000")+ 
      annotate("text", x = -log2(FC_limit), y = 0, label = paste0("FC = ", round(1/FC_limit,2)), hjust = 0.15, vjust = -0.25, angle=90, colour = "#990000")

    if(!is.na(pval_FDR0.20)) bp  <- bp + geom_hline(yintercept = (-1)*log(pval_FDR0.20), colour = "darkgreen", linetype = "dashed")
  }

  ## Resultados a devolver  
  list(limma %>% mutate_if(is.numeric, round, digits=5),
       ddCt  %>% mutate_if(is.numeric, round, digits=3), 
       bp,
       eBa)
}