library(limma)
library(DT)
base::source("https://raw.githubusercontent.com/idbenitez/Code_lab/main/Expresion.diferencial.R")


## Read data 

## Execute differential expression

### Crear vectores de nombres 
mi_names <- mi_names
adjusted_names <- NULL

res_de <- expresion_diferencial(aux2,# Base de datos 
                                "Group", # nombre del la variable grupo
                                NULL, # Vector de texto con los nombres de las variables de ajuste 
                                mi_names , # Vector de texto con los nombres de las features
                                "record_id") #  Nombre de la variable identificador de paciente


## Extraer tabla de resultados

x <- res_de[[2]]
datatable(x, 
          extensions = 'Buttons', 
          filter = 'none', 
          options = list(columnDefs = list(list(className = 'dt-center', targets = 1:2)),
                         pageLength = nrow(x[-c(1:2), ]), 
                         dom = 'Btr', scrollX = TRUE, scrollY = "300px", 
                         deferRender = TRUE, scroller = TRUE, 
                         buttons = list(
                           list(extend = 'csv', filename = "descriptive"),
                           list(extend = 'excel', filename = "descriptive"),
                           list(extend = 'pdf', filename = "descriptive"),
                           'print'), 
                         ordering = F, language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/English.json')), 
          escape = F,  
          colnames = colnames(x), caption = "Differential Expression with Limma")

## Volcano plot 

res_de[[3]]



