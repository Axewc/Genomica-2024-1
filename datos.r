# Instalar phytools
# install.packages("remotes")
# remotes::install_github("liamrevell/phytools")

# Instalar ape
# install.packages("ape")
# phytools y ape son paqueterías que nos permiten construir árboles filogenéticos con R.

# phytools contiene la base de datos que vamos a utilizar en ésta práctica:

# Llamar librerias
library(phytools)
library(ape)

# Cargar dataset
data(sunfish.data)
# Cargar arbol
data(sunfish.tree)

# Algunos ajustes extra
# Cambiar los nombres de hábito alimenticio (que es el caracter)
levels(sunfish.data$feeding.mode) <- c("no-piscivoro", "piscivoro")

# Colores bonitos, ¡CAMBIA LOS COLORES! 
cols <- setNames(c("#E4BE9E", "#71697A"), levels(sunfish.data[[1]]))

# Extraer variable con el caracter de interés
habito_alimenticio <- setNames(sunfish.data$feeding.mode, 
                               rownames(sunfish.data))