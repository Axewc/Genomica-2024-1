# Instalar phytools
install.packages("remotes")
remotes::install_github("liamrevell/phytools")

# Instalar ape
install.packages("ape")
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

# Para hacer una reconstrucción filogenética existen diversos métodos como el de parsimonia y de máxima verosimilitud.

# ¿Cuáles son las principales caracterísiticas y diferencias de éstos métodos?
#   El método de parsimonia es un método que busca la reconstrucción de un árbol filogenético que tenga el menor número de cambios evolutivos, es decir, que sea el más simple. 
#   El método de máxima verosimilitud busca la reconstrucción de un árbol filogenético que tenga la mayor probabilidad de haber generado los datos observados.

# Modelo de ajuste suponiendo tasas iguales
modelo_ER <- fitMk(sunfish.tree, habito_alimenticio, model = "ER")
#print(modelo_ER)

# Modelo de ajuste suponiendo tasas diferentes
modelo_ARD <- fitMk(sunfish.tree, habito_alimenticio, model = "ARD")
#print(modelo_ARD)

# ahora guardaremos la información de los dos modelos en un archivo de texto
sink("resultados.txt")
print(modelo_ER)
print(modelo_ARD)
sink()

# comparar modelos utilizando criterios AIC
AIC(modelo_ER, modelo_ARD)


# Según el criterio de AIC ¿Cuál es el mejor modelo ER o ARD?
#   El mejor modelo es el ARD, ya que tiene un menor valor de AIC.


# Semilla de aleatoriedad para tener resultados similares
set.seed(10)

# Mapeo estocástico con el modelo de elección, numero de escenarios = 100
# Recuerda cambiar el valor  de X
mapeo_arboles_100 <- make.simmap(sunfish.tree, habito_alimenticio,
                                 model = "ARD", nsim = 100)

# Plotear
# Rejilla para cada arbol
par(mfrow=c(10,10))
# Plotear cada simulacion
sapply(mapeo_arboles_100, plotSimmap, colors = cols, lwd = 2, ftype = "i",
       fsize = .1)
legend("bottomright", c("no piscivoro", "piscivoro"), pch = 21, pt.bg = cols,
       pt.cex = 2)


# ¿Qué ocurrió con el valor de Q? y ¿Por qué ocurrió esto?
#   El valor de Q disminuyó, esto ocurrió porque al aumentar el
#   número de escenarios, aumenta la probabilidad de que ocurra
#   un cambio evolutivo en el árbol filogenético.

# Exportar imagen a png
png("pez_luna_100X.png",
    res = 300,
    units = "in",
    height = 7,
    width = 5)
plot(summary(mapeo_arboles_100), colors = cols, ftype = "i")
legend("topleft", c("no piscivoro", "piscivoro"), pch = 21, pt.bg = cols,
       pt.cex = 2)
dev.off()

# Exportar imagen a png
# Esta vez con 1000 escenarios
mapeo_arboles_1000 <- make.simmap(sunfish.tree, habito_alimenticio,
                                  model = "ARD", nsim = 1000)
png("pez_luna_1000X.png",
    res = 300,
    units = "in",
    height = 7,
    width = 5)
plot(summary(mapeo_arboles_1000), colors = cols, ftype = "i")
legend("topleft", c("no piscivoro", "piscivoro"), pch = 21, pt.bg = cols,
       pt.cex = 2)
dev.off()


# ¿Cuál era el hábito alimentario ancestral más probable en este grupo de peces?
# Observando la gráfica, el hábito alimentario ancestral más probable
# Según el gráfico, el hábito alimentario ancestral más probable en este grupo
# de peces era el no piscívoro, es decir, que no se alimentaban de otros peces.
# La probabilidad de que el ancestro común fuera no piscívoro es del 0.86,
# mientras que la probabilidad de que fuera piscívoro es solo del 0.14

# ¿Crees que un estado de caracter es mejor que otro? ¿Cuál? ¿Por qué?
# El no piscívoro es mejor, finalmente es el que tiene una mayor probabilidad
# de ser el estado ancestral.
# Aunque no creo que haya un estado mejor que otro, ya que dependerá de las
# condiciones ambientales y de la disponibilidad de alimento.
# Pero ancestralmente, el no piscívoro es mejor para la especie.


# En los últimos años, el derretimiento de los polos se ha acelerado y se ha
# visto que algunos osos polares “están desarrollando” membranas interdigitales.
# Al parecer, los individuos con éste caracter nadan mucho más rápido y aumentan
# la probabilidad de atrapar a sus presas.
# Plantea tres hipótesis del surigmiento de éste nuevo caracter como si fueras:
# 1) Lamarck, 2) Darwin, y 3) miembro del la teoría sintética de la evolución.

# 1) Lamarck: Los osos polares desarrollaron membranas interdigitales porque
# nadaban mucho más rápido y aumentaban la probabilidad de atrapar a sus presas.

# 2) Darwin: Los osos polares desarrollaron membranas interdigitales porque
# nadaban mucho más rápido y aumentaban la probabilidad de atrapar a sus presas,
# por lo que tenían una mayor probabilidad de sobrevivir y reproducirse.

# 3) Teoría sintética de la evolución: Los osos polares desarrollaron membranas
# interdigitales porque nadaban mucho más rápido y aumentaban la probabilidad
# de atrapar a sus presas, por lo que tenían una mayor probabilidad de
# sobrevivir y reproducirse, y por lo tanto, de transmitir sus genes a la
# siguiente generación.
