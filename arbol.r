# Semilla de aleatoriedad para tener resultados similares
set.seed(10)

# Mapeo estocástico con el modelo de elección, numero de escenarios = 100
# Recuerda cambiar el valor  de X
mapeo_arboles_100 <-make.simmap(sunfish.tree,habito_alimenticio,model="ARD",nsim=100)

# Plotear
# Rejilla para cada arbol
par(mfrow=c(10,10))
# Plotear cada simulacion
sapply(mapeo_arboles_100, plotSimmap, colors=cols, lwd=2, ftype="i", fsize=.1)
legend("bottomright",c("no piscivoro","piscivoro"),pch=21,pt.bg=cols,pt.cex=2)

# Guardar cada arbol en un archivo png el nombre de cada archivo debe ser diferente
# Recuerda cambiar el valor  de X
png("pez_luna_100X.png",
    res = 300,
    units = "in",
    height = 7,
    width = 5)
plot(summary(mapeo_arboles_100),colors=cols,ftype="i")
legend("topleft",c("no piscivoro","piscivoro"),pch=21,pt.bg=cols,pt.cex=2)
dev.off()


# ¿Qué ocurrió con el valor de Q? y ¿Por qué ocurrió esto?
#   El valor de Q disminuyó, esto ocurrió porque al aumentar el número de escenarios, 
#   aumenta la probabilidad de que ocurra un cambio evolutivo en el árbol filogenético.

# Exportar imagen a png
png("pez_luna_100X.png",
    res = 300,
    units = "in",
    height = 7,
    width = 5)
plot(summary(mapeo_arboles_100),colors=cols,ftype="i")
legend("topleft",c("no piscivoro","piscivoro"),pch=21,pt.bg=cols,pt.cex=2)
dev.off()

# ¿Cuál era el hábito alimentario ancestral más probable en este grupo de peces?
# Observando la gráfica, el hábito alimentario ancestral más probable es el de piscivoro.
# y su probabilidad es de 0.9999.

# ¿Crees que un estado de caracter es mejor que otro? ¿Cuál? ¿Por qué?.
#   Creo que el estado de caracter mejor es el de piscivoro
#   porque tiene una probabilidad de 0.9999, mientras que el de no piscivoro tiene una probabilidad de 0.0001.


# En los últimos años, el derretimiento de los polos se ha acelerado y 
# se ha visto que algunos osos polares “están desarrollando” membranas interdigitales. 
# Al parecer, los individuos con éste caracter nadan mucho más rápido y aumentan la probabilidad de atrapar a sus presas. 
# Plantea tres hipótesis del surigmiento de éste nuevo caracter como si fueras: 1) Lamarck, 2) Darwin, y 3) un miembro del equipo que desarrolló la teoría sintética de la evolución.
# 1) Lamarck: Los osos polares desarrollaron membranas interdigitales porque nadaban mucho más rápido y aumentaban la probabilidad de atrapar a sus presas.
# 2) Darwin: Los osos polares desarrollaron membranas interdigitales porque nadaban mucho más rápido y aumentaban la probabilidad de atrapar a sus presas, por lo que tenían una mayor probabilidad de sobrevivir y reproducirse.
# 3) Teoría sintética de la evolución: Los osos polares desarrollaron membranas interdigitales porque nadaban mucho más rápido y aumentaban la probabilidad de atrapar a sus presas, por lo que tenían una mayor probabilidad de sobrevivir y reproducirse, y por lo tanto, de transmitir sus genes a la siguiente generación.
