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
