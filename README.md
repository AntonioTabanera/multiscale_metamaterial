# multiscale_metamaterial
Multiscale model for Metamaterial calculations

El fichero MAIN_metamaterial_MEF: corre el caso de una probeta de metamaterial sin ninguna simplificación hasta alcanzar el desplazamiento especificado. Permite comparar 

El fichero MAIN_solido_continuo_hexaedros_actualizacion_youngs_tangentes_MEF resuelve el solido homogeneo con elemento tridimensionales en el rango lineal y no lineal. La curva se obtiene de la simulacion del metamaterial completo

El fichero MAIN_solido_continuo_hexaedros_solo_young_inicial_MEF resuelve el solido homogeneo con elemento tridimensionales en el rango lineal

El fichero MAIN_multiescala_hexaedros_mallado_fino_grueso_MEF: Calcula el solido homogeneo mediante mallado fino y luego con los hexaedros gruesos (los de mayor tamaño) es los que usa para acoplar con las regiones de metamaterial
				Si en este fichero se introduce en "HIPERPARAMETROS": gradiente=true, en cada step de desplazamientos, aplica el metodo del gradiente para corregir que las fuerzas deben ser continuas en dos hexaedros adyacentes (no solo los desplazamientos)--> la funcion objetivo a optimizar es la que se consigue sumando las fuerzas en los nodos comunes que estan en dos hexaedros adyacentes 
