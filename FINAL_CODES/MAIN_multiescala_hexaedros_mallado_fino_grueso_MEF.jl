include(raw"celdillas_metamaterial.jl"); using .celdillas_metamaterial 
include(raw"funciones_MEF_metamaterial.jl"); using .funciones_MEF_metamaterial
include(raw"representar_paraview.jl"); using .representar_paraview
include(raw"representar_paraview_sucesion.jl"); using .representar_paraview_sucesion
include(raw"mallado_solido_tridimensional_hexaedros.jl"); using .mallado_solido_tridimensional_hexaedros
include(raw"curva_sigma_epsilon_conjunto.jl"); using .curva_sigma_epsilon_conjunto
include(raw"funciones_MEF_solido_hexaedros.jl"); using .funciones_MEF_solido_hexaedros
include(raw"funciones_acoplamiento_solido_metamaterial_hexaedros.jl"); using .funciones_acoplamiento_solido_metamaterial_hexaedros
include(raw"curva_sigma_epsilon_material.jl"); using .curva_sigma_epsilon_material
Young_material()
using Plots
using LinearAlgebra


### GENERAR LOS FICHEROS DE RESULTADOS
numero_fichero="prueba2023" # "FLEX_25"

#generar_tabla_tension_deformacion_Osgood()
fichero_ensayo=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_energ_element_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_2=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_inc_energ_element_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_3=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_energ_total_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_4=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_inc_energ_total_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_5=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_fuerza_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_6=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_energ_zona_metam_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_7=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_inc_energ_zona_metam_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_8=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_energ_metam_total_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_9=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_inc_energ_metam_total_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_10=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_curva_tension_deformacion_hexaedros_gruesos_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_11=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_deformacion_barra_max_hexaedros_gruesos_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo



### DIMENSIONES PROBETA RECTANGULAR Y CELDILLAS METAMATERIAL 
o=8
Nx=o#10 #numero de celdillas en x
Ny=o#6 #numero de celdillas en y
Nz=o#2   #numero de celdillas en z

X=10*Nx    #dimension x probeta
Y=10*Ny #dimension y probeta
Z=10*Nz #dimension z probeta
radio=5  #radio de las barras de metamaterial


### NUMERO DE DIVISIONES EN LOS DIFERENTES EJES PARA GENERAR EL MALLADO DEL SOLIDO TRIDIMENSIONAL GRUESO
p=4
N_3D_x_grueso=p
N_3D_y_grueso=p
N_3D_z_grueso=p


### NUMERO DE DIVISIONES EN LOS DIFERENTES EJES PARA GENERAR EL MALLADO DEL SOLIDO TRIDIMENSIONAL FINO
q=8
N_3D_x_fino=q
N_3D_y_fino=q
N_3D_z_fino=q



### HIPERPARAMETROS 
iter_desplaz=2000;      #Numero de iteraciones en desplazamiento hasta alcanzar el valor del desplazamiento final
iter_max=3;             #Numero de iteraciones para converger en cada desplazamiento
tolerancia=1e-5;        #Diferencia entre el valor SSR de dos iteraciones sucesivas para salir de las iteraciones
tolerancia_2=1e-6;      #Diferencia entre el valor SSR de dos iteraciones no sucesivas para salir de las iteraciones
diferencial=1e-8        #Diferencial de longitud para calcular derivadas
eps_limite=1e-12        #Valor de deformacion umbral que si supera un elemento hexaedrico, entoces la region pasa a analizarse con metamaterial
modelo_calculo="modelo_daño_hexaedros_gruesos"  
                                                #"modelo_energia_hexaedros_gruesos": Este modelo actualiza a todos los hexaedros finos con el mismo valor que el del hexaedro grueso. Los actualiza ajustando el Young para que la energía se iguale con la de la zona calculada por el metamaterial
                                                #"modelo_daño_hexaedros_finos": Este modelo actualiza a todos los hexaedros finos con valores distintos a partir de el del hexaedro grueso. Los actualiza ajustando el Young para que la energía se iguale con la de la zona calculada por el metamaterial
                                                #"modelo_daño_hexaedros_gruesos": Este modelo actualiza a todos los hexaedros finos con el mismo valor que el del hexaedro grueso. Los actualiza ajustando el Young en funcion del parametro de daño

exponente_reparto_young=0   #Si el exponente es cero, entonces, se les pone a todos los hexaedros finos el mismo Young que al hexaedro grueso

gradiente=true#false#       #Para aplicar o no el método del gradiente para redcir errores al hacer el acoplamiento multiescala

# HIPERPARAMETROS ADICIONALES PARA EL METODO DEL GRADIENTE
step_gradiente_desplaz=0.02*10/(194000*radio^2)         #Constante para definir el nuevo valor actualizado de desplazamiento en el metodo de gradiente
step_gradiente_giros=0.02/(10^2)*10^3/(194000*radio^4)  #Constante para definir el nuevo valor actualizado de giro en el metodo de gradiente
iter_max_solido_grad=1      #Numero de iteraciones para converger en cada desplazamiento (sustituye a iter_max)
iter_max_gradiente=5        #Numero de iteraciones maximas para cada vez que se intenta corregir con el metodo de gradiente



hiperparametros=Array{Any}(undef,1,7)
hiperparametros[2]=iter_max;  hiperparametros[3]=tolerancia; hiperparametros[4]=tolerancia_2; hiperparametros[5]=diferencial
hiperparametros[1]=eps_limite
hiperparametros[6]=modelo_calculo
hiperparametros[7]=exponente_reparto_young
hiperparametros_gradiente=hcat(hiperparametros, [step_gradiente_desplaz step_gradiente_giros iter_max_gradiente])
hiperparametros_gradiente[2]=iter_max_solido_grad



########################### DISTINTAS RESTRICCIONES

### RESTRICCIONES PARA TRACCION, TORSION, FLEXION PURA
# restriccion 1
punto_restriccion=[0 "all" "all"]#[0 "all" "all"]      # Coordenadas del nodo [x,y,z]
gdl_restriccion=[1 1 1]#[1 0 0]#                       # grados de libertad: [ux, uy, uz]


### RESTRICCIONES PARA TRACCION SIN EFECTOS DE LA BASE
# # restriccion 1
# punto_restriccion=[0 "all" "all"]
# gdl_restriccion=[1 0 0]
# # restriccion 2
# punto=[0 Y/2 "all"]
# gdl=[0 1 0]#[0 0 1]
# punto_restriccion=vcat(punto_restriccion,punto)
# gdl_restriccion=vcat(gdl_restriccion,gdl)
# # restriccion 3
# punto=[0 "all" Z/2]
# gdl=[0 0 1]#[0 1 0]
# punto_restriccion=vcat(punto_restriccion,punto)
# gdl_restriccion=vcat(gdl_restriccion,gdl)


### RESTRICCIONES PARA FLEXION A 3 PUNTOS
# # restriccion 1
# punto_restriccion=[0 "all" "all"]
# gdl_restriccion=[1 1 1]
# # restriccion 2
# punto= [X "all" "all"]
# gdl=[1 1 1]
# punto_restriccion=vcat(punto_restriccion,punto)
# gdl_restriccion=vcat(gdl_restriccion,gdl)


########################### DISTINTOS TIPOS DE CARGA

### TRACCION 
# # desplazamiento 1
# punto_desplazamiento=[X "all" "all"]#
# valor_desplazamiento=[1.0e1 0  0]


### CARGADO EN UN EXTREMO Y EMPOTRADO EN OTRO
# desplazamiento 1
punto_desplazamiento=[X "all" "all"]#
valor_desplazamiento=[ 0  0 10e0 ]


### FLEXION PURA (UN GIRO)
# # desplazamiento 1
# punto_desplazamiento=[X "all" Z]# [X "all" "all"]# 
# valor_desplazamiento=[1.5e-6 0  0 ]#[ 0  0 1.5e-6]#
# # desplazamiento 2
# punto=[X "all" 0]
# valor=[-1.5e-6 0 0 ]
# punto_desplazamiento=vcat(punto_desplazamiento,punto)
# valor_desplazamiento=vcat(valor_desplazamiento,valor)

# ## Para el ensayo de flexion pura, siempre que N_3D_z sea par
# for i in 1:N_3D_z/2-1
#     global punto, valor, punto_desplazamiento, valor_desplazamiento
#     punto=[X "all" Z/N_3D_z*i+Z/2]
#     valor=[1.5e-6*2*i/N_3D_z  0  0 ]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
#     punto=[X "all" -Z/N_3D_z*i+Z/2]
#     valor=[-1.5e-6*2*i/N_3D_z  0  0]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
# end


### TORSION
# # desplazamiento 1
# punto_desplazamiento=[X "all" Z]#[X "all" "all"]#
# valor_desplazamiento=[0 -1.5e-6  0 ]#[0 0 0  0 1e-6 0]#
# # desplazamiento 2
# punto=[X "all" 0]
# valor=[0 1.5e-6  0]#
# punto_desplazamiento=vcat(punto_desplazamiento,punto)
# valor_desplazamiento=vcat(valor_desplazamiento,valor)
# # desplazamiento 3
# punto=[X 0 "all"]
# valor=[ 0 0 -1.5e-6]#
# punto_desplazamiento=vcat(punto_desplazamiento,punto)
# valor_desplazamiento=vcat(valor_desplazamiento,valor)
# # desplazamiento 4
# punto=[X Y "all"]
# valor=[ 0 0 1.5e-6]#
# punto_desplazamiento=vcat(punto_desplazamiento,punto)
# valor_desplazamiento=vcat(valor_desplazamiento,valor)

# # para el ensayo de torsion, siempre que N_3D_z y N_3D_y sean par
# for i in 1:N_3D_z_fino/2-1
#     global punto, valor, punto_desplazamiento, valor_desplazamiento
#     punto=[X "all" Z/N_3D_z_fino*i+Z/2]
#     valor=[ 0 -1.5e-6*2*i/N_3D_z_fino 0 ]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
#     punto=[X "all" -Z/N_3D_z_fino*i+Z/2]
#     valor=[ 0 1.5e-6*2*i/N_3D_z_fino 0]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
# end
# for i in 1:N_3D_y_fino/2-1
#     global punto, valor, punto_desplazamiento, valor_desplazamiento
#     punto=[X Y/N_3D_y_fino*i+Y/2 "all" ]
#     valor=[ 0 0 1.5e-6*2*i/N_3D_y_fino ]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
#     punto=[X -Y/N_3D_y_fino*i+Y/2 "all"]
#     valor=[0  0 -1.5e-6*2*i/N_3D_y_fino]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
# end


### FLEXION A 3 PUNTOS (con N_3D_x par)
# # desplazamiento 1
# punto_desplazamiento=[X/2 "all" "all"]
# valor_desplazamiento=[0 0  -1.5e-6 ]


### RESOLUCION PARA INCREMENTOS DE DESPLAZAMIENTOS
valor_desplazamiento=valor_desplazamiento/iter_desplaz


### COMPROBAR QUE HAY UN NUMERO ENTERO DE HEXAEDROS FINOS EN CADA LADO POR CADA HEXAEDRO GRUESO
if N_3D_x_fino % N_3D_x_grueso==0
    if N_3D_y_fino % N_3D_y_grueso==0
        if N_3D_z_fino % N_3D_z_grueso==0    
        else
            println("No hay un numero entero de hexaedros finos en z por cada uno grueso")
        end
    else
        println("No hay un numero entero de hexaedros finos en y por cada uno grueso")
    end
else
    println("No hay un numero entero de hexaedros finos en x por cada uno grueso")
end


######################################### GENERAR LA GEOMETRIA ############################
println("##### INICIO DEL CALCULO")
### MALLADO GRUESO 3D DEL SOLIDO CONTINUO
coords_3D_grueso, conectividad_3D_grueso=mallado_3D_hexaedros(X,Y,Z,N_3D_x_grueso,N_3D_y_grueso,N_3D_z_grueso)
volumenes_hexaedros_grueso= volumenes_elementos_hexaedros(coords_3D_grueso, conectividad_3D_grueso)
vertices_caras_grueso, caras_hexaedros_grueso =numeracion_caras_hexaedros(conectividad_3D_grueso)
caras_contorno_solido_grueso =determinar_caras_contorno_solido(X,Y,Z, vertices_caras_grueso, coords_3D_grueso)
println("##### MALLADO GRUESO 3D DEL SOLIDO CONTINUO")

### MALLADO FINO 3D DEL SOLIDO CONTINUO
coords_3D_fino, conectividad_3D_fino=mallado_3D_hexaedros(X,Y,Z,N_3D_x_fino,N_3D_y_fino,N_3D_z_fino)
volumenes_hexaedros_fino= volumenes_elementos_hexaedros(coords_3D_fino, conectividad_3D_fino)
vertices_caras_fino, caras_hexaedros_fino =numeracion_caras_hexaedros(conectividad_3D_fino)
caras_contorno_solido_fino =determinar_caras_contorno_solido(X,Y,Z, vertices_caras_fino, coords_3D_fino)
println("##### MALLADO FINO 3D DEL SOLIDO CONTINUO")

### GENERAR METAMATERIAL
coords_metam,conectividad_metam=celdilla_FCC(X,Y,Z,Nx,Ny,Nz)
desplazamientos_metam=zeros(length(coords_metam[:,1]),3);giros_metam=zeros(length(coords_metam[:,1]),3);fuerzas_metam=zeros(length(coords_metam[:,1]),3);momentos_metam=zeros(length(coords_metam[:,1]),3)
plot_paraview_estructura_metamaterial(coords_metam,conectividad_metam,desplazamientos_metam,giros_metam,fuerzas_metam,momentos_metam)
println("##### GENERADO METAMATERIAL")

### DETERMINAR QUE NODOS DEL METAMATERIAL ESTAN DENTRO DE CADA ELEMENTO HEXAEDRICO GRUESO
coords_metam, conectividad_metam, conectiv_metam_hexaedros_exacto, conectiv_metam_hexaedros_extendido, conectiv_metam_hexaedros_reducido, porcentajes_barras_dentro_hexaedros,coords_metam_hexaedros_exacto, coords_metam_hexaedros_extendido, coords_metam_hexaedros_reducido, coords_metam_fronteras_exactas_hexaedros, tensiones_pandeo_elementos = nodos_metamaterial_en_elementos( X, Y, Z, coords_metam, conectividad_metam, coords_3D_grueso, conectividad_3D_grueso, volumenes_hexaedros_grueso, radio) 
println("##### DETERMINADOS QUE NODOS DEL METAMATERIAL ESTAN DENTRO DE CADA ELEMENTO HEXAEDRICO GRUESO")

### DETERMINAR QUE HEXAEDROS DEL MALLADO FINO ESTAN DENTRO DE CADA HEXAEDRO DEL GRUESO
hexaedros_mallado_fino_en_grueso = determinar_hexaedros_mallado_fino_en_grueso(N_3D_x_fino, N_3D_y_fino, N_3D_z_fino, N_3D_x_grueso, N_3D_y_grueso, N_3D_z_grueso)
println("##### DETERMINADOS QUE HEXAEDROS DEL MALLADO FINO ESTAN DENTRO DE CADA HEXAEDRO DEL GRUESO")

### DETERMINAR QUE NODOS DEL METAMATERIAL ESTAN DENTRO DE CADA HEXAEDRO DEL MALLADO FINO
nodos_metam_en_hexaedros_finos = clasificar_nodos_metam_en_hexaedros(coords_metam, coords_3D_fino, conectividad_3D_fino, volumenes_hexaedros_fino)
println("##### DETERMINADOS QUE NODOS DEL METAMATERIAL ESTAN DENTRO DE CADA HEXAEDRO DEL MALLADO FINO")

### FACTOR BARRAS EN LAS CARAS (GRADIENTE) --> SU CONTRIBUCION SE MULTIPLICA POR 2
factor_barras_caras_solido=factor_barras_caras_solido_gradiente(X,Y,Z,conectividad_metam, coords_metam)


#$$$$
println("##### geometria generada")

# for i in 1:64
#     println(length(findall(in(0.5),collect(values(porcentajes_barras_dentro_hexaedros[i])))))
# end

# for elemento in 1:length(conectiv_metam_hexaedros_exacto)
#     for i in 1:length(conectiv_metam_hexaedros_exacto[elemento][:,3])
#         barra=conectiv_metam_hexaedros_exacto[elemento][i,3]
#         println(elemento," ",conectiv_metam_hexaedros_exacto[elemento][i,1]," ",conectiv_metam_hexaedros_exacto[elemento][i,2]," ",porcentajes_barras_dentro_hexaedros[elemento][barra])
#     end
# end


# REPRESENTAR EN PARAVIEW EL MALLADO Y EL METAMATERIAL
iters="geometria_grueso";u=zeros(3*length(coords_3D_grueso[:,1]));f=zeros(3*length(coords_3D_grueso[:,1]))
plot_paraview_solido_sucesion_hexaedros(coords_3D_grueso, conectividad_3D_grueso, u, f, iters)

iters="geometria_fino";u=zeros(3*length(coords_3D_fino[:,1]));f=zeros(3*length(coords_3D_fino[:,1]))
plot_paraview_solido_sucesion_hexaedros(coords_3D_fino, conectividad_3D_fino, u, f, iters)


# FRONTERAS INTERIORES
dict_cero=Dict()
for elemento in 1:length(conectividad_3D_grueso[:,1])
    global dict_cero
    for nodo_metam in coords_metam_hexaedros_reducido[elemento]
        dict_cero[nodo_metam]=[0,0,0]
    end
    iter="geometria_interior"
    plot_paraview_estructura_metamaterial_multiescala_sucesion_hexaedros(coords_metam_hexaedros_reducido[elemento],conectiv_metam_hexaedros_reducido[elemento],dict_cero,dict_cero,dict_cero,dict_cero,elemento, iter)
end

# FRONTERAS EXTERIORES
dict_cero=Dict()
for elemento in 1:length(conectividad_3D_grueso[:,1])
    global dict_cero
    for nodo_metam in coords_metam_hexaedros_extendido[elemento]
        dict_cero[nodo_metam]=[0,0,0]
    end
    iter="geometria_exterior"
    plot_paraview_estructura_metamaterial_multiescala_sucesion_hexaedros(coords_metam_hexaedros_extendido[elemento],conectiv_metam_hexaedros_extendido[elemento],dict_cero,dict_cero,dict_cero,dict_cero,elemento, iter)
end

# FRONTERA EXACTA
dict_cero=Dict()
for elemento in 1:length(conectividad_3D_grueso[:,1])
    global dict_cero
    for nodo_metam in coords_metam_hexaedros_exacto[elemento]
        dict_cero[nodo_metam]=[0,0,0]
    end
    iter="geometria_exacta"
    plot_paraview_estructura_metamaterial_multiescala_sucesion_hexaedros(coords_metam_hexaedros_exacto[elemento],conectiv_metam_hexaedros_exacto[elemento],dict_cero,dict_cero,dict_cero,dict_cero,elemento, iter)
end


### AJUSTAR DE MANERA MAS PRECISA LOS ELEMENTOS QUE ESTAN EN LAS FRONTERAS (PARA CONOCER MEJOR CUANTA PARTE DE LAS BARRAS ESTA FUERA DEL hexaedro Y CUANTA DENTRO)
#println(volumenes_hexaedros)
#coeficientes_ajuste_geometria = ajuste_coeficientes_geometria_multiescala(X,Y,Z, coords_3D, conectividad_3D, volumenes_hexaedros, coords_metam_hexaedros_extendido, coords_metam_hexaedros_reducido, conectiv_metam_hexaedros_extendido, conectiv_metam_hexaedros_reducido, Nodos_metamaterial_hexaedros_semi_ext, radio)
#println();println("AJUSTE DE COEFICIENTES");println()
#println(coeficientes_ajuste_geometria);println()
#readline()


### INICIALIZAR VARIABLES SOLIDO (SE ESTUDIA EL QUE TIENE MALLADO FINO) 
N_hexaedros_finos=length(conectividad_3D_fino[:,1]); N_nodos_solido=length(coords_3D_fino[:,1]); N_hexaedros_gruesos=length(conectividad_3D_grueso[:,1]);
Youngs_tangentes_hexaedros_finos=ones(N_hexaedros_finos)*curva_tension_deformacion_conjunto(10) #modulo de young inicial igual al de referencia en todos los hexaedros
u_solido_acumulado=zeros(3*N_nodos_solido); f_solido_acumulado=zeros(3*N_nodos_solido); epsilon_voigt_acumulado=zeros(6,8,N_hexaedros_finos)
energias_hexaedros_acumulado=zeros(N_hexaedros_finos); energia_acumulada_solido=0
fuerzas_hexaedros_acumulado=zeros(24, N_hexaedros_finos)
hexaedros_analizados_metam=[]
Youngs_tangentes_hexaedros_gruesos=ones(N_hexaedros_gruesos)*curva_tension_deformacion_conjunto(10)


### INICIALIZAR VARIABLES METAMATERIAL
gdl_restring_globales_metam=[]; inc_valor_restricciones_fijos_global_metam=[]
gdl_total_zona_hexaedros=Dict(); gdl_fijos_zona_hexaedros=Dict()    
gdl_total_hexaedro=Dict(); gdl_fijos_hexaedro=Dict(); gdl_restring_totales_hexaedro_local=Dict();gdl_restring_fijos_hexaedro_local=Dict(); gdl_restring_no_fijos_hexaedro_local=Dict()
Youngs_tangentes_memoria_metam=Dict();  #este vector se ira llenando a medida que se vayan dañando los hexaedros
eps_x_ptos_int_acumulado_memoria=0; eps_xy_ptos_int_acumulado_memoria=0; eps_xz_ptos_int_acumulado_memoria=0;
u_metam_hexaedros_acumulado=Dict(); f_metam_hexaedros_acumulado=Dict() # diccionarios que se llenaran con los desplazamientos y fuerzas de las distintas zonas de metamaterial a medida que se dañen los hexaedros
f_acumulada_elementos_memoria=Dict(); energias_acumuladas_elementos_memoria=Dict(); invariante_def_memoria_metam=Dict()
gdl_contorno_solido_completo_no_fronteras=[]


### INICIALIZAR VARIABLES GRADIENTE
gdl_no_fijos_total=[]; gdl_no_fijos_zona_hexaedros=Dict(); inc_valor_restricciones_no_fijas_total=[] ; K_step=[]; gradiente_fuerzas_nodos=[]


inc_energias_solidos_gruesos=[]; energias_solidos_gruesos=[]; energia_elemento_solido_grueso_metam=[]; inc_energia_elemento_solido_grueso_metam=[]
fuerza_reaccion_x_iters=[]; fuerza_reaccion_y_iters=[]; fuerza_reaccion_z_iters=[]; desplaz_iters=[]
deform_solidos_gruesos=zeros(N_hexaedros_gruesos); tension_solidos_gruesos=zeros(N_hexaedros_gruesos);
datos_contribucion_barras_gdl=Dict();

for iter_desplazamiento in 1:iter_desplaz

    println();println("##############  ITERACION EN DESPLAZAMIENTO ",iter_desplazamiento,"  ##############");println()

    global coords_zonas_elementos, conectividad_zonas_elementos, hiperparametros, hexaedros_mallado_fino_en_grueso, nodos_metam_en_hexaedros_finos, hexaedros_analizados_metam, f_solido_acumulado, fuerzas_hexaedros_acumulado ,energias_hexaedros_acumulado, Youngs_tangentes_hexaedros_finos, u_metam_hexaedros_acumulado, f_metam_hexaedros_acumulado, Youngs_tangentes_memoria_metam, energia_acumulada_solido, f_solido_acumulado, u_solido_acumulado, f_elementos_acumulado, punto_restriccion, gdl_restriccion, punto_desplazamiento, valor_desplazamiento, eps_x_ptos_int_acumulado_memoria, eps_xy_ptos_int_acumulado_memoria, eps_xz_ptos_int_acumulado_memoria, epsilon_voigt_acumulado, coords_metam_hexaedros_exacto, coords_metam_hexaedros_reducido, conectiv_metam_hexaedros_reducido, Nodos_metamaterial_hexaedros_semi_ext, gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_no_fijos_hexaedro_local, gdl_restring_fijos_hexaedro_local, iter_desplaz, f_acumulada_elementos_memoria, energias_acumuladas_elementos_memoria, fuerza_reaccion_x_iters, fuerza_reaccion_y_iters, fuerza_reaccion_z_iters, desplaz_iters, inc_energias_solidos_gruesos, energias_solidos_gruesos, gdl_contorno_solido_completo_no_fronteras, iter_max , energia_elemento_solido_grueso_metam, inc_energia_elemento_solido_grueso_metam, Youngs_tangentes_hexaedros_gruesos, deform_solidos_gruesos, tension_solidos_gruesos, N_hexaedros_gruesos, fichero_ensayo, fichero_ensayo_2, fichero_ensayo_3, fichero_ensayo_4, fichero_ensayo_5, fichero_ensayo_6, fichero_ensayo_7, fichero_ensayo_8, fichero_ensayo_9, fichero_ensayo_10, fichero_ensayo_11, inc_u_solido, invariante_def_memoria_metam, gdl_no_fijos_total, gdl_no_fijos_zona_hexaedros, gdl_restring_totales_hexaedro_local, hiperparametros_gradiente, factor_barras_caras_solido, inc_valor_restricciones_no_fijas_total, K_step, gradiente_fuerzas_nodos, datos_contribucion_barras_gdl

    #println(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[:,1]])
    #para que inicialmente converja mejor
    if iter_desplazamiento == 1
        hiperparametros[2]=iter_max
    else
        hiperparametros[2]=iter_max;
    end


    ### RESOLUCION ITERATIVA
    # Con la zona exacta
    if gradiente == false
        hexaedros_analizados_metam, inc_u_solido, inc_f_solido, u_solido_acumulado,  f_solido_acumulado, fuerzas_hexaedros_acumulado, energias_hexaedros_acumulado, inc_energia_solido, inc_energias_solidos_gruesos, energias_solidos_gruesos, inc_deform_equiv_solidos_gruesos, inc_u_metam_hexaedros, inc_f_metam_hexaedros, f_acumulada_elementos_memoria, Youngs_tangentes_hexaedros_finos, Youngs_tangentes_hexaedros_gruesos, Youngs_tangentes_memoria_metam, epsilon_voigt_acumulado, energias_acumuladas_elementos_memoria, eps_x_ptos_int_acumulado_memoria, eps_xy_ptos_int_acumulado_memoria, eps_xz_ptos_int_acumulado_memoria , gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_totales_hexaedro_local, gdl_contorno_solido_completo_no_fronteras, energia_elemento_solido_grueso_metam, inc_energia_elemento_solido_grueso_metam, invariante_def_memoria_metam = resolucion_iterativa_acoplamiento_zonas_mallado_fino_grueso(hiperparametros, X,Y,Z, N_3D_x_grueso, N_3D_y_grueso, N_3D_z_grueso, hexaedros_analizados_metam, hexaedros_mallado_fino_en_grueso, nodos_metam_en_hexaedros_finos, u_solido_acumulado, f_solido_acumulado, fuerzas_hexaedros_acumulado, energias_hexaedros_acumulado, epsilon_voigt_acumulado, f_acumulada_elementos_memoria, energias_acumuladas_elementos_memoria, eps_x_ptos_int_acumulado_memoria, eps_xy_ptos_int_acumulado_memoria, eps_xz_ptos_int_acumulado_memoria, coords_metam, conectividad_metam ,coords_3D_fino, conectividad_3D_fino, conectividad_3D_grueso, volumenes_hexaedros_grueso, Youngs_tangentes_hexaedros_finos, Youngs_tangentes_hexaedros_gruesos, Youngs_tangentes_memoria_metam, tensiones_pandeo_elementos, punto_restriccion, gdl_restriccion, punto_desplazamiento, valor_desplazamiento, coords_metam_hexaedros_exacto, conectiv_metam_hexaedros_exacto, porcentajes_barras_dentro_hexaedros, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_totales_hexaedro_local, gdl_contorno_solido_completo_no_fronteras, invariante_def_memoria_metam, iter_desplazamiento, radio)
    
    elseif (gradiente == true) && (length(hexaedros_analizados_metam)==N_hexaedros_gruesos)
        # if iter_desplazamiento==2
        #     hiperparametros_gradiente[10]=20
        # else
        #     hiperparametros_gradiente[10]=20
        # end
        hexaedros_analizados_metam, inc_u_solido, inc_f_solido, u_solido_acumulado,  f_solido_acumulado, fuerzas_hexaedros_acumulado, energias_hexaedros_acumulado, inc_energia_solido, inc_energias_solidos_gruesos, energias_solidos_gruesos, inc_deform_equiv_solidos_gruesos, inc_u_metam_hexaedros, inc_f_metam_hexaedros, f_acumulada_elementos_memoria, Youngs_tangentes_hexaedros_finos, Youngs_tangentes_hexaedros_gruesos, Youngs_tangentes_memoria_metam, epsilon_voigt_acumulado, energias_acumuladas_elementos_memoria, eps_x_ptos_int_acumulado_memoria, eps_xy_ptos_int_acumulado_memoria, eps_xz_ptos_int_acumulado_memoria , gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_no_fijos_total, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_no_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_fijos_hexaedro_local, gdl_restring_no_fijos_hexaedro_local, gdl_contorno_solido_completo_no_fronteras, energia_elemento_solido_grueso_metam, inc_energia_elemento_solido_grueso_metam, invariante_def_memoria_metam, inc_valor_restricciones_no_fijas_total, K_step, gradiente_fuerzas_nodos, datos_contribucion_barras_gdl = resolucion_iterativa_acoplamiento_zonas_mallado_fino_grueso_gradiente_fuerzas(hiperparametros_gradiente, X,Y,Z, N_3D_x_grueso, N_3D_y_grueso, N_3D_z_grueso, hexaedros_analizados_metam, hexaedros_mallado_fino_en_grueso, nodos_metam_en_hexaedros_finos, u_solido_acumulado, f_solido_acumulado, fuerzas_hexaedros_acumulado, energias_hexaedros_acumulado, epsilon_voigt_acumulado, f_acumulada_elementos_memoria, energias_acumuladas_elementos_memoria, eps_x_ptos_int_acumulado_memoria, eps_xy_ptos_int_acumulado_memoria, eps_xz_ptos_int_acumulado_memoria, coords_metam, conectividad_metam ,coords_3D_fino, conectividad_3D_fino, conectividad_3D_grueso, volumenes_hexaedros_grueso, Youngs_tangentes_hexaedros_finos, Youngs_tangentes_hexaedros_gruesos, Youngs_tangentes_memoria_metam, tensiones_pandeo_elementos, punto_restriccion, gdl_restriccion, punto_desplazamiento, valor_desplazamiento, coords_metam_hexaedros_exacto, conectiv_metam_hexaedros_exacto, porcentajes_barras_dentro_hexaedros, factor_barras_caras_solido, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_no_fijos_total, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_no_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_fijos_hexaedro_local, gdl_restring_no_fijos_hexaedro_local, gdl_contorno_solido_completo_no_fronteras, invariante_def_memoria_metam, inc_valor_restricciones_no_fijas_total, K_step, gradiente_fuerzas_nodos, datos_contribucion_barras_gdl, iter_desplazamiento, radio)
         #K_step=[]                                                                                                                                                                                                                                                                                                                                                           
    elseif (gradiente == true) && (length(hexaedros_analizados_metam)<N_hexaedros_gruesos)
        hexaedros_analizados_metam, inc_u_solido, inc_f_solido, u_solido_acumulado,  f_solido_acumulado, fuerzas_hexaedros_acumulado, energias_hexaedros_acumulado, inc_energia_solido, inc_energias_solidos_gruesos, energias_solidos_gruesos, inc_deform_equiv_solidos_gruesos, inc_u_metam_hexaedros, inc_f_metam_hexaedros, f_acumulada_elementos_memoria, Youngs_tangentes_hexaedros_finos, Youngs_tangentes_hexaedros_gruesos, Youngs_tangentes_memoria_metam, epsilon_voigt_acumulado, energias_acumuladas_elementos_memoria, eps_x_ptos_int_acumulado_memoria, eps_xy_ptos_int_acumulado_memoria, eps_xz_ptos_int_acumulado_memoria , gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_totales_hexaedro_local, gdl_contorno_solido_completo_no_fronteras, energia_elemento_solido_grueso_metam, inc_energia_elemento_solido_grueso_metam, invariante_def_memoria_metam = resolucion_iterativa_acoplamiento_zonas_mallado_fino_grueso(hiperparametros, X,Y,Z, N_3D_x_grueso, N_3D_y_grueso, N_3D_z_grueso, hexaedros_analizados_metam, hexaedros_mallado_fino_en_grueso, nodos_metam_en_hexaedros_finos, u_solido_acumulado, f_solido_acumulado, fuerzas_hexaedros_acumulado, energias_hexaedros_acumulado, epsilon_voigt_acumulado, f_acumulada_elementos_memoria, energias_acumuladas_elementos_memoria, eps_x_ptos_int_acumulado_memoria, eps_xy_ptos_int_acumulado_memoria, eps_xz_ptos_int_acumulado_memoria, coords_metam, conectividad_metam ,coords_3D_fino, conectividad_3D_fino, conectividad_3D_grueso, volumenes_hexaedros_grueso, Youngs_tangentes_hexaedros_finos, Youngs_tangentes_hexaedros_gruesos, Youngs_tangentes_memoria_metam, tensiones_pandeo_elementos, punto_restriccion, gdl_restriccion, punto_desplazamiento, valor_desplazamiento, coords_metam_hexaedros_exacto, conectiv_metam_hexaedros_exacto, porcentajes_barras_dentro_hexaedros, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_totales_hexaedro_local, gdl_contorno_solido_completo_no_fronteras, invariante_def_memoria_metam, iter_desplazamiento, radio)
    else
    end

    ### ENERGIA TOTAL DEL SOLIDO
    energia_acumulada_solido = energia_acumulada_solido + inc_u_solido'*f_solido_acumulado - inc_energia_solido
    println()
    println(energia_acumulada_solido)
    println("ENERGIA ACUMULADA SOLIDO= ",sum(energias_hexaedros_acumulado));println();
    #u_solido_acumulado=u_solido_acumulado+inc_u_solido
    #f_solido_acumulado=f_solido_acumulado+inc_f_solido
    #readline()

    ### REPRESENTAR EN PARAVIEW EL SOLIDO
    plot_paraview_solido_sucesion_hexaedros(coords_3D_fino, conectividad_3D_fino, u_solido_acumulado, f_solido_acumulado, iter_desplazamiento)


    ### REPRESENTAR EN PARAVIEW LAS ZONA DE METAMATERIAL DE DENTRO DE LOS hexaedros DAÑADOS
    # obtener los desplazamientos y fuerzas totales a partir de los incrementos
    for elemento in hexaedros_analizados_metam

        # si todavia no se habia registrado ningun desplazamiento
        if haskey(u_metam_hexaedros_acumulado, elemento)==false
            u_metam_hexaedros_acumulado[elemento]=inc_u_metam_hexaedros[elemento]
            f_metam_hexaedros_acumulado[elemento]=inc_f_metam_hexaedros[elemento]
        else
            u_metam_hexaedros_acumulado[elemento]=u_metam_hexaedros_acumulado[elemento]+inc_u_metam_hexaedros[elemento]
            f_metam_hexaedros_acumulado[elemento]=f_metam_hexaedros_acumulado[elemento]+inc_f_metam_hexaedros[elemento]             
        end

        desplazamientos, giros, fuerzas, momentos = generar_diccionarios(u_metam_hexaedros_acumulado[elemento], f_metam_hexaedros_acumulado[elemento], coords_metam_hexaedros_exacto[elemento])

        plot_paraview_estructura_metamaterial_multiescala_sucesion_hexaedros(coords_metam_hexaedros_exacto[elemento], conectiv_metam_hexaedros_exacto[elemento], desplazamientos,giros,fuerzas,momentos,string(elemento),iter_desplazamiento)
        #readline() coords_metam_hexaedros_exacto[elemento],conectiv_metam_hexaedros_exacto[elemento]
 
    end       

    

    ### ALMACENAR VARIABLES PARA PLOTEAR FUERZA VS DESPLAZAMIENTO EN X, Y y Z
    gdl_restring_globales_solido, inc_valor_restricciones_solido=condic_contorno_solido(punto_restriccion, gdl_restriccion, punto_desplazamiento, valor_desplazamiento, coords_3D_fino)
    indices=findall(!in(0),inc_valor_restricciones_solido)
    gdl_desplaz_nulo=gdl_restring_globales_solido[indices]
    
    fuerza_reaccion_x=0; fuerza_reaccion_y=0; fuerza_reaccion_z=0; desplaz=0
    for gdl in gdl_desplaz_nulo

        # Si es un gdl en X
        if gdl%3==1
            fuerza_reaccion_x=fuerza_reaccion_x+f_solido_acumulado[gdl]
        # Si es un gdl en Y
        elseif gdl%3==2
            fuerza_reaccion_y=fuerza_reaccion_y+f_solido_acumulado[gdl]
        # Si es un gdl en Z
        else
            fuerza_reaccion_z=fuerza_reaccion_z+f_solido_acumulado[gdl]
        end
    end

    desplaz=maximum(inc_valor_restricciones_solido)*iter_desplazamiento

    
    push!(fuerza_reaccion_x_iters,abs(fuerza_reaccion_x))
    push!(fuerza_reaccion_y_iters,abs(fuerza_reaccion_y))
    push!(fuerza_reaccion_z_iters,abs(fuerza_reaccion_z))
    push!(desplaz_iters, desplaz)

    # PLOT DE FUERZA EN X
    plot(desplaz_iters, fuerza_reaccion_x_iters)
    xlabel!("desplazamiento")
    display(ylabel!("fuerza reaccion X"))

    # PLOT DE FUERZA EN Y
    plot(desplaz_iters, fuerza_reaccion_y_iters)
    xlabel!("desplazamiento")
    display(ylabel!("fuerza reaccion Y"))

    # PLOT DE FUERZA EN Z
    plot(desplaz_iters, fuerza_reaccion_z_iters)
    xlabel!("desplazamiento")
    display(ylabel!("fuerza reaccion Z"))


    ### AÑADIR INFORMACION A LOS FICHEROS DE DATOS
    fichero_ensayo=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_energ_element_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_2=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_inc_energ_element_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_3=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_energ_total_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_4=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_inc_energ_total_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_6=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_energ_zona_metam_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_5=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_fuerza_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_7=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_inc_energ_zona_metam_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_8=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_energ_metam_total_$(lpad(numero_fichero,2,"0")).txt","a")
    fichero_ensayo_9=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_inc_energ_metam_total_$(lpad(numero_fichero,2,"0")).txt","a")
    fichero_ensayo_10=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_curva_tension_deformacion_hexaedros_gruesos_$(lpad(numero_fichero,2,"0")).txt","a")
    fichero_ensayo_11=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_multiescala_deformacion_barra_max_hexaedros_gruesos_$(lpad(numero_fichero,2,"0")).txt","a")  

    for i in 1:N_hexaedros_gruesos
        println(fichero_ensayo, energias_solidos_gruesos[i])
        println(fichero_ensayo_2, inc_energias_solidos_gruesos[i])
        println(fichero_ensayo_6, energia_elemento_solido_grueso_metam[i])
        println(fichero_ensayo_7, inc_energia_elemento_solido_grueso_metam[i])       
        println(fichero_ensayo_11, maximum(values(invariante_def_memoria_metam[i])))
    end
    println(fichero_ensayo_3, sum(energias_solidos_gruesos))
    println(fichero_ensayo_4, sum(inc_energias_solidos_gruesos))
    println(fichero_ensayo_8, sum(values(energia_elemento_solido_grueso_metam)))
    println(fichero_ensayo_9, sum(values(inc_energia_elemento_solido_grueso_metam)))

    println(fichero_ensayo_5, fuerza_reaccion_x_iters[end]," ",desplaz_iters[end])




    ### VER COMO VARIA LA CURVA SIGMA-EPSILON EN CADA TETRAEDRO (A NIVEL LOCAL)
    deform_solidos_gruesos=deform_solidos_gruesos+inc_deform_equiv_solidos_gruesos
    inc_tension_solidos_gruesos=Youngs_tangentes_hexaedros_gruesos.*inc_deform_equiv_solidos_gruesos
    tension_solidos_gruesos=tension_solidos_gruesos+inc_tension_solidos_gruesos

    for i in 1:N_hexaedros_gruesos
        println(fichero_ensayo_10, deform_solidos_gruesos[i]," ", tension_solidos_gruesos[i]," ", Youngs_tangentes_hexaedros_gruesos[i])
    end


    close(fichero_ensayo)
    close(fichero_ensayo_2)
    close(fichero_ensayo_3)
    close(fichero_ensayo_4)
    close(fichero_ensayo_5)
    close(fichero_ensayo_6)
    close(fichero_ensayo_7)
    close(fichero_ensayo_8)
    close(fichero_ensayo_9)
    close(fichero_ensayo_10)
    close(fichero_ensayo_11)

end


############

# for i in 1:iter_desplaz
#     println(fichero_ensayo_5, fuerza_reaccion_z_iters[i]," ",desplaz_iters[i] )
# end

    ### COMPROBAR SUMA DE FUERZAS 
    # indice=findall(x->(x<X+0.0001)&&(x>X-0.001),coords_3D_fino[:,1])
    # gdls=3*indice.-2
    # F=sum(f_solido_acumulado[gdls])
    # Young_equiv=(F/(Y*Z))/(u_solido_acumulado[gdls[1]]/X)
    # println();println("Fuerza= ",F,"    Young_equiv= ",Young_equiv)
    # push!(fuerza_plot,F);push!(valor_desplaz_plot,valor_desplazamiento[1]*iter_desplazamiento)
