include(raw"celdillas_metamaterial.jl"); using .celdillas_metamaterial 
include(raw"funciones_MEF_metamaterial.jl"); using .funciones_MEF_metamaterial
include(raw"representar_paraview.jl"); using .representar_paraview
include(raw"representar_paraview_sucesion.jl"); using .representar_paraview_sucesion
include(raw"curva_sigma_epsilon_material.jl"); using .curva_sigma_epsilon_material
include(raw"mallado_solido_tridimensional_hexaedros.jl"); using .mallado_solido_tridimensional_hexaedros


using Plots
using LinearAlgebra


comparar_metamaterial_multiescala="off" #"on" #para poder comparar con el multiescala, antes se debe haber ejecutado el main multiecala para identificar las distintas regiones en que se divide el solido

### GENERAR LOS FICHEROS DE RESULTADOS (si comparar_metamaterial_multiescala="on")
numero_fichero="prueba2023" #"TRACC_25" # "FLEX_25"

#generar_tabla_tension_deformacion_Osgood()
fichero_ensayo=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_metamaterial_energ_element_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_2=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_metamaterial_inc_energ_element_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_3=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_metamaterial_energ_total_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_4=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_metamaterial_inc_energ_total_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo
fichero_ensayo_5=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_metamaterial_fuerza_$(lpad(numero_fichero,2,"0")).txt","w") #crear el archivo para recoger datos del ensayo


### DIMENSIONES PROBETA RECTANGULAR Y CELDILLAS METAMATERIAL 
p=8
Nx=p#10 #numero de celdillas en x
Ny=p#6 #numero de celdillas en y
Nz=p#2   #numero de celdillas en z

X=10*Nx    #dimension x probeta
Y=10*Ny #dimension y probeta
Z=10*Nz #dimension z probeta

radio=5  #radio de las barras de metamaterial


### HIPERPARAMETROS
iter_desplaz=10   #Numero de iteraciones hasta alcanzar el valor del desplazamiento final


### CONDICIONES DE CONTORNO DE LA PROBETA (Puntos de la geometria donde restringimos o imponemos un desplazamiento): 
    # Para cada restriccion hay que indicar el punto donde se aplica y los GDL restringidos
        # coordenada --> restringir el movimiento de los puntos con esa coordenada en el eje correspondiente; 
        # "all" --> restringir el movimiento a todos los puntos con esas coordenadas
        # GDL restringidos en coordenadas globales: (ux,uy,uz,thetax,thetay,thetaz): 0 implica que no se restringe y 1 que si se restringe
    # Las nuevas restricciones se añaden con su nombre con [punto_restriccion;punto]+[gdl_restriccion;gdl]
    # (Igual para los desplazamientos)


    ########################### DISTINTAS RESTRICCIONES

    ### RESTRICCIONES PARA TRACCION, TORSION, FLEXION PURA
    # restriccion 1
    punto_restriccion=[0 "all" "all"]#[0 "all" 0]# Coordenadas del nodo [x,y,z]
    gdl_restriccion=[1 1 1 0 0 0]#[1 1 1 1 1 1]# grados de libertad: [ux, uy, uz, thetax, thetay, thetaz]


    ### RESTRICCIONES PARA TRACCION SIN EFECTOS DE LA BASE
    # # restriccion 1
    # punto_restriccion=[0 "all" "all"]
    # gdl_restriccion=[1 0 0 0 0 0]#[1 1 1]#
    # # restriccion 2
    # punto=[0 Y/2 "all"]#["all" "all" Z]#[X "all" "all"]#
    # gdl=[0 1 0 0 0 0]#[0 0 1]#[1 1 1]#
    # punto_restriccion=vcat(punto_restriccion,punto)
    # gdl_restriccion=vcat(gdl_restriccion,gdl)
    # # restriccion 3
    # punto=[0 "all" Z/2]#["all" Y "all"]#
    # gdl=[0 0 1 0 0 0]#[0 1 0]#
    # punto_restriccion=vcat(punto_restriccion,punto)
    # gdl_restriccion=vcat(gdl_restriccion,gdl)


    ### RESTRICCIONES PARA FLEXION A 3 PUNTOS
    # # restriccion 1
    # punto_restriccion=[0 "all" "all"]#[0 "all" 0]#
    # gdl_restriccion=[1 1 1 0 0 0]#[1 1 1 1 1 1]#
    # # restriccion 2
    # punto= [X "all" "all"] #["all" Y "all"]#[0 Y/2 Z/2]#
    # gdl=[1 1 1 0 0 0]#[1 1 1 0 0 0]#
    # punto_restriccion=vcat(punto_restriccion,punto)
    # gdl_restriccion=vcat(gdl_restriccion,gdl)


    ########################### DISTINTOS TIPOS DE CARGA

    ### TRACCION 
    #desplazamiento 1
    punto_desplazamiento=[X "all" "all"]#
    valor_desplazamiento=[1.5e-0 0  0   0   0  0 ]#[1.5e-6 0 0 0  0  0]#


    ### CARGADO EN UN EXTREMO Y EMPOTRADO EN OTRO
    # # desplazamiento 1
    # punto_desplazamiento=[X "all" "all"]#
    # valor_desplazamiento=[ 0  0 1.5e-0  0   0  0 ]


    ### FLEXION PURA (UN GIRO) (con Nz par)
    # # desplazamiento 1
    # punto_desplazamiento=[X "all" Z]#[X "all" "all"]#
    # valor_desplazamiento=[1.5e-6 0  0   0   0  0 ]#[0 0 0  0 1e-6 0]#
    # # desplazamiento 2
    # punto=[X "all" 0]
    # valor=[-1.5e-6  0  0  0   0  0 ]#
    # punto_desplazamiento=vcat(punto_desplazamiento,punto)
    # valor_desplazamiento=vcat(valor_desplazamiento,valor)
    
    # ## Para el ensayo de flexion pura, siempre que Nz sea par
    # for i in 1:Nz/2-1
    #     global punto, valor, punto_desplazamiento, valor_desplazamiento
    #     punto=[X "all" Z/Nz*i+Z/2]
    #     valor=[1.5e-6*2*i/Nz  0  0  0   0  0 ]#
    #     punto_desplazamiento=vcat(punto_desplazamiento,punto)
    #     valor_desplazamiento=vcat(valor_desplazamiento,valor)
    #     punto=[X "all" -Z/Nz*i+Z/2]
    #     valor=[-1.5e-6*2*i/Nz  0  0  0   0  0 ]#
    #     punto_desplazamiento=vcat(punto_desplazamiento,punto)
    #     valor_desplazamiento=vcat(valor_desplazamiento,valor)
    # end


    ### TORSION
    # # desplazamiento 1
    # punto_desplazamiento=[X "all" Z]#[X "all" "all"]#
    # valor_desplazamiento=[0 -1.5e-6  0   0   0  0 ]#[0 0 0  0 1e-6 0]#
    # # desplazamiento 2
    # punto=[X "all" 0]
    # valor=[0 1.5e-6  0  0   0  0 ]#
    # punto_desplazamiento=vcat(punto_desplazamiento,punto)
    # valor_desplazamiento=vcat(valor_desplazamiento,valor)
    # # desplazamiento 3
    # punto=[X 0 "all"]
    # valor=[ 0 0 -1.5e-6 0  0  0 ]#
    # punto_desplazamiento=vcat(punto_desplazamiento,punto)
    # valor_desplazamiento=vcat(valor_desplazamiento,valor)
    # # desplazamiento 4
    # punto=[X Y "all"]
    # valor=[ 0 0 1.5e-6 0  0  0 ]#
    # punto_desplazamiento=vcat(punto_desplazamiento,punto)
    # valor_desplazamiento=vcat(valor_desplazamiento,valor)

    # # para el ensayo de torsion, siempre que Nz y Ny sean par
    # for i in 1:Nz/2-1
    #     global punto, valor, punto_desplazamiento, valor_desplazamiento
    #     punto=[X "all" Z/Nz*i+Z/2]
    #     valor=[ 0 -1.5e-6*2*i/Nz 0 0 0 0]#
    #     punto_desplazamiento=vcat(punto_desplazamiento,punto)
    #     valor_desplazamiento=vcat(valor_desplazamiento,valor)
    #     punto=[X "all" -Z/Nz*i+Z/2]
    #     valor=[ 0 1.5e-6*2*i/Nz 0 0 0 0]#
    #     punto_desplazamiento=vcat(punto_desplazamiento,punto)
    #     valor_desplazamiento=vcat(valor_desplazamiento,valor)
    # end
    # for i in 1:Ny/2-1
    #     global punto, valor, punto_desplazamiento, valor_desplazamiento
    #     punto=[X Y/Ny*i+Y/2 "all" ]
    #     valor=[ 0 0 1.5e-6*2*i/Ny 0 0 0]#
    #     punto_desplazamiento=vcat(punto_desplazamiento,punto)
    #     valor_desplazamiento=vcat(valor_desplazamiento,valor)
    #     punto=[X -Y/Ny*i+Y/2 "all"]
    #     valor=[0  0 -1.5e-6*2*i/Ny 0 0 0]#
    #     punto_desplazamiento=vcat(punto_desplazamiento,punto)
    #     valor_desplazamiento=vcat(valor_desplazamiento,valor)
    # end


    ### FLEXION A 3 PUNTOS (con Nx par)
    # # desplazamiento 1
    # punto_desplazamiento=[X/2 "all" "all"]#[X "all" "all"]#
    # valor_desplazamiento=[0 0  -1.5e-6   0   0  0 ]#[0 0 0  0 1e-6 0]#


    
# IMPONER EL VALOR DEL DESPLAZAMIENTO COMO INCREMENTOS
valor_desplazamiento=valor_desplazamiento/iter_desplaz


#########################################################################################################
# CREAR LAS celdillas FCC (esto depende del tipo de celdillas)
#########################################################################################################
coords,conectividad=celdilla_FCC(X,Y,Z,Nx,Ny,Nz)#_nodos_itermedios
N_nodos=length(coords[:,1]); N_elementos=length(conectividad[:,3])
tensiones_pandeo=tensiones_criticas_pandeo(coords, conectividad, Young_material(), radio)

println("nº nodos= ", N_nodos)
println("nº elementos metamaterial= ", N_elementos)


######################################### DIFERENTES CARGAS APLICADAS ############################################
desplaz_x=[];desplaz_y=[];desplaz_z=[];fuerza_x_total_iter=[]


# inicializar los modulos de Young de cada elemento
dict_ceros=Dict(); pandeo_cero=Dict(); eps_x_ptos_int_acumulado=Dict(); eps_xy_ptos_int_acumulado=Dict(); eps_xz_ptos_int_acumulado=Dict(); f_acumulada_elementos=Dict() ;energias_acumuladas_elementos=Dict(); invariante_def_metam=Dict()
for i in 1:length(conectividad[:,3])
    element_metam=conectividad[i,3]
    dict_ceros[element_metam]=0 #diccionario generico de ceros
    eps_x_ptos_int_acumulado[element_metam]=[0,0]      
    eps_xy_ptos_int_acumulado[element_metam]=[0,0]
    eps_xz_ptos_int_acumulado[element_metam]=[0,0]
    energias_acumuladas_elementos[element_metam]=0 
    invariante_def_metam[element_metam]=0 
    f_acumulada_elementos[element_metam] =zeros(12)      
end
Youngs_tangentes_metam=curva_tension_deformacion_tabla_Osgood(dict_ceros, "deform_equiv",dict_ceros)  #si no se tienen registros de modulos de young anteriores, se empieza por el valor base del material
young_material=Young_material()

desplazamientos=Dict(); giros=Dict(); fuerzas=Dict(); momentos=Dict()
# inicializar desplazamientos y fuerzas
for nodo in coords[:,4]
    desplazamientos[nodo] = [0, 0, 0]
    giros[nodo] = [0, 0, 0]
    fuerzas[nodo] = [0, 0, 0]
    momentos[nodo] = [0, 0, 0]
end
f_acumulado=zeros(6*N_nodos); u_acumulado=zeros(6*N_nodos); energia_total_zona=0; energia_lineal_2=0


inc_u=0; fuerza_reaccion_x_iters=[]; fuerza_reaccion_y_iters=[]; fuerza_reaccion_z_iters=[]; fuerza_reaccion_x_porcentaje_iters=[]; fuerza_reaccion_y_porcentaje_iters=[]; fuerza_reaccion_z_porcentaje_iters=[]; desplaz_iters=[]

for iter in 1:iter_desplaz

    global inc_u, f_acumulado, u_acumulado, inc_energia_elementos, deform_equiv,energias_unitarias_antig,energia_total_zona,punto_desplazamiento,punto_restriccion,valor_desplazamiento,desplaz_x,desplaz_y,desplaz_z, fuerza_x_total_iter, coords,conectividad, gdl_restring_globales, gdl_desplaz_globales,u, desplazamientos,fuerzas, Youngs_tangentes_metam, eps_x_ptos_int_acumulado, eps_xy_ptos_int_acumulado, eps_xz_ptos_int_acumulado, f_acumulada_elementos, energias_acumuladas_elementos, fuerza_reaccion_x_iters, fuerza_reaccion_y_iters, fuerza_reaccion_z_iters, fuerza_reaccion_x_porcentaje_iters, fuerza_reaccion_y_porcentaje_iters, fuerza_reaccion_z_porcentaje_iters, desplaz_iters, invariante_def_metam, fichero_ensayo, fichero_ensayo_2, fichero_ensayo_3, fichero_ensayo_4, fichero_ensayo_5
    println(); println("iter", iter)


    # IMPONER CONDICIONES DE CONTORNO
    gdl_restring_globales, valor_restricciones_globales=condic_contorno(punto_restriccion,gdl_restriccion,punto_desplazamiento,valor_desplazamiento,coords)

    # RESOLVER EL PROBLEMA MEDIANTE INCREMENTOS
    inc_u, inc_f, f_acumulada_elementos, inc_energia, inc_energia_elementos, energias_acumuladas_elementos, Youngs_tangentes_metam, eps_x_ptos_int_acumulado, eps_xy_ptos_int_acumulado, eps_xz_ptos_int_acumulado, invariante_def_metam = resolucion_iterativa_metamaterial(coords, conectividad, gdl_restring_globales, valor_restricciones_globales, Youngs_tangentes_metam, tensiones_pandeo, eps_x_ptos_int_acumulado, eps_xy_ptos_int_acumulado, eps_xz_ptos_int_acumulado, energias_acumuladas_elementos, f_acumulada_elementos, invariante_def_metam, radio)
    
    
    # OBTENER LA ENERGIA DEL CONJUNTO A PARTIR DE LAS SUMAS DE LOS DISTINTOS INCREMENTOS LINEALES
    energia_total_zona=energia_total_zona + inc_energia + f_acumulado'*inc_u 
    println(sum(values(energias_acumuladas_elementos)))
    println("energia total de la zona= ",energia_total_zona)

    
    # OBTENER LOS VALORES DE LAS VARIABLES ACTUALES A PARTIR DE LOS INCREMENTOS
    inc_desplazamientos,inc_giros,inc_fuerzas,inc_momentos=generar_diccionarios(inc_u, inc_f, coords)

    for nodo in coords[:,4]
        desplazamientos[nodo] = desplazamientos[nodo] + inc_desplazamientos[nodo]
        giros[nodo] = giros[nodo] + inc_giros[nodo]
        fuerzas[nodo] = fuerzas[nodo] + inc_fuerzas[nodo]
        momentos[nodo] = momentos[nodo] + inc_momentos[nodo]
    end
    f_acumulado=f_acumulado + inc_f
    u_acumulado=u_acumulado + inc_u


    # REPRESENTAR EN PARAVIEW
    plot_paraview_estructura_metamaterial_sucesion(coords,conectividad,desplazamientos,giros,fuerzas,momentos,iter)


    ### ALMACENAR VARIABLES PARA PLOTEAR FUERZA VS DESPLAZAMIENTO EN X, Y y Z
    indices=findall(!in(0),valor_restricciones_globales)
    gdl_desplaz_nulo=gdl_restring_globales[indices]
    
    fuerza_reaccion_x=0; fuerza_reaccion_y=0; fuerza_reaccion_z=0; desplaz=0
    for gdl in gdl_desplaz_nulo

        # Si es un gdl en X
        if gdl%3==1
            fuerza_reaccion_x=fuerza_reaccion_x+f_acumulado[gdl]
        # Si es un gdl en Y
        elseif gdl%3==2
            fuerza_reaccion_y=fuerza_reaccion_y+f_acumulado[gdl]
        # Si es un gdl en Z
        elseif gdl%3==0
            fuerza_reaccion_z=fuerza_reaccion_z+f_acumulado[gdl]
        else
        end
    end
    println(fuerza_reaccion_x)
    desplaz=maximum(valor_restricciones_globales)*iter

    
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



    ### PARA LA CURVA TENSION-DEFORMACION
    indice=findall(x->(x<X+0.0001)&&(x>X-0.001),coords[:,1])
    push!(desplaz_x, desplazamientos[coords[indice[1],4]][1])

    nodo_mitad_exterior_y=Int64((Ny+1)*floor(Nx/2)+1+(Nx+1)*(Ny+1)*Nz)  #(Nx+1)*(Ny+1)*floor((Nz+1)/2)+(Ny+1)*floor((Nx+1)/2)+1) #punto de la zona media en la cara exterior para calcular su desplazamiento en Y para obtener el poisson
    push!(desplaz_y,desplazamientos[nodo_mitad_exterior_y][2])

    nodo_mitad_exterior_z=Int64((Ny+1)*(floor(Nx/2)+1))  #Int64((Nx+1)*floor((Ny+1)/2)+floor((Nx+1)/2)+1) #punto de la zona media en la cara exterior para calcular su desplazamiento en Y para obtener el poisson
    push!(desplaz_z,desplazamientos[nodo_mitad_exterior_z][3])    


    ### OBTENER LA FUERZA APLICADA EN LOS PUNTOS CON DESPLAZAMIENTO IMPUESTO
    fuerza_x=[0.0]; fuerza_y=[0.0]; fuerza_z=[0.0]; momento_x=[0.0]; momento_y=[0.0]; momento_z=[0.0];

    global coords_vertices_hexaedro, volumen_hexaedro
    coords_vertices_hexaedro=[0 0 0; 0 Y 0; X 0 0; X Y 0; 0 0 Z; 0 Y Z; X 0 Z; X Y Z] 
    volumen_hexaedro=X*Y*Z

    indices=findall(!in(0),valor_restricciones_globales)
    for gdl in gdl_restring_globales[indices]
        if gdl%6==1
            push!(fuerza_x,f_acumulado[gdl])
        elseif gdl%6==2
            push!(fuerza_y,f_acumulado[gdl])
        elseif gdl%6==3
            push!(fuerza_z,f_acumulado[gdl])
        elseif gdl%6==4
            push!(momento_x,f_acumulado[gdl])
        elseif gdl%6==5
            push!(momento_y,f_acumulado[gdl])
        elseif gdl%6==0
            push!(momento_z,f_acumulado[gdl])
        else 
        end
    end

    fuerza_x_total=sum(fuerza_x);fuerza_y_total=sum(fuerza_y);fuerza_z_total=sum(fuerza_z)
    momento_x_total=sum(momento_x);momento_y_total=sum(momento_y);momento_z_total=sum(momento_z)
    println()
    println("f_aplicada_x= ", fuerza_x_total,"  f_aplicada_y= ", fuerza_y_total,"  f_aplicada_z= ", fuerza_z_total)
    println("M_aplicado_x= ", momento_x_total,"  M_aplicado_y= ", momento_y_total,"  M_aplicado_z= ", momento_z_total)
    #println("Young equiv= ", fuerza_x_total/(Y*Z)/(desplazamientos[nodos[1]][1]/X))
    push!(fuerza_x_total_iter, fuerza_x_total)


    ### COMPROBAR QUE EL SUMATORIO DE FUERZAS Y MOMENTOS ES IGUAL A CERO    
    momentos_fuerzas=Dict() # momento generado por las fuerzas en los nodos respecto al origen
    suma_fuerzas=[0,0,0]; suma_momentos=[0,0,0]; suma_momentos_fuerzas=[0,0,0]
    for nodo in coords[:,4]
        momentos_fuerzas[nodo]=cross(coords[nodo,1:3],fuerzas[nodo])
        suma_fuerzas=suma_fuerzas+fuerzas[nodo]
        suma_momentos=suma_momentos+momentos[nodo]
        suma_momentos_fuerzas=suma_momentos_fuerzas+momentos_fuerzas[nodo]
    end
    println()
    println("sum_F= ",suma_fuerzas[1]," ",suma_fuerzas[2]," ",suma_fuerzas[3])
    println("sum_M= ",suma_momentos[1]," ",suma_momentos[2]," ",suma_momentos[3])
    println("sum_M_F= ",suma_momentos_fuerzas[1]," ",suma_momentos_fuerzas[2]," ",suma_momentos_fuerzas[3])
    println()







    # coords_vertices_hexaedro=[0 0 0; 0 Y 0; X 0 0; X Y 0; 0 0 Z; 0 Y Z; X 0 Z; X Y Z] 
    # volumen_hexaedro=X*Y*Z
    porcentaje_barra=zeros(N_elementos); inc_energia_solo_dentro_barras=zeros(N_elementos); inc_energia_solo_dentro_barras_no_porcentajes=zeros(N_elementos); energia_solo_dentro_barras=zeros(N_elementos);
    for barra in 1:N_elementos
        nodo_1=conectividad[barra,1]; coords_nodo_1=coords[nodo_1,1:3]
        nodo_2=conectividad[barra,2]; coords_nodo_2=coords[nodo_2,1:3]
        valor_1, caras_pertenece_nodo_1=is_in_hexaedro(coords_nodo_1, coords_vertices_hexaedro, volumen_hexaedro)
        valor_2, caras_pertenece_nodo_2=is_in_hexaedro(coords_nodo_2, coords_vertices_hexaedro, volumen_hexaedro)

        # Hay que determinar cuantas caras tienen en comun y cuales son
        caras_comunes=[] 
        for cara in caras_pertenece_nodo_1 
            if findall(in(cara), caras_pertenece_nodo_2) != []
                push!(caras_comunes, cara)
            else
            end
        end
        N_caras_comunes=length(caras_comunes)

        if N_caras_comunes==0
            porcentaje_barra[barra]=1
        elseif N_caras_comunes==1
            porcentaje_barra[barra]=0.5
        else
            porcentaje_barra[barra]=0.25
        end

        energia_solo_dentro_barras[barra]=energias_acumuladas_elementos[barra]*porcentaje_barra[barra]
        inc_energia_solo_dentro_barras[barra]=inc_energia_elementos[barra]*porcentaje_barra[barra]
        inc_energia_solo_dentro_barras_no_porcentajes[barra]=inc_energia_elementos[barra]
        #println(nodo_1," ",nodo_2," ",porcentaje_barra[barra], " ",inc_energia_solo_dentro_barras[barra])
    end

    inc_energia_solo_dentro_total=sum(inc_energia_solo_dentro_barras)
    #Young_equiv=inc_energia_solo_dentro_total/(1/2*volumen_hexaedro*(1.5e-6/X)^2)
    println("inc_energia= ", inc_energia_solo_dentro_total)
    #println("a ",Young_equiv)

    energia_solo_dentro_total=sum(energia_solo_dentro_barras)
    #Young_equiv=energia_solo_dentro_total/(1/2*volumen_hexaedro*(1.5e-6/X)^2)
    println("energia= ", energia_solo_dentro_total)
    #println("a ",Young_equiv)


    ### ALMACENAR VARIABLES PARA PLOTEAR FUERZA VS DESPLAZAMIENTO EN X, Y y Z
    indices=findall(!in(0),valor_restricciones_globales)
    gdl_desplaz_nulo=gdl_restring_globales[indices]
    
    fuerza_reaccion_x_porcentaje=0; fuerza_reaccion_y_porcentaje=0; fuerza_reaccion_z_porcentaje=0; desplaz=0
    for gdl in gdl_desplaz_nulo
        nodo=ceil(gdl/6)
        indices=findall(in(nodo),conectividad[:,1:2])

        for indice in indices
            if indice[2]==1
                k=0
            else
                k=6
            end

            # Si es un gdl en X
            if gdl%3==1
                fuerza_reaccion_x_porcentaje=fuerza_reaccion_x_porcentaje+f_acumulada_elementos[indice[1]][1+k]*porcentaje_barra[indice[1]]
            # Si es un gdl en Y
            elseif gdl%3==2
                fuerza_reaccion_y_porcentaje=fuerza_reaccion_y_porcentaje+f_acumulada_elementos[indice[1]][2+k]*porcentaje_barra[indice[1]]
            # Si es un gdl en Z
            elseif gdl%3==0
                fuerza_reaccion_z_porcentaje=fuerza_reaccion_z_porcentaje+f_acumulada_elementos[indice[1]][3+k]*porcentaje_barra[indice[1]]
            else
            end
        end
    end
    println(fuerza_reaccion_x)
    desplaz=maximum(valor_restricciones_globales)*iter

    
    push!(fuerza_reaccion_x_porcentaje_iters,abs(fuerza_reaccion_x_porcentaje))
    push!(fuerza_reaccion_y_porcentaje_iters,abs(fuerza_reaccion_y_porcentaje))
    push!(fuerza_reaccion_z_porcentaje_iters,abs(fuerza_reaccion_z_porcentaje))
    #push!(desplaz_iters, desplaz)

    # PLOT DE FUERZA EN X
    plot(desplaz_iters, fuerza_reaccion_x_porcentaje_iters)
    xlabel!("desplazamiento")
    display(ylabel!("fuerza reaccion X"))

    # PLOT DE FUERZA EN Y
    plot(desplaz_iters, fuerza_reaccion_y_porcentaje_iters)
    xlabel!("desplazamiento")
    display(ylabel!("fuerza reaccion Y"))

    # PLOT DE FUERZA EN Z
    plot(desplaz_iters, fuerza_reaccion_z_porcentaje_iters)
    xlabel!("desplazamiento")
    display(ylabel!("fuerza reaccion Z"))


    ########
    # energia_solo_dentro_total=sum(inc_energia_solo_dentro_barras_no_porcentajes)
    # Young_equiv=energia_solo_dentro_total/(1/2*volumen_hexaedro*(1.5e-6/X)^2)
    # println("inc_energia_no_porcentajes= ", energia_solo_dentro_total)
    # println("a ",Young_equiv)

    # ### MODULO DE YOUNG DEL SOLIDO CONTINUO EQUIVALENTE (solo la zona lineal)
    # # display(plot(desplaz_x,fuerza_x_total_iter))
    # # display(xlabel!("desplazamiento"))
    # # display(ylabel!("fuerza"))

    # deformacion_x=desplaz_x/X; tension_x=fuerza_x_total_iter/(Y*Z)
    # # display(plot(deformacion_x,tension_x))
    # # display(xlabel!("deformacion_x"))
    # # display(ylabel!("tension_x"))

    # Young_equivalente_iter=tension_x./deformacion_x   #en MPa
    # Young_equivalente=sum(Young_equivalente_iter)/length(Young_equivalente_iter)


    # ### COEFICIENTE DE POISSON XY DEL SOLIDO CONTINUO EQUIVALENTE
    # deformacion_y=(desplaz_y)/(Y/2)
    # Poisson_equivalente_xy_iter=deformacion_y./deformacion_x
    # Poisson_equivalente_xy=sum(Poisson_equivalente_xy_iter)/length(Poisson_equivalente_xy_iter)


    # ### COEFICIENTE DE POISSON XZ DEL SOLIDO CONTINUO EQUIVALENTE
    # deformacion_z=(desplaz_z)/(Z/2)
    # Poisson_equivalente_xz_iter=deformacion_z./deformacion_x
    # Poisson_equivalente_xz=sum(Poisson_equivalente_xz_iter)/length(Poisson_equivalente_xz_iter)

    # Poisson_equivalente=(Poisson_equivalente_xy+Poisson_equivalente_xz)/2


    # println("E= ",Young_equivalente);println("nu_xy= ",Poisson_equivalente_xy);println("nu_xz= ",Poisson_equivalente_xz)


    ### AÑADIR INFORMACION A LOS FICHEROS DE DATOS
    fichero_ensayo=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_metamaterial_energ_element_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_2=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_metamaterial_inc_energ_element_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_3=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_metamaterial_energ_total_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_4=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_metamaterial_inc_energ_total_$(lpad(numero_fichero,2,"0")).txt","a") 
    fichero_ensayo_5=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\COMPARACIONES\\comparaciones_metamaterial_fuerza_$(lpad(numero_fichero,2,"0")).txt","a") 
    

    ########################### COMPARAR CON EL MULTIESCALA
    if comparar_metamaterial_multiescala == "on"
        N=4
        N_elementos_solidos=N*N*N; 
        inc_energia_zona_metam_int=zeros(N_elementos_solidos);inc_energia_zona_metam_ext=zeros(N_elementos_solidos);inc_energia_zona_metam_exacta=zeros(N_elementos_solidos) ;energia_zona_metam_exacta=zeros(N_elementos_solidos)
        for elemento in 1:N_elementos_solidos
            
            #global  inc_energia_zona_metam_int, inc_energia_zona_metam_ext, inc_energia_zona_metam_exacta

            # for elemento_metam in conectiv_metam_hexaedros_reducido[elemento][:,3]       
            #     inc_energia_zona_metam_int[elemento] = inc_energia_zona_metam_int[elemento] + inc_energia_elementos[elemento_metam]
            # end
            for elemento_metam in conectiv_metam_hexaedros_extendido[elemento][:,3]       
                #inc_energia_zona_metam_ext[elemento] = inc_energia_zona_metam_ext[elemento] + inc_energia_elementos[elemento_metam]
                inc_energia_zona_metam_exacta[elemento] = inc_energia_zona_metam_exacta[elemento] + inc_energia_elementos[elemento_metam]*porcentajes_barras_dentro_hexaedros[elemento][elemento_metam]
                energia_zona_metam_exacta[elemento] = energia_zona_metam_exacta[elemento] + energias_acumuladas_elementos[elemento_metam]*porcentajes_barras_dentro_hexaedros[elemento][elemento_metam]
            end

            println(elemento, " ",inc_energia_zona_metam_int[elemento], " ",inc_energia_zona_metam_ext[elemento], " ", inc_energia_zona_metam_exacta[elemento], " ", energia_zona_metam_exacta[elemento])
        
            println(fichero_ensayo, energia_zona_metam_exacta[elemento])
            println(fichero_ensayo_2, inc_energia_zona_metam_exacta[elemento])
        end

        println(fichero_ensayo_3, sum(energia_zona_metam_exacta))
        println(fichero_ensayo_4, sum(inc_energia_zona_metam_exacta))

        println()
        println("total ", sum(inc_energia_zona_metam_int), " ", sum(inc_energia_zona_metam_ext), " ", sum(inc_energia_zona_metam_exacta), " ", sum(energia_zona_metam_exacta))
        println()


        println(fichero_ensayo_5, fuerza_reaccion_x_iters[end]," ", fuerza_reaccion_x_porcentaje_iters[end]," ",desplaz_iters[end] )
    else
    end

    close(fichero_ensayo)
    close(fichero_ensayo_2)
    close(fichero_ensayo_3)
    close(fichero_ensayo_4)
    close(fichero_ensayo_5)
end






##################
# fichero_ensayo=open(raw"C:\__UNIVERSIDAD__\AEROESPACIAL\4º GIA\TFG\comparaciones.txt","w") #crear el archivo para recoger datos del ensayo

# for i in 1:50
#     println(fichero_ensayo, fuerza_reaccion_x_iters[i]," ",desplaz_iters[i] )
# end

# close(fichero_ensayo)